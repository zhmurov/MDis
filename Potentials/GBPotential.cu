/*
 * GBPotential.cu
 *
 *  Created on: 14.05.2012
 *      Author: zhmurov
 */

#include "GBPotential.cuh"

#include "../Core/global.h"
#include "../Util/Log.h"

namespace gb_potential {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<GB Potential> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void create(){
	if(getYesNoParameter(PARAMETER_GB_ON, 0)){
		potential.compute = compute;
		potential.destroy = destroyPotential;
		sprintf(potential.name, "GB");
		potentials[potentialsCount++] = &potential;

		energyOutputGB.computeValues = &computeEnergyGB;
		allocateCPU((void**)&energyOutputGB.values, parameters.Ntr*sizeof(float));
		strcpy(energyOutputGB.name, ENERGY_OUTPUT_NAME_GB);
		energyOutputs[energyOutputsCount] = &energyOutputGB;
		energyOutputsCount ++;

		if(getYesNoParameter(PARAMETER_SA_ON, DEFAULT_SA_ON)){
			energyOutputSA.computeValues = &computeEnergySA;
			allocateCPU((void**)&energyOutputSA.values, parameters.Ntr*sizeof(float));
			strcpy(energyOutputSA.name, ENERGY_OUTPUT_NAME_SA);
			energyOutputs[energyOutputsCount] = &energyOutputSA;
			energyOutputsCount ++;
		}

		updater.update = update;
		updater.destroy = destroyUpdater;
		updater.frequency = getIntegerParameter(PARAMETER_PAIRS_FREQ);
		sprintf(updater.name, "GB Pairs");
		updaters[updatersCount++] = &updater;
		init();
	}
}

void init(){

	gbBlockSize = BLOCK_SIZE;
	gbBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;

	int i, j;
	int traj, itot, jtot;

	char gbModel[256];
	getMaskedParameter(gbModel, PARAMETER_GB_MODEL, DEFAULT_GB_MODEL);
	if(strncmp(gbModel, GB_MODEL_HCT_STRING, 3) == 0){
		LOG << "Using Hawkins-Cramer-Truhlar (HCT) GB model.";
		gbData.model = GB_MODEL_HCT;
	} else if(strncmp(gbModel, GB_MODEL_OBC_STRING, 3) == 0){
		LOG << "Using Onufriev-Bashford-Case (OBC) GB model.";
		gbData.model = GB_MODEL_OBC;
	} else {
		DIE("Unknown GB model. Only '%s' and '%s' are allowed", GB_MODEL_HCT_STRING, GB_MODEL_OBC_STRING);
	}

	gbData.ein = getFloatParameter(PARAMETER_GENBORN_EIN, DEFAULT_GENBORN_EIN);
	gbData.eout = getFloatParameter(PARAMETER_GENBORN_EOUT, DEFAULT_GENBORN_EOUT);
	float epsilon = 1.0f/gbData.ein - 1.0f/gbData.eout;
	gbData.halfCoulomb = COULOMB_CONSTANT*0.5f;
	gbData.epsilon = epsilon;
	gbData.epsilonTimesHalfCoulomb = epsilon*gbData.halfCoulomb;
	float ions = getFloatParameter(PARAMETER_GB_IONS, 0.0f); // Monovalent ion concentration, [M]
	gbData.kkappa = 0.73 * 3.16 * sqrtf(ions); // 1/nm

	allocateCPU((void**)&gbData.h_rho, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&gbData.d_rho, gsystem.Ntot*sizeof(float));

	allocateCPU((void**)&gbData.h_S, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&gbData.d_S, gsystem.Ntot*sizeof(float));

	allocateCPU((void**)&gbData.h_Srho, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&gbData.d_Srho, gsystem.Ntot*sizeof(float));

	allocateCPU((void**)&gbData.h_saMult, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&gbData.d_saMult, gsystem.Ntot*sizeof(float));

	char gbParFilename[256];
	getMaskedParameter(gbParFilename, PARAMETER_GB_PARAMETERS_FILE);
	FILE* file = fopen(gbParFilename, "r");
	gbParCount = 0;
	char buffer[BUF_SIZE];
	if(file != NULL){
		while(fgets(buffer, BUF_SIZE, file) != NULL){
			if(buffer[0] != '[' && buffer[0] != '\n' && buffer[0] != ';'){
				//printf("%s", buffer);
				gbParCount++;
			}
		}
		gbPar = (GBParameters*)calloc(gbParCount, sizeof(GBParameters));
		rewind(file);
		gbParCount = 0;
//		printf("# Atom type   R_min(A)   R_i(A)     p_i        sigma_i(kcal/molA)   Description\n");
		while(fgets(buffer, BUF_SIZE, file) != NULL){
			if(buffer[0] != '[' && buffer[0] != '\n' && buffer[0] != ';'){
				char* pch = strtok(buffer, " \t");
				strcpy(gbPar[gbParCount].atomType, pch);
				pch = strtok(NULL, " \t");
				gbPar[gbParCount].sar = atof(pch);
				pch = strtok(NULL, " \t");
				gbPar[gbParCount].st = atof(pch);
				pch = strtok(NULL, " \t");
				gbPar[gbParCount].pi = atof(pch);
				pch = strtok(NULL, " \t");
				gbPar[gbParCount].gbr = atof(pch);
				pch = strtok(NULL, " \t;\n\r");
				gbPar[gbParCount].hct = atof(pch);
				pch = strtok(NULL, ";\n\r");
/*				printf("%-*s%-*f%-*f%-*f%-*f# %s\n",
						14, gbPar[gbParCount].atomType,
						11, findNonbondedType(gbPar[gbParCount].atomType, &ff).RminOver2*10.0f,
						11, gbPar[gbParCount].gbr*10.0f,
						11, gbPar[gbParCount].pi,
						21, 0.003f,
						pch);*/
				gbParCount++;
			}
		}
		LOG << "GB Parameters read:";
		for(i = 0; i < gbParCount; i++){
			printf("%s\t%f\t%f\t%f\t%f\t%f\n",
					gbPar[i].atomType,
					gbPar[i].sar,
					gbPar[i].st,
					gbPar[i].pi,
					gbPar[i].gbr,
					gbPar[i].hct);
		}
		fclose(file);
	} else {
		DIE("GB parameters file \"%s\" not found.", gbParFilename);
	}
	gbData.dielOffset = getFloatParameter(PARAMETER_GB_DIELECTRIC_OFFSET, DEFAULT_GB_DIELECTRIC_OFFSET);

	gbData.saOn = getYesNoParameter(PARAMETER_SA_ON, DEFAULT_SA_ON);
	gbData.sigma = getFloatParameter(PARAMETER_SA_SIGMA, DEFAULT_SA_SIGMA);
	gbData.Rprobe = getFloatParameter(PARAMETER_SA_RPROBE, DEFAULT_SA_RPROBE);

	for(i = 0; i < gsystem.Ntot; i++){
		int found = 0;
		int ii = i % gsystem.N; // Atom ID in topology
		for(j = 0; j < gbParCount; j++){
			if(strcmp(topology.atoms[ii].type, gbPar[j].atomType) == 0){
				gbData.h_rho[i] = gbPar[j].gbr-gbData.dielOffset;
				gbData.h_S[i] = gbPar[j].hct;
				gbData.h_Srho[i] = gbData.h_S[i]*gbData.h_rho[i];
				float r = gbData.h_rho[i] + gbData.dielOffset;
				float r6 = r*r*r;
				r6 = r6*r6;
				gbData.h_saMult[i] = 6.0f*4.0f*M_PI*gbData.sigma*(r + gbData.Rprobe)*(r + gbData.Rprobe)*r6;
				/*
				printf("%d\t%s\t%d\t%s\t%f\t%f\t%f\n",
						i, topology.atoms[ii].resName, topology.atoms[ii].resid, topology.atoms[ii].name,
						gbData.h_rho[i], gbData.h_S[i], gbData.h_Srho[i]);
				*/
				found = 1;
			}
		}
		if(!found){
			DIE("GB parameters for atom type %s have not been found.", topology.atoms[ii].type);
		}
	}

	cudaMemcpy(gbData.d_rho, gbData.h_rho, gsystem.Ntot*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(gbData.d_S, gbData.h_S, gsystem.Ntot*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(gbData.d_Srho, gbData.h_Srho, gsystem.Ntot*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(gbData.d_saMult, gbData.h_saMult, gsystem.Ntot*sizeof(float), cudaMemcpyHostToDevice);

	allocateCPU((void**)&gbData.h_alpha, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&gbData.d_alpha, gsystem.Ntot*sizeof(float));

	for(i = 0; i < gsystem.Ntot; i++){
		gbData.h_alpha[i] = 0.0f;
	}
	cudaMemcpy(gbData.d_alpha, gbData.h_alpha, gsystem.Ntot*sizeof(float), cudaMemcpyHostToDevice);

	allocateCPU((void**)&gbData.h_gbEnergies, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&gbData.d_gbEnergies, gsystem.Ntot*sizeof(float));

	allocateCPU((void**)&gbData.h_saEnergies, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&gbData.d_saEnergies, gsystem.Ntot*sizeof(float));

	allocateCPU((void**)&gbData.h_dGda, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&gbData.d_dGda, gsystem.Ntot*sizeof(float));

	allocateCPU((void**)&gbData.h_dGdr, gsystem.Ntot*sizeof(float4));
	allocateGPU((void**)&gbData.d_dGdr, gsystem.Ntot*sizeof(float4));

	cudaBindTexture(0, t_gbalpha, gbData.d_alpha,
					gsystem.Ntot*sizeof(float));

	if(gbData.model == GB_MODEL_OBC){
		gbData.alpha = getFloatParameter(PARAMETER_GB_OBC_ALPHA, DEFAULT_GB_OBC_ALPHA);
		gbData.beta = getFloatParameter(PARAMETER_GB_OBC_BETA, DEFAULT_GB_OBC_BETA);
		gbData.gamma = getFloatParameter(PARAMETER_GB_OBC_GAMMA, DEFAULT_GB_OBC_GAMMA);
		allocateCPU((void**)&gbData.h_dHdS, gsystem.Ntot*sizeof(float));
		allocateGPU((void**)&gbData.d_dHdS, gsystem.Ntot*sizeof(float));
	}


	allocateCPU((void**)&gbData.h_pairsCount, gsystem.Ntot*sizeof(int));
	allocateGPU((void**)&gbData.d_pairsCount, gsystem.Ntot*sizeof(int));

	gbData.pairsCutoff = getFloatParameter(PARAMETER_PAIRS_CUTOFF_NB, DEFAULT_PAIRS_CUTOFF_NB);
	gbData.pairsCutoff2 = gbData.pairsCutoff*gbData.pairsCutoff;

	pairlistsExtension = getIntegerParameter(PARAMETER_PAIRSLISTS_EXTENSION, DEFAULT_PAIRSLISTS_EXTENSION);

	int maxPairs = getIntegerParameter(PARAMETER_MAX_PAIRS_GB, 0);
	float4 ri, rj, dr;
	if (maxPairs == 0){
		LOG << "Counting pairs for GB potential...";
		for(i = 0; i < gsystem.Ntot; i++){
			gbData.h_pairsCount[i] = 0;
		}
		for(traj = 0; traj < parameters.Ntr; traj++){
			for(i = 0; i < gsystem.N; i++){
				itot = traj*gsystem.N + i;
				ri = gsystem.h_coord[itot];
				for(j = 0; j < gsystem.N; j++){
					if(i != j){
						jtot = traj*gsystem.N + j;
						rj = gsystem.h_coord[jtot];
						dr.x = rj.x - ri.x;
						dr.y = rj.y - ri.y;
						dr.z = rj.z - ri.z;
						float dr2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
						if(dr2 < gbData.pairsCutoff2){
							gbData.h_pairsCount[itot]++;
						}
					}
				}
			}
		}
		for(i = 0; i < gsystem.Ntot; i++){
			if(maxPairs < gbData.h_pairsCount[i]){
				maxPairs = gbData.h_pairsCount[i];
			}
		}
		LOG << "Maximum GB pairs found " << maxPairs;
	}else{
		LOG << "Using explicitly specified GB pairs count = " << maxPairs;
	}
	LOG << "Adding " << pairlistsExtension << " to the length of GB lists to avoid memory conflicts.";
	maxPairs += pairlistsExtension;

	allocateCPU((void**)&gbData.h_pairs, maxPairs*gsystem.widthTot*sizeof(int));
	allocateGPU((void**)&gbData.d_pairs, maxPairs*gsystem.widthTot*sizeof(int));

	cudaMemcpyToSymbol(c_gbData, &gbData, sizeof(GBData), 0, cudaMemcpyHostToDevice);
}


__device__ float Lij(float rij, float Srhoj, float rhoi){
	if(rij+Srhoj <= rhoi){
		return 1;
	} else if(rij-Srhoj <= rhoi){
		return rhoi;
	} else {
		return rij-Srhoj;
	}
}

__device__ float Uij(float rij, float Srhoj, float rhoi){
	if(rij+Srhoj <= rhoi){
		return 1;
	} else {
		return rij+Srhoj;
	}
}

__device__ float Hij(float rij, float Srhoj, float rhoi){
	if(rij+Srhoj <= rhoi){
		return 0;
	} else {
		float lij;
		float uij = rij+Srhoj;
		if(rij-Srhoj <= rhoi){
			lij = rhoi;
		} else {
			lij = rij-Srhoj;
		}
		float rijinv = 1.0f/rij;
		float lijinv = 1.0f/lij;
		float uijinv = 1.0f/uij;
		float lminu = lijinv - uijinv;
		float hij = lijinv + uijinv;
		hij *=  Srhoj*Srhoj*rijinv - rij;
		hij *= 0.25f;
		hij += 1.0f;
		hij *= lminu;
		hij += 0.5f*rijinv*logf(lij*uijinv);
		return hij;
	}
}

__device__ float dHij(float rij, float Srhoj, float rhoi){
	if(rij+Srhoj <= rhoi){
		return 0;
	} else {
		float lij;
		float uij = rij+Srhoj;
		float dlij;
		float duij = 1.0f;
		if(rij-Srhoj <= rhoi){
			lij = rhoi;
			dlij = 0.0f;
		} else {
			lij = rij-Srhoj;
			dlij = 1.0f;
		}
		float rijinv = 1.0f/rij;
		float lijinv = 1.0f/lij;
		float uijinv = 1.0f/uij;

		float lijinv2 = lijinv*lijinv;
		float uijinv2 = uijinv*uijinv;


		float ul0 = uijinv2 - lijinv2;
		float ul1 = rijinv*logf(lijinv*uij);
		ul1 += dlij*lijinv - duij*uijinv;
		ul1 *= 0.5f*rijinv;

		float ul2 = -dlij*lijinv2 + duij*uijinv2;
		float ul3 = dlij*lijinv2*lijinv - duij*uijinv2*uijinv;
		float dhij = ul1 + ul2;
		float Srhoj2rinv = Srhoj*Srhoj*rijinv;
		dhij += 0.25f*ul0*(1.0f + Srhoj2rinv*rijinv);
		dhij += 0.5f*ul3*(rij - Srhoj2rinv);
		return dhij;
	}
}

__global__ void gbRadii_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		int j;
		float rhoi = c_gbData.d_rho[d_i];
		float4 r1 = tex1Dfetch(t_coord, d_i);
		float alpha = 0.0f;
		int p;
		for(p = 0; p < c_gbData.d_pairsCount[d_i]; p++){
			j = c_gbData.d_pairs[p*c_gsystem.widthTot + d_i];
			float4 r2 = tex1Dfetch(t_coord, j);
			float Sjrhoj = c_gbData.d_Srho[j];
			r2.x -= r1.x;
			r2.y -= r1.y;
			r2.z -= r1.z;
			float rij2 = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			float rij = sqrtf(rij2);
			float hij = Hij(rij, Sjrhoj, rhoi);
			alpha += hij;
		}
		if(c_gbData.model == GB_MODEL_HCT){
			alpha = 1.0f/rhoi - 0.5f*alpha;
			alpha = 1.0f/alpha;
			float minalpha = rhoi + c_gbData.dielOffset;
			if(alpha < minalpha){
				alpha = minalpha;
			}
		} else if(c_gbData.model == GB_MODEL_OBC){
			float I = 0.5f*alpha*rhoi;
			float I2 = I*I;
			float arg = (c_gbData.alpha - c_gbData.beta*I + c_gbData.gamma*I2)*I;
			float tangent = tanhf(arg);
			float dHdS = 1.0f - tangent*tangent;
			dHdS *= rhoi / (rhoi + c_gbData.dielOffset);
			dHdS *= 0.5f*c_gbData.alpha - c_gbData.beta*I + 1.5f*c_gbData.gamma*I2;
			alpha = tangent / (rhoi + c_gbData.dielOffset);
			alpha = 1.0f/rhoi - alpha;
			alpha = 1.0f/alpha;
			dHdS *= alpha*alpha;
			c_gbData.d_dHdS[d_i] = dHdS;
		}
		c_gbData.d_alpha[d_i] = alpha;
	}
}

__global__ void gbdGda_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 r = c_gsystem.d_coord[d_i];//tex1Dfetch(t_coord, d_i);
		float4 f = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
		float4 r2;
		int at2;
		float ai = c_gbData.d_alpha[d_i];//tex1Dfetch(t_alpha, d_i);
		int j;
		float dGda = 0.0f;
		float mult, arg, mult2;
		int p;
		for(p = 0; p < c_gbData.d_pairsCount[d_i]; p++){
			j = c_gbData.d_pairs[p*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, j);
			at2 = (int)r2.w;
			r2.x -= r.x;
			r2.y -= r.y;
			r2.z -= r.z;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;

			arg = 0.25f*r2.w;
			float aj = tex1Dfetch(t_gbalpha, j);
			arg /= ai*aj;
			float expon = __expf(-arg);

			mult = r2.w + ai*aj*expon;
			mult = 1.0f/mult;
			mult = sqrtf(mult);
			mult = mult*mult*mult;
			mult *= tex1Dfetch(t_charges, at2);

			mult *= c_gbData.epsilonTimesHalfCoulomb;
			mult2 = mult;

			mult *= 1.0f+arg;
			mult *= expon;

			dGda += mult*aj;

			mult2 *= 2.0f*(1.0f-0.25f*expon);

			f.x += mult2*r2.x;
			f.y += mult2*r2.y;
			f.z += mult2*r2.z;
		}
		float qi = tex1Dfetch(t_charges, (int)r.w);
		dGda += c_gbData.epsilonTimesHalfCoulomb*qi/(ai*ai);
		dGda *= qi;
		c_gbData.d_dGda[d_i] = dGda;
		f.x *= qi;
		f.y *= qi;
		f.z *= qi;
		c_gbData.d_dGdr[d_i] = f;
	}
}

__global__ void gbdGda_ion_kernel(){ // On {nvcc 4.2 arch sm_20} uses 27 register, non-ionic kernel uses 23
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 r = c_gsystem.d_coord[d_i];//tex1Dfetch(t_coord, d_i);
		float4 f = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
		float4 r2;
		int at2;
		float ai = c_gbData.d_alpha[d_i];//tex1Dfetch(t_alpha, d_i);
		int j;
		float dGda = 0.0f;
		float mult, arg, mult2;
		int p;
		for(p = 0; p < c_gbData.d_pairsCount[d_i]; p++){
			j = c_gbData.d_pairs[p*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, j);
			at2 = (int)r2.w;
			r2.x -= r.x;
			r2.y -= r.y;
			r2.z -= r.z;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;

			arg = 0.25f*r2.w;
			float aj = tex1Dfetch(t_gbalpha, j);
			arg /= ai*aj;
			float expon = __expf(-arg);

			mult = r2.w + ai*aj*expon;
			mult = sqrtf(mult); // mult = f_{GB}
			float kappaf = - c_gbData.kkappa * mult;
			float dpref = (1.0f + kappaf)*__expf(kappaf);
			
			mult = 1.0f/mult;
			mult = mult*mult*mult;
			mult *= tex1Dfetch(t_charges, at2);
			// mult = q2 / f_{GB}^3
			const float pref = 1.0f/c_gbData.ein - dpref/c_gbData.eout; // Prefactor
			mult *= pref*c_gbData.halfCoulomb;
			
			mult2 = mult;

			mult *= 1.0f+arg;
			mult *= expon;
			// mult = (1+gamma*r^2/ai/aj)*exp(-gamma*r^2/ai/aj) * q2 / f_{GB}^3 = (q2/aj/f_{GB}^2) * df_{GB}/dai

			dGda += mult*aj;
			mult2 *= 2.0f*(1.0f-0.25f*expon);
			// mult2 = 2*(1-gamma*exp(-gamma*r^2/ai/aj)) * q2 / f_{GB}^3

			f.x += mult2*r2.x;
			f.y += mult2*r2.y;
			f.z += mult2*r2.z;
		}
		float qi = tex1Dfetch(t_charges, (int)r.w);
		float kappaf = - ai * c_gbData.kkappa;
		float pref = 1.0f / c_gbData.ein - (1.0f + kappaf)*__expf(kappaf) / c_gbData.eout;
		dGda += c_gbData.halfCoulomb*pref*qi/(ai*ai);
		dGda *= qi;
		c_gbData.d_dGda[d_i] = dGda;
		f.x *= qi;
		f.y *= qi;
		f.z *= qi;
		c_gbData.d_dGdr[d_i] = f;
	}
}

__global__ void gb_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float rhoi = c_gbData.d_rho[d_i];
		float Sirhoi = c_gbData.d_Srho[d_i];
		float4 r1 = c_gsystem.d_coord[d_i];//tex1Dfetch(t_coord, d_i);
		float3 f;
		float4 r2;
		float ai = tex1Dfetch(t_gbalpha, d_i);
		float dGdai = c_gbData.d_dGda[d_i];
		float dHdS;
		if(c_gbData.model == GB_MODEL_HCT){
			dHdS = 0.5f*ai*ai;
		} else if(c_gbData.model == GB_MODEL_OBC){
			dHdS = c_gbData.d_dHdS[d_i];
		}
		dGdai *= dHdS;
		float dGdaj;
		f.x = 0.0f;
		f.y = 0.0f;
		f.z = 0.0f;
		int j;
		int p;
		float dGSAdai=0.0f;
		float invai7;
		if(c_gbData.saOn){
			dGSAdai = c_gbData.d_saMult[d_i];
			invai7 = powf(ai, -7);
			/*invai7 = 1.0f/ai;
			invai7 = invai7*invai7;
			invai7 = invai7*invai7*invai7;
			invai7 = invai7/ai;*/
			dGSAdai *= invai7;
			dGSAdai *= dHdS;
		}
		for(p = 0; p < c_gbData.d_pairsCount[d_i]; p++){
			j = c_gbData.d_pairs[p*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, j);
			float Sjrhoj = c_gbData.d_Srho[j];
			float rhoj = c_gbData.d_rho[j];
			r2.x -= r1.x;
			r2.y -= r1.y;
			r2.z -= r1.z;
			float rij2 = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			float rij = sqrtf(rij2);

			float aj = tex1Dfetch(t_gbalpha, j);
			//float dai = (Hij(rij+0.0001f, Sjrhoj, rhoi)-Hij(rij-0.0001f, Sjrhoj, rhoi))/0.0002f;
			//float daj = (Hij(rij+0.0001f, Sirhoi, rhoj)-Hij(rij-0.0001f, Sirhoi, rhoj))/0.0002f;
			float dai = dHij(rij, Sjrhoj, rhoi);
			float daj = dHij(rij, Sirhoi, rhoj);

			dGdaj = c_gbData.d_dGda[j];
			float dHdS;
			if(c_gbData.model == GB_MODEL_HCT){
				dHdS = 0.5f*aj*aj;
			} else if(c_gbData.model == GB_MODEL_OBC){
				dHdS = c_gbData.d_dHdS[j];
			}

			dGdaj *= dHdS;

			float mult;
			float dGSAdaj = 0.0f;
			if(c_gbData.saOn){
				dGSAdaj = c_gbData.d_saMult[j];
				float invaj7 = powf(aj, -7);
				/*invai7 = 1.0f/aj;
				invaj7 = invaj7*invaj7;
				invaj7 = invaj7*invaj7*invaj7;
				invaj7 = invaj7/aj;*/
				dGSAdaj *= invaj7;
				dGSAdaj *= dHdS;
				mult = ((dGdai-dGSAdai)*dai + daj*(dGdaj-dGSAdaj))/rij;
			} else {
				mult = (dGdai*dai + daj*dGdaj)/rij;
			}

			f.x += mult*r2.x;
			f.y += mult*r2.y;
			f.z += mult*r2.z;
		}
		float4 dGdr = c_gbData.d_dGdr[d_i];

		f.x += dGdr.x;
		f.y += dGdr.y;
		f.z += dGdr.z;

		r2 = c_gsystem.d_forces[d_i];
		//r2.x = 0; r2.y = 0; r2.z = 0;
		r2.x += f.x;
		r2.y += f.y;
		r2.z += f.z;
		c_gsystem.d_forces[d_i] = r2;
	}
}

void compute(){
	gbRadii_kernel<<<gbBlockCount, gbBlockSize>>>();
	checkCUDAError("GB Radii");
	if (gbData.kkappa == 0.0f) {
		gbdGda_kernel<<<gbBlockCount, gbBlockSize>>>();
		checkCUDAError("GB dGda");
	} else {
		gbdGda_ion_kernel<<<gbBlockCount, gbBlockSize>>>();
		checkCUDAError("GB dGda (ion)");
	}
	gb_kernel<<<gbBlockCount, gbBlockSize>>>();
	checkCUDAError("GB Forces");
	//cudaMemcpy(gbData.h_alpha, gbData.d_alpha, gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	/*if(gbData.model == GB_MODEL_OBC){
		cudaMemcpy(gbData.h_dfdS, gbData.d_dfdS, gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	}
	cudaMemcpy(gbData.h_gbEnergies, gbData.d_gbEnergies, gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);*/
	//cudaMemcpy(gsystem.h_forces, gsystem.d_forces, gsystem.Ntot*sizeof(float4), cudaMemcpyDeviceToHost);
	/*float pot = 0.0f;
	*/
	/*int traj, i;
	FILE* file = fopen("forces_mdis.dat", "w");
	for(traj = 0; traj < parameters.Ntr; traj++){
		for(i = 0; i < gsystem.N; i++){*/
			/*printf("%d\t%s\t%d\t%s\t%f\t%f\n",
					i, topology.atoms[i].resName, topology.atoms[i].resid, topology.atoms[i].name,
					gbData.h_alpha[i], gbData.h_gbEnergies[i]);*/
		//	if(gbData.model == GB_MODEL_HCT){
				//fprintf(file, "%d\t%f\t%f\t%f\t%f\t%f\n",
					//	i, gbData.h_alpha[i], gbData.h_gbEnergies[i], gsystem.h_forces[i].x, gsystem.h_forces[i].y, gsystem.h_forces[i].z);
				/*fprintf(file, "%d\t%f\t%f\t%f\t%f\n",
					i,  gbData.h_alpha[i], gsystem.h_forces[i].x, gsystem.h_forces[i].y, gsystem.h_forces[i].z);*/
		//	}
			/*if(gbData.model == GB_MODEL_OBC){
				fprintf(file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n",
						i, gbData.h_alpha[i], gbData.h_gbEnergies[i], gsystem.h_forces[i].x, gsystem.h_forces[i].y, gsystem.h_forces[i].z, gbData.h_dfdS[i]);
			}*/
			//pot += gbData.h_gbEnergies[i];
/*		}
	}
	fclose(file);
	//pot = -COULOMB_CONSTANT*0.5f*(1.0f - 1.0f/80.0f)*pot;//*KCALL_PER_KJ;
	//printf("Total energy: %f\n", pot);

	exit(0);*/
	/*cudaMemcpy(gbData.h_alpha, gbData.d_alpha, gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	float totalSASA = 0.0f;
	float sasaEnergy = 0.0f;
	FILE* file = fopen("gbsasa.dat", "w");
	int i;
	for(i = 0; i < gsystem.Ntot; i++){
		//if(topology.atoms[i].name[0] != 'H'){
			float r = gbData.h_rho[i] + gbData.dielOffset;
			float mult = r/gbData.h_alpha[i];
			mult = mult*mult;
			mult = mult*mult*mult;
			float sasa = 4.0f*M_PI*(r + 0.14f)*(r + 0.14f)*mult;
			fprintf(file, "%d\t%f\n", i, sasa);
			totalSASA += sasa;
			sasaEnergy += sasa*0.0054f*100.0f/KCALL_PER_KJ;
		//}
	}
	fclose(file);
	printf("Total SASA (GB): %f; Energy: %f\n", totalSASA, sasaEnergy);
	exit(0);*/
}

__global__ void gbEnergy_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float pot = 0.0f;
		float4 r1 = tex1Dfetch(t_coord, d_i);
		float4 r2;
		int at1, at2;
		at1 = (int)r1.w;
		float ai = c_gbData.d_alpha[d_i];
		int j;
		int p;
		for(p = 0; p < c_gbData.d_pairsCount[d_i]; p++){
			j = c_gbData.d_pairs[p*c_gsystem.widthTot + d_i];
		//for(j = 0; j < c_gsystem.Ntot; j++){
			r2 =  tex1Dfetch(t_coord, j);
			at2 = (int)r2.w;
			r2.x -= r1.x;
			r2.y -= r1.y;
			r2.z -= r1.z;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			float arg = -0.25f*r2.w;
			float aj = c_gbData.d_alpha[j];
			arg /= ai*aj;
			float expon = expf(arg);
			float mult;
			mult = r2.w + ai*aj*expon;
			mult = sqrtf(mult);

			mult = 1.0f/mult;
			if (c_gbData.kkappa != 0.0f) {
				float pref = __expf(- c_gbData.kkappa / mult);
				mult *= 1.0f / c_gbData.ein - pref / c_gbData.eout;
			}
			mult *= tex1Dfetch(t_charges, at2);
			pot += mult;
		}
		float qi = tex1Dfetch(t_charges, at1);
		if (c_gbData.kkappa != 0.0f) {
			float pref = __expf(- c_gbData.kkappa * ai);
			qi *= 1.0f / c_gbData.ein - pref / c_gbData.eout;
		}
		pot += qi/ai;
		pot *= qi;
		c_gbData.d_gbEnergies[d_i] = pot;
	}
}

inline void computeEnergyGB(){
	gbRadii_kernel<<<gbBlockCount, gbBlockSize>>>();
	checkCUDAError("GB Radii");
	gbEnergy_kernel<<<gbBlockCount, gbBlockSize>>>();
	checkCUDAError("GB Energy");

	cudaMemcpy(gbData.h_gbEnergies, gbData.d_gbEnergies,
			gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	checkCUDAError("GB Energy - MemCpy");
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		float pot = 0.0f;
		for(i = 0; i < gsystem.N; i++){
			/*if(nonBondedData.h_ljEnergies[i] > 0.0f){
				printf("Atom %d has another one close to it (LJ Energy = %f)\n", i, nonBondedData.h_ljEnergies[i]);
			}*/
			//if(topology.atoms[i].type[0] == 'H'){
			pot += gbData.h_gbEnergies[i + traj*gsystem.N];
			//}
		}
		//pot /= 2.0f;
		if (gbData.kkappa != 0.0f)
			energyOutputGB.values[traj] = -pot*gbData.halfCoulomb;
		else
			energyOutputGB.values[traj] = -pot*gbData.epsilonTimesHalfCoulomb;
	}
}

__global__ void gbEnergySA_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		//float r = c_gbData.d_rho[d_i] + c_gbData.dielOffset;
		float mult = 1.0/c_gbData.d_alpha[d_i];
		mult = mult*mult;
		mult = mult*mult*mult;
		c_gbData.d_saEnergies[d_i] = c_gbData.d_saMult[d_i]*mult/6.0f;
	}
}

inline void computeEnergySA(){
	gbEnergySA_kernel<<<gbBlockCount, gbBlockSize>>>();
	checkCUDAError("SA Energy");

	cudaMemcpy(gbData.h_saEnergies, gbData.d_saEnergies,
			gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	checkCUDAError("SA Energy - MemCpy");
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		float pot = 0.0f;
		for(i = 0; i < gsystem.N; i++){
			/*if(nonBondedData.h_ljEnergies[i] > 0.0f){
				printf("Atom %d has another one close to it (LJ Energy = %f)\n", i, nonBondedData.h_ljEnergies[i]);
			}*/
			//if(topology.atoms[i].type[0] == 'H'){
			pot += gbData.h_saEnergies[i + traj*gsystem.N];
			//}
		}
		//pot /= 2.0f;
		energyOutputSA.values[traj] = pot;
	}
}

__global__ void updatePairsLists_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		int i;
		float4 r1 = tex1Dfetch(t_coord, d_i);
		float4 r2;
		int pairsCount = 0;
		int traj = d_i/c_gsystem.N;
		for(i = traj*c_gsystem.N; i < (traj + 1)*c_gsystem.N; i++){
			if(i != d_i){
				r2 = tex1Dfetch(t_coord, i);
				r2 -= r1;
				DO_PBC(r2);
				r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
				if(r2.w < c_gbData.pairsCutoff2){
					c_gbData.d_pairs[pairsCount*c_gsystem.widthTot + d_i] = i;
					pairsCount++;
				}
			}
		}
		c_gbData.d_pairsCount[d_i] = pairsCount;
	}
}

void update(){
	updatePairsLists_kernel<<<gbBlockCount, gbBlockSize>>>();
}

void destroyPotential(){

}

void destroyUpdater(){

}

#undef LOG

}
