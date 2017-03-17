/*
 * GenBornPotential.cu
 *
 *  Created on: Jan 10, 2011
 *      Author: zhmurov
 */

#include "../Core/global.h"
#include "../Core/md.cuh"
#include "../Util/Log.h"
#include "GenBornPotential.cuh"
#include "../Updaters/PairsListsUpdater.cuh"
#include "../Updaters/RestartOutputManager.cuh"

#define BUF_SIZE 256

namespace genborn_potential {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<genborn_potential> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void initGenBornParameter(char* name, char* value);

void create(){
	if (!getYesNoParameter(PARAMETER_GENBORN_ON, 0))
		return;
	genBornPotential.compute = &compute;
	genBornPotential.destroy = &destroy;
	sprintf(genBornPotential.name, "GenBorn potential");
	potentials[potentialsCount] = &genBornPotential;
	potentialsCount ++;
	genBornEnergyOutput.computeValues = &computeEnergy;
	allocateCPU((void**)&genBornEnergyOutput.values, parameters.Ntr*sizeof(float));
	strcpy(genBornEnergyOutput.name, ENERGY_OUTPUT_NAME_GEN_BORN);
	energyOutputs[energyOutputsCount] = &genBornEnergyOutput;
	energyOutputsCount ++;

	init();
}

void init(){
	LOG << "Initializing Generalized Born 'GenBorn' potential...";

	genBornBlockSize = BLOCK_SIZE;
	genBornBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;

//	initNonBondedPotential();
	int i;

	allocateCPU((void**)&genBornData.h_genBornEnergies, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&genBornData.d_genBornEnergies, gsystem.Ntot*sizeof(float));

	allocateCPU((void**)&genBornData.h_alpha, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&genBornData.d_alpha, gsystem.Ntot*sizeof(float));

	allocateCPU((void**)&genBornData.h_dGda, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&genBornData.d_dGda, gsystem.Ntot*sizeof(float));

	allocateCPU((void**)&genBornData.h_dGdr, gsystem.Ntot*sizeof(float4));
	allocateGPU((void**)&genBornData.d_dGdr, gsystem.Ntot*sizeof(float4));

	allocateCPU((void**)&genBornData.h_RVdW, atomTypesCount*sizeof(float));
	allocateGPU((void**)&genBornData.d_RVdW, atomTypesCount*sizeof(float));

	allocateCPU((void**)&genBornData.h_VVdW, atomTypesCount*sizeof(float));
	allocateGPU((void**)&genBornData.d_VVdW, atomTypesCount*sizeof(float));

	allocateCPU((void**)&genBornData.h_LP1, atomTypesCount*sizeof(float));
	allocateGPU((void**)&genBornData.d_LP1, atomTypesCount*sizeof(float));

	char filename[100];
	getMaskedParameter(filename, PARAMETER_GENBORN_PARAMETERS_FILE);
	LOG << "Reading GenBorn parameters from '" << filename << "'";
	FILE* gbFile = safe_fopen(filename, "r");
	char buffer[BUF_SIZE];
	DPRINTF("GenBorn parameters:\n");
	while(fgets(buffer, BUF_SIZE, gbFile) != NULL){
		if(buffer[0] != '#' && buffer[0] != ' ' && buffer[0] != '\n' && buffer[0] != '\t'){
			char* name = strtok(buffer, " \t");
			char* value = strtok(NULL, " \t");
			initGenBornParameter(name, value);
		}
	}
	fclose(gbFile);

	float ein = getFloatParameter(PARAMETER_GENBORN_EIN, DEFAULT_GENBORN_EIN);
	float eout = getFloatParameter(PARAMETER_GENBORN_EOUT, DEFAULT_GENBORN_EOUT);
	genBornData.epsilon = 1.0f/ein - 1.0f/eout;
	genBornData.halfCoulomb = COULOMB_CONSTANT*0.5f;
	genBornData.epsilonTimesHalfCoulomb = genBornData.epsilon*genBornData.halfCoulomb;
	genBornData.cutAlpha = getFloatParameter(PARAMETER_GENBORN_CUTALPHA, DEFAULT_GENBORN_CUTALPHA);

	//Converting to MD units
	genBornData.P1 *= ANGSTROM_TO_NM;
	genBornData.P2 *= ANGSTROM_TO_NM*KCAL_TO_KJOULE;
	genBornData.P3 *= ANGSTROM_TO_NM*KCAL_TO_KJOULE;
	genBornData.P4 *= ANGSTROM_TO_NM*KCAL_TO_KJOULE;
	//Computing values for Lambda' and P1' - P5'
	if(genBornData.Lambda != 0){
		LOG << "Using original Still definition of GB self-energy";
		genBornData.Lambda = -COULOMB_CONSTANT/(2.0f*genBornData.Lambda);
		genBornData.P1 = COULOMB_CONSTANT*genBornData.P1/2.0f;
	} else {
		LOG << "Using linearized definition of Still GB self-energy.";
		genBornData.Phi *= ANGSTROM_TO_NM;
	}

	LOG << "Computing VdW radii/volumes for all atoms.";
	DPRINTF("Atom type:\tR_VdW:\t\tV_VdW:\t\tL+P:\n");
	for(i = 0; i < atomTypesCount; i++){
		genBornData.h_RVdW[i] = atomTypes[i].RminOver2;
		if(genBornData.h_RVdW[i] < 0.08f){
			genBornData.h_RVdW[i] = 0.08f;
		}
		genBornData.h_VVdW[i] = (4.0f/3.0f)*M_PI*genBornData.h_RVdW[i]*genBornData.h_RVdW[i]*genBornData.h_RVdW[i];
		if(genBornData.Lambda != 0){
			genBornData.h_LP1[i] = (genBornData.Lambda/genBornData.h_RVdW[i]) +
					(genBornData.P1/(genBornData.h_RVdW[i]*genBornData.h_RVdW[i]));
		} else {
			genBornData.h_LP1[i] = -0.5f*COULOMB_CONSTANT/(genBornData.h_RVdW[i] + genBornData.Phi + genBornData.P1);
		}
		DPRINTF("%d.%s\t\t%5.4f\t\t%5.4f\t\t%5.4f\n", i, atomTypes[i].name,
				genBornData.h_RVdW[i], genBornData.h_VVdW[i], genBornData.h_LP1[i]);
	}

	float roff = getFloatParameter(PARAMETER_NB_CUTOFF);
	float ron = getFloatParameter(PARAMETER_NB_SWITCH);
	float roff2 = roff*roff;
	float ron2 = ron*ron;
	float oneOverRoffMinRon = 1.0f/(roff2-ron2);
	float oneOverRoffMinRonCub = oneOverRoffMinRon*oneOverRoffMinRon*oneOverRoffMinRon;

	genBornData.roff = roff;
	genBornData.ron = ron;
	genBornData.roff2 = roff2;
	genBornData.ron2 = ron2;

	genBornData.C1sw = 2.0f*oneOverRoffMinRonCub;
	genBornData.C2sw = roff2;
	genBornData.C3sw = (roff2-3.0f*ron2)/2.0f;

	genBornData.C1dsw = 12.0*oneOverRoffMinRonCub;

	cudaMemcpy(genBornData.d_RVdW, genBornData.h_RVdW,
			atomTypesCount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(genBornData.d_VVdW, genBornData.h_VVdW,
			atomTypesCount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(genBornData.d_LP1, genBornData.h_LP1,
			atomTypesCount*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpyToSymbol(c_RVdW, genBornData.h_RVdW,
			atomTypesCount*sizeof(float), 0);
	cudaMemcpyToSymbol(c_VVdW, genBornData.h_VVdW,
				atomTypesCount*sizeof(float), 0);
	cudaMemcpyToSymbol(c_LP1, genBornData.h_LP1,
				atomTypesCount*sizeof(float), 0);

	cudaBindTexture(0, t_RVdW, genBornData.d_RVdW,
			atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_VVdW, genBornData.d_VVdW,
				atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_LP1, genBornData.d_LP1,
				atomTypesCount*sizeof(float));

	cudaBindTexture(0, t_alpha, genBornData.d_alpha,
				gsystem.Ntot*sizeof(float));


	cudaMemcpyToSymbol(c_genBornData, &genBornData,
				sizeof(GenBornData), 0, cudaMemcpyHostToDevice);

	LOG << "Done initializing 'GenBorn' potential.";
}

__device__ float swfunct(float rij2){
	if(rij2 > c_genBornData.roff2){
		return 0.0f;
	} else
	if(rij2 < c_genBornData.ron2){
		return 1.0f;
	} else {
		float result = c_genBornData.C1sw*(c_genBornData.C2sw - rij2)*(c_genBornData.C2sw - rij2)*(c_genBornData.C3sw + rij2);
		return result;
	}
}

__device__ float dswfunct(float rij2){
	if(rij2 > c_genBornData.roff2){
		return 0.0f;
	} else
	if(rij2 < c_genBornData.ron2){
		return 0.0f;
	} else {
		return c_genBornData.C1dsw*(rij2-c_genBornData.ron2)*(rij2-c_genBornData.roff2);
	}
}

__global__ void genBornRadii_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		s_Ri[threadIdx.x] = c_gsystem.d_coord[d_i];//tex1Dfetch(t_coord, d_i);
		float4 r2;
		int at1, at2;
		at1 = (int)s_Ri[threadIdx.x].w;
		/*s_f[threadIdx.x].x = 0.0f;
		s_f[threadIdx.x].y = 0.0f;
		s_f[threadIdx.x].z = 0.0f;
		s_f[threadIdx.x].w = c_LP1[at1];*/
		float alpha = tex1Dfetch(t_LP1, at1);
		float mult;
		float arg;
		float CCF;
		int i;
		//int pair = 0;
		for(i = 0; i < c_pairsListsData.d_pairs12ListCount[d_i]; i++){
			r2 =  tex1Dfetch(t_coord, c_pairsListsData.d_pairs12List[i*c_gsystem.widthTot + d_i]);
			at2 = (int)r2.w;
			r2.x -= s_Ri[threadIdx.x].x;
			r2.y -= s_Ri[threadIdx.x].y;
			r2.z -= s_Ri[threadIdx.x].z;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			r2.w = 1.0f/r2.w;
			mult = c_genBornData.P2;
			mult *= r2.w;
			mult *= r2.w;
			//s_f[threadIdx.x].w += mult*c_VVdW[at2];
			alpha += mult*tex1Dfetch(t_VVdW, at2);//c_VVdW[at2];
			//mult *= -4.0f*r2.w;
			//c_genBornData.d_sigma[pair*c_gsystem.widthTot + d_i] = mult;
			/*mult *= c_VVdW[at2];
			s_f[threadIdx.x].x += mult*r2.x;
			s_f[threadIdx.x].y += mult*r2.y;
			s_f[threadIdx.x].z += mult*r2.z;*/
			//pair++;
		}
		for(i = 0; i < c_pairsListsData.d_pairs13ListCount[d_i]; i++){
			r2 =  tex1Dfetch(t_coord, c_pairsListsData.d_pairs13List[i*c_gsystem.widthTot + d_i]);
			at2 = (int)r2.w;
			r2.x -= s_Ri[threadIdx.x].x;
			r2.y -= s_Ri[threadIdx.x].y;
			r2.z -= s_Ri[threadIdx.x].z;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			r2.w = 1.0f/r2.w;
			mult = c_genBornData.P3;
			mult *= r2.w;
			mult *= r2.w;
			//s_f[threadIdx.x].w += mult*c_VVdW[at2];
			alpha += mult*tex1Dfetch(t_VVdW, at2);//c_VVdW[at2];
			//mult *= -4.0f*r2.w;
			//c_genBornData.d_sigma[pair*c_gsystem.widthTot + d_i] = mult;
			/*mult *= c_VVdW[at2];
			s_f[threadIdx.x].x += mult*r2.x;
			s_f[threadIdx.x].y += mult*r2.y;
			s_f[threadIdx.x].z += mult*r2.z;*/
			//pair++;
		}

		for(i = 0; i < c_pairsListsData.d_pairs14ListCount[d_i]; i++){
			r2 =  tex1Dfetch(t_coord, c_pairsListsData.d_pairs14List[i*c_gsystem.widthTot + d_i]);
			at2 = (int)r2.w;
			r2.x -= s_Ri[threadIdx.x].x;
			r2.y -= s_Ri[threadIdx.x].y;
			r2.z -= s_Ri[threadIdx.x].z;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			arg = tex1Dfetch(t_RVdW, at1) + tex1Dfetch(t_RVdW, at2);//c_RVdW[at1] + c_RVdW[at2];
			arg = arg*arg;
			arg = r2.w/arg;
			arg *= c_genBornData.P5;
			r2.w = 1.0f/r2.w;
			mult = c_genBornData.P4;
			mult *= r2.w;
			mult *= r2.w;
			if(arg > 1.0f){
				//s_f[threadIdx.x].w += mult*c_VVdW[at2];
				alpha += mult*tex1Dfetch(t_VVdW, at2);//c_VVdW[at2];
				//mult *= -4.0f;
				//mult *= r2.w;
			} else {
				arg *= M_PI;
				CCF = __cosf(arg);
				CCF = 1.0f - CCF;
				mult *= CCF;
				//s_f[threadIdx.x].w += mult*0.25f*CCF*c_VVdW[at2];
				alpha += mult*0.25f*CCF*tex1Dfetch(t_VVdW, at2);//c_VVdW[at2];
				//mult *= r2.w;
				//mult *= -CCF+__sinf(arg)*arg;
			}
			//c_genBornData.d_sigma[pair*c_gsystem.widthTot + d_i] = mult;
			/*mult *= c_VVdW[at2];
			s_f[threadIdx.x].x += mult*r2.x;
			s_f[threadIdx.x].y += mult*r2.y;
			s_f[threadIdx.x].z += mult*r2.z;*/
			//pair++;
		}
		for(i = 0; i < c_pairsListsData.d_pairsLJListCount[d_i]; i++){
			r2 =  tex1Dfetch(t_coord, c_pairsListsData.d_pairsLJList[i*c_gsystem.widthTot + d_i]);
			at2 = (int)r2.w;
			r2.x -= s_Ri[threadIdx.x].x;
			r2.y -= s_Ri[threadIdx.x].y;
			r2.z -= s_Ri[threadIdx.x].z;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			if(r2.w < c_nonBondedData.ljCutoff2){
				arg = tex1Dfetch(t_RVdW, at1) + tex1Dfetch(t_RVdW, at2);//c_RVdW[at1] + c_RVdW[at2];
				arg = arg*arg;
				arg = r2.w/arg;
				arg *= c_genBornData.P5;
				r2.w = 1.0f/r2.w;
				mult = c_genBornData.P4;
				mult *= r2.w;
				mult *= r2.w;
				if(arg > 1.0f){
					//s_f[threadIdx.x].w += mult*c_VVdW[at2];
					alpha += mult*tex1Dfetch(t_VVdW, at2);//c_VVdW[at2];
					//mult *= -4.0f;
					//mult *= r2.w;
				} else {
					arg *= M_PI;
					CCF = __cosf(arg);
					CCF = 1.0f - CCF;
					mult *= CCF;
					//s_f[threadIdx.x].w += mult*0.25f*CCF*c_VVdW[at2];
					alpha += mult*0.25f*CCF*tex1Dfetch(t_VVdW, at2);//c_VVdW[at2];
					//mult *= r2.w;
					//mult *= -CCF+__sinf(arg)*arg;
				}
			//	c_genBornData.d_sigma[pair*c_gsystem.widthTot + d_i] = mult;
				/*mult *= c_VVdW[at2];
				s_f[threadIdx.x].x += mult*r2.x;
				s_f[threadIdx.x].y += mult*r2.y;
				s_f[threadIdx.x].z += mult*r2.z;*/
				//pair++;
			}
		}
		/*s_f[threadIdx.x].w = -0.5f*COULOMB_CONSTANT/s_f[threadIdx.x].w;
		mult = 2.0f*s_f[threadIdx.x].w*s_f[threadIdx.x].w/COULOMB_CONSTANT;
		s_f[threadIdx.x].x *= mult;
		s_f[threadIdx.x].y *= mult;
		s_f[threadIdx.x].z *= mult;
		c_genBornData.d_alphaG[d_i] = s_f[threadIdx.x];
		c_genBornData.d_alpha[d_i] = s_f[threadIdx.x].w;*/
		//alpha *= c_genBornData.epsilon;
		alpha = -c_genBornData.halfCoulomb/alpha;
		if(alpha <= 0 || alpha > c_genBornData.cutAlpha){
			alpha = c_genBornData.cutAlpha;
		}
		c_genBornData.d_alpha[d_i] = alpha;
	}
}

__global__ void genBorndGda_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 r = c_gsystem.d_coord[d_i];//tex1Dfetch(t_coord, d_i);
		float4 f = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
		float4 r2;
		//int at1 = (int)s_Ri[threadIdx.x].w;
		int at2;
		float ai = c_genBornData.d_alpha[d_i];//tex1Dfetch(t_alpha, d_i);
		int i;
		float dGda = 0.0f;//ai*ai;
		//dGda = c_charges[at1]/dGda;
		float mult, arg, mult2;
		float pot;
		for(i = 0; i < c_pairsListsData.d_pairsExclusionListCount[d_i]; i++){
			int j = c_pairsListsData.d_pairsExclusionList[i*c_gsystem.widthTot + d_i];
			if(j != d_i){
				r2 =  tex1Dfetch(t_coord, j);
				at2 = (int)r2.w;
				r2.x -= r.x;
				r2.y -= r.y;
				r2.z -= r.z;
				DO_PBC(r2);
				r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;

				arg = 0.25f*r2.w;
				float aj = c_genBornData.d_alpha[j];//tex1Dfetch(t_alpha, j);
				arg /= ai*aj;
				float expon = __expf(-arg);

				mult = r2.w + ai*aj*expon;
				mult = 1.0f/mult;
				mult = sqrtf(mult);
				mult = mult*mult*mult;
				mult *= tex1Dfetch(t_charges, at2);//c_charges[at2];

				mult2 = mult;

				mult *= 1.0f+arg;
				mult *= expon;

				dGda += mult*aj;

				mult2 *= 2.0f*(1.0f-0.25f*expon);

				/*mult *= c_genBornData.d_sigma[i*c_gsystem.widthTot + d_i];
				mult *= aj*aj*2.0f/COULOMB_CONSTANT;
				mult *= c_VVdW[at1];
				mult *= ai;

				mult += mult2;*/

				f.x += mult2*r2.x;
				f.y += mult2*r2.y;
				f.z += mult2*r2.z;
			}
		}

		for(i = 0; i < c_pairsListsData.d_pairsLJListCount[d_i]; i++){
			int j = c_pairsListsData.d_pairsLJList[i*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, j);
			at2 = (int)r2.w;
			r2.x -= r.x;
			r2.y -= r.y;
			r2.z -= r.z;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;

			if(r2.w < c_nonBondedData.ljCutoff2){

				arg = 0.25f*r2.w;
				float aj = c_genBornData.d_alpha[j];//tex1Dfetch(t_alpha, j);
				arg /= ai*aj;
				float expon = __expf(-arg);

				mult = r2.w + ai*aj*expon;
				mult = 1.0f/mult;

				mult = sqrtf(mult);
				pot = mult;
				mult = mult*mult*mult;

				mult *= tex1Dfetch(t_charges, at2);//c_charges[at2];

				pot *= tex1Dfetch(t_charges, at2);

				mult2 = mult;

				mult *= 1.0f+arg;
				mult *= expon;
				mult *= swfunct(r2.w);

				dGda += mult*aj;

				mult2 *=2.0f*(1.0f-0.25f*expon);

				/*mult *= c_genBornData.d_sigma[i*c_gsystem.widthTot + d_i];
				mult *= aj*aj*2.0f/COULOMB_CONSTANT;
				mult *= c_VVdW[at1];
				mult *= ai;

				mult += mult2;*/

				mult2 *= swfunct(r2.w);
				mult2 -= pot*dswfunct(r2.w);

				f.x += mult2*r2.x;
				f.y += mult2*r2.y;
				f.z += mult2*r2.z;
			}

		}

		//f = make_float4(0.0f, 0.0f, 0.0f, 0.0f);

		//dGda = 0.0f;
		c_genBornData.d_dGda[d_i] = dGda;
		c_genBornData.d_dGdr[d_i] = f;
	}
}

__global__ void genBorn_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 r = c_gsystem.d_coord[d_i];//tex1Dfetch(t_coord, d_i);
		//s_f[threadIdx.x] = c_genBornData.d_alphaG[d_i];
		float3 f;
		float4 r2;
		int at1, at2;
		at1 = (int)r.w;
		float ai = tex1Dfetch(t_alpha, d_i);
		int i, j;
		float CCF;
		float dGdai = c_genBornData.d_dGda[d_i]*ai*ai;
		float dGdaj;
		float mult = tex1Dfetch(t_charges, at1);
		dGdai += mult;//c_charges[at1];
		dGdai *= mult;//c_charges[at1];
		//s_f[threadIdx.x].x *= mult;
		//s_f[threadIdx.x].y *= mult;
		//s_f[threadIdx.x].z *= mult;
		f.x = 0.0f;
		f.y = 0.0f;
		f.z = 0.0f;
		for(i = 0; i < c_pairsListsData.d_pairs12ListCount[d_i]; i++){
			j = c_pairsListsData.d_pairs12List[i*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, j);
			at2 = (int)r2.w;
			r2.x -= r.x;
			r2.y -= r.y;
			r2.z -= r.z;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			float aj = tex1Dfetch(t_alpha, j);
			r2.w = 1.0f/r2.w;
			mult = c_genBornData.P2;
			mult *= r2.w;
			mult *= r2.w;
			mult *= -4.0f*r2.w;

			float dai = tex1Dfetch(t_VVdW, at2);
			float daj = tex1Dfetch(t_VVdW, at1);
			mult *= TWO_OVER_COULOMB_CONSTANT;//2.0f/COULOMB_CONSTANT;
			//mult *= c_genBornData.epsilon;

			dai *= mult;
			daj *= mult;

			mult = tex1Dfetch(t_charges, at2);
			dGdaj = c_genBornData.d_dGda[j];
			dGdaj *= aj;
			dGdaj *= aj;
			dGdaj += mult;//c_charges[at2];
			dGdaj *= mult;//c_charges[at2];

			mult = dGdai*dai + daj*dGdaj;

			f.x += mult*r2.x;
			f.y += mult*r2.y;
			f.z += mult*r2.z;
		}

		for(i = 0; i < c_pairsListsData.d_pairs13ListCount[d_i]; i++){
			j = c_pairsListsData.d_pairs13List[i*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, j);
			at2 = (int)r2.w;
			r2.x -= r.x;
			r2.y -= r.y;
			r2.z -= r.z;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			float aj = tex1Dfetch(t_alpha, j);
			r2.w = 1.0f/r2.w;
			mult = c_genBornData.P3;
			mult *= r2.w;
			mult *= r2.w;
			mult *= -4.0f*r2.w;

			float dai = tex1Dfetch(t_VVdW, at2);
			float daj = tex1Dfetch(t_VVdW, at1);
			mult *= TWO_OVER_COULOMB_CONSTANT;//2.0f/COULOMB_CONSTANT;
			//mult *= c_genBornData.epsilon;

			dai *= mult;
			daj *= mult;

			mult = tex1Dfetch(t_charges, at2);
			dGdaj = c_genBornData.d_dGda[j];
			dGdaj *= aj;
			dGdaj *= aj;
			dGdaj += mult;//c_charges[at2];
			dGdaj *= mult;//c_charges[at2];

			mult = dGdai*dai + daj*dGdaj;

			f.x += mult*r2.x;
			f.y += mult*r2.y;
			f.z += mult*r2.z;
		}

		for(i = 0; i < c_pairsListsData.d_pairs14ListCount[d_i]; i++){
			j = c_pairsListsData.d_pairs14List[i*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, j);
			at2 = (int)r2.w;
			r2.x -= r.x;
			r2.y -= r.y;
			r2.z -= r.z;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			float aj = tex1Dfetch(t_alpha, j);

			float arg = tex1Dfetch(t_RVdW, at1) + tex1Dfetch(t_RVdW, at2);
			arg = arg*arg;
			arg = r2.w/arg;
			arg *= c_genBornData.P5;
			r2.w = 1.0f/r2.w;
			mult = c_genBornData.P4;
			mult *= r2.w;
			mult *= r2.w;
			if(arg > 1.0f){
				mult *= -4.0f;
				mult *= r2.w;
			} else {
				arg *= M_PI;
				CCF = __cosf(arg);
				CCF = 1.0f - CCF;
				mult *= CCF;
				mult *= r2.w;
				mult *= -CCF+__sinf(arg)*arg;
			}

			float dai = tex1Dfetch(t_VVdW, at2);
			float daj = tex1Dfetch(t_VVdW, at1);
			mult *= TWO_OVER_COULOMB_CONSTANT;//2.0f/COULOMB_CONSTANT;
			//mult *= c_genBornData.epsilon;

			dai *= mult;//*tex1Dfetch(t_VVdW, at2);//*TWO_OVER_COULOMB_CONSTANT;///*c_VVdW[at2]*/*2.0f/COULOMB_CONSTANT;
			daj *= mult;//*tex1Dfetch(t_VVdW, at1);//*TWO_OVER_COULOMB_CONSTANT;///*c_VVdW[at1]*/*2.0f/COULOMB_CONSTANT;

			mult = tex1Dfetch(t_charges, at2);
			dGdaj = c_genBornData.d_dGda[j];
			dGdaj *= aj;
			dGdaj *= aj;
			dGdaj += mult;//c_charges[at2];
			dGdaj *= mult;//c_charges[at2];

			mult = dGdai*dai + daj*dGdaj;

			f.x += mult*r2.x;
			f.y += mult*r2.y;
			f.z += mult*r2.z;
		}

		for(i = 0; i < c_pairsListsData.d_pairsLJListCount[d_i]; i++){
			j = c_pairsListsData.d_pairsLJList[i*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, j);
			at2 = (int)r2.w;
			r2.x -= r.x;
			r2.y -= r.y;
			r2.z -= r.z;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			if(r2.w < c_nonBondedData.ljCutoff2){
				float aj = tex1Dfetch(t_alpha, j);

				float arg = tex1Dfetch(t_RVdW, at1) + tex1Dfetch(t_RVdW, at2);
				arg = arg*arg;
				arg = r2.w/arg;
				arg *= c_genBornData.P5;
				r2.w = 1.0f/r2.w;
				mult = c_genBornData.P4;
				mult *= r2.w;
				mult *= r2.w;
				if(arg > 1.0f){
					mult *= -4.0f;
					mult *= r2.w;
				} else {
					arg *= M_PI;
					CCF = __cosf(arg);
					CCF = 1.0f - CCF;
					mult *= CCF;
					mult *= r2.w;
					mult *= -CCF+__sinf(arg)*arg;
				}

				float dai = tex1Dfetch(t_VVdW, at2);
				float daj = tex1Dfetch(t_VVdW, at1);
				mult *= TWO_OVER_COULOMB_CONSTANT;//2.0f/COULOMB_CONSTANT;
				//mult *= c_genBornData.epsilon;

				dai *= mult;//*tex1Dfetch(t_VVdW, at2);//*TWO_OVER_COULOMB_CONSTANT;///*c_VVdW[at2]*/*2.0f/COULOMB_CONSTANT;
				daj *= mult;//*tex1Dfetch(t_VVdW, at1);//*TWO_OVER_COULOMB_CONSTANT;///*c_VVdW[at1]*/*2.0f/COULOMB_CONSTANT;

				mult = tex1Dfetch(t_charges, at2);
				dGdaj = c_genBornData.d_dGda[j];
				dGdaj *= aj;
				dGdaj *= aj;
				dGdaj += mult;//c_charges[at2];
				dGdaj *= mult;//c_charges[at2];

				mult = dGdai*dai + daj*dGdaj;

				f.x += mult*r2.x;
				f.y += mult*r2.y;
				f.z += mult*r2.z;
			}
		}

		float4 dGdr = c_genBornData.d_dGdr[d_i];
		mult = tex1Dfetch(t_charges, at1);

		f.x += mult*dGdr.x;
		f.y += mult*dGdr.y;
		f.z += mult*dGdr.z;

		mult = c_genBornData.epsilonTimesHalfCoulomb;//0.5f*COULOMB_CONSTANT*c_genBornData.epsilon;//2.0f/COULOMB_CONSTANT;
		r2 = c_gsystem.d_forces[d_i];
		r2.x += f.x*mult;
		r2.y += f.y*mult;
		r2.z += f.z*mult;
		c_gsystem.d_forces[d_i] = r2;
	}
}

inline void compute(){
	genBornRadii_kernel<<<genBornBlockCount, genBornBlockSize>>>();
	/*cudaThreadSynchronize();
	checkCUDAError("GenBorn Radii");
	cudaMemcpy(genBornData.h_alpha, genBornData.d_alpha,
					gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	int traj, i;
	char filename[100];
	sprintf(filename, "gbalpha_%d.txt", step);
	FILE* gbalpha = fopen(filename, "w");
	for(traj = 0; traj < parameters.Ntr; traj++){
		for(i = 0; i < gsystem.N; i++){
			fprintf(gbalpha, "%d\t%s\t%d\t%s\t%f\n",
					i, topology.atoms[i].resName, topology.atoms[i].resid, topology.atoms[i].name, genBornData.h_alpha[i]);
			if(genBornData.h_alpha[i] < 0.0f){
				printf("GB Radii is negative for %d\t%s\t%d\t%s\t%f\n",
									i, topology.atoms[i].resName, topology.atoms[i].resid, topology.atoms[i].name, genBornData.h_alpha[i]);
				exit(0);
			}
			if(genBornData.h_alpha[i] > genBornData.cutAlpha){
				printf("GB Radii is large for %d\t%s\t%d\t%s\t%f\n",
									i, topology.atoms[i].resName, topology.atoms[i].resid, topology.atoms[i].name, genBornData.h_alpha[i]);
			}
			pdbOutputData.atoms[i].occupancy = genBornData.h_alpha[i];
			pdbOutputData.atoms[i].x = topology.atoms[i].x*10.0f;
			pdbOutputData.atoms[i].y = topology.atoms[i].y*10.0f;
			pdbOutputData.atoms[i].z = topology.atoms[i].z*10.0f;
		}
	}
	fclose(gbalpha);
	writePDB("bornRadii.pdb", &pdbOutputData);*/
	//if(step > 0){exit(0);}
/*	float pot = 0.0f;
	int j;
	float4 ri, rj;
	int ati, atj;
	for(i = 0; i < atomCount; i++){
		ri = gsystem.h_coord[i];
		ati = (int)ri.w;
		for(j = 0; j < atomCount; j++){
			rj = gsystem.h_coord[j];
			atj = (int)rj.w;
			rj.x -= ri.x;
			rj.y -= ri.y;
			rj.z -= ri.z;
			float r2 = rj.x*rj.x + rj.y*rj.y + rj.z*rj.z;
			float alphaij = genBornData.h_alphaG[i].w*genBornData.h_alphaG[j].w;
			float arg = -0.25f*r2/alphaij;
			float expon = expf(arg);
			pot += atomTypes[ati].charge*atomTypes[atj].charge/(sqrtf(r2+alphaij*expon));
		}
	}
	pot *= -0.5f*genBornData.epsilon*COULOMB_CONSTANT/KCAL_TO_KJOULE;
	printf("Solvation energy: %f.\n", pot);
	printf("\n\nGrad alpha:\n");
	for(i = 0; i < atomCount; i++){
		printf("%d: (%f, %f, %f) %f\n", i, genBornData.h_alphaG[i].x, genBornData.h_alphaG[i].y, genBornData.h_alphaG[i].z,
					sqrtf(genBornData.h_alphaG[i].x*genBornData.h_alphaG[i].x + genBornData.h_alphaG[i].y*genBornData.h_alphaG[i].y + genBornData.h_alphaG[i].z*genBornData.h_alphaG[i].z));
	}*/
	genBorndGda_kernel<<<genBornBlockCount, genBornBlockSize>>>();
	//cudaThreadSynchronize();
/*	cudaMemcpy(genBornData.h_dGda, genBornData.d_dGda, gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(genBornData.h_dGdr, genBornData.d_dGdr, gsystem.Ntot*sizeof(float4), cudaMemcpyDeviceToHost);
	FILE* gb_dGdadr = fopen("gb_dGdadr.txt", "w");
	int naned = 0;
	for(traj = 0; traj < parameters.Ntr; traj++){
		for(i = 0; i < gsystem.N; i++){
			fprintf(gb_dGdadr, "%d\t%s\t%d\t%s\t%f\t%f\t%f\t%f\t%f\n",
					i, topology.atoms[i].resName, topology.atoms[i].resid, topology.atoms[i].name,
					genBornData.h_dGda[i],
					genBornData.h_dGdr[i].x,
					genBornData.h_dGdr[i].y,
					genBornData.h_dGdr[i].z,
					genBornData.h_dGdr[i].w);
			if(isnan(genBornData.h_dGda[i]) ||
					isnan(genBornData.h_dGdr[i].x) ||
					isnan(genBornData.h_dGdr[i].y) ||
					isnan(genBornData.h_dGdr[i].z) ||
					isnan(genBornData.h_dGdr[i].w)){
				printf("%d\t%s\t%d\t%s\t%f\t%f\t%f\t%f\t%f\n",
									i, topology.atoms[i].resName, topology.atoms[i].resid, topology.atoms[i].name,
									genBornData.h_dGda[i],
									genBornData.h_dGdr[i].x,
									genBornData.h_dGdr[i].y,
									genBornData.h_dGdr[i].z,
									genBornData.h_dGdr[i].w);
				naned = 1;
			}
		}
	}
	if(naned){
		exit(0);
	}
	fclose(gb_dGdadr);*/
#ifdef CUDA_USE_L1
	cudaFuncSetCacheConfig(genBorn_kernel, cudaFuncCachePreferL1);
#endif
	genBorn_kernel<<<genBornBlockCount, genBornBlockSize>>>();
/*	cudaMemcpy(gsystem.h_forces, gsystem.d_forces, gsystem.Ntot*sizeof(float4), cudaMemcpyDeviceToHost);
	float3 force = make_float3(0.0f, 0.0f, 0.0f);
	for(i = 0; i < gsystem.Ntot; i++){
		force.x += gsystem.h_forces[i].x;
		force.y += gsystem.h_forces[i].y;
		force.z += gsystem.h_forces[i].z;*/
		/*printf("%d: (%f, %f, %f) %f\n", i, gsystem.h_forces[i].x, gsystem.h_forces[i].y, gsystem.h_forces[i].z,
					sqrtf(gsystem.h_forces[i].x*gsystem.h_forces[i].x + gsystem.h_forces[i].y*gsystem.h_forces[i].y + gsystem.h_forces[i].z*gsystem.h_forces[i].z));*/
	/*}
	printf("Net force (GenBorn): (%f, %f, %f) %f\n", force.x, force.y, force.z,
			sqrtf(force.x*force.x + force.y*force.y + force.z*force.z));*/
	//exit(0);
	checkCUDAError("GenBorn");
}

__global__ void genBornEnergy_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		s_Ri[threadIdx.x] = c_gsystem.d_coord[d_i];//tex1Dfetch(t_coord, d_i);
		float4 r2;
		int at1, at2;
		at1 = (int)s_Ri[threadIdx.x].w;
		float ai = tex1Dfetch(t_alpha, d_i);
		int i, j;
		float pot = 0.0f;
		for(i = 0; i < c_pairsListsData.d_pairsExclusionListCount[d_i]; i++){
			j = c_pairsListsData.d_pairsExclusionList[i*c_gsystem.widthTot + d_i];
			//int j = c_pairsListsData.d_pairsLJList[i*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, j);
			at2 = (int)r2.w;
			r2.x -= s_Ri[threadIdx.x].x;
			r2.y -= s_Ri[threadIdx.x].y;
			r2.z -= s_Ri[threadIdx.x].z;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			float arg = -0.25f*r2.w;
			float aj = tex1Dfetch(t_alpha, j);
			arg /= ai*aj;
			float expon = expf(arg);
			float mult;
			mult = r2.w + ai*aj*expon;
			mult = sqrtf(mult);
			mult = 1.0f/mult;
			mult *= tex1Dfetch(t_charges, at1)*tex1Dfetch(t_charges, at2);
			pot += mult;
		}
		for(i = 0; i < c_pairsListsData.d_pairsLJListCount[d_i]; i++){
			j = c_pairsListsData.d_pairsLJList[i*c_gsystem.widthTot + d_i];
			//int j = c_pairsListsData.d_pairsLJList[i*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, j);
			at2 = (int)r2.w;
			r2.x -= s_Ri[threadIdx.x].x;
			r2.y -= s_Ri[threadIdx.x].y;
			r2.z -= s_Ri[threadIdx.x].z;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			if(r2.w < c_genBornData.roff2){
				float arg = -0.25f*r2.w;
				float aj = tex1Dfetch(t_alpha, j);
				arg /= ai*aj;
				float expon = expf(arg);
				float mult;
				mult = r2.w + ai*aj*expon;
				mult = sqrtf(mult);
				mult = 1.0f/mult;
				mult *= tex1Dfetch(t_charges, at1);
				mult *= tex1Dfetch(t_charges, at2);
				mult *= swfunct(r2.w);
				pot += mult;
			}
		}
		c_genBornData.d_genBornEnergies[d_i] = pot;
	}
}


inline void computeEnergy(){
	genBornRadii_kernel<<<genBornBlockCount, genBornBlockSize>>>();
	checkCUDAError("GenBorn Radii");
	genBornEnergy_kernel<<<genBornBlockCount, genBornBlockSize>>>();
	checkCUDAError("GenBorn Energy");

	cudaMemcpy(genBornData.h_genBornEnergies, genBornData.d_genBornEnergies,
					gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		float pot = 0.0f;
		for(i = 0; i < gsystem.N; i++){
			/*if(nonBondedData.h_ljEnergies[i] > 0.0f){
				printf("Atom %d has another one close to it (LJ Energy = %f)\n", i, nonBondedData.h_ljEnergies[i]);
			}*/
			//if(topology.atoms[i].type[0] == 'H'){
			pot += genBornData.h_genBornEnergies[i + traj*gsystem.N];
			//}
		}
		//pot /= 2.0f;
		genBornEnergyOutput.values[traj] = -pot*genBornData.epsilonTimesHalfCoulomb;
	}
	checkCUDAError("GenBorn Energy - MemCpy");

}

void destroy(){

}

void initGenBornParameter(char* name, char* value){
	if(name[0] == 'P'){
		char i = name[1];
		switch (i) {
			case '1':
				genBornData.P1 = atof(value);
				DPRINTF("P1 = %f\n", genBornData.P1);
				break;
			case '2':
				genBornData.P2 = atof(value);
				DPRINTF("P2 = %f\n", genBornData.P2);
				break;
			case '3':
				genBornData.P3 = atof(value);
				DPRINTF("P3 = %f\n", genBornData.P3);
				break;
			case '4':
				genBornData.P4 = atof(value);
				DPRINTF("P4 = %f\n", genBornData.P4);
				break;
			case '5':
				genBornData.P5 = atof(value);
				DPRINTF("P5 = %f\n", genBornData.P5);
				break;
		}
	}
	if(strncmp("Lambda", name, 6) == 0){
		genBornData.Lambda = atof(value);
		DPRINTF("Lambda = %f\n", genBornData.Lambda);
	}
	if(strncmp("Phi", name, 6) == 0){
		genBornData.Phi = atof(value);
		DPRINTF("Phi = %f\n", genBornData.Phi);
	}
}

#undef LOG

} // namespace genborn_potential
