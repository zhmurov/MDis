/*
 * GBSWPotential.cu
 *
 *  Created on: Sep 22, 2011
 *      Author: zhmurov
 */

#include "GBSWPotential.cuh"

#include "../Core/global.h"
#include "../Util/Log.h"
#include "../IO/quadratureio.h"
#include "../IO/gbswpbradiiio.h"

namespace gbsw_potential {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<GBSW Potential> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void create(){
	if(getYesNoParameter(PARAMETER_GBSW_ON, 0)){
		potential.compute = compute;
		potential.destroy = destroy;
		sprintf(potential.name, "GBSW");
		potentials[potentialsCount++] = &potential;
		init();
	}
}

void init(){

	gbswBlockSize = BLOCK_SIZE;
	gbswBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;

	int i, j;

	gbswData.cutoff = getFloatParameter(PARAMETER_GBSW_CUTOFF, getFloatParameter(PARAMETER_NB_CUTOFF));
	gbswData.w = getFloatParameter(PARAMETER_GBSW_SMOOTH_LENGTH, DEFAULT_GBSW_SMOOTH_LENGTH);
	gbswData.a0 = getFloatParameter(PARAMETER_GBSW_A0, DEFAULT_GBSW_A0);
	gbswData.a1 = getFloatParameter(PARAMETER_GBSW_A1, DEFAULT_GBSW_A1);

	gbswData.threeOverFourW = 3.0f/(4.0f*gbswData.w);
	gbswData.oneOverFourW3 = 1.0f/(4.0f*gbswData.w*gbswData.w*gbswData.w);
	gbswData.threeOverFourW3 = 3.0f/(4.0f*gbswData.w*gbswData.w*gbswData.w);

	allocateCPU((void**)&gbswData.h_alpha, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&gbswData.d_alpha, gsystem.Ntot*sizeof(float));

	allocateCPU((void**)&gbswData.h_dG1, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&gbswData.d_dG1, gsystem.Ntot*sizeof(float));

	for(i = 0; i < gsystem.Ntot; i++){
		gbswData.h_alpha[i] = 0.0f;
	}
	cudaMemcpy(gbswData.d_alpha, gbswData.h_alpha, gsystem.Ntot*sizeof(float), cudaMemcpyHostToDevice);

	allocateCPU((void**)&gbswData.h_RPB, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&gbswData.d_RPB, gsystem.Ntot*sizeof(float));

	char pbradiifilename[PARAMETER_LENGTH];
	getMaskedParameter(pbradiifilename, PARAMETER_GBSW_PBRADII_FILENAME);
	PBRFileData pbRadii;
	readPBRadiiFile(pbradiifilename, &pbRadii);

	float s = -1.2964f*gbswData.w + 0.9914f; // Nina et. al., Biophys. Chem. 1999, v.78, p.89.
	                                         // Coefficients were obtained using linear regression of the data in Table 2
	                                         // of Nina et. al. paper since the formula given after the table
	                                         // has only two significant numbers when data itself have three.
	s = 0.979f; // Temporary. For better agreement with CHARMM
	printf("Scaling factor for PB radii is %f.\n", s);
	for(i = 0; i < gsystem.Ntot; i++){
		for(j = 0; j < pbRadii.rCount; j++){
			if(comparePBRadiiWithWC(topology.atoms[i].name, pbRadii.r[j].name) &&
					comparePBRadiiWithWC(topology.atoms[i].resName, pbRadii.r[j].resName)){
				gbswData.h_RPB[i] = pbRadii.r[j].value/10.0f;
			}
		}

		if(gbswData.h_RPB[i] != 0){
			gbswData.h_RPB[i] = s*(gbswData.h_RPB[i] + gbswData.w);
		}
		printf("%d:  %s (%s%d):   %f\n", i,
				topology.atoms[i].name,
				topology.atoms[i].resName,
				topology.atoms[i].resid,
				gbswData.h_RPB[i]);
	}

	gbswData.nabla = FLT_MAX;//gbswData.h_RPB[0];
	for(i = 1; i < gsystem.Ntot; i++){
		if(gbswData.h_RPB[i] < gbswData.nabla){// && gbswData.h_RPB[i] != 0){
			gbswData.nabla = gbswData.h_RPB[i];
		}
		//gbswData.h_RPB[i] -= gbswData.w;
		//gbswData.h_RPB[i] *= 1000.0f;
	}
	if(gbswData.nabla == 0.0f){
		gbswData.nabla = 0.05f;
	}

	cudaMemcpy(gbswData.d_RPB, gbswData.h_RPB, gsystem.Ntot*sizeof(float), cudaMemcpyHostToDevice);

	AngularQuadrature aquadr;
	char angularQuadrFilename[PARAMETER_LENGTH];
	getMaskedParameter(angularQuadrFilename, PARAMETER_GBSW_ANG_QUADR_FILENAME);
	readAngularQuadrature(angularQuadrFilename, &aquadr);

	allocateCPU((void**)&gbswData.h_angularQuadr, aquadr.pointCount*sizeof(float4));
	allocateGPU((void**)&gbswData.d_angularQuadr, aquadr.pointCount*sizeof(float4));

	gbswData.angularPointsCount = aquadr.pointCount;
	for(i = 0; i < aquadr.pointCount; i++){
		gbswData.h_angularQuadr[i].x = aquadr.points[i].x;
		gbswData.h_angularQuadr[i].y = aquadr.points[i].y;
		gbswData.h_angularQuadr[i].z = aquadr.points[i].z;
		gbswData.h_angularQuadr[i].w = aquadr.points[i].w;
	}

	cudaMemcpy(gbswData.d_angularQuadr, gbswData.h_angularQuadr,
				gbswData.angularPointsCount*sizeof(float4), cudaMemcpyHostToDevice);

	int nrquadrtmp = getIntegerParameter(PARAMETER_GBSW_RAD_QUADR_POINTCOUNT, DEFAULT_GBSW_RAD_QUADR_POINTCOUNT);

	if(nrquadrtmp == 0){   // Two-step definition of quadrature with total of 24 points if
		                   // the number of points is not defined in the configuration file
		                   // See: CHARMM gbsw.src file for details (line 491 in ver. c35b5)

		gbswData.radialPointsCount = 24;

		RadialQuadrature rquadr1;
		RadialQuadrature rquadr2;
		allocateCPU((void**)&gbswData.h_radialQuadr, gbswData.radialPointsCount*sizeof(float2));
		allocateGPU((void**)&gbswData.d_radialQuadr, gbswData.radialPointsCount*sizeof(float2));
		if(gbswData.nabla - gbswData.w > 0.05f){
			gbswData.nabla = gbswData.nabla - gbswData.w;
		}
		float rmidint = gbswData.nabla + 0.05f;
		createGaussianLegendreRadialQuadrature(gbswData.nabla, rmidint, 5, &rquadr1);
		createGaussianLegendreRadialQuadrature(rmidint, gbswData.cutoff, gbswData.radialPointsCount-5, &rquadr2);

		allocateCPU((void**)&gbswData.h_radialQuadr, gbswData.radialPointsCount*sizeof(float2));

		for(i = 0; i < rquadr1.pointCount; i++){
			gbswData.h_radialQuadr[i].x = rquadr1.points[i].r;
			gbswData.h_radialQuadr[i].y = rquadr1.points[i].w;
		}
		for(i = 0; i < rquadr2.pointCount; i++){
			gbswData.h_radialQuadr[rquadr1.pointCount + i].x = rquadr2.points[i].r;
			gbswData.h_radialQuadr[rquadr1.pointCount + i].y = rquadr2.points[i].w;
		}

	} else {
		RadialQuadrature rquadr;
		createGaussianLegendreRadialQuadrature(gbswData.nabla, gbswData.cutoff,
				getIntegerParameter(PARAMETER_GBSW_RAD_QUADR_POINTCOUNT, DEFAULT_GBSW_RAD_QUADR_POINTCOUNT), &rquadr);

		gbswData.radialPointsCount = rquadr.pointCount;

		allocateCPU((void**)&gbswData.h_radialQuadr, gbswData.radialPointsCount*sizeof(float2));

		for(i = 0; i < rquadr.pointCount; i++){
			gbswData.h_radialQuadr[i].x = rquadr.points[i].r;
			gbswData.h_radialQuadr[i].y = rquadr.points[i].w;
		}
	}

	allocateGPU((void**)&gbswData.d_radialQuadr, gbswData.radialPointsCount*sizeof(float2));

	cudaMemcpy(gbswData.d_radialQuadr, gbswData.h_radialQuadr,
				gbswData.radialPointsCount*sizeof(float2), cudaMemcpyHostToDevice);

	PDB quadrPDB;
	quadrPDB.ssCount = 0;
	quadrPDB.atomCount = gbswData.angularPointsCount*gbswData.radialPointsCount;
	quadrPDB.atoms = (PDBAtom*)calloc(quadrPDB.atomCount, sizeof(PDBAtom));
	for(i = 0; i < gbswData.angularPointsCount; i++){
		float4 raq = gbswData.h_angularQuadr[i];
		for(j = 0; j < gbswData.radialPointsCount; j++){
			float2 rrq = gbswData.h_radialQuadr[j];
			float4 rmn;
			rmn.x = raq.x*rrq.x;
			rmn.y = raq.y*rrq.x;
			rmn.z = raq.z*rrq.x;

			quadrPDB.atoms[gbswData.radialPointsCount*i + j].id = gbswData.radialPointsCount*i + j;
			sprintf(quadrPDB.atoms[gbswData.radialPointsCount*i + j].name, "CA");
			sprintf(quadrPDB.atoms[gbswData.radialPointsCount*i + j].resName, "QUA");
			quadrPDB.atoms[gbswData.radialPointsCount*i + j].chain = 'Q';
			quadrPDB.atoms[gbswData.radialPointsCount*i + j].altLoc = ' ';
			quadrPDB.atoms[gbswData.radialPointsCount*i + j].resid = i;
			quadrPDB.atoms[gbswData.radialPointsCount*i + j].x = rmn.x;
			quadrPDB.atoms[gbswData.radialPointsCount*i + j].y = rmn.y;
			quadrPDB.atoms[gbswData.radialPointsCount*i + j].z = rmn.z;
		}
	}
	writePDB("quadr.pdb", &quadrPDB);

	gbswData.oneOverNabla = 1.0f/gbswData.nabla;
	gbswData.oneOverFourNabla4 = 1.0f/(4.0f*gbswData.nabla*gbswData.nabla*gbswData.nabla*gbswData.nabla);
	gbswData.oneOverFourPi = 1.0f;//1.0f/(4.0f*M_PI);

	cudaBindTexture(0, t_angularQuadr, gbswData.d_angularQuadr, gbswData.angularPointsCount*sizeof(float4));
	cudaBindTexture(0, t_radialQuadr, gbswData.d_radialQuadr, gbswData.radialPointsCount*sizeof(float2));

	allocateCPU((void**)&gbswData.h_gbswEnergies, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&gbswData.d_gbswEnergies, gsystem.Ntot*sizeof(float));

	allocateCPU((void**)&gbswData.h_isH, gsystem.Ntot*sizeof(int));
	allocateGPU((void**)&gbswData.d_isH, gsystem.Ntot*sizeof(int));
	for(i = 0; i < gsystem.Ntot; i++){
		if(topology.atoms[i].name[0] == 'H'){
			gbswData.h_isH[i] = 1;
		} else {
			gbswData.h_isH[i] = 0;
		}
	}
	cudaMemcpy(gbswData.d_isH, gbswData.h_isH, gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);

	allocateCPU((void**)&gbswData.h_dGda, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&gbswData.d_dGda, gsystem.Ntot*sizeof(float));

	allocateCPU((void**)&gbswData.h_dGdr, gsystem.Ntot*sizeof(float4));
	allocateGPU((void**)&gbswData.d_dGdr, gsystem.Ntot*sizeof(float4));

	cudaBindTexture(0, t_gbswalpha, gbswData.d_alpha,
					gsystem.Ntot*sizeof(float));

	cudaMemcpyToSymbol(c_gbswData, &gbswData, sizeof(GBSWData), 0, cudaMemcpyHostToDevice);
}

__device__ float gbswHij(float r, int j){
	float rPB = c_gbswData.d_RPB[j];
	float Rc1 = rPB - c_gbswData.w;
	float Rc2 = rPB + c_gbswData.w;
	if(r < Rc1){
		return 0;
	} else
	if(r < Rc2){
		float arg = r - rPB;
		return 0.5f + c_gbswData.threeOverFourW*arg - c_gbswData.oneOverFourW3*arg*arg*arg;
	} else {
		return 1;
	}
}

__device__ float gbswdHij(float r, int j){
	float rPB = c_gbswData.d_RPB[j];
	float Rc1 = rPB - c_gbswData.w;
	float Rc2 = rPB + c_gbswData.w;
	if(r < Rc1){
		return 0;
	} else
	if(r < Rc2){
		float arg = r - rPB;
		return c_gbswData.threeOverFourW - c_gbswData.threeOverFourW3*arg*arg;
	} else {
		return 0;
	}
}

__global__ void gbswRadii_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		//if(c_gbswData.d_isH[d_i] != 1){
		float4 r1 = tex1Dfetch(t_coord, d_i);
		float4 r2;
		int m, n, j;
		float dG0 = 0.0f;
		float dG1 = 0.0f;
		for(m = 0; m < c_gbswData.angularPointsCount; m++){
			float4 raq = tex1Dfetch(t_angularQuadr, m);
			for(n = 0; n < c_gbswData.radialPointsCount; n++){
				float4 rmn;
				float2 rrq = tex1Dfetch(t_radialQuadr, n);
				rmn.x = raq.x*rrq.x;
				rmn.y = raq.y*rrq.x;
				rmn.z = raq.z*rrq.x;
				rmn.w = rmn.x*rmn.x + rmn.y*rmn.y + rmn.z*rmn.z; //== rrq.x^2
				float w = raq.w*rrq.y;
				float h = 1.0f;
				for(j = 0; j < c_gsystem.N; j++){
					//if(c_gbswData.d_isH[j] != 1){
					//if(j != d_i){
					r2 = tex1Dfetch(t_coord, j);
					float4 r;
					r.x = r1.x + rmn.x - r2.x;
					r.y = r1.y + rmn.y - r2.y;
					r.z = r1.z + rmn.z - r2.z;
					r.w = r.x*r.x + r.y*r.y + r.z*r.z;
					r.w = sqrtf(r.w);
					h *= gbswHij(r.w, j);
					//}
					//}
				}
				h = 1.0f - h;
				w *= h;
				rmn.w = 1.0f/rmn.w;
				w *= rmn.w;
				dG0 += w;
				w *= rmn.w;
				w *= sqrtf(rmn.w);
				dG1 += w;
			}
		}
		dG0 = c_gbswData.oneOverNabla - c_gbswData.oneOverFourPi*dG0;
		dG1 = c_gbswData.oneOverFourNabla4 - c_gbswData.oneOverFourPi*dG1;
		dG1 = powf(dG1, 0.25f);
		c_gbswData.d_dG1[d_i] = dG1;
		dG0 *= c_gbswData.a0;
		dG1 *= c_gbswData.a1;
		dG0 += dG1;
		c_gbswData.d_alpha[d_i] = 1.0f/dG0;
		//}
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
		float ai = tex1Dfetch(t_gbswalpha, d_i);
		int j;
		float dGda = 0.0f;//ai*ai;
		//dGda = c_charges[at1]/dGda;
		float mult, arg, mult2;
		float pot;
		for(j = 0; j < c_gsystem.N; j++){
			if(j != d_i){
				r2 =  tex1Dfetch(t_coord, j);
				at2 = (int)r2.w;
				r2.x -= r.x;
				r2.y -= r.y;
				r2.z -= r.z;
				DO_PBC(r2);
				r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;

				arg = 0.25f*r2.w;
				float aj = tex1Dfetch(t_gbswalpha, j);
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

		c_gbswData.d_dGda[d_i] = dGda;
		c_gbswData.d_dGdr[d_i] = f;
	}
}

__global__ void gbsw_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){

	}
}

__global__ void gbswEnergy_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float pot = 0.0f;
		//if(c_gbswData.d_isH[d_i] != 1){
		float4 r1 = tex1Dfetch(t_coord, d_i);
		float4 r2;
		int at1, at2;
		at1 = (int)r1.w;
		float ai = c_gbswData.d_alpha[d_i];
		int j;
		for(j = 0; j < c_gsystem.Ntot; j++){
			//if(c_gbswData.d_isH[j] != 1){
			//if(j != d_i){
			r2 =  tex1Dfetch(t_coord, j);
			at2 = (int)r2.w;
			r2.x -= r1.x;
			r2.y -= r1.y;
			r2.z -= r1.z;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			float arg = -0.25f*r2.w;
			float aj = c_gbswData.d_alpha[j];
			arg /= ai*aj;
			float expon = expf(arg);
			float mult;
			mult = r2.w + ai*aj*expon;
			mult = sqrtf(mult);
			mult = 1.0f/mult;
			mult *= tex1Dfetch(t_charges, at1)*tex1Dfetch(t_charges, at2);
			pot += mult;
			//}
			//}
		//}
		}
		c_gbswData.d_gbswEnergies[d_i] = pot;
	}
}

void compute(){
	gbswRadii_kernel<<<gbswBlockCount, gbswBlockSize>>>();
	gbswEnergy_kernel<<<gbswBlockCount, gbswBlockSize>>>();
	cudaMemcpy(gbswData.h_alpha, gbswData.d_alpha, gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(gbswData.h_gbswEnergies, gbswData.d_gbswEnergies, gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	float pot = 0.0f;
	int traj, i;
	for(traj = 0; traj < parameters.Ntr; traj++){
		for(i = 0; i < gsystem.N; i++){
			printf("%d\t%s\t%d\t%s\t%f\t%f\n",
					i, topology.atoms[i].resName, topology.atoms[i].resid, topology.atoms[i].name,
					gbswData.h_alpha[i], gbswData.h_gbswEnergies[i]);
			pot += gbswData.h_gbswEnergies[i];
		}
	}
	pot = -COULOMB_CONSTANT*0.5f*(1.0f - 1.0f/80.0f)*pot*KCALL_PER_KJ;
	printf("Total energy: %f\n", pot);
	exit(0);
}


void destroy(){

}

#undef LOG

} // namespace gbsw_potential
