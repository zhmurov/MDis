/*
 * HarmonicConstraints.cu
 *
 *  Created on: May 25, 2011
 *	  Author: zhmurov
 */
#include "../Core/global.h"
#include "../Core/md.cuh"
#include "../Util/Log.h"
#include "HarmonicConstraints.cuh"

namespace harmonic_constraints {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<harmonic_constraints> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void create(){
	if (!getYesNoParameter(PARAMETER_CONSTRAINTS_ENABLED, 0))
		return;
	potential.compute = &compute;
	potential.destroy = &destroy;
	sprintf(potential.name, "Harmonic constraints");
	potentials[potentialsCount] = &potential;
	potentialsCount ++;
	init();
}

void init(){
	LOG << "Initializing harmonic constraints...";
	harmonicConstraintsBlockSize = BLOCK_SIZE;
	harmonicConstraintsBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;

	harmonicConstraintsData.Ks = getFloatParameter(PARAMETER_CONSTRAINTS_KS);

	if(getYesNoParameter(PARAMETER_CONSTRAINTS_LIMIT_FORCE, 0)){
		harmonicConstraintsData.limitForce = true;
		harmonicConstraintsData.maxF = getFloatParameter(PARAMETER_CONSTRAINTS_MAX_FORCE);
	}

	allocateCPU((void**)&harmonicConstraintsData.h_fixedConstraints, gsystem.Ntot*sizeof(float4));
	allocateGPU((void**)&harmonicConstraintsData.d_fixedConstraints, gsystem.Ntot*sizeof(float4));

	char refPDBFilename[PARAMETER_LENGTH];
	getMaskedParameter(refPDBFilename, PARAMETER_CONSTRAINTS_REF_FILENAME);
	PDB refPDB;
	readPDB(refPDBFilename, &refPDB);
	int i, j;
	harmonicConstraintsData.maxRelativeConstraintsPerAtom = 0;
	int fixedConstraintsCount = 0;
	int relativeConstraintsCount = 0;
	char group[PARAMETER_LENGTH];
	getMaskedParameter(group, PARAMETER_CONSTRAINTS_GROUP, PARAMETER_VALUE_CONSTRAINTS_GROUP_NONE);
	int groupOn;
	int groupIn;
	if(strcmp(group, PARAMETER_VALUE_CONSTRAINTS_GROUP_NONE) == 0){
		groupOn = 0;
		groupIn = -1;
	} else {
		groupOn = 1;
		if(strcmp(group, PARAMETER_VALUE_CONSTRAINTS_GROUP_IN) == 0){
			groupIn = 1;
		} else if(strcmp(group, PARAMETER_VALUE_CONSTRAINTS_GROUP_BW) == 0){
			groupIn = 0;
		} else {
			DIE("Constraints group should be %s, %s, or %s", PARAMETER_VALUE_CONSTRAINTS_GROUP_NONE,
					PARAMETER_VALUE_CONSTRAINTS_GROUP_IN, PARAMETER_VALUE_CONSTRAINTS_GROUP_BW);
		}
	}
	float cutoff = getFloatParameter(PARAMETER_CONSTRAINTS_CUTOFF, 0.0f);

	allocateCPU((void**)&harmonicConstraintsData.h_relativeConstraintsCount, gsystem.Ntot*sizeof(int));
	allocateGPU((void**)&harmonicConstraintsData.d_relativeConstraintsCount, gsystem.Ntot*sizeof(int));

	for(i = 0; i < refPDB.atomCount; i++){
		if(refPDB.atoms[i].occupancy == 1.0){
			harmonicConstraintsData.h_fixedConstraints[i].x = refPDB.atoms[i].x/10.0f;
			harmonicConstraintsData.h_fixedConstraints[i].y = refPDB.atoms[i].y/10.0f;
			harmonicConstraintsData.h_fixedConstraints[i].z = refPDB.atoms[i].z/10.0f;
			if(refPDB.atoms[i].beta != 0){
				harmonicConstraintsData.h_fixedConstraints[i].w = refPDB.atoms[i].beta;
			} else {
				harmonicConstraintsData.h_fixedConstraints[i].w = harmonicConstraintsData.Ks;
			}
			fixedConstraintsCount++;
		} else {
			harmonicConstraintsData.h_fixedConstraints[i].w = 0.0f;
		}
	}

	float4 dr;

	for(i = 0; i < refPDB.atomCount; i++){
		harmonicConstraintsData.h_relativeConstraintsCount[i] = 0;
	}
	for(i = 0; i < refPDB.atomCount; i++){
		if(refPDB.atoms[i].occupancy == 2.0){
			for(j = i + 1; j < refPDB.atomCount; j++){
				if(refPDB.atoms[j].occupancy == 2.0){
					if((groupOn == 0) ||
						(groupOn == 1 && groupIn == 1 && refPDB.atoms[i].beta == refPDB.atoms[j].beta) ||
						(groupOn == 1 && groupIn == 0 && refPDB.atoms[i].beta != refPDB.atoms[j].beta)){

							dr.x = (refPDB.atoms[i].x - refPDB.atoms[j].x)/10.0f;
							dr.y = (refPDB.atoms[i].y - refPDB.atoms[j].y)/10.0f;
							dr.z = (refPDB.atoms[i].z - refPDB.atoms[j].z)/10.0f;
							dr.w = abs(dr);

							if(cutoff == 0.0f || dr.w < cutoff){
								harmonicConstraintsData.h_relativeConstraintsCount[i] ++;
								harmonicConstraintsData.h_relativeConstraintsCount[j] ++;
							}
					}
				}
			}
		}
	}


/*	for(i = 0; i < refPDB.atomCount; i++){
		harmonicConstraintsData.h_relativeConstraintsCount[i] = 0;
		if(refPDB.atoms[i].occupancy == 2.0){
			if(groupOn == 0){
				harmonicConstraintsData.h_relativeConstraintsCount[i] ++;
				relativeConstraintsCount ++;
			} else {
				for(j = 0; j < refPDB.atomCount; j++){
					if(j != i && refPDB.atoms[j].occupancy == 2.0){
						if(groupIn == 1 && refPDB.atoms[i].beta == refPDB.atoms[j].beta){
							harmonicConstraintsData.h_relativeConstraintsCount[i]++;
							relativeConstraintsCount ++;
						} else
						if(groupIn == 0 && refPDB.atoms[i].beta != refPDB.atoms[j].beta){
							harmonicConstraintsData.h_relativeConstraintsCount[i]++;
							relativeConstraintsCount ++;
						}

					}
				}
			}
		}
	}*/

	for(i = 0; i < gsystem.N; i++){
		if(harmonicConstraintsData.maxRelativeConstraintsPerAtom < harmonicConstraintsData.h_relativeConstraintsCount[i]){
			harmonicConstraintsData.maxRelativeConstraintsPerAtom = harmonicConstraintsData.h_relativeConstraintsCount[i];
		}
	}

	LOG << "Found " << fixedConstraintsCount << " fixed constraints and " << relativeConstraintsCount << " relative.";

	allocateCPU((void**)&harmonicConstraintsData.h_relativeConstraints,
			harmonicConstraintsData.maxRelativeConstraintsPerAtom*gsystem.widthTot*sizeof(int));
	allocateGPU((void**)&harmonicConstraintsData.d_relativeConstraints,
			harmonicConstraintsData.maxRelativeConstraintsPerAtom*gsystem.widthTot*sizeof(int));

	allocateCPU((void**)&harmonicConstraintsData.h_relativeConstraintsData,
			harmonicConstraintsData.maxRelativeConstraintsPerAtom*gsystem.widthTot*sizeof(int));
	allocateGPU((void**)&harmonicConstraintsData.d_relativeConstraintsData,
			harmonicConstraintsData.maxRelativeConstraintsPerAtom*gsystem.widthTot*sizeof(int));


	LOG << "Relative constraints:";

	for(i = 0; i < refPDB.atomCount; i++){
		harmonicConstraintsData.h_relativeConstraintsCount[i] = 0;
	}
	for(i = 0; i < refPDB.atomCount; i++){
		if(refPDB.atoms[i].occupancy == 2.0){
			for(j = i + 1; j < refPDB.atomCount; j++){
				if(refPDB.atoms[j].occupancy == 2.0){
					if((groupOn == 0) ||
						(groupOn == 1 && groupIn == 1 && refPDB.atoms[i].beta == refPDB.atoms[j].beta) ||
						(groupOn == 1 && groupIn == 0 && refPDB.atoms[i].beta != refPDB.atoms[j].beta)){

							dr.x = (refPDB.atoms[i].x - refPDB.atoms[j].x)/10.0f;
							dr.y = (refPDB.atoms[i].y - refPDB.atoms[j].y)/10.0f;
							dr.z = (refPDB.atoms[i].z - refPDB.atoms[j].z)/10.0f;
							dr.w = abs(dr);

							if(cutoff == 0.0f || dr.w < cutoff){

								LOG << refPDB.atoms[i].id << " " << refPDB.atoms[i].resName << refPDB.atoms[i].resid << "(" << refPDB.atoms[i].name << ")" <<
										" - " << refPDB.atoms[j].id << " " << refPDB.atoms[j].resName << refPDB.atoms[j].resid << "(" << refPDB.atoms[j].name << ")";

								harmonicConstraintsData.h_relativeConstraints[harmonicConstraintsData.h_relativeConstraintsCount[i]*gsystem.widthTot + i] = j;
								harmonicConstraintsData.h_relativeConstraints[harmonicConstraintsData.h_relativeConstraintsCount[j]*gsystem.widthTot + j] = i;

								harmonicConstraintsData.h_relativeConstraintsData[harmonicConstraintsData.h_relativeConstraintsCount[i]*gsystem.widthTot + i] = dr.w;
								harmonicConstraintsData.h_relativeConstraintsData[harmonicConstraintsData.h_relativeConstraintsCount[j]*gsystem.widthTot + j] = dr.w;

								harmonicConstraintsData.h_relativeConstraintsCount[i] ++;
								harmonicConstraintsData.h_relativeConstraintsCount[j] ++;

							}
					}
				}
			}
		}
	}

	cudaMemcpy(harmonicConstraintsData.d_fixedConstraints, harmonicConstraintsData.h_fixedConstraints,
			gsystem.Ntot*sizeof(float4), cudaMemcpyHostToDevice);

	cudaMemcpy(harmonicConstraintsData.d_relativeConstraintsCount, harmonicConstraintsData.h_relativeConstraintsCount,
			gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(harmonicConstraintsData.d_relativeConstraints, harmonicConstraintsData.h_relativeConstraints,
			harmonicConstraintsData.maxRelativeConstraintsPerAtom*gsystem.widthTot*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(harmonicConstraintsData.d_relativeConstraintsData, harmonicConstraintsData.h_relativeConstraintsData,
			harmonicConstraintsData.maxRelativeConstraintsPerAtom*gsystem.widthTot*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpyToSymbol(c_harmonicConstraintsData, &harmonicConstraintsData,
			sizeof(HarmonicConstraintsData), 0, cudaMemcpyHostToDevice);

	LOG << "Done initializing harmonic constraints.";
}

__global__ void harmonicConstraintsPotential_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 coord = tex1Dfetch(t_coord, d_i);
		float4 f = c_gsystem.d_forces[d_i];
		float4 r2 = c_harmonicConstraintsData.d_fixedConstraints[d_i];
		float r;
		if(r2.w != 0.0f){
			r2.x -= coord.x;
			r2.y -= coord.y;
			r2.z -= coord.z;
			DO_PBC(r2);
			r = abs(r2);
			f.x += r2.w*r2.x;
			f.y += r2.w*r2.y;
			f.z += r2.w*r2.z;
		}
		int i;
		for(i = 0; i < c_harmonicConstraintsData.d_relativeConstraintsCount[d_i]; i++){
			int j = c_harmonicConstraintsData.d_relativeConstraints[i*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, j);//c_gsystem.d_coord[bond.j];
			r2.w = c_harmonicConstraintsData.d_relativeConstraintsData[i*c_gsystem.widthTot + d_i];
			r2.x -= coord.x;
			r2.y -= coord.y;
			r2.z -= coord.z;
			DO_PBC(r2);
			r = abs(r2);
			float mult = c_harmonicConstraintsData.Ks*(r - r2.w) / r;
			f.x += mult*r2.x;
			f.y += mult*r2.y;
			f.z += mult*r2.z;
		}
		c_gsystem.d_forces[d_i] = f;
	}
}

__global__ void harmonicConstraintsLimitedForcePotential_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 coord = tex1Dfetch(t_coord, d_i);
		float4 f = c_gsystem.d_forces[d_i];
		float4 r2 = c_harmonicConstraintsData.d_fixedConstraints[d_i];
		float r;
		if(r2.w != 0.0f){
			r2.x -= coord.x;
			r2.y -= coord.y;
			r2.z -= coord.z;
			DO_PBC(r2);
			r = abs(r2);
			f.x += r2.w*r2.x;
			f.y += r2.w*r2.y;
			f.z += r2.w*r2.z;
		}
		int i;
		for(i = 0; i < c_harmonicConstraintsData.d_relativeConstraintsCount[d_i]; i++){
			int j = c_harmonicConstraintsData.d_relativeConstraints[i*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, j);//c_gsystem.d_coord[bond.j];
			r2.w = c_harmonicConstraintsData.d_relativeConstraintsData[i*c_gsystem.widthTot + d_i];
			r2.x -= coord.x;
			r2.y -= coord.y;
			r2.z -= coord.z;
			DO_PBC(r2);
			r = abs(r2);
			float mult = c_harmonicConstraintsData.Ks*(r - r2.w);
			if(mult > c_harmonicConstraintsData.maxF){
				mult = c_harmonicConstraintsData.maxF;
			}
			mult /= r;
			f.x += mult*r2.x;
			f.y += mult*r2.y;
			f.z += mult*r2.z;
		}
		c_gsystem.d_forces[d_i] = f;
	}
}

inline void compute(){
	if(harmonicConstraintsData.limitForce){
		harmonicConstraintsLimitedForcePotential_kernel<<<harmonicConstraintsBlockCount, harmonicConstraintsBlockSize>>>();
	} else {
		harmonicConstraintsPotential_kernel<<<harmonicConstraintsBlockCount, harmonicConstraintsBlockSize>>>();
	}
}

void destroy(){

}

#undef LOG

} // namespace harmonic_constraints
