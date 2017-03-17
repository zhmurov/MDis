/*
 * NonBondedPotential.cu
 *
 *  Created on: Aug 5, 2010
 *      Author: zhmurov
 */

#include "../Core/global.h"
#include "../Core/md.cuh"
#include "../Util/Log.h"
#include "NonBondedPotential.cuh"
#include "../Updaters/PairsListsUpdater.cuh"

// This optimization is now completely disabled, but why throw away good code?
#ifdef CUDA_USE_L1
# define CUDA_NB_USE_L1
#endif

namespace non_bonded_potential {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<non_bonded_potential> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void create(){
	char nbType[100];
	getMaskedParameter(nbType, PARAMETER_NB_TYPE);
	if(strcmp(nbType, PARAMETER_VALUE_NB_TYPE_RDIESHIFT) == 0){
		nbPotentialType = NB_POTENTIAL_TYPE_RDIE_SHIFT;
	} else
	if(strcmp(nbType, PARAMETER_VALUE_NB_TYPE_CDIESWITCH) == 0){
		nbPotentialType = NB_POTENTIAL_TYPE_CDIE_SWITCH;
	} else {
		DIE("Only RDie/shift and CDie/switch non-bonded potentials are supported.");
	}
	if(nbPotentialType == NB_POTENTIAL_TYPE_RDIE_SHIFT){
		nonBondedPotential.compute = &computeRDieShift;
	} else
	if(nbPotentialType == NB_POTENTIAL_TYPE_CDIE_SWITCH){
		nonBondedPotential.compute = &computeCDieSwitch;
	}
	nonBondedPotential.destroy = &destroy;
	sprintf(nonBondedPotential.name, "Non-bonded potential");
	potentials[potentialsCount] = &nonBondedPotential;
	potentialsCount ++;
	if(nbPotentialType == NB_POTENTIAL_TYPE_RDIE_SHIFT){
		ljEnergyOutput.computeValues = &computeLJPotentialEnergyRDieShift;
	} else
	if(nbPotentialType == NB_POTENTIAL_TYPE_CDIE_SWITCH){
		ljEnergyOutput.computeValues = &computeLJPotentialEnergyCDieSwitch;
	}
	allocateCPU((void**)&ljEnergyOutput.values, parameters.Ntr*sizeof(float));
	strcpy(ljEnergyOutput.name, ENERGY_OUTPUT_NAME_LENNARD_JONES);
	energyOutputs[energyOutputsCount] = &ljEnergyOutput;
	energyOutputsCount ++;
	if(nbPotentialType == NB_POTENTIAL_TYPE_RDIE_SHIFT){
		coulombEnergyOutput.computeValues = &computeCoulombPotentialEnergyRDieShift;
	} else
	if(nbPotentialType == NB_POTENTIAL_TYPE_CDIE_SWITCH){
		coulombEnergyOutput.computeValues = &computeCoulombPotentialEnergyCDieSwitch;
	}

	allocateCPU((void**)&coulombEnergyOutput.values, parameters.Ntr*sizeof(float));
	strcpy(coulombEnergyOutput.name, ENERGY_OUTPUT_NAME_ELECTROSTATIC);
	energyOutputs[energyOutputsCount] = &coulombEnergyOutput;
	energyOutputsCount ++;
	init();
}


void init(){

	LOG << "Initializing non-bonded potential...";

	nonBondedBlockSize = BLOCK_SIZE;
	nonBondedBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;

	allocateCPU((void**)&nonBondedData.h_ljEnergies, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&nonBondedData.d_ljEnergies, gsystem.Ntot*sizeof(float));
	allocateCPU((void**)&nonBondedData.h_coulombEnergies, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&nonBondedData.d_coulombEnergies, gsystem.Ntot*sizeof(float));

	allocateCPU((void**)&nonBondedData.h_charges, atomTypesCount*sizeof(float));
	allocateCPU((void**)&nonBondedData.h_ljparameters, atomTypesCount*sizeof(float2));
	allocateCPU((void**)&nonBondedData.h_ljEpsilon, atomTypesCount*sizeof(float));
	allocateCPU((void**)&nonBondedData.h_ljR0, atomTypesCount*sizeof(float));
	allocateGPU((void**)&nonBondedData.d_charges, atomTypesCount*sizeof(float));
	allocateGPU((void**)&nonBondedData.d_ljparameters, atomTypesCount*sizeof(float2));
	allocateGPU((void**)&nonBondedData.d_ljEpsilon, atomTypesCount*sizeof(float));
	allocateGPU((void**)&nonBondedData.d_ljR0, atomTypesCount*sizeof(float));

	allocateCPU((void**)&nonBondedData.h_charges14, atomTypesCount*sizeof(float));
	allocateCPU((void**)&nonBondedData.h_ljparameters14, atomTypesCount*sizeof(float2));
	allocateCPU((void**)&nonBondedData.h_ljEpsilon14, atomTypesCount*sizeof(float));
	allocateCPU((void**)&nonBondedData.h_ljR014, atomTypesCount*sizeof(float));
	allocateGPU((void**)&nonBondedData.d_charges14, atomTypesCount*sizeof(float));
	allocateGPU((void**)&nonBondedData.d_ljparameters14, atomTypesCount*sizeof(float2));
	allocateGPU((void**)&nonBondedData.d_ljEpsilon14, atomTypesCount*sizeof(float));
	allocateGPU((void**)&nonBondedData.d_ljR014, atomTypesCount*sizeof(float));

	nonBondedData.ljCutoff = getFloatParameter(PARAMETER_LJ_CUTOFF, -1.0f);
	if(nonBondedData.ljCutoff == -1.0f){
		nonBondedData.ljCutoff = getFloatParameter(PARAMETER_NB_CUTOFF, DEFAULT_LJ_CUTOFF);
	}
	nonBondedData.ljSwitch = getFloatParameter(PARAMETER_LJ_SWITCH, -1.0f);
	if(nonBondedData.ljSwitch == -1.0f){
		nonBondedData.ljSwitch = getFloatParameter(PARAMETER_NB_SWITCH, DEFAULT_LJ_SWITCH);
	}
	nonBondedData.coulombCutoff = getFloatParameter(PARAMETER_COULOMB_CUTOFF, -1.0f);
	if(nonBondedData.coulombCutoff == -1.0f){
		nonBondedData.coulombCutoff = getFloatParameter(PARAMETER_NB_CUTOFF, DEFAULT_COULOMB_CUTOFF);
	}

	nonBondedData.ljCutoff2 = nonBondedData.ljCutoff*nonBondedData.ljCutoff;
	nonBondedData.ljSwitch2 = nonBondedData.ljSwitch*nonBondedData.ljSwitch;
	nonBondedData.coulombCutoff2 = nonBondedData.coulombCutoff*nonBondedData.coulombCutoff;

	nonBondedData.oneOverLjCutoff6 = 1.0f/(nonBondedData.ljCutoff2*nonBondedData.ljCutoff2*nonBondedData.ljCutoff2);
	nonBondedData.oneOverCoulombCutoff4 = 1.0f/(nonBondedData.coulombCutoff2*nonBondedData.coulombCutoff2);

	nonBondedData.oneOverCutoff2MinSwitch2Cub = 1.0f/(nonBondedData.ljCutoff2 - nonBondedData.ljSwitch2);
	nonBondedData.oneOverCutoff2MinSwitch2Cub = nonBondedData.oneOverCutoff2MinSwitch2Cub *
			nonBondedData.oneOverCutoff2MinSwitch2Cub *
			nonBondedData.oneOverCutoff2MinSwitch2Cub;

	nonBondedData.C1 = (-nonBondedData.ljCutoff2 + 3.0f*nonBondedData.ljSwitch2)*nonBondedData.ljCutoff2;
	nonBondedData.C2 = -nonBondedData.ljCutoff2 + 9.0f*nonBondedData.ljSwitch2;
	nonBondedData.C3 = COULOMB_CONSTANT*nonBondedData.oneOverCutoff2MinSwitch2Cub;

	int i;
	float scaleE14 = getFloatParameter(PARAMETER_SCALE_14_FACTOR, 1.0);

	DPRINTF("Non-bonded Parameters: \n");
	DPRINTF("Atom type:\tCharge:\t\tEpsilon:\tR0/2:\t\tCharge(1-4):\tEpsilon(1-4):\tR0(1-4)/2:\n");
	for(i = 0; i < atomTypesCount; i++){
		nonBondedData.h_charges[i] = atomTypes[i].charge;
		nonBondedData.h_ljparameters[i].x = sqrtf(12.0f*abs(atomTypes[i].epsilon));
		nonBondedData.h_ljparameters[i].y = atomTypes[i].RminOver2;
		nonBondedData.h_charges14[i] = atomTypes[i].charge*sqrtf(scaleE14);
		nonBondedData.h_ljparameters14[i].x = sqrtf(12.0f*abs(atomTypes[i].epsilon_14));
		nonBondedData.h_ljparameters14[i].y = atomTypes[i].RminOver2_14;

		nonBondedData.h_ljEpsilon[i] = nonBondedData.h_ljparameters[i].x;
		nonBondedData.h_ljR0[i] = nonBondedData.h_ljparameters[i].y;
		nonBondedData.h_ljEpsilon14[i] = nonBondedData.h_ljparameters14[i].x;
		nonBondedData.h_ljR014[i] = nonBondedData.h_ljparameters14[i].y;

		/*nonBondedData.h_charges14[i] = atomTypes[i].charge;
		nonBondedData.h_ljparameters14[i].x = atomTypes[i].epsilon;
		nonBondedData.h_ljparameters14[i].y = atomTypes[i].RminOver2;*/

		DPRINTF("%d.%s\t\t%5.4f\t\t%5.4f\t\t%5.4f\t\t%5.4f\t\t%5.4f\t\t%5.4f\n", i, atomTypes[i].name,
				nonBondedData.h_charges[i], nonBondedData.h_ljparameters[i].x, nonBondedData.h_ljparameters[i].y,
				nonBondedData.h_charges14[i], nonBondedData.h_ljparameters14[i].x, nonBondedData.h_ljparameters14[i].y);
	}

	cudaMemcpy(nonBondedData.d_charges, nonBondedData.h_charges, atomTypesCount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(nonBondedData.d_ljparameters, nonBondedData.h_ljparameters, atomTypesCount*sizeof(float2), cudaMemcpyHostToDevice);
	cudaMemcpy(nonBondedData.d_ljEpsilon, nonBondedData.h_ljEpsilon, atomTypesCount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(nonBondedData.d_ljR0, nonBondedData.h_ljR0, atomTypesCount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(nonBondedData.d_charges14, nonBondedData.h_charges14, atomTypesCount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(nonBondedData.d_ljparameters14, nonBondedData.h_ljparameters14, atomTypesCount*sizeof(float2), cudaMemcpyHostToDevice);
	cudaMemcpy(nonBondedData.d_ljEpsilon14, nonBondedData.h_ljEpsilon14, atomTypesCount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(nonBondedData.d_ljR014, nonBondedData.h_ljR014, atomTypesCount*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpyToSymbol(c_charges, nonBondedData.h_charges, atomTypesCount*sizeof(float), 0);
	cudaMemcpyToSymbol(c_ljparameters, nonBondedData.h_ljparameters, atomTypesCount*sizeof(float2), 0);
	cudaMemcpyToSymbol(c_charges14, nonBondedData.h_charges14, atomTypesCount*sizeof(float), 0);
	cudaMemcpyToSymbol(c_ljparameters14, nonBondedData.h_ljparameters14, atomTypesCount*sizeof(float2), 0);

	cudaBindTexture(0, t_charges, nonBondedData.d_charges, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_charges14, nonBondedData.d_charges14, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_ljparameters, nonBondedData.d_ljparameters, atomTypesCount*sizeof(float2));
	cudaBindTexture(0, t_ljparameters14, nonBondedData.d_ljparameters14, atomTypesCount*sizeof(float2));

	cudaBindTexture(0, t_ljEpsilon, nonBondedData.d_ljEpsilon, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_ljR0, nonBondedData.d_ljR0, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_ljEpsilon14, nonBondedData.d_ljEpsilon14, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_ljR014, nonBondedData.d_ljR014, atomTypesCount*sizeof(float));

	/*allocateGPU((void**)&nonBondedData.d_atoms, atomCount*sizeof(Atom));
	cudaMemcpy(nonBondedData.d_atoms, atoms, atomCount*sizeof(Atom), cudaMemcpyHostToDevice);*/

	cudaMemcpyToSymbol(c_nonBondedData, &nonBondedData,
				sizeof(NonBondedData), 0, cudaMemcpyHostToDevice);

	checkCUDAError("Copying data to GPU");
#ifdef CUDA_NB_USE_L1
	cudaFuncSetCacheConfig(nonBondedPotentialRDieShift_kernel, cudaFuncCachePreferL1);
	checkCUDAError("configuring nonBondedPotentialRDieShift_kernel cache policy");
#endif

	LOG << "Done initializing non-bonded potential.";
}

__global__ void nonBondedPotentialRDieShift_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
#ifndef CUDA_NB_USE_L1
		float4 r1 = tex1Dfetch(t_coord, d_i);
#else
		float4 r1 = c_gsystem.d_coord[d_i];
#endif
		float4 f = c_gsystem.d_forces[d_i];
		float4 r2;
		float mult;
		int i;
		int at1, at2;
		at1 = (int)r1.w;
	//	s_x1[threadIdx.x] = tex1Dfetch(t_ljEpsilon14, at1);
	//	s_x2[threadIdx.x] = tex1Dfetch(t_ljR014, at1);
	//	s_x3[threadIdx.x] = tex1Dfetch(t_charges14, at1);
		//float2 par1 = tex1Dfetch(t_ljparameters14, at1);//c_ljparameters14[at1];
		//float2 par1;
		//par1.x = s_x1[threadIdx.x];
		//par1.y = s_x2[threadIdx.x];
		//float2 par2;

		float r;
		for(i = 0; i < c_pairsListsData.d_pairs14ListCount[d_i]; i++){
#ifndef CUDA_NB_USE_L1
			r2 =  tex1Dfetch(t_coord, c_pairsListsData.d_pairs14List[i*c_gsystem.widthTot + d_i]);
#else
			r2 =  c_gsystem.d_coord[c_pairsListsData.d_pairs14List[i*c_gsystem.widthTot + d_i]];
#endif
			at2 = (int)r2.w;
			r2 -= r1;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			mult = 0.0f;
			if(r2.w < c_nonBondedData.ljCutoff2){
				//par2 = tex1Dfetch(t_ljparameters14, at2);//c_ljparameters14[at2];
				r = tex1Dfetch(t_ljR014, at1);//par1.y + par2.y;
				r1.w = tex1Dfetch(t_ljR014, at2);
				r += r1.w;
				r = r*r;
				mult = r/r2.w;
				r1.w = mult*mult*mult;
				mult = (r1.w - r1.w*r1.w)/r2.w;

				r1.w = r/c_nonBondedData.ljCutoff2;
				r1.w = r1.w*r1.w*r1.w;
				mult -= (r1.w - r1.w*r1.w)*r2.w*r2.w*c_nonBondedData.oneOverLjCutoff6;

				r = tex1Dfetch(t_ljEpsilon14, at1);
				r1.w = tex1Dfetch(t_ljEpsilon14, at2);
				r *= r1.w;

				mult *= r;//par1.x*par2.x;
			}

			if(r2.w < c_nonBondedData.coulombCutoff2){

				//par2.x = tex1Dfetch(t_charges14, at1);
				//par2.y = tex1Dfetch(t_charges14, at2);
				r = tex1Dfetch(t_charges14, at1);

				r1.w = 1.0f/r2.w;

				r2.w = tex1Dfetch(t_charges14, at2);

				r1.w = r1.w*r1.w;
				r1.w = -r1.w;
				r1.w += c_nonBondedData.oneOverCoulombCutoff4;

				r *= r2.w;
				r *= 2.0f;
				r *= COULOMB_CONSTANT_RDIE;
				r *= r1.w;

				mult += r;
				//mult += 2.0f*COULOMB_CONSTANT_RDIE*par2.x*par2.y *
				//		(c_nonBondedData.oneOverCoulombCutoff4 - r1.w);

			}

			f.x += mult*r2.x;
			f.y += mult*r2.y;
			f.z += mult*r2.z;
		}
		//par1 = tex1Dfetch(t_ljparameters, at1);//c_ljparameters[at1];
		for(i = 0; i < c_pairsListsData.d_pairsLJListCount[d_i]; i++){
#ifndef CUDA_NB_USE_L1
			r2 =  tex1Dfetch(t_coord, c_pairsListsData.d_pairsLJList[i*c_gsystem.widthTot + d_i]);
#else
			r2 =  c_gsystem.d_coord[c_pairsListsData.d_pairsLJList[i*c_gsystem.widthTot + d_i]];
#endif
			at2 = (int)r2.w;
			//par2 = tex1Dfetch(t_ljparameters, at2);//c_ljparameters[at2];
			r2 -= r1;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			mult = 0.0f;
			if(r2.w < c_nonBondedData.ljCutoff2){
				r = tex1Dfetch(t_ljR0, at1);//par1.y + par2.y;
				r1.w = tex1Dfetch(t_ljR0, at2);
				r += r1.w;
				r = r*r;
				mult = r/r2.w;
				r1.w = mult*mult*mult;
				mult = (r1.w - r1.w*r1.w)/r2.w;

				r1.w = r/c_nonBondedData.ljCutoff2;
				r1.w = r1.w*r1.w*r1.w;
				mult -= (r1.w - r1.w*r1.w)*r2.w*r2.w*c_nonBondedData.oneOverLjCutoff6;

				r = tex1Dfetch(t_ljEpsilon, at1);
				r1.w = tex1Dfetch(t_ljEpsilon, at2);
				r *= r1.w;

				mult *= r;
			}

			if(r2.w < c_nonBondedData.coulombCutoff2){
				/*r1.w = *c_charges[at1]*c_charges[at2]/(r2.w*r2.w);

				r = r2.w / c_nonBondedData.coulombCutoff2;
				r1.w = r1.w * (2.0f - 8.0f*r + 6.0f*r*r2.w);

				mult -= r1.w;*/

			//	par2.x = tex1Dfetch(t_charges, at1);
				//par2.y = tex1Dfetch(t_charges, at2);

				r = tex1Dfetch(t_charges, at1);

				r1.w = 1.0f/r2.w;

				r2.w = tex1Dfetch(t_charges, at2);

				r1.w = r1.w*r1.w;
				r1.w = -r1.w;
				r1.w += c_nonBondedData.oneOverCoulombCutoff4;

				r *= r2.w;
				r *= 2.0f;
				r *= COULOMB_CONSTANT_RDIE;
				r *= r1.w;

				mult += r;

				//mult += 2.0f*COULOMB_CONSTANT_RDIE*par2.x*par2.y*
				//		(c_nonBondedData.oneOverCoulombCutoff4 - r1.w);

			}

			f.x += mult*r2.x;
			f.y += mult*r2.y;
			f.z += mult*r2.z;
		}
		for(i = 0; i < c_pairsListsData.d_pairsCoulombListCount[d_i]; i++){
#ifndef CUDA_NB_USE_L1
			r2 =  tex1Dfetch(t_coord, c_pairsListsData.d_pairsCoulombList[i*c_gsystem.widthTot + d_i]);
#else
			r2 =  c_gsystem.d_coord[c_pairsListsData.d_pairsCoulombList[i*c_gsystem.widthTot + d_i]];
#endif
			at2 = (int)r2.w;
			r2 -= r1;
			DO_PBC(r2);
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;

			if(r2.w < c_nonBondedData.coulombCutoff2){

				r = tex1Dfetch(t_charges, at1);

				r1.w = 1.0f/r2.w;

				r2.w = tex1Dfetch(t_charges, at2);

				r1.w = r1.w*r1.w;
				r1.w = -r1.w;
				r1.w += c_nonBondedData.oneOverCoulombCutoff4;

				r *= r2.w;
				r *= 2.0f;
				r *= COULOMB_CONSTANT_RDIE;
				r *= r1.w;

				mult += r;

				//mult = 2.0f*COULOMB_CONSTANT_RDIE*par2.x*par2.y*
				//		(c_nonBondedData.oneOverCoulombCutoff4 - r1.w);

			}

			f.x += mult*r2.x;
			f.y += mult*r2.y;
			f.z += mult*r2.z;
		}
		c_gsystem.d_forces[d_i] = f;
	}
}

inline void computeRDieShift(){
//!!
/*
	cudaBindTexture(0, t_charges, nonBondedData.d_charges, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_charges14, nonBondedData.d_charges14, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_ljparameters, nonBondedData.d_ljparameters, atomTypesCount*sizeof(float2));
	cudaBindTexture(0, t_ljparameters14, nonBondedData.d_ljparameters14, atomTypesCount*sizeof(float2));

	cudaBindTexture(0, t_ljEpsilon, nonBondedData.d_ljEpsilon, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_ljR0, nonBondedData.d_ljR0, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_ljEpsilon14, nonBondedData.d_ljEpsilon14, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_ljR014, nonBondedData.d_ljR014, atomTypesCount*sizeof(float));
*/
	nonBondedPotentialRDieShift_kernel<<<nonBondedBlockCount, nonBondedBlockSize>>>();
/*
	cudaUnbindTexture(t_charges);
	cudaUnbindTexture(t_charges14);
	cudaUnbindTexture(t_ljparameters);
	cudaUnbindTexture(t_ljparameters14);

	cudaUnbindTexture(t_ljEpsilon);
	cudaUnbindTexture(t_ljR0);
	cudaUnbindTexture(t_ljEpsilon14);
	cudaUnbindTexture(t_ljR014);
*/
//!!
	//checkCUDAError("lj14");
	//nonBondedPotential_kernel<<<nonBondedBlockCount, nonBondedBlockSize>>>();
	//checkCUDAError("lj");
	//nonBondedCoulomb_kernel<<<nonBondedBlockCount, nonBondedBlockSize>>>();
	//checkCUDAError("coulomb");
	/*cudaMemcpy(gsystem.h_forces, gsystem.d_forces, atomCount*sizeof(float4), cudaMemcpyDeviceToHost);
	int i;
	float3 force = make_float3(0.0f, 0.0f, 0.0f);
	for(i = 0; i < atomCount; i++){
		force.x += gsystem.h_forces[i].x;
		force.y += gsystem.h_forces[i].y;
		force.z += gsystem.h_forces[i].z;
		printf("%d: (%f, %f, %f) %f\n", i, gsystem.h_forces[i].x, gsystem.h_forces[i].y, gsystem.h_forces[i].z,
					sqrtf(gsystem.h_forces[i].x*gsystem.h_forces[i].x + gsystem.h_forces[i].y*gsystem.h_forces[i].y + gsystem.h_forces[i].z*gsystem.h_forces[i].z));
	}
	printf("Net force (non-bonded): (%f, %f, %f) %f\n", force.x, force.y, force.z,
			sqrtf(force.x*force.x + force.y*force.y + force.z*force.z));*/
}

__global__ void nonBonded14PotentialEnergyRDieShift_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 r1 = tex1Dfetch(t_coord, d_i);
		float4 r2;
		int at1, at2;
		at1 = (int)r1.w;
		float2 par1 = tex1Dfetch(t_ljparameters14, at1);
		float2 par2;
		int i;
		float r;
		float potLJ = 0.0f;
		float potCoulomb = 0.0f;
		for(i = 0; i < c_pairsListsData.d_pairs14ListCount[d_i]; i++){
			r2 =  tex1Dfetch(t_coord, c_pairsListsData.d_pairs14List[i*c_gsystem.widthTot + d_i]);
			at2 = (int)r2.w;
			r2 -= r1;
			DO_PBC(r2);
			r = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			if(r < c_nonBondedData.ljCutoff2){
				par2 = tex1Dfetch(t_ljparameters14, at2);
				if(par1.x != 0.0f && par2.x != 0.0f){
					r2.w = par1.y + par2.y;
					r2.w = r2.w*r2.w;

					r1.w = r2.w/r;
					r1.w = r1.w*r1.w*r1.w;
					//potLJ += (r1.w*r1.w - 2.0f*r1.w)*par1.x*par2.x/12.0f;
					par2.y = (r1.w*r1.w - 2.0f*r1.w);

					r1.w = r2.w/c_nonBondedData.ljCutoff2;
					r1.w = r1.w*r1.w*r1.w;

					par2.y += (r1.w*r1.w - r1.w)*2.0f*c_nonBondedData.oneOverLjCutoff6*r*r*r;

					par2.y += (4.0f*r1.w - 3.0f*r1.w*r1.w);

					par2.y *= par1.x*par2.x/12.0f;
					/*if(r > c_nonBondedData.ljSwitch2){
						r2.w *= (c_nonBondedData.ljCutoff2 - r)*(c_nonBondedData.ljCutoff2 - r);
						r2.w *= (c_nonBondedData.ljCutoff2 + 2.0f*r - 3.0f*c_nonBondedData.ljSwitch2);
						r2.w /= (c_nonBondedData.ljCutoff2 - c_nonBondedData.ljSwitch2);
						r2.w /= (c_nonBondedData.ljCutoff2 - c_nonBondedData.ljSwitch2);
						r2.w /= (c_nonBondedData.ljCutoff2 - c_nonBondedData.ljSwitch2);
					}*/
					potLJ += par2.y;

					/*r1.w = r2.w/c_nonBondedData.ljCutoff2;
					r1.w = r1.w*r1.w*r1.w;
					potLJ -= (r1.w*r1.w - 2.0f*r1.w)*par1.x*par2.x/12.0f;*/
				}
			}
			if(r < c_nonBondedData.coulombCutoff2){
				r1.w = COULOMB_CONSTANT_RDIE*tex1Dfetch(t_charges14, at1)*tex1Dfetch(t_charges14, at2)/r;
				r1.w *= (1.0f - r/c_nonBondedData.coulombCutoff2);
				r1.w *= (1.0f - r/c_nonBondedData.coulombCutoff2);
				potCoulomb += r1.w;
			}

		}
		c_nonBondedData.d_ljEnergies[d_i] = potLJ;
		c_nonBondedData.d_coulombEnergies[d_i] = potCoulomb;
	}
}

__global__ void nonBondedPotentialEnergyRDieShift_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 r1 = tex1Dfetch(t_coord, d_i);
		float4 r2;
		int at1, at2;
		at1 = (int)r1.w;
		float2 par1 = tex1Dfetch(t_ljparameters, at1);
		float2 par2;
		int i;
		float r;
		float potLJ = c_nonBondedData.d_ljEnergies[d_i];
		float potCoulomb = c_nonBondedData.d_coulombEnergies[d_i];
		for(i = 0; i < c_pairsListsData.d_pairsLJListCount[d_i]; i++){
			r2 =  tex1Dfetch(t_coord, c_pairsListsData.d_pairsLJList[i*c_gsystem.widthTot + d_i]);
			at2 = (int)r2.w;
			r2 -= r1;
			DO_PBC(r2);
			r = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			if(r < c_nonBondedData.ljCutoff2){
				par2 = tex1Dfetch(t_ljparameters, at2);
				r2.w = par1.y + par2.y;
				r2.w = r2.w*r2.w;
				r1.w = r2.w/r;
				r1.w = r1.w*r1.w*r1.w;
				//potLJ += (r1.w*r1.w - 2.0f*r1.w)*par1.x*par2.x/12.0f;
				par2.y = (r1.w*r1.w - 2.0f*r1.w);

				r1.w = r2.w/c_nonBondedData.ljCutoff2;
				r1.w = r1.w*r1.w*r1.w;

				par2.y += (r1.w*r1.w - r1.w)*2.0f*c_nonBondedData.oneOverLjCutoff6*r*r*r;

				par2.y += (4.0f*r1.w - 3.0f*r1.w*r1.w);

				par2.y *= par1.x*par2.x/12.0f;

				//r2.w = (r1.w*r1.w - 2.0f*r1.w)*par1.x*par2.x/12.0f;
				/*if(r > c_nonBondedData.ljSwitch2){
					r2.w *= (c_nonBondedData.ljCutoff2 - r)*(c_nonBondedData.ljCutoff2 - r);
					r2.w *= (c_nonBondedData.ljCutoff2 + 2.0f*r - 3.0f*c_nonBondedData.ljSwitch2);
					r2.w /= (c_nonBondedData.ljCutoff2 - c_nonBondedData.ljSwitch2);
					r2.w /= (c_nonBondedData.ljCutoff2 - c_nonBondedData.ljSwitch2);
					r2.w /= (c_nonBondedData.ljCutoff2 - c_nonBondedData.ljSwitch2);
				}*/
				potLJ += par2.y;
				/*r1.w = r2.w/c_nonBondedData.ljCutoff2;
				r1.w = r1.w*r1.w*r1.w;
				potLJ -= (r1.w*r1.w - 2.0f*r1.w)*par1.x*par2.x/12.0f;*/
			}
			if(r < c_nonBondedData.coulombCutoff2){
				r1.w = COULOMB_CONSTANT_RDIE*tex1Dfetch(t_charges, at1)*tex1Dfetch(t_charges, at2)/r;
				r1.w *= (1.0f - r/c_nonBondedData.coulombCutoff2);
				r1.w *= (1.0f - r/c_nonBondedData.coulombCutoff2);
				potCoulomb += r1.w;
			}

		}
		c_nonBondedData.d_ljEnergies[d_i] = potLJ;
		c_nonBondedData.d_coulombEnergies[d_i] = potCoulomb;
	}
}

__global__ void nonBondedCoulombEnergyRDieShift_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 r1 = tex1Dfetch(t_coord, d_i);
		float4 r2;
		int at1, at2;
		at1 = (int)r1.w;
		int i;
		float r;
		float potCoulomb = c_nonBondedData.d_coulombEnergies[d_i];
		for(i = 0; i < c_pairsListsData.d_pairsCoulombListCount[d_i]; i++){
			r2 =  tex1Dfetch(t_coord, c_pairsListsData.d_pairsCoulombList[i*c_gsystem.widthTot + d_i]);
			at2 = (int)r2.w;
			r2 -= r1;
			DO_PBC(r2);
			r = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;

			if(r < c_nonBondedData.coulombCutoff2){
				r1.w = COULOMB_CONSTANT_RDIE*tex1Dfetch(t_charges, at1)*tex1Dfetch(t_charges, at2)/r;
				r1.w *= (1.0f - r/c_nonBondedData.coulombCutoff2);
				r1.w *= (1.0f - r/c_nonBondedData.coulombCutoff2);
				potCoulomb += r1.w;
			}

		}
		c_nonBondedData.d_coulombEnergies[d_i] = potCoulomb;
	}
}

inline void computeLJPotentialEnergyRDieShift(){
//!!
/*
	cudaBindTexture(0, t_charges, nonBondedData.d_charges, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_charges14, nonBondedData.d_charges14, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_ljparameters, nonBondedData.d_ljparameters, atomTypesCount*sizeof(float2));
	cudaBindTexture(0, t_ljparameters14, nonBondedData.d_ljparameters14, atomTypesCount*sizeof(float2));

	cudaBindTexture(0, t_ljEpsilon, nonBondedData.d_ljEpsilon, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_ljR0, nonBondedData.d_ljR0, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_ljEpsilon14, nonBondedData.d_ljEpsilon14, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_ljR014, nonBondedData.d_ljR014, atomTypesCount*sizeof(float));
*/
	nonBonded14PotentialEnergyRDieShift_kernel<<<nonBondedBlockCount, nonBondedBlockSize>>>();
	checkCUDAError("LJ 1-4 energy");
	nonBondedPotentialEnergyRDieShift_kernel<<<nonBondedBlockCount, nonBondedBlockSize>>>();
	checkCUDAError("LJ energy");
	nonBondedCoulombEnergyRDieShift_kernel<<<nonBondedBlockCount, nonBondedBlockSize>>>();
	checkCUDAError("Coulomb energy");
/*
	cudaUnbindTexture(t_charges);
	cudaUnbindTexture(t_charges14);
	cudaUnbindTexture(t_ljparameters);
	cudaUnbindTexture(t_ljparameters14);

	cudaUnbindTexture(t_ljEpsilon);
	cudaUnbindTexture(t_ljR0);
	cudaUnbindTexture(t_ljEpsilon14);
	cudaUnbindTexture(t_ljR014);
*/
//!!

	cudaMemcpy(nonBondedData.h_ljEnergies, nonBondedData.d_ljEnergies,
				gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(nonBondedData.h_coulombEnergies, nonBondedData.d_coulombEnergies,
				gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		float pot = 0.0f;
		for(i = 0; i < gsystem.N; i++){
			/*if(nonBondedData.h_ljEnergies[i] > 0.0f){
				printf("Atom %d has another one close to it (LJ Energy = %f)\n", i, nonBondedData.h_ljEnergies[i]);
			}*/
			pot += nonBondedData.h_ljEnergies[i + traj*gsystem.N];
		}
		pot /= 2.0f;
		ljEnergyOutput.values[traj] = pot;
	}
}

inline void computeCoulombPotentialEnergyRDieShift(){
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		float pot = 0.0f;
		for(i = 0; i < gsystem.N; i++){
			pot += nonBondedData.h_coulombEnergies[i + traj*gsystem.N];
		}
		pot /= 2.0f;
		coulombEnergyOutput.values[traj] = pot;
	}
}

__global__ void nonBondedPotentialCDieSwitch_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 r1 = tex1Dfetch(t_coord, d_i);
		float4 f = c_gsystem.d_forces[d_i];
		float4 r2;
		float mult;
		int at1, at2;
		at1 = (int)r1.w;
		float2 par1 = tex1Dfetch(t_ljparameters14, at1);//c_ljparameters14[at1];
		float2 par2;
		int i;
		float r;
		for(i = 0; i < c_pairsListsData.d_pairs14ListCount[d_i]; i++){
			r2 =  tex1Dfetch(t_coord, c_pairsListsData.d_pairs14List[i*c_gsystem.widthTot + d_i]);
			at2 = (int)r2.w;
			r2 -= r1;
			DO_PBC(r2);
			r = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			mult = 0.0f;
			float m1, m2, m3;
			if(r < c_nonBondedData.ljCutoff2){
				if(r > c_nonBondedData.ljSwitch2){
					m2 = c_nonBondedData.ljCutoff2 - r;
					m3 = r - c_nonBondedData.ljSwitch2;
					m1 = m2 + 3.0f*m3;
					m2 *= m1;
					m1 = (c_nonBondedData.ljCutoff2 - r)*c_nonBondedData.oneOverCutoff2MinSwitch2Cub;
				} else {
					m1 = 1.0f;
					m2 = 1.0f;
					m3 = 0.0f;
				}
				par2 = tex1Dfetch(t_ljparameters14, at2);//c_ljparameters14[at2];
				if(par1.x != 0.0f && par2.x != 0.0f){
					r2.w = par1.y + par2.y;
					r2.w = r2.w*r2.w;
					mult = r2.w/r;
					r1.w = mult*mult*mult;
					mult = m2*(r1.w - r1.w*r1.w)/r + m3*(2.0f*r1.w - r1.w*r1.w);

					mult *= par1.x*par2.x*m1;
				} else {
					mult = 0.0f;
				}

				par2.x = tex1Dfetch(t_charges14, at1);
				par2.y = tex1Dfetch(t_charges14, at2);
				if(par2.x != 0.0f && par2.y != 0.0f){
					if(r > c_nonBondedData.ljSwitch2){
						r1.w = c_nonBondedData.C2-10.0f*r;
						r1.w *= r;
						r1.w += c_nonBondedData.C1;
						r2.w = c_nonBondedData.ljCutoff2 - r;
						r1.w *= r2.w;
						r2.w = 1.0f/r;
						mult += par2.x*par2.y*r2.w*sqrtf(r2.w)*r1.w*c_nonBondedData.C3;
					} else {
						r2.w = 1.0f/r;
						mult -= COULOMB_CONSTANT*par2.x*par2.y*r2.w*sqrtf(r2.w);
					}
					//r1.w = 1.0f/r;
					//r1.w = (m2*r1.w + 12.0f*m3)*sqrtf(r1.w);

					//mult -= COULOMB_CONSTANT*par2.x*par2.y*r1.w*m1;
				}
				/*if(c_charges14[at1] != 0.0f && c_charges14[at2] != 0.0f){
					r1.w = 1.0f/r;
					r1.w = r1.w*sqrtf(r1.w);

					mult -= COULOMB_CONSTANT*c_charges14[at1]*c_charges14[at2]*r1.w;
				}*/

				f.x += mult*r2.x;
				f.y += mult*r2.y;
				f.z += mult*r2.z;
			}
		}
		par1 = tex1Dfetch(t_ljparameters, at1);//c_ljparameters[at1];
		for(i = 0; i < c_pairsListsData.d_pairsLJListCount[d_i]; i++){
			r2 =  tex1Dfetch(t_coord, c_pairsListsData.d_pairsLJList[i*c_gsystem.widthTot + d_i]);
			at2 = (int)r2.w;
			r2 -= r1;
			DO_PBC(r2);
			r = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			mult = 0.0f;
			float m1, m2, m3;
			if(r < c_nonBondedData.ljCutoff2){
				if(r > c_nonBondedData.ljSwitch2){
					m2 = c_nonBondedData.ljCutoff2 - r;
					m3 = r - c_nonBondedData.ljSwitch2;
					m1 = m2 + 3.0f*m3;
					m2 *= m1;
					m1 = (c_nonBondedData.ljCutoff2 - r)*c_nonBondedData.oneOverCutoff2MinSwitch2Cub;
				} else {
					m1 = 1.0f;
					m2 = 1.0f;
					m3 = 0.0f;
				}
				par2 = tex1Dfetch(t_ljparameters, at2);//c_ljparameters[at2];
				if(par1.x != 0.0f && par2.x != 0.0f){
					r2.w = par1.y + par2.y;
					r2.w = r2.w*r2.w;
					mult = r2.w/r;
					r1.w = mult*mult*mult;
					mult = m2*(r1.w - r1.w*r1.w)/r + m3*(2.0f*r1.w - r1.w*r1.w);

					mult *= par1.x*par2.x*m1;
				} else {
					mult = 0.0f;
				}

				par2.x = tex1Dfetch(t_charges, at1);
				par2.y = tex1Dfetch(t_charges, at2);
				if(par2.x != 0.0f && par2.y != 0.0f){
					if(r > c_nonBondedData.ljSwitch2){
						r1.w = c_nonBondedData.C2-10.0f*r;
						r1.w *= r;
						r1.w += c_nonBondedData.C1;
						r2.w = c_nonBondedData.ljCutoff2 - r;
						r1.w *= r2.w;
						r2.w = 1.0f/r;
						mult += par2.x*par2.y*r2.w*sqrtf(r2.w)*r1.w*c_nonBondedData.C3;
					} else {
						r2.w = 1.0f/r;
						mult -= COULOMB_CONSTANT*par2.x*par2.y*r2.w*sqrtf(r2.w);
					}
					//r1.w = 1.0f/r;
					//r1.w = (m2*r1.w + 12.0f*m3)*sqrtf(r1.w);
					//mult -= COULOMB_CONSTANT*par2.x*par2.y*r1.w*m1;
				}
				/*if(c_charges[at1] != 0.0f && c_charges[at2] != 0.0f){
					r1.w = 1.0f/r;
					r1.w = r1.w*sqrtf(r1.w);

					mult -= COULOMB_CONSTANT*c_charges[at1]*c_charges[at2]*r1.w;
				}*/
				f.x += mult*r2.x;
				f.y += mult*r2.y;
				f.z += mult*r2.z;
			}
		}
		/*for(i = 0; i < c_pairsListsData.d_pairsCoulombListCount[d_i]; i++){
			r2 =  tex1Dfetch(t_coord, c_pairsListsData.d_pairsCoulombList[i*c_gsystem.widthTot + d_i]);
			at2 = (int)r2.w;
			par2 = c_ljparameters[at2];
			r2.x -= r1.x;
			r2.y -= r1.y;
			r2.z -= r1.z;
			r = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			mult = 0.0f;
			if(r < c_nonBondedData.coulombCutoff2){
				if(c_charges[at1] != 0.0f && c_charges[at2] != 0.0f){
					r1.w = 1.0f/r;
					r1.w = r1.w*sqrtf(r1.w);

					mult -= COULOMB_CONSTANT*c_charges[at1]*c_charges[at2]*r1.w;
				}

				f.x += mult*r2.x;
				f.y += mult*r2.y;
				f.z += mult*r2.z;
			}
		}*/
		c_gsystem.d_forces[d_i] = f;
	}
}

inline void computeCDieSwitch(){
	nonBondedPotentialCDieSwitch_kernel<<<nonBondedBlockCount, nonBondedBlockSize>>>();
}

__global__ void nonBondedPotentialEnergyCDieSwitch_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	float potLJ = 0.0f;
	float potCoulomb = 0.0f;
	if(d_i < c_gsystem.Ntot){
		float4 r1 = tex1Dfetch(t_coord, d_i);
		float4 r2;
		float mult;
		int at1, at2;
		at1 = (int)r1.w;
		float2 par1 = c_ljparameters14[at1];
		float2 par2;
		int i;
		float r;
		for(i = 0; i < c_pairsListsData.d_pairs14ListCount[d_i]; i++){
			r2 =  tex1Dfetch(t_coord, c_pairsListsData.d_pairs14List[i*c_gsystem.widthTot + d_i]);
			at2 = (int)r2.w;
			r2 -= r1;
			DO_PBC(r2);
			r = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			mult = 0.0f;
			float m1, m2, m3;
			if(r < c_nonBondedData.ljCutoff2){
				if(r > c_nonBondedData.ljSwitch2){
					m2 = c_nonBondedData.ljCutoff2 - r;
					m3 = r - c_nonBondedData.ljSwitch2;
					m1 = m2 + 3.0f*m3;
					m2 *= m1;
					m1 = (c_nonBondedData.ljCutoff2 - r)*c_nonBondedData.oneOverCutoff2MinSwitch2Cub;
				} else {
					m1 = 1.0f;
					m2 = 1.0f;
					m3 = 0.0f;
				}
				par2 = c_ljparameters14[at2];
				if(par1.x != 0.0f && par2.x != 0.0f){
					r2.w = par1.y + par2.y;
					r2.w = r2.w*r2.w;
					mult = r2.w/r;
					r1.w = mult*mult*mult;
					mult = (r1.w*r1.w - 2.0f*r1.w);

					mult *= par1.x*par2.x*m1*m2/12.0f;
					//mult *= par1.x*par2.x/12.0f;
				} else {
					mult = 0.0f;
				}
				potLJ += mult;

				if(c_charges14[at1] != 0.0f && c_charges14[at2] != 0.0f){
					r1.w = 1.0f/r;
					potCoulomb += COULOMB_CONSTANT*c_charges14[at1]*c_charges14[at2]*sqrtf(r1.w)*m1*m2;
					//potCoulomb += COULOMB_CONSTANT*c_charges14[at1]*c_charges14[at2]*sqrtf(r1.w);
				}

			}

		}
		par1 = c_ljparameters[at1];
		for(i = 0; i < c_pairsListsData.d_pairsLJListCount[d_i]; i++){
			r2 =  tex1Dfetch(t_coord, c_pairsListsData.d_pairsLJList[i*c_gsystem.widthTot + d_i]);
			at2 = (int)r2.w;
			par2 = c_ljparameters[at2];
			r2 -= r1;
			DO_PBC(r2);
			r = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			mult = 0.0f;
			float m1, m2, m3;
			if(r < c_nonBondedData.ljCutoff2){
				if(r > c_nonBondedData.ljSwitch2){
					m2 = c_nonBondedData.ljCutoff2 - r;
					m3 = r - c_nonBondedData.ljSwitch2;
					m1 = m2 + 3.0f*m3;
					m2 *= m1;
					m1 = (c_nonBondedData.ljCutoff2 - r)*c_nonBondedData.oneOverCutoff2MinSwitch2Cub;
				} else {
					m1 = 1.0f;
					m2 = 1.0f;
					m3 = 0.0f;
				}
				par2 = c_ljparameters[at2];
				if(par1.x != 0.0f && par2.x != 0.0f){
					r2.w = par1.y + par2.y;
					r2.w = r2.w*r2.w;
					mult = r2.w/r;
					r1.w = mult*mult*mult;
					mult = (r1.w*r1.w - 2.0f*r1.w);

					mult *= par1.x*par2.x*m1*m2/12.0f;
					//mult *= par1.x*par2.x/12.0f;
				} else {
					mult = 0.0f;
				}
				potLJ += mult;

				if(c_charges[at1] != 0.0f && c_charges[at2] != 0.0f){
					r1.w = 1.0f/r;
					potCoulomb += COULOMB_CONSTANT*c_charges[at1]*c_charges[at2]*sqrtf(r1.w)*m1*m2;
					//potCoulomb += COULOMB_CONSTANT*c_charges[at1]*c_charges[at2]*sqrtf(r1.w);
				}

			}
		}
		c_nonBondedData.d_ljEnergies[d_i] = potLJ;
		c_nonBondedData.d_coulombEnergies[d_i] = potCoulomb;
	}
}

inline void computeLJPotentialEnergyCDieSwitch(){
	nonBondedPotentialEnergyCDieSwitch_kernel<<<nonBondedBlockCount, nonBondedBlockSize>>>();
	cudaMemcpy(nonBondedData.h_ljEnergies, nonBondedData.d_ljEnergies,
				gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(nonBondedData.h_coulombEnergies, nonBondedData.d_coulombEnergies,
				gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		float pot = 0.0f;
		for(i = 0; i < gsystem.N; i++){
			/*if(nonBondedData.h_ljEnergies[i] > 0.0f){
				printf("Atom %d has another one close to it (LJ Energy = %f)\n", i, nonBondedData.h_ljEnergies[i]);
			}*/
			pot += nonBondedData.h_ljEnergies[i + traj*gsystem.N];
		}
		pot /= 2.0f;
		ljEnergyOutput.values[traj] = pot;
	}
}

inline void computeCoulombPotentialEnergyCDieSwitch(){
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		float pot = 0.0f;
		for(i = 0; i < gsystem.N; i++){
			pot += nonBondedData.h_coulombEnergies[i + traj*gsystem.N];
		}
		pot /= 2.0f;
		coulombEnergyOutput.values[traj] = pot;
	}
}

void destroy(){

}

float2 computeNBBetweenAtoms(int i, int j){
	float2 pot;
	float4 r1, r2;
	r1 =  gsystem.h_coord[i];
	r2 =  gsystem.h_coord[j];
	int at1 = (int)r1.w;
	int at2 = (int)r2.w;
	r2 -= r1;
	DO_PBC(r2);
	float r = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
	if(r < nonBondedData.ljCutoff2){
		int is14 = 0;
		int isExcluded = 0;
		int k;
		for(k = 0; k < pairsListsData.h_pairs14ListCount[i]; k++){
			if(pairsListsData.h_pairs14List[k*gsystem.widthTot + i] == j){
				is14 = 1;
			}
		}
		for(k = 0; k < pairsListsData.h_pairsExclusionListCount[i]; k++){
			if(pairsListsData.h_pairsExclusionList[k*gsystem.widthTot + i] == j){
				isExcluded = 1;
			}
		}
		if(is14 || !isExcluded){
			float2 par1, par2;
			float q1, q2;
			if(is14){
				par1 = nonBondedData.h_ljparameters14[at1];
				par2 = nonBondedData.h_ljparameters14[at2];
				q1 = nonBondedData.h_charges14[at1];
				q2 = nonBondedData.h_charges14[at2];
			} else {
				par1 = nonBondedData.h_ljparameters[at1];
				par2 = nonBondedData.h_ljparameters[at2];
				q1 = nonBondedData.h_charges[at1];
				q2 = nonBondedData.h_charges[at2]	;
			}

			r2.w = par1.y + par2.y;
			r2.w = r2.w*r2.w;
			r1.w = r2.w/r;
			r1.w = r1.w*r1.w*r1.w;

			pot.x = (r1.w*r1.w - 2.0f*r1.w);

			/*r1.w = r2.w/nonBondedData.ljCutoff2;
			r1.w = r1.w*r1.w*r1.w;

			pot.x += (r1.w*r1.w - r1.w)*2.0f*nonBondedData.oneOverLjCutoff6*r*r*r;

			pot.x += (4.0f*r1.w - 3.0f*r1.w*r1.w);*/

			pot.x *= par1.x*par2.x/12.0f;

			r1.w = COULOMB_CONSTANT_RDIE*q1*q2/r;
			r1.w *= (1.0f - r/nonBondedData.coulombCutoff2);
			r1.w *= (1.0f - r/nonBondedData.coulombCutoff2);
			pot.y = r1.w;

			return pot;
		} else {
			return make_float2(0.0f, 0.0f);
		}
	} else {
		return make_float2(0.0f, 0.0f);
	}

}

#undef LOG

} // namespace non_bonded_potential
