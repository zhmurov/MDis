/*
 * GenBornPotential.cuh
 *
 *  Created on: Jan 10, 2011
 *      Author: zhmurov
 */

#pragma once

namespace genborn_potential {

int genBornBlockSize;
int genBornBlockCount;

typedef struct {
	int maxPairs;
	float epsilon;
	float halfCoulomb;
	float epsilonTimesHalfCoulomb;
	float P1;
	float P2;
	float P3;
	float P4;
	float P5;
	float Lambda;
	float Phi;
	float* h_alpha;
	float* d_alpha;
	float* h_dGda;
	float* d_dGda;
	float4* h_dGdr;
	float4* d_dGdr;
	float* h_RVdW;
	float* d_RVdW;
	float* h_VVdW;
	float* d_VVdW;
	float* h_LP1;
	float* d_LP1;

	float roff, ron;
	float roff2, ron2;
	float C1sw, C2sw, C3sw;
	float C1dsw, C2dsw;

	float cutAlpha;

	float* h_genBornEnergies;
	float* d_genBornEnergies;

} GenBornData;

GenBornData genBornData;
__device__ __constant__ GenBornData c_genBornData;

texture<float, 1, cudaReadModeElementType> t_alpha;

__device__ __constant__ float c_RVdW[MAX_ATOM_TYPES];
texture<float, 1, cudaReadModeElementType> t_RVdW;

__device__ __constant__ float c_VVdW[MAX_ATOM_TYPES];
texture<float, 1, cudaReadModeElementType> t_VVdW;

__device__ __constant__ float c_LP1[MAX_ATOM_TYPES];
texture<float, 1, cudaReadModeElementType> t_LP1;

Potential genBornPotential;
EnergyOutput genBornEnergyOutput;

void create();
void init();
inline void compute();
inline void computeEnergy();
void destroy();

} // namespace genborn_potential
