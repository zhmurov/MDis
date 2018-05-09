/*
 * FragmemPotential.cuh
 *
 *  Created on: Mar 28, 2018
 *      Author: zhmurov
 */

#pragma once

namespace fragmem_potential {
	
int fragmemBlockSize;
int fragmemBlockCount;

typedef struct __align__(16) {
	int j;
	float r0;
	float oneOverSigma2;
	float weight;
} GFragmemPair;

typedef struct {
	int maxFragmemsPerAtom;
	int* h_fragmemCount;
	int* d_fragmemCount;
	GFragmemPair* h_fragmems;
	GFragmemPair* d_fragmems;
	float* d_fragmemEnergy;
	float* h_fragmemEnergy;
	float strength;
} GFragmemData;

GFragmemData fragmemData;
__device__ __constant__ GFragmemData c_fragmemData;

Potential potential;
EnergyOutput energyOutput;

void create();
void init();
inline void compute();
inline void computeEnergy();
void destroy();

} // namespace fragmem_potential
