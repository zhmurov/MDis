/*
 * ImproperPotential.cuh
 *
 *  Created on: Aug 5, 2010
 *      Author: zhmurov
 */

#pragma once

namespace improper_potential {

int improperBlockSize;
int improperBlockCount;
int improperSummBlockSize;
int improperSummBlockCount;

typedef struct __align__(16) {
	float kpsi;
	float psi0;
	float n;
	int multiplicityRef;
} GImproperParameters;

typedef struct {
	int I;
	int Itot;
	int maxImpropersPerAtom;
	int* h_improperCount;
	int* d_improperCount;
	int4* h_impropers;
	int4* d_impropers;
	int4* h_improperRefs;
	int4* d_improperRefs;
	GImproperParameters* h_improperParameters;
	GImproperParameters* d_improperParameters;
	GImproperParameters* h_multiplicityParameters;
	GImproperParameters* d_multiplicityParameters;
	float4* h_improperForces;
	float4* d_improperForces;
	float* h_improperEnergies;
	float* d_improperEnergies;
} GImproperData;

GImproperData improperData;
__device__ __constant__ GImproperData c_improperData;

Potential potential;
EnergyOutput energyOutput;

void create();
void init();
inline void compute();
inline void computeEnergy();
void destroy();

} // namespace improper_potential