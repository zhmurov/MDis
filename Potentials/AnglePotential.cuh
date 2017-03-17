/*
 * AnglePotential.cuh
 *
 *  Created on: Aug 4, 2010
 *      Author: zhmurov
 */

#pragma once

namespace angle_potential {

int angleBlockSize;
int angleBlockCount;
int angleSummBlockSize;
int angleSummBlockCount;

typedef struct {
	int A;
	int Atot;
	int maxAnglesPerAtom;
	int* h_angleCount;
	int* d_angleCount;
	int4* h_angles;
	int4* d_angles;
	int4* h_angleRefs;
	int4* d_angleRefs;
	float4* h_angleForces;
	float4* d_angleForces;
	float2* h_angleTypes;
	float2* d_angleTypes;
	float* h_angleEnergies;
	float* d_angleEnergies;
} GAngleData;

GAngleData angleData;
__device__ __constant__ GAngleData c_angleData;
texture<float2, 1, cudaReadModeElementType> t_angleTypes;

Potential potential;
EnergyOutput energyOutput;

void create();
void init();
inline void compute();
inline void computeEnergy();
void destroy();


} // namespace angle_potential