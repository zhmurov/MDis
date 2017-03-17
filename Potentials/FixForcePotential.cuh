/*
 * FixForcePotential.cuh
 *
 *  Created on: Nov 14, 2010
 *      Author: zhmurov
 */

#pragma once

#define FIXFORCE_ATOM_FREE		0
#define FIXFORCE_ATOM_FIXED		1
#define FIXFORCE_ATOM_PULLED	2

namespace fixforce_potential {

int fixForceBlockSize;
int fixForceBlockCount;

char** fixForceOutputFilename;
FILE* fixForceOutputFile;

typedef struct {
	int* h_atomMasks;
	int* d_atomMasks;
	float4* h_extForces;
	float4* d_extForces;
} FixForceData;

typedef struct {
	float4* h_initialTipPositions;
	float4* d_initialTipPositions;
	float4* h_tipPositions;
	float4* d_tipPositions;
	float* h_ks;
	float* d_ks;
	float4* h_pullVelocities;
	float4* d_pullVelocities;
} SMDData;

int fixForceTimesMoved;
__device__ __constant__ int c_fixForceTimesMoved;

FixForceData fixForceData;
__device__ __constant__ FixForceData c_fixForceData;

SMDData smdData;
__device__ __constant__ SMDData c_smdData;

Potential fixForcePotential;
Updater smdForceUpdater;

void create();
void init();
void initFConst();
void initSMD();
inline void computeFConst();
inline void computeSMD();
void destroy();

void updatePullingForceSMD();
void destroySMDForceUpdater();

} // namespace fixforce_potential
