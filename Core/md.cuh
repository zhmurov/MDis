/*
 * md.cuh
 *
 *  Created on: Jul 28, 2010
 *      Author: zhmurov
 */

#pragma once

#include <cuda.h> // To retrieve CUDA version
#include "gsystem.h"

#define MAX_POTENTIALS_COUNT	10
#define MAX_UPDATERS_COUNT		10
#define MAX_CONSTR_ALGS_COUNT	 3
#define MAX_ENERGY_OUTPUT_COUNT	20
#define MAX_RESTARTERS_COUNT	10

#define BLOCK_SIZE	256

/*typedef struct {
	float x;
	float y;
	float z;
	int type;
} Coord;*/

/*typedef struct __align__ (16) {
	float x;
	float y;
	float z;
	float T;
} Vel;

typedef struct __align__ (16) {
	float x;
	float y;
	float z;
	float m;
} Force;*/

/*
typedef struct {

	int N;
	int width;
	int Nsim;
	int Ntot;
	int widthTot;

	float4* h_coord;
	float4* h_vel;
	float4* h_forces;

	float4* d_coord;
	float4* d_vel;
	float4* d_forces;

	int* h_atomTypes;
	int* d_atomTypes;

	float* h_m;
	float* d_m;

} GSystem;

typedef struct {
	char name[100];
	void (*compute)();
	void (*destroy)();
} Potential;

typedef struct {
	char name[100];
	int frequency;
	void (*update)();
	void (*destroy)();
} Updater;

typedef struct {
	char name[100];
	float h;
	void (*integrate)();
	void (*destroy)();
} Integrator;

typedef struct {
	char name[100];
	float* values;
	void (*computeValues)();
} EnergyOutput;

*/

long long int step;
long long int firststep;

double* trajectoryTime;

GSystem gsystem;
__device__ __constant__ GSystem c_gsystem;
texture<float4, 1, cudaReadModeElementType> t_coord; // Coordinates
//texture<float4, 1, cudaReadModeElementType> t_midcoord;
texture<float, 1, cudaReadModeElementType> t_m;
texture<int, 1, cudaReadModeElementType> t_atomTypes;

__device__ __shared__ float4 s_Ri[BLOCK_SIZE];

cudaDeviceProp deviceProps;

extern Potential** potentials;
extern Updater** updaters;
extern Integrator* integrator;
extern ConstrAlg** constrAlgs;
extern EnergyOutput** energyOutputs;
extern Restarter** restarters;

extern int potentialsCount;
extern int updatersCount;
extern int constrAlgsCount;
extern int energyOutputsCount;
extern int restartersCount;

void copyCoordinatesFromGPU(bool force = false);
void copyCoordinatesToGPU(int traj, int force);
void copyVelocitiesFromGPU();
void copyVelocitiesToGPU(int traj, int force);
void addRestarter(const char* name, void (*save)(FILE*), void (*load)(FILE*));
