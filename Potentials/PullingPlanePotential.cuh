/*
 * PullingPlanePotential.cuh
 *
 *  Created on: Sep 19, 2013
 *      Author: zip
 */

#pragma once

#define PLANE_ATOM_FREE 0
#define PLANE_ATOM_FIXED 1
#define PLANE_ATOM_PULLED 2

namespace pulling_plane_potential
{

int planePullingGridSize;
int planePullingBlockSize;

typedef struct
{
	int atomID;
	int planeID;
	float bDistance;
} WorkList;

typedef struct
{
	int Npulled;
	int Nfixed;
	int Ntotal;
	int Nplane;
	int* prc;
	float pullForce;
	float pullSpeed;
	float pullSpring;
	float planeMass;
	float fixSpring;
	float* h_force;
	float* d_force;
	float* h_planeDisplacement;
	float* d_planeDisplacement;
	float4 planeNorm;
	float4 planePoint;
	WorkList* h_workList;
	WorkList* d_workList;
} PotentialData;

PotentialData potentialData;
__device__ __constant__ PotentialData c_potentialData;

Potential potential;
Updater planeLocationUpdater;

FILE* outputFile;
char** outputFilename;
int outputFreq;
int logOutputFreq;

void create();
void init();
void destroy();

void initWorkList();

void updatePlaneLocation();
void destroyPlaneLocationUpdater();

inline void computePlanePulling();
inline void computeEnergy();

} //namespace pulling_plane_potential
