/*
 * PushingPlanePotential.cuh
 *
 *  Created on: Apr 16, 2016
 *      Author: kir_min
 */

#define FIX_PLANE 0
#define PUSHING_PLANE 1

namespace pushing_plane_potential
{

int blockSize;
int blockCount;

typedef struct
{
	int planeCount;
	float4* h_planeNormal;
	float4* d_planeNormal;
	float4* h_planePosition0;
	float4* d_planePosition0;
	float4* h_planePosition;
	float4* d_planePosition;
	float sigma;
	float epsilon;
	float4* h_forces;
	float4* d_forces;
	float pushSpeed;
	float springConstant;
	float dampingConstant;
	long long int updateFreq;
	long long int outputFreq;
	float* r;

} PotentialData;

PotentialData potentialData;
__device__ __constant__ PotentialData c_potentialData;

Potential potential;
Updater planeLocationUpdater;
PDB refpdb;

//char outputfilename[1024];
char** outputfilename;
char reffilename[1024];

void create();
void init();
void destroy();

void updatePlaneLocation();
void destroyPlaneLocationUpdater();

inline void computePlanePushing();

} //namespace pushing_plane_potential
