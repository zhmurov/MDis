/*
 * PushingPlanePotential.cuh
 *
 *  Created on: Apr 20, 2016
 *      Author: kir_min
 */

namespace drum_potential
{

int blockSize;
int blockCount;

typedef struct
{
	float4 Norm;
	float4 Point;
	float sigma;
	float epsilon;
	float R;

} PotentialData;

PotentialData potentialData;
__device__ __constant__ PotentialData c_potentialData;

Potential potential;
char outputfilename[1024];

void create();
void init();
void destroy();

inline void computeDrumPotential();

} //namespace drum_potential
