/*
 * RepulsiveBoundaryPotential.cuh
 *
 *  Created on: Feb 28, 2011
 *      Author: zhmurov
 */

#pragma once

namespace repulsive_boundary {

int repulsiveBoundaryBlockSize;
int repulsiveBoundaryBlockCount;

typedef struct {
	float4 geometry;
	float epsilon;
	float sigma;
	float sixEpsilonSigma6;
} RepulsiveBoundary;


Potential potential;

void create();
void init();
inline void compute();
void destroy();

} // namespace repulsive_boundary

repulsive_boundary::RepulsiveBoundary repulsiveBoundaryData;
__device__ __constant__ repulsive_boundary::RepulsiveBoundary c_repulsiveBoundaryData;
