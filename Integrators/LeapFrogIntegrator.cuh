/*
 * LeapFrogIntegrator.cuh
 *
 *  Created on: Jul 29, 2010
 *      Author: zhmurov
 */

#pragma once

#include "../Core/md.cuh"

namespace leapfrog_integrator {

int blockSize;
int blockCount;

typedef struct {
	float h; // timestep, ps
	float gamma; // damping coefficient, 1/ps
} LeapFrog;

LeapFrog leapFrog;
__device__ __constant__ LeapFrog c_leapFrog;

Integrator leapFrogIntegrator;

void create();
void init();
void inline computeLeapFrogIntegrator();
void inline finalizeLeapFrogIntegrator();



} // namespace leapfrog_integrator
