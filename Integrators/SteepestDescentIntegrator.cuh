/*
 * SteepestDescentIntegrator.cuh
 *
 *  Created on: Nov 13, 2010
 *      Author: zhmurov
 */

#pragma once

namespace sd_integrator {

int blockSize;
int blockCount;

typedef struct {
	float h;
	float maxForce;
	float timestepChangeVel;
	float timestepFactor;
	float finalTimestep;
} SteepestDescent;

SteepestDescent steepestDescent;
__device__ __constant__ SteepestDescent c_steepestDescent;

Integrator steepestDescentIntegrator;
Updater steepestDescentTimestepUpdater;

void create();
void init();
void inline computeSteepestDescentIntegrator();
void inline finalizeSteepestDescentIntegrator();

void inline updateTimestep();
void destroyUpdater();

} // namespace sd_integrator
