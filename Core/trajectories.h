/*
 * trajectories.h
 *
 *  Created on: Apr 18, 2011
 *      Author: zhmurov
 */

#ifndef TRAJECTORIES_H_
#define TRAJECTORIES_H_

#include <vector_types.h>

typedef struct {
	int run;
	char runString[64];
	float* E;
	float Etot;
	float T_MB;
} Trajectory;

int* trajMap;

void initTrajectories();

void computeEnergies();

#endif /* TRAJECTORIES_H_ */
