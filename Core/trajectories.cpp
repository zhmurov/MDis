/*
 * trajectories.cpp
 *
 *  Created on: Apr 18, 2011
 *      Author: zhmurov
 */
#include "trajectories.h"
#include "parameters.h"
#include "gsystem.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Trajectory* trajectories;
long int lastStepComputed;

void initTrajectories(){
	trajectories = (Trajectory*)calloc(parameters.Ntr, sizeof(Trajectory));
	int i, traj;
	for(i = 0; i < parameters.Ntr; i++){
		traj = parameters.firstrun + i;
		trajectories[i].run = traj;
		sprintf(trajectories[i].runString, "%d\0", traj);
	}
	lastStepComputed = -1;
	trajMap = (int*)calloc(parameters.Ntr, sizeof(int));
	for(i = 0; i < parameters.Ntr; i++){
		trajMap[i] = i;
	}
}

void computeEnergies(){
	if(lastStepComputed != step){
		int i;
		for(i = 0; i < energyOutputsCount; i++){
			energyOutputs[i]->computeValues();
		}
		lastStepComputed = step;
	}
}

void swapTrajectories(int t1, int t2){

}
