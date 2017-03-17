/*
 * CoordinatesOutputManagerDCD.cuh
 *
 *  Created on: Aug 2, 2010
 *      Author: zhmurov
 */

#pragma once

namespace coordinates_output_dcd {

typedef struct {
	float* X;
	float* Y;
	float* Z;
	FILE* file;
} DCDOutput;

Updater updater;
DCDOutput dcdOutput;

void create();
void init();
void save();
void saveTrajectory(int traj);
void destroy();

} // namespace coordinates_output_dcd