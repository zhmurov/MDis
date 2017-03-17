/*
 * EnergyOutputManager.cuh
 *
 *  Created on: Aug 2, 2010
 *      Author: zhmurov
 */

#pragma once

namespace energy_output {

#define KCALL_PER_KJ	0.239006f


FILE* energyOutputFile;
typedef struct {
	float* T; // Temperature for a trajectory [K]
	float* E; // Sum of all energies for a trajectory [kJ]
} EnergyOutputData;

Updater updater;

void create();
void init();
void save();
void saveRigidBody(int traj);
void destroy();

} // namespace energy_output

extern energy_output::EnergyOutputData energyOutputData;
extern int constrCount;
