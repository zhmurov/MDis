/*
 * RestartOutputManager.cuh
 *
 *  Created on: Nov 13, 2010
 *      Author: zhmurov
 */

#pragma once

namespace restart_output {

Updater updater;
PDB pdbOutputData;
XYZ xyzOutputData;

void create();
void init();
void savePDBs();
void saveXYZs();
void destroy();

} // namespace restart_output