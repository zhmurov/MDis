/*
 * ffreader.h
 *
 *  Created on: Jul 4, 2009
 *      Author: zhmurov
 */

#pragma once

#include "../Core/forcefield.h"
#include "../Core/topology.h"

#define FF_SECTION_BONDS		0
#define FF_SECTION_ANGLES		1
#define FF_SECTION_DIHEDRALS	2
#define FF_SECTION_IMPROPERS	3
#define FF_SECTION_CMAP			4
#define FF_SECTION_NONBONDED	5

/*
 * Public methods
 */
void readCHARMMFF(const char* filename, ForceField* ff, char* fftype);
void convertAtomTypesCHARMM19(const char* topFilename, Atom* atoms, int atomCount);
