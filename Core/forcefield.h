/*
 * forcefield.h
 *
 *  Created on: Jul 4, 2009
 *      Author: zhmurov
 */

#pragma once

#include "topology.h"

#define FF_TYPE_CHARMM19			"CHARMM19"
#define FF_TYPE_CHARMM22			"CHARMM22"

typedef struct{
	char name[10];
} FFAtomType;

typedef struct{
	FFAtomType atomType1;
	FFAtomType atomType2;
	float kb;
	float b0;
} FFBondType;

typedef struct{
	FFAtomType atomType1;
	FFAtomType atomType2;
	FFAtomType atomType3;
	float ktheta;
	float theta0;
	float kub;
	float s0;
} FFAngleType;

typedef struct{
	FFAtomType atomType1;
	FFAtomType atomType2;
	FFAtomType atomType3;
	FFAtomType atomType4;
	float kchi;
	float delta;
	int n;
} FFDihedralType;

typedef struct{
	FFAtomType atomType1;
	FFAtomType atomType2;
	FFAtomType atomType3;
	FFAtomType atomType4;
	float kpsi;
	float psi0;
	int n;
} FFImproperType;

typedef struct{
	FFAtomType atomType1;
	FFAtomType atomType2;
	FFAtomType atomType3;
	FFAtomType atomType4;
	FFAtomType atomType5;
	FFAtomType atomType6;
	FFAtomType atomType7;
	FFAtomType atomType8;
	int gridSize;
	float* data;
} FFCMAPType;

typedef struct{
	FFAtomType atomType;
	float ignored;
	float epsilon;
	float RminOver2;
	float ignored_14;
	float epsilon_14;
	float RminOver2_14;
} FFNonbondedType;

typedef struct{
	int bondTypesCount;
	int angleTypesCount;
	int dihedralTypesCount;
	int improperTypesCount;
	int cmapTypesCount;
	int nonbondedTypesCount;
	FFBondType* bondTypes;
	FFAngleType* angleTypes;
	FFDihedralType* dihedralTypes;
	FFImproperType* improperTypes;
	FFNonbondedType* nonbondedTypes;
	FFCMAPType* cmapTypes;
} ForceField;

void readCHARMMFF(const char* filename, ForceField* ff, char* fftype);
extern void convertAtomTypesCHARMM19(const char* topFilename, Atom* atoms, int atomCount);
