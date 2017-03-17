/*
 * atomTypes.h
 *
 *  Created on: Jul 29, 2010
 *      Author: zhmurov
 */

#pragma once

#define MAX_ATOM_TYPES	100

typedef struct {
	char name[ATOM_NAME_LENGTH];
	float charge;
	float mass;
	float ignored;
	float epsilon;
	float RminOver2;
	float ignored_14;
	float epsilon_14;
	float RminOver2_14;
} AtomType;

extern AtomType* atomTypes;
extern int atomTypesCount;

extern void initAtomTypes();
extern void printAtomType(AtomType atomType);
