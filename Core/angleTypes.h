/*
 * angleTypes.h
 *
 *  Created on: Jul 29, 2010
 *      Author: zhmurov
 */

#pragma once

#define ATOM_NAME_LENGTH	10

typedef struct {
	char atomType1[ATOM_NAME_LENGTH];
	char atomType2[ATOM_NAME_LENGTH];
	char atomType3[ATOM_NAME_LENGTH];
	float ktheta;
	float theta0;
} AngleType;

extern AngleType* angleTypes;
extern int angleTypesCount;

extern void initAngleTypes();
extern void printAngleType(AngleType angleType);
