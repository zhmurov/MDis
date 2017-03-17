/*
 * angleTypes.c
 *
 *  Created on: Jul 29, 2010
 *      Author: zhmurov
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../Util/wrapper.h"
#include "topology.h"
#include "angleTypes.h"

AngleType* angleTypes;
int angleTypesCount;

int compareAngles(Angle a1, Angle a2);
int getAngleTypeId(Angle angle);

void initAngleTypes(){
	printf("Looking for angle types...\n");
	int i, j;
	angleTypesCount = 0;
	int typeExists = 0;
	for(i = 0; i < topology.angleCount; i++){
		//for(j = 0; j < i; j++){
		j = 0;
		while(!typeExists && j < i){
			if(compareAngles(topology.angles[i], topology.angles[j])){
				typeExists = 1;
			}
			j++;
		}
		if(!typeExists){
			angleTypesCount ++;
		} else {
			typeExists = 0;
		}
	}
	printf("Found %d unique angle types.\n", angleTypesCount);
	angleTypes = (AngleType*)calloc(angleTypesCount, sizeof(AngleType));
	angleTypesCount = 0;
	for(i = 0; i < topology.angleCount; i++){
		//for(j = 0; j < i; j++){
		j = 0;
		while(!typeExists && j < i){
			if(compareAngles(topology.angles[i], topology.angles[j])){
				typeExists = 1;
			}
			j++;
		}
		if(!typeExists){

			strcpy(angleTypes[angleTypesCount].atomType1, topology.atoms[topology.angles[i].i].type);
			strcpy(angleTypes[angleTypesCount].atomType2, topology.atoms[topology.angles[i].j].type);
			strcpy(angleTypes[angleTypesCount].atomType3, topology.atoms[topology.angles[i].k].type);

			angleTypes[angleTypesCount].ktheta = topology.angles[i].ktheta;
			angleTypes[angleTypesCount].theta0 = topology.angles[i].theta0;

			//printf("%d: ", angleTypesCount);
			//printAngleType(angleTypes[angleTypesCount]);

			angleTypesCount ++;
		} else {
			typeExists = 0;
		}
	}

	for(i = 0; i < topology.angleCount; i++){
		topology.angles[i].type = getAngleTypeId(topology.angles[i]);
	}
	printf("Done looking for angle types.\n");
}

int compareAngles(Angle a1, Angle a2){
	if(a1.ktheta == a2.ktheta && a1.theta0 == a2.theta0){
		return 1;
	} else {
		return 0;
	}
}

int getAngleTypeId(Angle angle){
	int j;
	for(j = 0; j < angleTypesCount; j++){
		AngleType type = angleTypes[j];
		if(angle.ktheta == type.ktheta && angle.theta0 == type.theta0){
			return j;
		}
	}
	DIE("Can't find angle type for the following angle: (%s-%s-%s)",
			topology.atoms[angle.i].type,
			topology.atoms[angle.j].type,
			topology.atoms[angle.k].type);
    return -1;
}

void printAngleType(AngleType angleType){
	printf("%6s%6s%6s|\t\tktheta = %5.3f;\ttheta0 = %5.3f.\n",
		angleType.atomType1,
		angleType.atomType2,
		angleType.atomType3,
		angleType.ktheta,
		angleType.theta0);
}
