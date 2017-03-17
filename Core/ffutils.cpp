/*
 * ffutils.c
 *
 *  Created on: Jul 5, 2009
 *      Author: zhmurov
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "../Core/forcefield.h"
#include "../Util/wrapper.h"

FFBondType findBondType(char atomType1[5], char atomType2[5], ForceField* forceField){
	int i;
	for(i = 0; i < forceField->bondTypesCount; i++){
		if((strncmp(atomType1, forceField->bondTypes[i].atomType1.name, 5) == 0) &&
				(strncmp(atomType2, forceField->bondTypes[i].atomType2.name, 5) == 0)){
			return forceField->bondTypes[i];
		}
		if((strncmp(atomType2, forceField->bondTypes[i].atomType1.name, 5) == 0) &&
				(strncmp(atomType1, forceField->bondTypes[i].atomType2.name, 5) == 0)){
			return forceField->bondTypes[i];
		}
	}
	for(i = 0; i < forceField->bondTypesCount; i++){
		if((strncmp("X", forceField->bondTypes[i].atomType1.name, 1) == 0) &&
				(strncmp(atomType2, forceField->bondTypes[i].atomType2.name, 5) == 0)){
			return forceField->bondTypes[i];
		}
		if((strncmp("X", forceField->bondTypes[i].atomType1.name, 1) == 0) &&
				(strncmp(atomType1, forceField->bondTypes[i].atomType2.name, 5) == 0)){
			return forceField->bondTypes[i];
		}
		if((strncmp(atomType1, forceField->bondTypes[i].atomType1.name, 5) == 0) &&
				(strncmp("X", forceField->bondTypes[i].atomType2.name, 1) == 0)){
			return forceField->bondTypes[i];
		}
		if((strncmp(atomType2, forceField->bondTypes[i].atomType1.name, 5) == 0) &&
				(strncmp("X", forceField->bondTypes[i].atomType2.name, 1) == 0)){
			return forceField->bondTypes[i];
		}
	}
	DIE("ERROR: Can't find bond type for atoms '%s' and '%s'.", atomType1, atomType2);
}

FFAngleType findAngleType(char atomType1[5], char atomType2[5], char atomType3[5], ForceField* forceField){
	int i;
	for(i = 0; i < forceField->angleTypesCount; i++){
		if((strncmp(atomType1, forceField->angleTypes[i].atomType1.name, 5) == 0 ||
				strncmp("X", forceField->angleTypes[i].atomType1.name, 1) == 0) &&
			(strncmp(atomType2, forceField->angleTypes[i].atomType2.name, 5) == 0 ||
				strncmp("X", forceField->angleTypes[i].atomType2.name, 1) == 0) &&
			(strncmp(atomType3, forceField->angleTypes[i].atomType3.name, 5) == 0 ||
				strncmp("X", forceField->angleTypes[i].atomType3.name, 1) == 0) ){
			return forceField->angleTypes[i];
		}
		if((strncmp(atomType3, forceField->angleTypes[i].atomType1.name, 5) == 0 ||
				strncmp("X", forceField->angleTypes[i].atomType1.name, 1) == 0) &&
			(strncmp(atomType2, forceField->angleTypes[i].atomType2.name, 5) == 0 ||
				strncmp("X", forceField->angleTypes[i].atomType2.name, 1) == 0) &&
			(strncmp(atomType1, forceField->angleTypes[i].atomType3.name, 5) == 0 ||
				strncmp("X", forceField->angleTypes[i].atomType3.name, 1) == 0) ){
			return forceField->angleTypes[i];
		}
	}
	DIE("ERROR: Can't find angle type for atoms '%s', '%s' and '%s'.", atomType1, atomType2, atomType3);
}

int findDihedralType(char atomType1[5], char atomType2[5], char atomType3[5], char atomType4[5],
		ForceField* forceField, FFDihedralType* dt){
	int i;
	int multiplicity = 0;
	int found = 0;
	for(i = 0; i < forceField->dihedralTypesCount; i++){
		while(i < forceField->dihedralTypesCount &&
				strncmp(atomType1, forceField->dihedralTypes[i].atomType1.name, 5) == 0 &&
				strncmp(atomType2, forceField->dihedralTypes[i].atomType2.name, 5) == 0 &&
				strncmp(atomType3, forceField->dihedralTypes[i].atomType3.name, 5) == 0 &&
				strncmp(atomType4, forceField->dihedralTypes[i].atomType4.name, 5) == 0){
			/*printf("Using dihedral between'%s', '%s', '%s' and '%s'. Curent mult: %d.\n",
					forceField->dihedralTypes[i].atomType1.name,
					forceField->dihedralTypes[i].atomType2.name,
					forceField->dihedralTypes[i].atomType3.name,
					forceField->dihedralTypes[i].atomType4.name,
					multiplicity);*/

			dt[multiplicity] = forceField->dihedralTypes[i];
			multiplicity ++;
			i++;
			found = 1;
		}
		if(found == 1){
			return multiplicity;
		}
	}
	for(i = 0; i < forceField->dihedralTypesCount; i++){
		while(i < forceField->dihedralTypesCount &&
				strncmp(atomType4, forceField->dihedralTypes[i].atomType1.name, 5) == 0 &&
				strncmp(atomType3, forceField->dihedralTypes[i].atomType2.name, 5) == 0 &&
				strncmp(atomType2, forceField->dihedralTypes[i].atomType3.name, 5) == 0 &&
				strncmp(atomType1, forceField->dihedralTypes[i].atomType4.name, 5) == 0){
			/*printf("Using dihedral between'%s', '%s', '%s' and '%s'. Curent mult: %d.\n",
					forceField->dihedralTypes[i].atomType1.name,
					forceField->dihedralTypes[i].atomType2.name,
					forceField->dihedralTypes[i].atomType3.name,
					forceField->dihedralTypes[i].atomType4.name,
					multiplicity);*/
			dt[multiplicity] = forceField->dihedralTypes[i];
			multiplicity ++;
			i++;
			found = 1;
		}
		if(found == 1){
			return multiplicity;
		}
	}
	if(multiplicity > 0){
		return multiplicity;
	}
	//printf("WARNING: Can't find parameters for dihedral '%s', '%s', '%s' and '%s'. "
	//		"Using parameters with undefined atoms ('X').\n", atomType1, atomType2, atomType3, atomType4);
	for(i = 0; i < forceField->dihedralTypesCount; i++){
		while(i < forceField->dihedralTypesCount &&
				(strncmp(atomType1, forceField->dihedralTypes[i].atomType1.name, 5) == 0 ||
				strncmp("X", forceField->dihedralTypes[i].atomType1.name, 1) == 0) &&
			(strncmp(atomType2, forceField->dihedralTypes[i].atomType2.name, 5) == 0) &&
			(strncmp(atomType3, forceField->dihedralTypes[i].atomType3.name, 5) == 0) &&
			(strncmp(atomType4, forceField->dihedralTypes[i].atomType4.name, 5) == 0 ||
				strncmp("X", forceField->dihedralTypes[i].atomType4.name, 1) == 0)){
			/*printf("Using dihedral between'%s', '%s', '%s' and '%s'. Curent mult: %d.\n",
					forceField->dihedralTypes[i].atomType1.name,
					forceField->dihedralTypes[i].atomType2.name,
					forceField->dihedralTypes[i].atomType3.name,
					forceField->dihedralTypes[i].atomType4.name,
					multiplicity);*/
			dt[multiplicity] = forceField->dihedralTypes[i];
			multiplicity ++;
			i++;
			found = 1;
		}
		if(found == 1){
			return multiplicity;
		}
	}
	for(i = 0; i < forceField->dihedralTypesCount; i++){
		while(i < forceField->dihedralTypesCount &&
				(strncmp(atomType4, forceField->dihedralTypes[i].atomType1.name, 5) == 0 ||
				strncmp("X", forceField->dihedralTypes[i].atomType1.name, 1) == 0) &&
			(strncmp(atomType3, forceField->dihedralTypes[i].atomType2.name, 5) == 0) &&
			(strncmp(atomType2, forceField->dihedralTypes[i].atomType3.name, 5) == 0) &&
			(strncmp(atomType1, forceField->dihedralTypes[i].atomType4.name, 5) == 0 ||
				strncmp("X", forceField->dihedralTypes[i].atomType4.name, 1) == 0)){
			/*printf("Using dihedral between'%s', '%s', '%s' and '%s'. Curent mult: %d.\n",
					forceField->dihedralTypes[i].atomType1.name,
					forceField->dihedralTypes[i].atomType2.name,
					forceField->dihedralTypes[i].atomType3.name,
					forceField->dihedralTypes[i].atomType4.name,
					multiplicity);*/
			dt[multiplicity] = forceField->dihedralTypes[i];
			multiplicity ++;
			i++;
			found = 1;
		}
		if(found == 1){
			return multiplicity;
		}
	}
	if(multiplicity > 0){
		return multiplicity;
	}
	DIE("ERROR: Can't find parameters for dihedral '%s', '%s', '%s' and '%s'.\n",
			atomType1, atomType2, atomType3, atomType4);

}

int findImproperType(char atomType1[5], char atomType2[5], char atomType3[5], char atomType4[5],
		ForceField* forceField, FFImproperType* it){
	int i;
	int multiplicity = 0;
	int found = 0;
	for(i = 0; i < forceField->improperTypesCount; i++){
		found = 0;
		if(strncmp(atomType1, forceField->improperTypes[i].atomType1.name, 5) == 0 &&
				strncmp(atomType2, forceField->improperTypes[i].atomType2.name, 5) == 0  &&
				strncmp(atomType3, forceField->improperTypes[i].atomType3.name, 5) == 0  &&
				strncmp(atomType4, forceField->improperTypes[i].atomType4.name, 5) == 0){
			it[multiplicity] = forceField->improperTypes[i];
			multiplicity ++;
			found = 1;
		}
		if(found == 0  &&
				strncmp(atomType4, forceField->improperTypes[i].atomType1.name, 5) == 0 &&
				strncmp(atomType3, forceField->improperTypes[i].atomType2.name, 5) == 0 &&
				strncmp(atomType2, forceField->improperTypes[i].atomType3.name, 5) == 0 &&
				strncmp(atomType1, forceField->improperTypes[i].atomType4.name, 5) == 0 ){
			it[multiplicity] = forceField->improperTypes[i];
			multiplicity ++;
		}
	}
	if(multiplicity > 0){
		return multiplicity;
	}

	//printf("WARNING: Can't find parameters for improper '%s', '%s', '%s' and '%s'. "
	//			"Using parameters with undefined atoms ('X').\n", atomType1, atomType2, atomType3, atomType4);


	for(i = 0; i < forceField->improperTypesCount; i++){
		found = 0;
		if((strncmp(atomType1, forceField->improperTypes[i].atomType1.name, 5) == 0 ||
				strncmp("X", forceField->improperTypes[i].atomType1.name, 1) == 0) &&
			(strncmp(atomType2, forceField->improperTypes[i].atomType2.name, 5) == 0 ||
				strncmp("X", forceField->improperTypes[i].atomType2.name, 1) == 0) &&
			(strncmp(atomType3, forceField->improperTypes[i].atomType3.name, 5) == 0 ||
				strncmp("X", forceField->improperTypes[i].atomType3.name, 1) == 0) &&
			(strncmp(atomType4, forceField->improperTypes[i].atomType4.name, 5) == 0 ||
				strncmp("X", forceField->improperTypes[i].atomType4.name, 1) == 0)){
			it[multiplicity] = forceField->improperTypes[i];
			multiplicity ++;
			return multiplicity;
		}
		if((strncmp(atomType4, forceField->improperTypes[i].atomType1.name, 5) == 0 ||
				strncmp("X", forceField->improperTypes[i].atomType1.name, 1) == 0) &&
			(strncmp(atomType3, forceField->improperTypes[i].atomType2.name, 5) == 0 ||
				strncmp("X", forceField->improperTypes[i].atomType2.name, 1) == 0) &&
			(strncmp(atomType2, forceField->improperTypes[i].atomType3.name, 5) == 0 ||
				strncmp("X", forceField->improperTypes[i].atomType3.name, 1) == 0) &&
			(strncmp(atomType1, forceField->improperTypes[i].atomType4.name, 5) == 0 ||
				strncmp("X", forceField->improperTypes[i].atomType4.name, 1) == 0)){
			it[multiplicity] = forceField->improperTypes[i];
			multiplicity ++;
			return multiplicity;
		}
	}
	DIE("ERROR: Can't find bond type for improper '%s', '%s', '%s' and '%s'.\n",
			atomType1, atomType2, atomType3, atomType4);
}

FFCMAPType findCMAPType(
		char atomType1[5], char atomType2[5], char atomType3[5], char atomType4[5],
		char atomType5[5], char atomType6[5], char atomType7[5], char atomType8[5],
		ForceField* forceField){
	int i;
	for(i = 0; i < forceField->cmapTypesCount; i++){
		if(strcmp(atomType1, forceField->cmapTypes[i].atomType1.name) == 0 &&
				strcmp(atomType2, forceField->cmapTypes[i].atomType2.name) == 0 &&
				strcmp(atomType3, forceField->cmapTypes[i].atomType3.name) == 0 &&
				strcmp(atomType4, forceField->cmapTypes[i].atomType4.name) == 0 &&
				strcmp(atomType5, forceField->cmapTypes[i].atomType5.name) == 0 &&
				strcmp(atomType6, forceField->cmapTypes[i].atomType6.name) == 0 &&
				strcmp(atomType7, forceField->cmapTypes[i].atomType7.name) == 0 &&
				strcmp(atomType8, forceField->cmapTypes[i].atomType8.name) == 0){
			return forceField->cmapTypes[i];
		}

	}
}

FFNonbondedType findNonbondedType(char atomType[5], ForceField* forceField){
	int i;
	//printf("nonbonded types count: %d\n", forceField->nonbondedTypesCount);
	for(i = 0; i < forceField->nonbondedTypesCount; i++){
		//printf("%s\n", forceField->nonbondedTypes[i].atomType.name);
		if(strcmp(atomType, forceField->nonbondedTypes[i].atomType.name) == 0){
			//printf("Found type for %s\n", atomType);
			return forceField->nonbondedTypes[i];
		}
	}
	for(i = 0; i < forceField->nonbondedTypesCount; i++){
		if(atomType[0] == forceField->nonbondedTypes[i].atomType.name[0] &&
				((forceField->nonbondedTypes[i].atomType.name[1] == '*') ||
				(forceField->nonbondedTypes[i].atomType.name[1] == '%'))){
			//printf("Non-bonded type: %s = %s \n", atomType, forceField->nonbondedTypes[i].atomType.name);
			return forceField->nonbondedTypes[i];
		}
	}
	DIE("Can't find non-bonded type for %s\n", atomType);
}
