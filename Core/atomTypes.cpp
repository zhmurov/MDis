/*
 * atomTypes.c
 *
 *  Created on: Jul 29, 2010
 *      Author: zhmurov
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "topology.h"
#include "atomTypes.h"
#include "../Util/wrapper.h"

AtomType* atomTypes;
int atomTypesCount;

int compareAtoms(Atom a1, Atom a2);
int getAtomTypeId(Atom atom);
void testAtomTypes();

void initAtomTypes(){
	printf("Looking for atom types...\n");
	int i, j;
	atomTypesCount = 0;
	int typeExists = 0;
	for(i = 0; i < topology.atomCount; i++){
		j = 0;
		while(j < i && !typeExists){
		//for(j = 0; j < i; j++){
			if(compareAtoms(topology.atoms[i], topology.atoms[j])){
				typeExists = 1;
			}
			j++;
		}
		if(!typeExists){
			atomTypesCount ++;
		} else {
			typeExists = 0;
		}
	}
	printf("Found %d unique atom types:\n", atomTypesCount);
	atomTypes = (AtomType*)calloc(atomTypesCount, sizeof(AtomType));
	atomTypesCount = 0;
	for(i = 0; i < topology.atomCount; i++){
		j = 0;
		while(j < i && !typeExists){
		//for(j = 0; j < i; j++){
			if(compareAtoms(topology.atoms[i], topology.atoms[j])){
				typeExists = 1;
			}
			j++;
		}
		if(!typeExists){
			strcpy(atomTypes[atomTypesCount].name, topology.atoms[i].type);
			atomTypes[atomTypesCount].charge = topology.atoms[i].charge;
			atomTypes[atomTypesCount].mass = topology.atoms[i].mass;
			atomTypes[atomTypesCount].ignored = topology.atoms[i].ignored;
			atomTypes[atomTypesCount].epsilon = topology.atoms[i].epsilon;
			atomTypes[atomTypesCount].RminOver2 = topology.atoms[i].RminOver2;
			atomTypes[atomTypesCount].ignored_14 = topology.atoms[i].ignored_14;
			atomTypes[atomTypesCount].epsilon_14 = topology.atoms[i].epsilon_14;
			atomTypes[atomTypesCount].RminOver2_14 = topology.atoms[i].RminOver2_14;

			printf("%d: ", atomTypesCount);
			printAtomType(atomTypes[atomTypesCount]);

			atomTypesCount ++;
		} else {
			typeExists = 0;
		}
	}

	for(i = 0; i < topology.atomCount; i++){
		topology.atoms[i].typeId = getAtomTypeId(topology.atoms[i]);
		/*printf("Atom %d(%s) has atom type %s(id: %d)\n",
				atoms[i].id, atoms[i].type, atomTypes[atoms[i].typeId].name, atoms[i].typeId);*/
	}
	testAtomTypes();
	printf("Done looking for atom types.\n");

}

int compareAtoms(Atom a1, Atom a2){
	if(a1.charge == a2.charge && a1.mass == a2.mass &&
			a1.ignored == a2.ignored && a1.epsilon == a2.epsilon && a1.RminOver2 == a2.RminOver2 &&
			a1.ignored_14 == a2.ignored_14 && a1.epsilon_14 == a2.epsilon_14 && a1.RminOver2_14 == a2.RminOver2_14){
		return 1;
	} else {
		return 0;
	}
}

int getAtomTypeId(Atom atom){
	int j;
	for(j = 0; j < atomTypesCount; j++){
		AtomType type = atomTypes[j];
		if(atom.charge == type.charge && atom.mass == type.mass &&
				atom.ignored == type.ignored && atom.epsilon == type.epsilon && atom.RminOver2 == type.RminOver2 &&
				atom.ignored_14 == type.ignored_14 && atom.epsilon_14 == type.epsilon_14 && atom.RminOver2_14 == type.RminOver2_14
		){
			return j;
		}
	}
	DIE("ERROR: Can't find atom type for the following atom %d (%s%d).",
			atom.id, atom.resName, atom.resid);
}

void printAtomType(AtomType atomType){
	printf("%6s;\tq = %5.3f;\tm = %5.3f;\te = %5.3f(%5.3f);\tR/2 = %5.3f(%5.3f).\n",
		atomType.name,
		atomType.charge,
		atomType.mass,
		atomType.epsilon,
		atomType.epsilon_14,
		atomType.RminOver2,
		atomType.RminOver2_14);
}

void testAtomTypes(){
	int i;
	for(i = 0; i < topology.atomCount; i++){
		Atom atom = topology.atoms[i];
		if(atomTypes[atom.typeId].epsilon != atom.epsilon ||
				atomTypes[atom.typeId].RminOver2 != atom.RminOver2 ||
				atomTypes[atom.typeId].charge != atom.charge ||
				atomTypes[atom.typeId].epsilon_14 != atom.epsilon_14 ||
				atomTypes[atom.typeId].RminOver2_14 != atom.RminOver2_14){
			DIE("Atom types error!");
		}
	}
}
