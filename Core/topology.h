/*
 * topology.h
 *
 *  Created on: Jul 29, 2010
 *      Author: zhmurov
 */

#pragma once

#define ATOM_NAME_LENGTH	10
#define ATOM_TYPE_LENGTH	10
#define RESID_NAME_LENGTH	10
#define ATOM_SEGMENT_LENGTH	10

#define MAX_DIHEDRAL_MULTIPLICITY 16
#define MAX_IMPROPER_MULTIPLICITY 4

typedef struct {
	int typeId;
	int id;
	char name[ATOM_NAME_LENGTH];
	char type[ATOM_TYPE_LENGTH];
	int resid;
	char resName[RESID_NAME_LENGTH];
	char chain;
	char segment[ATOM_SEGMENT_LENGTH];
	double x;
	double y;
	double z;
	float vx;
	float vy;
	float vz;
	float charge;
	float mass;
	float ignored;
	float epsilon;
	float RminOver2;
	float ignored_14;
	float epsilon_14;
	float RminOver2_14;
} Atom;

typedef struct {
	int i;
	int j;
	int type;
	float kb;
	float b0;
} Bond;

typedef struct {
	int i;
	int j;
	int k;
	int type;
	float ktheta;
	float theta0;
} Angle;

typedef struct {
	int i;
	int j;
	float kub;
	float s0;
} UreyBradley;

typedef struct {
	int i;
	int j;
	int k;
	int l;
	int type;
	int multiplicity;
	float kchi[MAX_DIHEDRAL_MULTIPLICITY];
	int n[MAX_DIHEDRAL_MULTIPLICITY];
	float delta[MAX_DIHEDRAL_MULTIPLICITY];
} Dihedral;

typedef struct {
	int i;
	int j;
	int k;
	int l;
	int type;
	int multiplicity;
	float kpsi[MAX_IMPROPER_MULTIPLICITY];
	int n[MAX_IMPROPER_MULTIPLICITY];
	float psi0[MAX_IMPROPER_MULTIPLICITY];
} Improper;

typedef struct {
	int i1;
	int j1;
	int k1;
	int l1;
	int i2;
	int j2;
	int k2;
	int l2;
	int type;
} CMAPCorrection;

typedef struct {
	int i;
	int j;
} ExplicitExclusion;

typedef struct {
	int atomCount;
	int bondCount;
	int angleCount;
	int ureyBradleyCount;
	int dihedralCount;
	int improperCount;
	int cmapCount;
	int exclusionsCount;
	Atom* atoms;
	Bond* bonds;
	Angle* angles;
	UreyBradley* ureyBradley;
	Dihedral* dihedrals;
	Improper* impropers;
	CMAPCorrection* cmaps;
	ExplicitExclusion* exclusions;
} Topology;

extern Topology topology;

extern void initTopology();

void readCoordinates(char* filename);
void readVelocities(char* filename);

void generateVelocities(float T);

int findAtom(int atomId);
