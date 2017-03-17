#pragma once

#include <vector_types.h>
#include <stdio.h>

/*typedef struct {
	float x;
	float y;
	float z;
	int type;
} Coord;*/

/*typedef struct __align__ (16) {
	float x;
	float y;
	float z;
	float T;
} Vel;

typedef struct __align__ (16) {
	float x;
	float y;
	float z;
	float m;
} Force;*/

typedef struct {

	int N;
	int Nsim;
	int widthSim;
	int Ntot;
	int widthTot;

	double3* h_realCoord; //redundant

	float4* h_coord;
	float4* h_vel;
	float4* h_forces;

	float4* d_coord;
	float4* d_midcoord;
	float4* d_vel;
	float4* d_forces;

	int* h_atomTypes;
	int* d_atomTypes;

	float* h_m;
	float* d_m;

} GSystem;

typedef struct {
	char name[100];
	void (*compute)();
	void (*destroy)();
} Potential;

typedef struct {
	char name[100];
	int frequency;
	void (*update)();
	void (*destroy)();
} Updater;

typedef struct {
	char name[100];
	float h;
	void (*integrate)();
	void (*finalize)();
	void (*destroy)();
} Integrator;

typedef struct {
	char name[100];
	int  Nconstr;
	void (*compute)();
	void (*destroy)();
} ConstrAlg;

typedef struct {
	char name[100];
	float* values;
	void (*computeValues)();
} EnergyOutput;

// Use addRestarter() function to create restarters
typedef struct {
	char name[100]; // WARNING: name should not conatain whitespace chars (scanf trouble)
	void (*save)(FILE*); // Implementation should write state to FILE
	void (*load)(FILE*); // Implementation should restore state from FILE
} Restarter;

