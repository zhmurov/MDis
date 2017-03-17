/*
 * HarmonicConstraints.cuh
 *
 *  Created on: May 25, 2011
 *      Author: zhmurov
 */

#pragma once

namespace harmonic_constraints {

int harmonicConstraintsBlockSize;
int harmonicConstraintsBlockCount;

typedef struct {
	float Ks; // K-spring for constraint, kJ/mol/nm^2
	float4* h_fixedConstraints; // (x,y,z) --- desired coordinates of atom, w --- spring constant
	float4* d_fixedConstraints;
	int maxRelativeConstraintsPerAtom;
	int* h_relativeConstraintsCount;
	int* d_relativeConstraintsCount;
	int* h_relativeConstraints;
	int* d_relativeConstraints;
	float* h_relativeConstraintsData;
	float* d_relativeConstraintsData;
} HarmonicConstraintsData;

HarmonicConstraintsData harmonicConstraintsData;
__device__ __constant__ HarmonicConstraintsData c_harmonicConstraintsData;

Potential potential;

void create();
void init();
inline void compute();
void destroy();

} // namespace harmonic_constraints
