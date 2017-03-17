/*
 * DihedralPotential.cuh
 *
 *  Created on: Aug 4, 2010
 *      Author: zhmurov
 */

#pragma once

namespace dihedral_potential {

int dihedralBlockSize;
int dihedralBlockCount;
int dihedralSummBlockSize;
int dihedralSummBlockCount;

typedef struct {
	float kchi;
	float delta;
	float n;
	int multiplicityRef;
} GDihedralParameters;

typedef struct {
	int D;
	int Dtot;
	int maxDihedralsPerAtom;
	int* h_dihedralCount;
	int* d_dihedralCount;
	int4* h_dihedrals;
	int4* d_dihedrals;
	int4* h_dihedralRefs;
	int4* d_dihedralRefs;
	GDihedralParameters* h_dihedralParameters;
	GDihedralParameters* d_dihedralParameters;
	GDihedralParameters* h_multiplicityParameters;
	GDihedralParameters* d_multiplicityParameters;
	float4* h_dihedralForces;
	float4* d_dihedralForces;
	float* h_dihedralEnergies;
	float* d_dihedralEnergies;
} GDihedralData;

GDihedralData dihedralData;
__device__ __constant__ GDihedralData c_dihedralData;

Potential potential;
EnergyOutput energyOutput;

void create();
void init();
inline void compute();
inline void computeEnergy();
void destroy();

} // namespace dihedral_potential;
