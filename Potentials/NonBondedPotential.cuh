/*
 * NonBondedPtential.cuh
 *
 *  Created on: Aug 5, 2010
 *      Author: zhmurov
 */

#pragma once

namespace non_bonded_potential {

int nonBondedBlockSize;
int nonBondedBlockCount;

#define NB_POTENTIAL_TYPE_RDIE_SHIFT	1
#define NB_POTENTIAL_TYPE_CDIE_SWITCH	2

int nbPotentialType;

typedef struct {
	float* h_ljEnergies;
	float* d_ljEnergies;
	float* h_coulombEnergies;
	float* d_coulombEnergies;
	float2* h_ljparameters;
	float2* d_ljparameters;
	float2* h_ljparameters14;
	float2* d_ljparameters14;
	float* h_ljEpsilon;
	float* d_ljEpsilon;
	float* h_ljR0;
	float* d_ljR0;
	float* h_ljEpsilon14;
	float* d_ljEpsilon14;
	float* h_ljR014;
	float* d_ljR014;
	float* h_charges;
	float* d_charges;
	float* h_charges14;
	float* d_charges14;
	float ljCutoff;
	float ljSwitch;
	float coulombCutoff;
	float ljCutoff2;
	float ljSwitch2;
	float coulombCutoff2;
	float oneOverLjCutoff6;
	float oneOverCoulombCutoff4;
	float oneOverCutoff2MinSwitch2Cub;
	float C1;
	float C2;
	float C3;
	//Atom* d_atoms;
} NonBondedData;

NonBondedData nonBondedData;

Potential nonBondedPotential;
EnergyOutput ljEnergyOutput;
EnergyOutput coulombEnergyOutput;

void create();
void init();
inline void computeRDieShift();
inline void computeCDieSwitch();
inline void computeLJPotentialEnergyRDieShift();
inline void computeCoulombPotentialEnergyRDieShift();
inline void computeLJPotentialEnergyCDieSwitch();
inline void computeCoulombPotentialEnergyCDieSwitch();
void destroy();

} // namespace non_bonded_potential

__device__ __constant__ non_bonded_potential::NonBondedData c_nonBondedData;

__device__ __constant__ float c_charges[MAX_ATOM_TYPES];
texture<float, 1, cudaReadModeElementType> t_charges;

__device__ __constant__ float c_charges14[MAX_ATOM_TYPES];
texture<float, 1, cudaReadModeElementType> t_charges14;

__device__ __constant__ float2 c_ljparameters[MAX_ATOM_TYPES];
texture<float2, 1, cudaReadModeElementType> t_ljparameters;

__device__ __constant__ float2 c_ljparameters14[MAX_ATOM_TYPES];
texture<float2, 1, cudaReadModeElementType> t_ljparameters14;

texture<float, 1, cudaReadModeElementType> t_ljEpsilon;
texture<float, 1, cudaReadModeElementType> t_ljR0;
texture<float, 1, cudaReadModeElementType> t_ljEpsilon14;
texture<float, 1, cudaReadModeElementType> t_ljR014;
