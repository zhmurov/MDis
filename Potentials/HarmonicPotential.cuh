/*
 * HarmonicPotential.cuh
 *
 *  Created on: Aug 2, 2010
 *      Author: zhmurov
 */

#pragma once

namespace harmonic_potential {
	
int harmonicBlockSize;
int harmonicBlockCount;

typedef struct __align__(16) {
	int j;
	float kb;
	float b0;
} GHarmonicPair;

typedef struct {
	int maxHarmonicPerAtom;
	int* h_harmonicBondCount;
	int* d_harmonicBondCount;
	int* h_harmonicCount;
	int* d_harmonicCount;
	GHarmonicPair* h_harmonic;
	GHarmonicPair* d_harmonic;
	float* d_harmonicBondEnergy;
	float* h_harmonicBondEnergy;
	float* d_ureyBradleyEnergy;
	float* h_ureyBradleyEnergy;
} GHarmonicData;

GHarmonicData harmonicData;
__device__ __constant__ GHarmonicData c_harmonicData;

Potential potential;
EnergyOutput harmonicBondEnergyOutput;
EnergyOutput ureyBradleyEnergyOutput;

void create();
void init();
inline void compute();
inline void computeHarmonicBondEnergy();
inline void computeUreyBradleyEnergy();
void destroy();

} // namespace harmonic_potential
