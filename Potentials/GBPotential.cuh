/*
 * GBPotential.cuh
 *
 *  Created on: 14.05.2012
 *      Author: zhmurov
 */


#pragma once

#define GB_MODEL_HCT_STRING		"hct"
#define GB_MODEL_OBC_STRING		"obc"
#define GB_MODEL_HCT			1
#define GB_MODEL_OBC			2

int gbBlockSize;
int gbBlockCount;

namespace gb_potential {

	typedef struct {

		float* h_rho;
		float* d_rho;
		float* h_S;
		float* d_S;
		float* h_Srho;
		float* d_Srho;

		float* h_alpha;
		float* d_alpha;

		float* h_gbEnergies;
		float* d_gbEnergies;
		float* h_saEnergies;
		float* d_saEnergies;

		float* h_dGda;
		float* d_dGda;
		float4* h_dGdr;
		float4* d_dGdr;

		float pairsCutoff;
		float pairsCutoff2;
		int* h_pairsCount;
		int* d_pairsCount;
		int* h_pairs;
		int* d_pairs;

		int model;
		int saOn;

		float dielOffset;

		float sigma;
		float Rprobe;
		float* h_saMult;
		float* d_saMult;

		float epsilonTimesHalfCoulomb;
		float epsilon;
		float ein, eout;
		float halfCoulomb;
		float kkappa;

		float alpha;
		float beta;
		float gamma;

		float* h_dHdS;
		float* d_dHdS;

	} GBData;

	typedef struct {
		char atomType[ATOM_TYPE_LENGTH];
		float sar;
		float st;
		float pi;
		float gbr;
		float hct;
	} GBParameters;

	GBData gbData;
	__device__ __constant__ GBData c_gbData;

	int gbParCount;
	GBParameters* gbPar;

	Potential potential;
	Updater updater;
	EnergyOutput energyOutputGB;
	EnergyOutput energyOutputSA;

	void create();
	void init();
	void compute();
	void update();
	void computeEnergyGB();
	void computeEnergySA();
	void destroyPotential();
	void destroyUpdater();

}

texture<float, 1, cudaReadModeElementType> t_gbalpha;

