/*
 * SASAPotential.cuh
 *
 *  Created on: Aug 6, 2010
 *      Author: zhmurov
 */

#pragma once

namespace sasa_potential {

#define MAX_SASA_PAIRSLIST_ITEMS_PER_ATOM		512
#define MAX_SASA_PAIRS_PER_ATOM					128

int sasaBlockSize;
int sasaBlockCount;
int sasaPairsListBlockSize;
int sasaPairsListBlockCount;

/*typedef struct {
	float Ripr;
	float Ripr2;
	float piOverSi;
	float sigmaSi;
} GSASAParameters;*/

typedef struct {

	int maxPairsListItemsPerAtom;
	int maxSASAPairsPerAtom;

	float sasaPairsListCutoff;
	float sasaCutoff;

	float Rprobe;
	float pij_cov;
	float pij_nb;

	int* h_threadAtom;
	int* d_threadAtom;
	int* h_atomThread;
	int* d_atomThread;
	int threadsCount;
	int threadsCountTot;
	int widthTot;

	int* h_pairsListCount;
	int* d_pairsListCount;
	int* h_pairsList;
	int* d_pairsList;

	int* h_pairs12Counts;
	int* d_pairs12Counts;

	int* h_pairs134Counts;
	int* d_pairs134Counts;

	int* h_pairsBondedCounts;
	int* d_pairsBondedCounts;

	int* h_sasaListCount;
	int* d_sasaListCount;
	float4* h_BGrad;
	float4* d_BGrad;
	float* h_B;
	float* d_B;

	float4* h_BGradT;
	float4* d_BGradT;
	float* h_BT;
	float* d_BT;

	int* h_sasaList;
	int* d_sasaList;

	//GSASAParameters* h_sasaParameters;
	//GSASAParameters* d_sasaParameters;

	float* h_sasaRipr;
	float* h_sasaRipr2;
	float* h_sasaPiOverSi;
	float* h_sasaSigmaSi;
	float* h_sasaSi;

	float* d_sasaRipr;
	float* d_sasaRipr2;
	float* d_sasaPiOverSi;
	float* d_sasaSigmaSi;
	float* d_sasaSi;

	float* h_sasaEnergies;
	float* d_sasaEnergies;

	float* h_sasa;
	float* d_sasa;

} GSASAData;

GSASAData sasaData;
__device__ __constant__ GSASAData c_sasaData;

//__device__ __constant__ GSASAParameters c_sasaParameters[MAX_ATOM_TYPES];

texture<float, 1, cudaReadModeElementType> t_sasaRipr;
texture<float, 1, cudaReadModeElementType> t_sasaRipr2;
texture<float, 1, cudaReadModeElementType> t_sasaPiOverSi;
texture<float, 1, cudaReadModeElementType> t_sasaSigmaSi;

//texture<float4, 1, cudaReadModeElementType> t_sasaBGrad;
//texture<float, 1, cudaReadModeElementType> t_sasaB;


Potential potential;
Updater sasaPairsListUpdater;
EnergyOutput energyOutput;

void create();
void init();
inline void compute();
inline void computeEnergy();
void destroy();

inline void updateSASAPairsList();
void destroySASAPairsListUpdater();

} // namespace sasa_potential

