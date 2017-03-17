#pragma once
#define MAX_EL_NUM 1

namespace ccma_constrAlg {

int ccmaBlockSize;
int ccmaBlockCount1;
int ccmaBlockCount2;

typedef struct {
	int 	ccmaConstrCount;
	int		ccmaAtomsCount;
	float 	tol;
	int 	ccmaConstrTot;
	int		ccmaAtomsTot;

	int2*	h_atomPair;
	float*	h_M;
	int*	h_attachedBonds;
	int*	h_matrElIds;
	float*	h_matrEls;
	int*	h_numValidEls;
	float4* h_constrDir;
	int* 	h_ccmaConvergedFlag;

	int2*	d_atomPair;
	float*	d_M;
	int4*	d_attachedBonds;
	int*	d_matrElIds;
	float*	d_matrEls;
	int*	d_numValidEls;
	float4* d_constrDir;
	int* 	d_ccmaConvergedFlag;

	float*	d_constrSigma;
	float*	d_constrLambda;

} CcmaConstrData;

CcmaConstrData ccmaConstrData;
__device__ __constant__ CcmaConstrData c_ccmaConstrData;

int maxRowElements=0;
ConstrAlg constrAlg;

void create();
float findAngle();
void init();
inline void compute();
void destroy();


} // namespace ccma_constrAlg
