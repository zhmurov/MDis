#pragma once

namespace shake_constrAlg {

int shakeBlockSize;
int shakeBlockCount;

typedef struct {
	int 	shakeClustersCount;
	float 	tol;
	int 	shakeClustersTot;
	int* 	h_shakeClusters; 		//(4*int each) 	{X , H1, H2  , H3}
	float*	h_shakeClustersParam;	//(4*float each){mx, mh, invM, d }
	int4*	d_shakeClusters;
	float4*	d_shakeClustersParam;
} ShakeConstrData;

ShakeConstrData shakeConstrData;
__device__ __constant__ ShakeConstrData c_shakeConstrData;


ConstrAlg constrAlg;

void create();
void init();
inline void compute();
void destroy();

} // namespace shake_constrAlg
