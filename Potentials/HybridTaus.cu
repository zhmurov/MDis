/*
 * ht.cu
 *
 *  Created on: Feb 17, 2010
 *      Author: zhmurov
 */

#include "../Util/ran2.h"
#include "../Util/Log.h"

namespace hybrid_taus {
	
#define jflone 0x3f800000
#define jflmsk 0x007fffff
#define EPS 1.0e-8f

#include "../Util/Log.h"

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<hybrid_taus> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

int hasGauss = 0;
float gauss;

struct HybTau{
	uint4* h_seeds;
	uint4* d_seeds;
	uint4 mseed;
};

HybTau ht;
__device__ __constant__ HybTau c_ht;

void generateSeeds(uint4* seeds, int seed, int Np);

int initRand(int seed, int Np){
	LOG << "Initializing Hybrid Taus PRNG...";
	allocateCPU((void**)&ht.h_seeds, Np*sizeof(uint4));
	generateSeeds(ht.h_seeds, seed, Np);
	allocateGPU((void**)&ht.d_seeds, Np*sizeof(uint4));
	cudaMemcpy(ht.d_seeds, ht.h_seeds, Np*sizeof(uint4), cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_ht, &ht, sizeof(HybTau), 0, cudaMemcpyHostToDevice);
	LOG << "done";
	return 1;
}

void generateSeeds(uint4* seeds, int seed, int Np){
	for(int i = 0; i < Np; i++){
		do{
			seeds[i].x = (unsigned)(ran2::ran2(&seed)*UINT_MAX);
		} while(seeds[i].x < 128);
		do{
			seeds[i].y = (unsigned)(ran2::ran2(&seed)*UINT_MAX);
		} while(seeds[i].y < 128);
		do{
			seeds[i].z = (unsigned)(ran2::ran2(&seed)*UINT_MAX);
		} while(seeds[i].z < 128);
		do{
			seeds[i].w = (unsigned)(ran2::ran2(&seed)*UINT_MAX);
		} while(seeds[i].w < 128);
	}
	ht.mseed = seeds[0];
}

// Random number generator.
// Taked from GPU Gems 3, Chapter 37

__device__ inline float uintToFloat(unsigned uint){
	unsigned itemp = jflone | (jflmsk & uint);
	float result = (*(float *)&itemp) - 1.0f;
	if(result == 0){
		return EPS;
	} else {
		return result;
	}
}

__device__ inline unsigned TausStep(unsigned &z, int S1, int S2, int S3, unsigned M){
	unsigned b = (((z << S1)^z) >> S2);
	return z = (((z & M) << S3) ^b);
}

__device__ inline unsigned LCGStep(unsigned &z, unsigned A, unsigned C){
	return z = (A * z + C);
}

__device__ inline unsigned HybridTaus(uint4 &seed){
	return 	TausStep(seed.x, 13, 19, 12, 4294967294UL) ^
		TausStep(seed.y,  2, 25,  4, 4294967288UL) ^
		TausStep(seed.z,  3, 11, 17, 4294967280UL) ^
		LCGStep(seed.w, 1664525, 1013904223UL);
	/*return 2.3283064365387e-10 * (
			TausStep(seed.x, 13, 19, 12, 4294967294UL) ^
			TausStep(seed.y,  2, 25,  4, 4294967288UL) ^
			TausStep(seed.z,  3, 11, 17, 4294967280UL) ^
			LCGStep(seed.w, 1664525, 1013904223UL));*/
}


__device__ inline float4 rforce(int d_i){
	uint4 seed = c_ht.d_seeds[d_i];
	float4 result;
	float r = sqrtf(-2.0f * logf(uintToFloat(HybridTaus(seed))));
	float theta = 2.0f*M_PI*uintToFloat(HybridTaus(seed));
	result.x = r*__sinf(theta);
	result.y = r*__cosf(theta);
	r = sqrtf(-2.0f * logf(uintToFloat(HybridTaus(seed))));
	theta = 2.0f*M_PI*uintToFloat(HybridTaus(seed));
	result.z = r*__sinf(theta);
	result.w = r*__cosf(theta);
	c_ht.d_seeds[d_i] = seed;
	return result;
}

#undef jflone
#undef jflmsk
#undef EPS 


#undef LOG

} // namespace hybrid_taus

