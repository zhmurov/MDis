/*
 * GBSWPotential.cuh
 *
 *  Created on: Sep 22, 2011
 *      Author: zhmurov
 */

#pragma once

int gbswBlockSize;
int gbswBlockCount;

namespace gbsw_potential {

	typedef struct {

		float* h_alpha;
		float* d_alpha;

		float* h_dG1;
		float* d_dG1;

		float* h_RPB;
		float* d_RPB;

		float cutoff;
		float w;
		float a0;
		float a1;
		float nabla;

		float threeOverFourW;
		float oneOverFourW3;
		float threeOverFourW3;

		float oneOverNabla;
		float oneOverFourNabla4;
		float oneOverFourPi;

		int angularPointsCount;
		float4* h_angularQuadr;
		float4* d_angularQuadr;

		int radialPointsCount;
		float2* h_radialQuadr;
		float2* d_radialQuadr;

		float* h_gbswEnergies;
		float* d_gbswEnergies;

		int* h_isH;
		int* d_isH;

		float* h_dGda;
		float* d_dGda;
		float4* h_dGdr;
		float4* d_dGdr;


	} GBSWData;

	GBSWData gbswData;
	__device__ __constant__ GBSWData c_gbswData;

	Potential potential;

	void create();
	void init();
	void compute();
	void destroy();

}

texture<float4, 1, cudaReadModeElementType> t_angularQuadr;
texture<float2, 1, cudaReadModeElementType> t_radialQuadr;

texture<float, 1, cudaReadModeElementType> t_gbswalpha;
