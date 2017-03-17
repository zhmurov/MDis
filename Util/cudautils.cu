/*
 * cudautils.cu
 *
 *  Created on: Apr 14, 2011
 *      Author: zhmurov
 */

#include "cudautils.cuh"
#include "../Core/global.h"

__global__ void reduce_kernel(float* d_data, int d_output){
	__shared__ float s_data[REDUCE_BLOCKSIZE];
	if(blockIdx.x*blockDim.x + threadIdx.x < c_gsystem.N){
		int s_i = threadIdx.x;
		int i = 2 * blockIdx.x*blockDim.x + threadIdx.x;

		s_data[s_i] = d_data[i] + d_data[i + blockDim.x];

		__syncthreads();

		int s;
		for(s = blockDim.x/2; s > 32; s>>=1){
			if(s_i < s){
				s_data[s_i] += d_data[s_i + s];
			}
			__syncthreads();
		}

		if(s_i < 32){
			s_data[s_i] += d_data[s_i + 32];
			s_data[s_i] += d_data[s_i + 16];
			s_data[s_i] += d_data[s_i + 8];
			s_data[s_i] += d_data[s_i + 4];
			s_data[s_i] += d_data[s_i + 2];
			s_data[s_i] += d_data[s_i + 1];
		}

		if(s_i == 0){
			d_output[blockIdx.x] = s_data[0];
		}
	}
}

float reduce(float* d_data, int N){
	int blockNum = N/REDUCE_BLOCKSIZE + 1;
	float result = 0;
	int i;
	if(d_sums == NULL){
		allocateGPU((void**)&d_sums, blockNum*sizeof(float));
		allocateCPU((void**)&d_sums, blockNum*sizeof(float));
		for(i = 0; i < blockNum; i++){
			d_sums[i] = 0.0f;
		}
	}

	reduce_kernel<<<blockNum, REDUCE_BLOCKSIZE>>>(d_data, d_sums);

	if(blockNum > REDUCE_BLOCKSIZE){
		result = reduce(d_sums, blockNum);
	} else {
		cudaMemcpy(h_sums, d_sums, blockNum*sizeof(float), cudaMemcpyDeviceToHost);

		for(i = 0; i < blockNum; i++){
			result += h_sums[i];
		}
	}

	return result;
}
