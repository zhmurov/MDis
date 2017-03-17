/*
 * memory.cu
 *
 *  Created on: Nov 9, 2010
 *      Author: zhmurov
 */
#include <stdio.h>
#include "memory.h"
#include "Cuda.h"

long int memoryGPU = 0;
long int memoryCPU = 0;

void allocateCPU(void** pointer, long size, int pinned){
    printf("CPU+: ");
	printFormatedMemory(size);
    printf("\n");
#ifndef NDEBUG
	pinned = 0;
#endif
    if (pinned) {
 	    cudaMallocHost(pointer, size);
	    checkCUDAError("allocating memory on host");
	} else {
	    *pointer = malloc(size);
    }
	memoryCPU += size;
	printf("CPU: ");
	printFormatedMemory(memoryCPU);
	printf("\n");
}

void allocateGPU(void** pointer, long size){
    printf("GPU+: ");
	printFormatedMemory(size);
    printf("\n");
	cudaMalloc(pointer, size);
	checkCUDAError("allocating memory on device");
	memoryGPU += size;
	printf("GPU: ");
	printFormatedMemory(memoryGPU);
	printf("\n");
}

void allocateGPUPitch(void** pointer, size_t* pitch, int width, int height){
	cudaMallocPitch(pointer, pitch, width, height);
	checkCUDAError("allocating memory on device");
	memoryGPU += (*pitch) * height;
	printf("GPU: ");
	printFormatedMemory(memoryGPU);
	printf("\n");
}


void freeCPU(void** pointer, long size){
    printf("CPU-: ");
	printFormatedMemory(size);
    printf("\n");
	cudaFreeHost(*pointer);
	checkCUDAError("releasing memory on host");
	memoryCPU -= size;
	printf("CPU: ");
	printFormatedMemory(memoryCPU);
	printf("\n");
}

void freeGPU(void** pointer, long size){
    printf("GPU-: ");
	printFormatedMemory(size);
    printf("\n");
	cudaFree(*pointer);
	checkCUDAError("releasing memory on device");
	memoryGPU -= size;
	printf("GPU: ");
	printFormatedMemory(memoryGPU);
	printf("\n");
}

void printMemoryUsed(){
	printf("Memory used on host: ");
	printFormatedMemory(memoryCPU);
	printf("\nMemory used on device: ");
	printFormatedMemory(memoryGPU);
	printf("\n");
}

void printFormatedMemory(long int bytes){
	long int kb = bytes/1024;
	long int mb = kb/1024;
	if(mb > 0){
		printf("%ld%s", mb, MEMORY_MB);
	} else
	if(kb > 0){
		printf("%ld%s", kb, MEMORY_KB);
	} else {
		printf("%ld%s", bytes, MEMORY_BYTES);
	}
}
