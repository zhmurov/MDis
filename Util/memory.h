/*
 * memory.cuh
 *
 *  Created on: Nov 9, 2010
 *      Author: zhmurov
 */

#pragma once

#define MEMORY_BYTES "bytes"
#define MEMORY_KB	"KB"
#define MEMORY_MB	"MB"

extern long int memoryGPU;
extern long int memoryCPU;

// In `allocateCPU', setting parameter `pinned' to 0 will force using 
// `malloc' instead of `cudaMallocHost'. It should be used if you allocate
// big chunks of host memory that are not supposed to be copied to/from GPU
void allocateCPU(void** pointer, long size, int pinned = 1); 
void allocateGPU(void** pointer, long size);
void allocateGPUPitch(void** pointer, size_t* pitch, int width, int height);

void freeCPU(void** pointer, long size);
void freeGPU(void** pointer, long size);

void printMemoryUsed();
void printFormatedMemory(long bytes);
