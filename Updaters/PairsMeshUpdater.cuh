/*
 * PairsMeshUpdater.cuh
 *
 *  Created on: Sep 28, 2011
 *      Author: zhmurov
 */

#pragma once

#include "../Core/gsystem.h"

namespace pairsmesh {

dim3 meshBlock;
dim3 meshGrid;

typedef struct {

	int minmaxBlockCount;

	float4 minCoord;
	float4 maxCoord;
	float4* h_minCoordArr;
	float4* h_maxCoordArr;
	float4* d_minCoordArr;
	float4* d_maxCoordArr;
	float4 origin;
	float step;
	float4 length;
	float margin;
	float margin2;
	int4 size;
	int* h_atoms;
	int* d_atoms;
} Mesh;

Mesh mesh;
__device__ __constant__ Mesh c_mesh;

Updater updater;

void create();
void init();
void destroy();
void update();

}

texture<int, 1, cudaReadModeElementType> t_mesh;
