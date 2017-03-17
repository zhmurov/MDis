/*
 * PairsMeshUpdater.cu
 *
 *  Created on: Sep 28, 2011
 *      Author: zhmurov
 */

#include "PairsMeshUpdater.cuh"
#include "../Util/Log.h"

namespace pairsmesh {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<mesh> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

#define MINMAX_BLOCK_SIZE 64
#define MESH_BLOCK_SIZE 8

void inline defineMesh();

void create(){

	if(getYesNoParameter(PARAMETER_GBSW_ON, 0)){
		// Initialize updater
		updater.update = update;
		updater.destroy = destroy;
		updater.frequency = getIntegerParameter(PARAMETER_MESH_FREQ);
		sprintf(updater.name, "Pairs Mesh");
		updaters[updatersCount++] = &updater;

		init();
	}
}

void init(){

	LOG << "Initializing...";

	if(gsystem.Ntot % MINMAX_BLOCK_SIZE == 0){
		mesh.minmaxBlockCount = gsystem.Ntot/MINMAX_BLOCK_SIZE;
	} else {
		mesh.minmaxBlockCount = gsystem.Ntot/MINMAX_BLOCK_SIZE + 1;
	}

	meshBlock = dim3(MESH_BLOCK_SIZE, MESH_BLOCK_SIZE);

	mesh.margin = getFloatParameter(PARAMETER_MESH_MARGIN, DEFAULT_MESH_MARGIN);
	mesh.margin2 = mesh.margin*mesh.margin;
	mesh.step = getFloatParameter(PARAMETER_MESH_STEP, DEFAULT_MESH_STEP);

	allocateCPU((void**)&mesh.h_minCoordArr, mesh.minmaxBlockCount*sizeof(float4));
	allocateCPU((void**)&mesh.h_maxCoordArr, mesh.minmaxBlockCount*sizeof(float4));

	allocateGPU((void**)&mesh.d_minCoordArr, mesh.minmaxBlockCount*sizeof(float4));
	allocateGPU((void**)&mesh.d_maxCoordArr, mesh.minmaxBlockCount*sizeof(float4));

	int i;
	float4 r;
	mesh.minCoord = gsystem.h_coord[0];
	mesh.maxCoord = gsystem.h_coord[0];
	for(i = 1; i < gsystem.N; i++){
		r = gsystem.h_coord[i];

		if(r.x < mesh.minCoord.x){
			mesh.minCoord.x = r.x;
		}
		if(r.y < mesh.minCoord.y){
			mesh.minCoord.y = r.y;
		}
		if(r.z < mesh.minCoord.z){
			mesh.minCoord.z = r.z;
		}

		if(r.x > mesh.maxCoord.x){
			mesh.maxCoord.x = r.x;
		}
		if(r.y > mesh.maxCoord.y){
			mesh.maxCoord.y = r.y;
		}
		if(r.z > mesh.maxCoord.z){
			mesh.maxCoord.z = r.z;
		}
	}

	defineMesh();

	LOG << "Mesh dimensions:";
	//printf("\t x\ty\tz\n");
	//printf("Min:%5.3f\t%5.3f\t%5.3f\n", mesh.minCoord.x, mesh.minCoord.y, mesh.minCoord.z);
	LOG << "MIN: \t" << mesh.minCoord;
	LOG << "MAX: \t" << mesh.maxCoord;
	//printf("Max:%5.3f\t%5.3f\t%5.3f\n", mesh.maxCoord.x, mesh.maxCoord.y, mesh.maxCoord.z);

	LOG << "Mesh size: " << mesh.size;

	int extension = getIntegerParameter(PARAMETER_MESH_EXTENSION, DEFAULT_MESH_EXTENSION);

	allocateGPU((void**)&mesh.d_atoms, (mesh.size.w + extension)*sizeof(int));
	allocateCPU((void**)&mesh.h_atoms, (mesh.size.w + extension)*sizeof(int));

	cudaMemcpyToSymbol(c_mesh, &mesh, sizeof(Mesh), 0, cudaMemcpyHostToDevice);

	LOG << "Initialized.";
}

__global__ void findMinMax_kernel(){
	__shared__ float4 minCoord[MINMAX_BLOCK_SIZE];
	__shared__ float4 maxCoord[MINMAX_BLOCK_SIZE];
	minCoord[threadIdx.x] = make_float4(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX);
	maxCoord[threadIdx.x] = make_float4(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX);
	float4 coord;
	int i = blockIdx.x*blockDim.x;
	i *= 2;
	i += threadIdx.x;
	if(i < c_gsystem.Ntot){
		minCoord[threadIdx.x] = tex1Dfetch(t_coord, i);
		maxCoord[threadIdx.x] = minCoord[threadIdx.x];
		if(i + blockDim.x < c_gsystem.Ntot){
			coord = tex1Dfetch(t_coord, i + blockDim.x);

			minCoord[threadIdx.x].x = fmin(coord.x, minCoord[threadIdx.x].x);
			minCoord[threadIdx.x].y = fmin(coord.y, minCoord[threadIdx.x].y);
			minCoord[threadIdx.x].z = fmin(coord.z, minCoord[threadIdx.x].z);

			maxCoord[threadIdx.x].x = fmax(coord.x, maxCoord[threadIdx.x].x);
			maxCoord[threadIdx.x].y = fmax(coord.y, maxCoord[threadIdx.x].y);
			maxCoord[threadIdx.x].z = fmax(coord.z, maxCoord[threadIdx.x].z);
		}

		__syncthreads();

		int s;
		for(s = blockDim.x/2; s > 32; s >>= 1){
			if(threadIdx.x < s){
				minCoord[threadIdx.x].x = fmin(minCoord[threadIdx.x + s].x, minCoord[threadIdx.x].x);
				minCoord[threadIdx.x].y = fmin(minCoord[threadIdx.x + s].y, minCoord[threadIdx.x].y);
				minCoord[threadIdx.x].z = fmin(minCoord[threadIdx.x + s].z, minCoord[threadIdx.x].z);

				maxCoord[threadIdx.x].x = fmax(maxCoord[threadIdx.x + s].x, maxCoord[threadIdx.x].x);
				maxCoord[threadIdx.x].y = fmax(maxCoord[threadIdx.x + s].y, maxCoord[threadIdx.x].y);
				maxCoord[threadIdx.x].z = fmax(maxCoord[threadIdx.x + s].z, maxCoord[threadIdx.x].z);

				__syncthreads();

			}
		}

		if(threadIdx.x < 32){

			s = threadIdx.x + 32;

			minCoord[threadIdx.x].x = fmin(minCoord[s].x, minCoord[threadIdx.x].x);
			minCoord[threadIdx.x].y = fmin(minCoord[s].y, minCoord[threadIdx.x].y);
			minCoord[threadIdx.x].z = fmin(minCoord[s].z, minCoord[threadIdx.x].z);

			maxCoord[threadIdx.x].x = fmax(maxCoord[s].x, maxCoord[threadIdx.x].x);
			maxCoord[threadIdx.x].y = fmax(maxCoord[s].y, maxCoord[threadIdx.x].y);
			maxCoord[threadIdx.x].z = fmax(maxCoord[s].z, maxCoord[threadIdx.x].z);

			s -= 16;

			minCoord[threadIdx.x].x = fmin(minCoord[s].x, minCoord[threadIdx.x].x);
			minCoord[threadIdx.x].y = fmin(minCoord[s].y, minCoord[threadIdx.x].y);
			minCoord[threadIdx.x].z = fmin(minCoord[s].z, minCoord[threadIdx.x].z);

			maxCoord[threadIdx.x].x = fmax(maxCoord[s].x, maxCoord[threadIdx.x].x);
			maxCoord[threadIdx.x].y = fmax(maxCoord[s].y, maxCoord[threadIdx.x].y);
			maxCoord[threadIdx.x].z = fmax(maxCoord[s].z, maxCoord[threadIdx.x].z);

			s -= 8;

			minCoord[threadIdx.x].x = fmin(minCoord[s].x, minCoord[threadIdx.x].x);
			minCoord[threadIdx.x].y = fmin(minCoord[s].y, minCoord[threadIdx.x].y);
			minCoord[threadIdx.x].z = fmin(minCoord[s].z, minCoord[threadIdx.x].z);

			maxCoord[threadIdx.x].x = fmax(maxCoord[s].x, maxCoord[threadIdx.x].x);
			maxCoord[threadIdx.x].y = fmax(maxCoord[s].y, maxCoord[threadIdx.x].y);
			maxCoord[threadIdx.x].z = fmax(maxCoord[s].z, maxCoord[threadIdx.x].z);

			s -= 4;

			minCoord[threadIdx.x].x = fmin(minCoord[s].x, minCoord[threadIdx.x].x);
			minCoord[threadIdx.x].y = fmin(minCoord[s].y, minCoord[threadIdx.x].y);
			minCoord[threadIdx.x].z = fmin(minCoord[s].z, minCoord[threadIdx.x].z);

			maxCoord[threadIdx.x].x = fmax(maxCoord[s].x, maxCoord[threadIdx.x].x);
			maxCoord[threadIdx.x].y = fmax(maxCoord[s].y, maxCoord[threadIdx.x].y);
			maxCoord[threadIdx.x].z = fmax(maxCoord[s].z, maxCoord[threadIdx.x].z);

			s -= 2;

			minCoord[threadIdx.x].x = fmin(minCoord[s].x, minCoord[threadIdx.x].x);
			minCoord[threadIdx.x].y = fmin(minCoord[s].y, minCoord[threadIdx.x].y);
			minCoord[threadIdx.x].z = fmin(minCoord[s].z, minCoord[threadIdx.x].z);

			maxCoord[threadIdx.x].x = fmax(maxCoord[s].x, maxCoord[threadIdx.x].x);
			maxCoord[threadIdx.x].y = fmax(maxCoord[s].y, maxCoord[threadIdx.x].y);
			maxCoord[threadIdx.x].z = fmax(maxCoord[s].z, maxCoord[threadIdx.x].z);

			s--;

			minCoord[threadIdx.x].x = fmin(minCoord[s].x, minCoord[threadIdx.x].x);
			minCoord[threadIdx.x].y = fmin(minCoord[s].y, minCoord[threadIdx.x].y);
			minCoord[threadIdx.x].z = fmin(minCoord[s].z, minCoord[threadIdx.x].z);

			maxCoord[threadIdx.x].x = fmax(maxCoord[s].x, maxCoord[threadIdx.x].x);
			maxCoord[threadIdx.x].y = fmax(maxCoord[s].y, maxCoord[threadIdx.x].y);
			maxCoord[threadIdx.x].z = fmax(maxCoord[s].z, maxCoord[threadIdx.x].z);
		}

		if(threadIdx.x == 0){
			c_mesh.d_minCoordArr[blockIdx.x] = minCoord[0];
			c_mesh.d_maxCoordArr[blockIdx.x] = maxCoord[0];
		}

	}
}

void findMinMax(int n){
	int numBlocks = n/(2*MINMAX_BLOCK_SIZE);
	findMinMax_kernel<<<mesh.minmaxBlockCount, MINMAX_BLOCK_SIZE>>>();
	if(numBlocks > MINMAX_BLOCK_SIZE){
		findMinMax(numBlocks);
	} else {
		cudaMemcpy(mesh.h_minCoordArr, mesh.d_minCoordArr, numBlocks*sizeof(float4), cudaMemcpyDeviceToHost);
		cudaMemcpy(mesh.h_maxCoordArr, mesh.d_maxCoordArr, numBlocks*sizeof(float4), cudaMemcpyDeviceToHost);
		int i;
		for(i = 1; i < numBlocks; i++){
			mesh.h_minCoordArr[0].x = fmin(mesh.h_minCoordArr[0].x, mesh.h_minCoordArr[i].x);
			mesh.h_minCoordArr[0].y = fmin(mesh.h_minCoordArr[0].y, mesh.h_minCoordArr[i].y);
			mesh.h_minCoordArr[0].z = fmin(mesh.h_minCoordArr[0].z, mesh.h_minCoordArr[i].z);

			mesh.h_maxCoordArr[0].x = fmax(mesh.h_maxCoordArr[0].x, mesh.h_maxCoordArr[i].x);
			mesh.h_maxCoordArr[0].y = fmax(mesh.h_maxCoordArr[0].y, mesh.h_maxCoordArr[i].y);
			mesh.h_maxCoordArr[0].z = fmax(mesh.h_maxCoordArr[0].z, mesh.h_maxCoordArr[i].z);
		}
	}
	mesh.minCoord = mesh.h_minCoordArr[0];
	mesh.maxCoord = mesh.h_maxCoordArr[0];
}

void inline defineMesh(){

	mesh.origin.x = mesh.minCoord.x - mesh.margin;
	mesh.origin.y = mesh.minCoord.y - mesh.margin;
	mesh.origin.z = mesh.minCoord.z - mesh.margin;

	mesh.length.x = mesh.maxCoord.x - mesh.origin.x + mesh.margin;
	mesh.length.y = mesh.maxCoord.y - mesh.origin.y + mesh.margin;
	mesh.length.z = mesh.maxCoord.z - mesh.origin.z + mesh.margin;

	mesh.size.x = mesh.length.x/mesh.step;
	mesh.size.y = mesh.length.y/mesh.step;
	mesh.size.z = mesh.length.z/mesh.step;

	mesh.size.w = mesh.size.x*mesh.size.y*mesh.size.z;

	meshGrid = dim3(mesh.size.x/MESH_BLOCK_SIZE + 1, mesh.size.y/MESH_BLOCK_SIZE + 1);

	//cudaFree((void**)&mesh.d_atoms);
	//allocateGPU((void**)&mesh.d_atoms, mesh.size.w*sizeof(int));

	cudaMemcpyToSymbol(c_mesh, &mesh, sizeof(Mesh), 0, cudaMemcpyHostToDevice);

	//cudaFreeHost((void**)&mesh.h_atoms);
	//allocateCPU((void**)&mesh.h_atoms, mesh.size.w*sizeof(int));

}

__global__ void fillMesh_kernel(float z, int* d_mesh2D){
	int ix = blockIdx.x*blockDim.x + threadIdx.x;
	int iy = blockIdx.y*blockDim.y + threadIdx.y;
	if(ix < c_mesh.size.x && iy < c_mesh.size.y){
		int i;
		float x = c_mesh.origin.x + ix*c_mesh.step;
		float y = c_mesh.origin.y + iy*c_mesh.step;
		float4 coord = tex1Dfetch(t_coord, 0);
		float dx = x - coord.x;
		float dy = y - coord.y;
		float dz = z - coord.z;
		float r = dx*dx + dy*dy + dz*dz;
		float minr = r;
		float mini = 0;
		for(i = 1; i < c_gsystem.Ntot; i++){
			coord = tex1Dfetch(t_coord, i);
			dx = x - coord.x;
			dy = y - coord.y;
			dz = z - coord.z;
			r = dx*dx + dy*dy + dz*dz;
			if(r < minr){
				minr = r;
				mini = i;
			}
		}
		if(minr > c_mesh.margin2){
			mini = -1;
		}
		d_mesh2D[iy*c_mesh.size.x + ix] = mini;
	}

}

void inline fillMesh(){
	int iz;
	float z;
	for(iz = 0; iz < mesh.size.z; iz++){
		z = mesh.origin.z + iz*mesh.step;
		fillMesh_kernel<<<meshGrid, meshBlock>>>(z, &mesh.d_atoms[iz*mesh.size.x*mesh.size.y]);
	}
}

void update(){
	findMinMax(gsystem.Ntot);
	/*printf("\t x\ty\tz\n");
	printf("Min:%5.3f\t%5.3f\t%5.3f\n", mesh.minCoord.x, mesh.minCoord.y, mesh.minCoord.z);
	printf("Max:%5.3f\t%5.3f\t%5.3f\n", mesh.maxCoord.x, mesh.maxCoord.y, mesh.maxCoord.z);*/
	defineMesh();
	fillMesh();
/*	cudaMemcpy(mesh.h_atoms, mesh.d_atoms,
			mesh.size.w*sizeof(int), cudaMemcpyDeviceToHost);


	PDB pdb;
	readPDB("1bi7_ab.pdb", &pdb);

	PDB meshPDB;
	meshPDB.ssCount = 0;
	meshPDB.atomCount = mesh.size.w;
	meshPDB.atoms = (PDBAtom*)calloc(meshPDB.atomCount, sizeof(PDBAtom));
#define XYZ(xx,yy,zz) ((zz)*mesh.size.x*mesh.size.y + (yy)*mesh.size.x + (xx))
	int i, j, k;
	for(i = 0; i < mesh.size.x; i++){
		for(j = 0; j < mesh.size.y; j++){
			for(k = 0; k < mesh.size.z; k++){
				if(mesh.h_atoms[XYZ(i,j,k)] != -1){
					memcpy(&meshPDB.atoms[XYZ(i,j,k)], &pdb.atoms[mesh.h_atoms[XYZ(i,j,k)]], sizeof(PDBAtom));
					meshPDB.atoms[XYZ(i,j,k)].x = (mesh.origin.x + mesh.step*i)*10.0f;
					meshPDB.atoms[XYZ(i,j,k)].y = (mesh.origin.y + mesh.step*j)*10.0f;
					meshPDB.atoms[XYZ(i,j,k)].z = (mesh.origin.z + mesh.step*k)*10.0f;
				} else {
					memcpy(&meshPDB.atoms[XYZ(i,j,k)], &pdb.atoms[0], sizeof(PDBAtom));
					sprintf(meshPDB.atoms[XYZ(i,j,k)].name, "XX");
					meshPDB.atoms[XYZ(i,j,k)].x = 0;
					meshPDB.atoms[XYZ(i,j,k)].y = 0;
					meshPDB.atoms[XYZ(i,j,k)].z = 0;
				}
			}
		}
	}
	writePDB("mesh.pdb", &meshPDB);

	exit(0);*/
}

void destroy(){

}

#undef LOG

}
