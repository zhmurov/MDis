/*
 * FragmemPotential.cu
 *
 *  Created on: Mar 28, 2018
 *      Author: zhmurov
 */

#include "../Core/global.h"
#include "../Core/md.cuh"
#include "../Util/Log.h"
#include "FragmemPotential.cuh"

namespace fragmem_potential {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<fragmem_potential> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void create(){
	if (!getYesNoParameter(PARAMETER_FRAGMEM_POTENTIAL_ON, 0))
		return;
	potential.compute = &compute;
	potential.destroy = &destroy;
	sprintf(potential.name, "Fragmem potential");
	potentials[potentialsCount] = &potential;
	potentialsCount ++;
	energyOutput.computeValues = &computeEnergy;
	allocateCPU((void**)&energyOutput.values, parameters.Ntr*sizeof(float));
	strcpy(energyOutput.name, ENERGY_OUTPUT_NAME_FRAGMEM);
	energyOutputs[energyOutputsCount] = &energyOutput;
	energyOutputsCount ++;
	init();
}

void init(){
	LOG << "Initializing fragmem potential...";
	fragmemBlockSize = BLOCK_SIZE;
	fragmemBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;

	int i, j;

	allocateCPU((void**)&fragmemData.h_fragmemCount, gsystem.Ntot*sizeof(int));
	allocateGPU((void**)&fragmemData.d_fragmemCount, gsystem.Ntot*sizeof(int));
	allocateCPU((void**)&fragmemData.h_fragmemEnergy, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&fragmemData.d_fragmemEnergy, gsystem.Ntot*sizeof(float));

	fragmemData.strength = getFloatParameter(PARAMETER_FRAGMEM_STRENGTH);

	for(i = 0; i < gsystem.N; i++){
		fragmemData.h_fragmemCount[i] = 0;
	}

	char* pch;
	char filename[1024];
	getMaskedParameter(filename, PARAMETER_FRAGMEM_FILE);
	char buffer[BUF_SIZE];
	FILE* file = fopen(filename, "r");
	if(file == NULL){
		LOG << "Fragment memory file '" << filename << "' not found";
		exit(0);
	}

	while(fgets(buffer, BUF_SIZE, file) != NULL){
		pch = strtok(buffer, " \t");
		i = atoi(pch);
		pch = strtok(NULL, " \t");
		j = atoi(pch);
		fragmemData.h_fragmemCount[i] ++;
		fragmemData.h_fragmemCount[j] ++;
	}


	fragmemData.maxFragmemsPerAtom = 0;
	for(i = 0; i < gsystem.N; i++){
		if(fragmemData.h_fragmemCount[i] > fragmemData.maxFragmemsPerAtom){
			fragmemData.maxFragmemsPerAtom = fragmemData.h_fragmemCount[i];
		}
	}

	LOG << "Maximum fragmem pairs per atom is " << fragmemData.maxFragmemsPerAtom;

	allocateCPU((void**)&fragmemData.h_fragmems,
			gsystem.widthTot*fragmemData.maxFragmemsPerAtom*sizeof(GFragmemPair));
	allocateGPU((void**)&fragmemData.d_fragmems,
			gsystem.widthTot*fragmemData.maxFragmemsPerAtom*sizeof(GFragmemPair));

	rewind(file);

	for(i = 0; i < gsystem.N; i++){
		fragmemData.h_fragmemCount[i] = 0;
	}

	while(fgets(buffer, BUF_SIZE, file) != NULL){
		pch = strtok(buffer, " \t");
		i = atoi(pch);
		pch = strtok(NULL, " \t");
		j = atoi(pch);
		pch = strtok(NULL, " \t");
		float r0 = atof(pch) / 10.0f;
		pch = strtok(NULL, " \t");
		float weight = atof(pch);
		pch = strtok(NULL, " \t");
		float sigma = atof(pch) / 10.0f;
		float oneOverSigma2 = 1.0f/(sigma*sigma);

		fragmemData.h_fragmems[fragmemData.h_fragmemCount[i]*gsystem.widthTot + i].j = j;
		fragmemData.h_fragmems[fragmemData.h_fragmemCount[i]*gsystem.widthTot + i].r0 = r0;
		fragmemData.h_fragmems[fragmemData.h_fragmemCount[i]*gsystem.widthTot + i].oneOverSigma2 = oneOverSigma2;
		fragmemData.h_fragmems[fragmemData.h_fragmemCount[i]*gsystem.widthTot + i].weight = weight;

		fragmemData.h_fragmems[fragmemData.h_fragmemCount[j]*gsystem.widthTot + j].j = i;
		fragmemData.h_fragmems[fragmemData.h_fragmemCount[j]*gsystem.widthTot + j].r0 = r0;
		fragmemData.h_fragmems[fragmemData.h_fragmemCount[j]*gsystem.widthTot + j].oneOverSigma2 = oneOverSigma2;
		fragmemData.h_fragmems[fragmemData.h_fragmemCount[j]*gsystem.widthTot + j].weight = weight;

		fragmemData.h_fragmemCount[i] ++;
		fragmemData.h_fragmemCount[j] ++;
	}

	fclose(file);

	for(i = 0; i < gsystem.Ntot; i++){
		fragmemData.h_fragmemEnergy[i] = 0.0f;
	}

	int traj, m, itot;
	for(traj = 1; traj < parameters.Ntr; traj++){
		for(i = 0; i < gsystem.N; i++){
			itot = traj*gsystem.N + i;
			fragmemData.h_fragmemCount[itot] = fragmemData.h_fragmemCount[i];
			for(m = 0; m < fragmemData.maxFragmemsPerAtom; m++){
				fragmemData.h_fragmems[m*gsystem.widthTot + itot].j =
						fragmemData.h_fragmems[m*gsystem.widthTot + i].j + traj*gsystem.N;
				fragmemData.h_fragmems[m*gsystem.widthTot + itot].r0 =
						fragmemData.h_fragmems[m*gsystem.widthTot + i].r0;
				fragmemData.h_fragmems[m*gsystem.widthTot + itot].oneOverSigma2 =
						fragmemData.h_fragmems[m*gsystem.widthTot + i].oneOverSigma2;
				fragmemData.h_fragmems[m*gsystem.widthTot + itot].weight =
						fragmemData.h_fragmems[m*gsystem.widthTot + i].weight;
			}
		}
	}


	cudaMemcpy(fragmemData.d_fragmemCount, fragmemData.h_fragmemCount,
			gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(fragmemData.d_fragmems, fragmemData.h_fragmems,
			gsystem.widthTot*fragmemData.maxFragmemsPerAtom*sizeof(GFragmemPair), cudaMemcpyHostToDevice);
	cudaMemcpy(fragmemData.d_fragmemEnergy, fragmemData.h_fragmemEnergy,
			gsystem.Ntot*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_fragmemData, &fragmemData,
			sizeof(GFragmemData), 0, cudaMemcpyHostToDevice);
	LOG << "Done initializing fragmem potential.";
}

__global__ void fragmemPotential_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 coord = tex1Dfetch(t_coord, d_i);
		float4 f = c_gsystem.d_forces[d_i];
		float4 r2;
		float r, mult;
		int i;
		GFragmemPair fragmem;
		for(i = 0; i < c_fragmemData.d_fragmemCount[d_i]; i++){
			fragmem = c_fragmemData.d_fragmems[i*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, fragmem.j);//c_gsystem.d_coord[bond.j];
			r2 -= coord;
			DO_PBC(r2);
			r = sqrtf(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);

			float dr = r - fragmem.r0;
			float mult1 = dr*fragmem.oneOverSigma2;
			float arg = -0.5f*dr*mult1;

			mult = c_fragmemData.strength*fragmem.weight*expf(arg)*mult1/r;

			f.x += mult*r2.x;
			f.y += mult*r2.y;
			f.z += mult*r2.z;
		}
		c_gsystem.d_forces[d_i] = f;
	}
}

inline void compute(){
	fragmemPotential_kernel<<<fragmemBlockCount, fragmemBlockSize>>>();
	/*cudaMemcpy(gsystem.h_forces, gsystem.d_forces, atomCount*sizeof(float4), cudaMemcpyDeviceToHost);
	int i;
	float3 force = make_float3(0.0f, 0.0f, 0.0f);
	for(i = 0; i < atomCount; i++){
		force.x += gsystem.h_forces[i].x;
		force.y += gsystem.h_forces[i].y;
		force.z += gsystem.h_forces[i].z;
		printf("%d: (%f, %f, %f) %f\n", i, gsystem.h_forces[i].x, gsystem.h_forces[i].y, gsystem.h_forces[i].z,
					sqrtf(gsystem.h_forces[i].x*gsystem.h_forces[i].x + gsystem.h_forces[i].y*gsystem.h_forces[i].y + gsystem.h_forces[i].z*gsystem.h_forces[i].z));
	}
	printf("Net force (harmonic): (%f, %f, %f) %f\n", force.x, force.y, force.z,
			sqrtf(force.x*force.x + force.y*force.y + force.z*force.z));*/
}

__global__ void fragmemPotentialEnergy_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 coord = tex1Dfetch(t_coord, d_i);
		float pot = 0.0f;
		float4 r2;
		float r;
		int i;
		GFragmemPair fragmem;
		for(i = 0; i < c_fragmemData.d_fragmemCount[d_i]; i++){
			fragmem = c_fragmemData.d_fragmems[i*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, fragmem.j);//c_gsystem.d_coord[bond.j];
			r2 -= coord;
			DO_PBC(r2);
			r = sqrtf(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);

			float dr = r - fragmem.r0;
			float arg = -0.5f*dr*dr*fragmem.oneOverSigma2;

			pot += -c_fragmemData.strength*fragmem.weight*expf(arg);
		}
		c_fragmemData.d_fragmemEnergy[d_i] = pot;
	}
}

inline void computeEnergy(){
	fragmemPotentialEnergy_kernel<<<fragmemBlockCount, fragmemBlockSize>>>();
	cudaMemcpy(fragmemData.h_fragmemEnergy,fragmemData.d_fragmemEnergy,
				gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		float pot = 0.0f;
		for(i = 0; i < gsystem.N; i++){
			pot += fragmemData.h_fragmemEnergy[i + traj*gsystem.N];
		}
		energyOutput.values[traj] = pot*0.5f;
	}
	checkCUDAError("fragmem energy");


	/*int traj, m, i, j, itot, jtot, count, totcount;
	double r, r0;
	for(traj = 0; traj < parameters.Ntr; traj++){
		count = 0;
		totcount = 0;
		for(i = 0; i < gsystem.N; i++){
			itot = traj*gsystem.N + i;
			for(m = 0; m < fragmemData.maxFragmemsPerAtom; m++){
				j = fragmemData.h_fragmems[m*gsystem.widthTot + itot].j;
				jtot = fragmemData.h_fragmems[m*gsystem.widthTot + i].j + traj*gsystem.N;
				r0 = fragmemData.h_fragmems[m*gsystem.widthTot + i].r0;
				float dx = gsystem.h_coord[itot].x - gsystem.h_coord[jtot].x;
				float dy = gsystem.h_coord[itot].y - gsystem.h_coord[jtot].y;
				float dz = gsystem.h_coord[itot].z - gsystem.h_coord[jtot].z;
				r = sqrtf(dx*dx + dy*dy + dz*dz);
				if(fabs(r - r0) < 0.2){
					count ++;
				}
				totcount ++;
			}
		}
		printf("%d\t%f\n", traj, double(count)/double(totcount));
	}*/
}

void destroy(){

}

#undef LOG

} // namespace harmonic_potential
