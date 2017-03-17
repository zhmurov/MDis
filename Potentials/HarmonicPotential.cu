/*
 * HarmonicPotential.cu
 *
 *  Created on: Aug 2, 2010
 *      Author: zhmurov
 */

#include "../Core/global.h"
#include "../Core/md.cuh"
#include "../Util/Log.h"
#include "HarmonicPotential.cuh"

namespace harmonic_potential {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<harmonic_potential> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void create(){
	potential.compute = &compute;
	potential.destroy = &destroy;
	sprintf(potential.name, "Harmonic potential");
	potentials[potentialsCount] = &potential;
	potentialsCount ++;
	harmonicBondEnergyOutput.computeValues = &computeHarmonicBondEnergy;
	allocateCPU((void**)&harmonicBondEnergyOutput.values, parameters.Ntr*sizeof(float));
	strcpy(harmonicBondEnergyOutput.name, ENERGY_OUTPUT_NAME_HARMONIC);
	energyOutputs[energyOutputsCount] = &harmonicBondEnergyOutput;
	energyOutputsCount ++;
	ureyBradleyEnergyOutput.computeValues = &computeUreyBradleyEnergy;
	allocateCPU((void**)&ureyBradleyEnergyOutput.values, parameters.Ntr*sizeof(float));
	strcpy(ureyBradleyEnergyOutput.name, ENERGY_OUTPUT_NAME_UREY_BRADLEY);
	energyOutputs[energyOutputsCount] = &ureyBradleyEnergyOutput;
	energyOutputsCount ++;
	init();
}

void init(){
	LOG << "Initializing harmonic potential...";
	harmonicBlockSize = BLOCK_SIZE;
	harmonicBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;

	int pp, i, j;

	allocateCPU((void**)&harmonicData.h_harmonicCount, gsystem.Ntot*sizeof(int));
	allocateGPU((void**)&harmonicData.d_harmonicCount, gsystem.Ntot*sizeof(int));

	allocateCPU((void**)&harmonicData.h_harmonicBondCount, gsystem.Ntot*sizeof(int));
	allocateGPU((void**)&harmonicData.d_harmonicBondCount, gsystem.Ntot*sizeof(int));

	for(i = 0; i < gsystem.N; i++){
		harmonicData.h_harmonicCount[i] = 0;
	}
	for(pp = 0; pp < topology.bondCount; pp++){
		i = topology.bonds[pp].i;
		j = topology.bonds[pp].j;
		harmonicData.h_harmonicCount[i] ++;
		harmonicData.h_harmonicCount[j] ++;
	}
	for(pp = 0; pp < topology.ureyBradleyCount; pp++){
		i = topology.ureyBradley[pp].i;
		j = topology.ureyBradley[pp].j;
		harmonicData.h_harmonicCount[i] ++;
		harmonicData.h_harmonicCount[j] ++;
	}
	harmonicData.maxHarmonicPerAtom = 0;//MAX_COVALENT_PER_ATOM;
	for(i = 0; i < gsystem.N; i++){
		if(harmonicData.h_harmonicCount[i] > harmonicData.maxHarmonicPerAtom){
			harmonicData.maxHarmonicPerAtom = harmonicData.h_harmonicCount[i];
		}
	}
	LOG << "Maximum harmonic potential pairs (covalent bonds and Urey-Bradley pairs) per atom is " << harmonicData.maxHarmonicPerAtom;

	allocateCPU((void**)&harmonicData.h_harmonic,
			gsystem.widthTot*harmonicData.maxHarmonicPerAtom*sizeof(GHarmonicPair));
	allocateGPU((void**)&harmonicData.d_harmonic,
			gsystem.widthTot*harmonicData.maxHarmonicPerAtom*sizeof(GHarmonicPair));

	allocateCPU((void**)&harmonicData.h_harmonicBondEnergy, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&harmonicData.d_harmonicBondEnergy, gsystem.Ntot*sizeof(float));
	allocateCPU((void**)&harmonicData.h_ureyBradleyEnergy, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&harmonicData.d_ureyBradleyEnergy, gsystem.Ntot*sizeof(float));

	for(i = 0; i < gsystem.Ntot; i++){
		harmonicData.h_harmonicBondEnergy[i] = 0.0f;
		harmonicData.h_ureyBradleyEnergy[i] = 0.0f;
	}

	for(i = 0; i < gsystem.N; i++){
		harmonicData.h_harmonicCount[i] = 0;
		harmonicData.h_harmonicBondCount[i] = 0;
	}
	int totalHarmonic = 0;
	for(pp = 0; pp < topology.bondCount; pp++){
		i = topology.bonds[pp].i;
		j = topology.bonds[pp].j;
		harmonicData.h_harmonic[harmonicData.h_harmonicCount[i]*gsystem.widthTot + i].j = j;
		harmonicData.h_harmonic[harmonicData.h_harmonicCount[j]*gsystem.widthTot + j].j = i;
		harmonicData.h_harmonic[harmonicData.h_harmonicCount[i]*gsystem.widthTot + i].kb = 2.0f*topology.bonds[pp].kb;
		harmonicData.h_harmonic[harmonicData.h_harmonicCount[j]*gsystem.widthTot + j].kb = 2.0f*topology.bonds[pp].kb;
		harmonicData.h_harmonic[harmonicData.h_harmonicCount[i]*gsystem.widthTot + i].b0 = topology.bonds[pp].b0;
		harmonicData.h_harmonic[harmonicData.h_harmonicCount[j]*gsystem.widthTot + j].b0 = topology.bonds[pp].b0;
		harmonicData.h_harmonicCount[i] ++;
		harmonicData.h_harmonicCount[j] ++;
		harmonicData.h_harmonicBondCount[i] ++;
		harmonicData.h_harmonicBondCount[j] ++;
		if(harmonicData.h_harmonicCount[i] > harmonicData.maxHarmonicPerAtom
				|| harmonicData.h_harmonicCount[j] > harmonicData.maxHarmonicPerAtom){
			DIE("Maximum number of harmonic bonds exceeded the limit of %d.", harmonicData.maxHarmonicPerAtom );
		}
		totalHarmonic++;
	}

	for(pp = 0; pp < topology.ureyBradleyCount; pp++){
		i = topology.ureyBradley[pp].i;
		j = topology.ureyBradley[pp].j;
		if(topology.ureyBradley[pp].kub != 0.0f){
			harmonicData.h_harmonic[harmonicData.h_harmonicCount[i]*gsystem.widthTot + i].j = j;
			harmonicData.h_harmonic[harmonicData.h_harmonicCount[j]*gsystem.widthTot + j].j = i;
			harmonicData.h_harmonic[harmonicData.h_harmonicCount[i]*gsystem.widthTot + i].kb = 2.0f*topology.ureyBradley[pp].kub;
			harmonicData.h_harmonic[harmonicData.h_harmonicCount[j]*gsystem.widthTot + j].kb = 2.0f*topology.ureyBradley[pp].kub;
			harmonicData.h_harmonic[harmonicData.h_harmonicCount[i]*gsystem.widthTot + i].b0 = topology.ureyBradley[pp].s0;
			harmonicData.h_harmonic[harmonicData.h_harmonicCount[j]*gsystem.widthTot + j].b0 = topology.ureyBradley[pp].s0;
			harmonicData.h_harmonicCount[i] ++;
			harmonicData.h_harmonicCount[j] ++;
			if(harmonicData.h_harmonicCount[i] > harmonicData.maxHarmonicPerAtom
					|| harmonicData.h_harmonicCount[j] > harmonicData.maxHarmonicPerAtom){
				DIE("Maximum number of harmonic bonds exceeded the limit of %d.", harmonicData.maxHarmonicPerAtom );
			}
			totalHarmonic++;
		}
	}

	int traj, h, itot;
	for(traj = 1; traj < parameters.Ntr; traj++){
		for(i = 0; i < gsystem.N; i++){
			itot = traj*gsystem.N + i;
			harmonicData.h_harmonicCount[itot] = harmonicData.h_harmonicCount[i];
			harmonicData.h_harmonicBondCount[itot] = harmonicData.h_harmonicBondCount[i];
			for(h = 0; h < harmonicData.maxHarmonicPerAtom; h++){
				harmonicData.h_harmonic[h*gsystem.widthTot + itot].j =
						harmonicData.h_harmonic[h*gsystem.widthTot + i].j + traj*gsystem.N;
				harmonicData.h_harmonic[h*gsystem.widthTot + itot].kb =
							harmonicData.h_harmonic[h*gsystem.widthTot + i].kb;
				harmonicData.h_harmonic[h*gsystem.widthTot + itot].b0 =
							harmonicData.h_harmonic[h*gsystem.widthTot + i].b0;
			}
		}
	}

	for(i = 0; i < gsystem.N; i++){
		for(j = 0; j < harmonicData.h_harmonicCount[i]; j++){
			if(harmonicData.h_harmonic[j*gsystem.widthTot + i].j > gsystem.N){
				DIE("Harmonic data is inconsistent on bond %d of atom %d (bonded to %d)",
						j, i, harmonicData.h_harmonic[j*gsystem.widthTot + i].j);
			}
		}
	}

	cudaMemcpy(harmonicData.d_harmonicCount, harmonicData.h_harmonicCount,
			gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(harmonicData.d_harmonicBondCount, harmonicData.h_harmonicBondCount,
			gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(harmonicData.d_harmonic, harmonicData.h_harmonic,
			gsystem.widthTot*harmonicData.maxHarmonicPerAtom*sizeof(GHarmonicPair), cudaMemcpyHostToDevice);
	cudaMemcpy(harmonicData.d_harmonicBondEnergy, harmonicData.h_harmonicBondEnergy,
			gsystem.Ntot*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(harmonicData.d_ureyBradleyEnergy, harmonicData.h_ureyBradleyEnergy,
			gsystem.Ntot*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_harmonicData, &harmonicData,
			sizeof(GHarmonicData), 0, cudaMemcpyHostToDevice);
	LOG << "Done initializing harmonic potential.";
}

__global__ void harmonicPotential_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 coord = tex1Dfetch(t_coord, d_i);
		float4 f = c_gsystem.d_forces[d_i];
		float4 r2;
		float r, mult;
		int i;
		GHarmonicPair harmonic;
		for(i = 0; i < c_harmonicData.d_harmonicCount[d_i]; i++){
			harmonic = c_harmonicData.d_harmonic[i*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, harmonic.j);//c_gsystem.d_coord[bond.j];
			r2 -= coord;
			DO_PBC(r2);
			r = sqrtf(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);
			mult = harmonic.kb*(r - harmonic.b0) / r;
			f.x += mult*r2.x;
			f.y += mult*r2.y;
			f.z += mult*r2.z;
		}
		c_gsystem.d_forces[d_i] = f;
	}
}

inline void compute(){
	harmonicPotential_kernel<<<harmonicBlockCount, harmonicBlockSize>>>();
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

__global__ void harmonicPotentialEnergy_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 coord = tex1Dfetch(t_coord, d_i);
		float pot = 0.0f;
		float4 r2;
		int i;
		GHarmonicPair harmonic;
		for(i = 0; i < c_harmonicData.d_harmonicBondCount[d_i]; i++){
			harmonic = c_harmonicData.d_harmonic[i*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, harmonic.j);//c_gsystem.d_coord[bond.j];
			r2 -= coord;
			DO_PBC(r2);
			r2.w = sqrtf(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);
			pot += harmonic.kb*(r2.w - harmonic.b0)*(r2.w - harmonic.b0)*0.5f;
		}
		c_harmonicData.d_harmonicBondEnergy[d_i] = pot;
	}
}

inline void computeHarmonicBondEnergy(){
	harmonicPotentialEnergy_kernel<<<harmonicBlockCount, harmonicBlockSize>>>();
	cudaMemcpy(harmonicData.h_harmonicBondEnergy, harmonicData.d_harmonicBondEnergy,
				gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		float pot = 0.0f;
		for(i = 0; i < gsystem.N; i++){
			pot += harmonicData.h_harmonicBondEnergy[i + traj*gsystem.N];
		}
		harmonicBondEnergyOutput.values[traj] = pot*0.5f;
	}
	checkCUDAError("harmonic energy");
}

__global__ void ureyBradleyPotentialEnergy_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 coord = tex1Dfetch(t_coord, d_i);
		float pot = 0.0f;
		float4 r2;
		int i;
		GHarmonicPair harmonic;
		for(i = c_harmonicData.d_harmonicBondCount[d_i]; i < c_harmonicData.d_harmonicCount[d_i]; i++){
			harmonic = c_harmonicData.d_harmonic[i*c_gsystem.widthTot + d_i];
			r2 =  tex1Dfetch(t_coord, harmonic.j);//c_gsystem.d_coord[bond.j];
			r2 -= coord;
			DO_PBC(r2);
			r2.w = sqrtf(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);
			pot += harmonic.kb*(r2.w - harmonic.b0)*(r2.w - harmonic.b0)*0.5f;
		}
		c_harmonicData.d_ureyBradleyEnergy[d_i] = pot;
	}
}

inline void computeUreyBradleyEnergy(){
	//checkCUDAError("before U-B energy");
	//cudaThreadSynchronize();
	//checkCUDAError("threas synch U-B energy");
	ureyBradleyPotentialEnergy_kernel<<<harmonicBlockCount, harmonicBlockSize>>>();
	//cudaThreadSynchronize();
	//checkCUDAError("before U-B energy memcpy");
	cudaMemcpy(harmonicData.h_ureyBradleyEnergy, harmonicData.d_ureyBradleyEnergy,
				gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		float pot = 0.0f;
		for(i = 0; i < gsystem.N; i++){
			pot += harmonicData.h_ureyBradleyEnergy[i + traj*gsystem.N];
		}
		ureyBradleyEnergyOutput.values[traj] = pot*0.5f;
	}
	checkCUDAError("U-B energy");
}

void destroy(){

}

#undef LOG

} // namespace harmonic_potential
