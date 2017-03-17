/*
 * AnglePotential.cu
 *
 *  Created on: Aug 4, 2010
 *      Author: zhmurov
 */

#include "../Core/global.h"
#include "../Core/md.cuh"
#include "../Util/Log.h"
#include "AnglePotential.cuh"

namespace angle_potential {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<angle_potential> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)


void create(){
	potential.compute = &compute;
	potential.destroy = &destroy;
	sprintf(potential.name, "Angle potential");
	potentials[potentialsCount] = &potential;
	potentialsCount ++;
	energyOutput.computeValues = &computeEnergy;
	allocateCPU((void**)&energyOutput.values, parameters.Ntr*sizeof(float));
	strcpy(energyOutput.name, ENERGY_OUTPUT_NAME_ANGLE);
	energyOutputs[energyOutputsCount] = &energyOutput;
	energyOutputsCount ++;
	init();
}

void init(){
	LOG << "Initializing angle potential...";
	angleData.A = topology.angleCount;
	angleData.Atot = topology.angleCount*parameters.Ntr;
	angleBlockSize = BLOCK_SIZE;
	angleBlockCount = angleData.Atot/BLOCK_SIZE + 1;
	angleSummBlockSize = BLOCK_SIZE;
	angleSummBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;
	if(angleData.Atot > 0){

		allocateCPU((void**)&angleData.h_angleCount, gsystem.Ntot*sizeof(int));
		allocateGPU((void**)&angleData.d_angleCount, gsystem.Ntot*sizeof(int));

		int i, a;

		//Testing
		/*int atomsInBlock = 0;
		int anglesInBlock = 0;
		int anglesInBlockBefore = 0;
		int currentBlock = 0;
		for(i = 0; i < gsystem.Nsim; i++){
			atomsInBlock ++;
			for(a = 0; a < topology.angleCount; a++){
				Angle angle = topology.angles[a];
				if(angle.i == i || angle.j == i || angle.k == i){
					anglesInBlock ++;
				}
			}
			if(anglesInBlock >= angleBlockSize){
				if(anglesInBlock == angleBlockSize){
					printf("Block %d: %d atoms, %d angles.\n", currentBlock, atomsInBlock, anglesInBlock);
				} else {
					printf("Block %d: %d atoms, %d angles.\n", currentBlock, atomsInBlock, anglesInBlockBefore);
					i--;
				}
				anglesInBlock = 0;
				atomsInBlock = 0;
				currentBlock ++;
			} else {
				anglesInBlockBefore = anglesInBlock;
			}

		}*/
		//exit(0);
		//Done Testing

		for(i = 0; i < gsystem.Ntot; i++){
			angleData.h_angleCount[i] = 0;
		}

		for(a = 0; a < angleData.A; a++){
			Angle angle = topology.angles[a];
			angleData.h_angleCount[angle.i]++;
			angleData.h_angleCount[angle.j]++;
			angleData.h_angleCount[angle.k]++;
		}

		angleData.maxAnglesPerAtom = 0;
		for(i = 0; i < gsystem.N; i++){
			if(angleData.h_angleCount[i] > angleData.maxAnglesPerAtom){
				angleData.maxAnglesPerAtom = angleData.h_angleCount[i];
			}
		}
		LOG << "Maximum angles per atom is " << angleData.maxAnglesPerAtom;

		allocateCPU((void**)&angleData.h_angles, angleData.Atot*sizeof(int4));
		allocateGPU((void**)&angleData.d_angles, angleData.Atot*sizeof(int4));
		allocateCPU((void**)&angleData.h_angleRefs, angleData.Atot*sizeof(int4));
		allocateGPU((void**)&angleData.d_angleRefs, angleData.Atot*sizeof(int4));
		allocateCPU((void**)&angleData.h_angleForces, gsystem.widthTot*angleData.maxAnglesPerAtom*sizeof(float4));
		allocateGPU((void**)&angleData.d_angleForces, gsystem.widthTot*angleData.maxAnglesPerAtom*sizeof(float4));

		allocateCPU((void**)&angleData.h_angleTypes, angleTypesCount*sizeof(float2));
		allocateGPU((void**)&angleData.d_angleTypes, angleTypesCount*sizeof(float2));

		allocateCPU((void**)&angleData.h_angleEnergies, angleData.Atot*sizeof(float));
		allocateGPU((void**)&angleData.d_angleEnergies, angleData.Atot*sizeof(float));

		for(i = 0; i < gsystem.Ntot; i++){
			angleData.h_angleCount[i] = 0;
		}

		for(a = 0; a < angleData.A; a++){
			Angle angle = topology.angles[a];
			angleData.h_angles[a].x = angle.i;
			angleData.h_angles[a].y = angle.j;
			angleData.h_angles[a].z = angle.k;
			angleData.h_angles[a].w = angle.type;
			angleData.h_angleRefs[a].x = angleData.h_angleCount[angle.i];
			angleData.h_angleRefs[a].y = angleData.h_angleCount[angle.j];
			angleData.h_angleRefs[a].z = angleData.h_angleCount[angle.k];
			angleData.h_angleCount[angle.i]++;
			angleData.h_angleCount[angle.j]++;
			angleData.h_angleCount[angle.k]++;
		}

		for(i = 0; i < gsystem.N; i++){
			if(angleData.h_angleCount[i] > angleData.maxAnglesPerAtom){
				DIE("Maximum angles per atom exceeded the limit of %d on atom %d",
						angleData.maxAnglesPerAtom, i);
			}
		}

		int traj, atot, itot;
		for(traj = 1; traj < parameters.Ntr; traj++){
			for(a = 0; a < angleData.A; a++){
				atot = angleData.A*traj + a;
				angleData.h_angles[atot].x = angleData.h_angles[a].x + gsystem.N*traj;
				angleData.h_angles[atot].y = angleData.h_angles[a].y + gsystem.N*traj;
				angleData.h_angles[atot].z = angleData.h_angles[a].z + gsystem.N*traj;
				angleData.h_angles[atot].w = angleData.h_angles[a].w;
				angleData.h_angleRefs[atot].x = angleData.h_angleRefs[a].x;
				angleData.h_angleRefs[atot].y = angleData.h_angleRefs[a].y;
				angleData.h_angleRefs[atot].z = angleData.h_angleRefs[a].z;
			}
			for(i = 0; i < gsystem.N; i++){
				itot = gsystem.N*traj + i;
				angleData.h_angleCount[itot] = angleData.h_angleCount[i];
			}
		}

		/*for(a = 0; a < topology.angleCount; a++){
			printf("%d: (%d-%d-%d, %d) \n", a,
				angleData.h_angles[a].x,
				angleData.h_angles[a].y,
				angleData.h_angles[a].z,
				angleData.h_angles[a].w);
		}*/

		for(i = 0; i < angleTypesCount; i++){
			angleData.h_angleTypes[i].x = angleTypes[i].ktheta;
			angleData.h_angleTypes[i].y = angleTypes[i].theta0;
		}

		cudaMemcpy(angleData.d_angleCount, angleData.h_angleCount,
				gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(angleData.d_angles, angleData.h_angles,
				angleData.Atot*sizeof(int4), cudaMemcpyHostToDevice);
		cudaMemcpy(angleData.d_angleRefs, angleData.h_angleRefs,
				angleData.Atot*sizeof(int4), cudaMemcpyHostToDevice);

		cudaMemcpy(angleData.d_angleTypes, angleData.h_angleTypes,
				angleTypesCount*sizeof(int2), cudaMemcpyHostToDevice);

		cudaBindTexture(0, t_angleTypes, angleData.d_angleTypes, angleTypesCount*sizeof(int2));
	}

	cudaMemcpyToSymbol(c_angleData, &angleData,
				sizeof(GAngleData), 0, cudaMemcpyHostToDevice);

	LOG << "Done initializing angle potential.";
}

__global__ void harmonicAnglePotential_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_angleData.Atot){
		int4 angle = c_angleData.d_angles[d_i];
		int4 ref = c_angleData.d_angleRefs[d_i];
		float4 r1 = tex1Dfetch(t_coord, angle.x);
		float4 r2 = tex1Dfetch(t_coord, angle.y);
		float4 r3 = tex1Dfetch(t_coord, angle.z);
		float2 par = tex1Dfetch(t_angleTypes, angle.w);
		float3 dr12, dr32;
		dr12.x = r1.x - r2.x;
		dr12.y = r1.y - r2.y;
		dr12.z = r1.z - r2.z;
		DO_PBC(dr12);
		dr32.x = r3.x - r2.x;
		dr32.y = r3.y - r2.y;
		dr32.z = r3.z - r2.z;
		DO_PBC(dr32);
		float r12inv = 1.0f/sqrtf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z);
		float r32inv = 1.0f/sqrtf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z);
		float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
		if(costheta > 1.0f){
			costheta = 1.0f;
		} else
		if(costheta < -1.0f){
			costheta = -1.0f;
		}
		float sintheta = sqrtf(1.0f - costheta*costheta);
		float theta = acos(costheta);
		float diff = theta - par.y;
		if(sintheta < 1.e-6){
			if(diff < 0){
				diff *= 2.0f*par.x;
			} else {
				diff *= -2.0f*par.x;
			}
		} else {
			diff *= (-2.0f*par.x) / sintheta;
		}
		float c1 = diff*r12inv;
		float c2 = diff*r32inv;

		float4 f1, f2, f3;
		f1.x = c1*(dr12.x*(r12inv*costheta) - dr32.x*r32inv);
		f1.y = c1*(dr12.y*(r12inv*costheta) - dr32.y*r32inv);
		f1.z = c1*(dr12.z*(r12inv*costheta) - dr32.z*r32inv);
		f2 = f1;
		f3.x = c2*(dr32.x*(r32inv*costheta) - dr12.x*r12inv);
		f3.y = c2*(dr32.y*(r32inv*costheta) - dr12.y*r12inv);
		f3.z = c2*(dr32.z*(r32inv*costheta) - dr12.z*r12inv);
		f2.x += f3.x;
		f2.y += f3.y;
		f2.z += f3.z;

		f2.x = -f2.x;
		f2.y = -f2.y;
		f2.z = -f2.z;

		c_angleData.d_angleForces[c_gsystem.widthTot*ref.x + angle.x] = f1;
		c_angleData.d_angleForces[c_gsystem.widthTot*ref.y + angle.y] = f2;
		c_angleData.d_angleForces[c_gsystem.widthTot*ref.z + angle.z] = f3;

	}
}

__global__ void summAngleForces_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 f = c_gsystem.d_forces[d_i];
		float4 df;
		int i;
		for(i = 0; i < c_angleData.d_angleCount[d_i]; i++){
			df = c_angleData.d_angleForces[c_gsystem.widthTot*i + d_i];
			f.x += df.x;
			f.y += df.y;
			f.z += df.z;
		}
		c_gsystem.d_forces[d_i] = f;
	}
}

inline void compute(){
	//checkCUDAError("before angle potential");
	harmonicAnglePotential_kernel<<<angleBlockCount, angleBlockSize>>>();
	//checkCUDAError("before sum angles");
	summAngleForces_kernel<<<angleSummBlockCount, angleSummBlockSize>>>();
	//checkCUDAError("after sum angles");
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
	printf("Net force (angles): (%f, %f, %f) %f\n", force.x, force.y, force.z,
			sqrtf(force.x*force.x + force.y*force.y + force.z*force.z));*/
}

__global__ void harmonicAnglePotentialEnergy_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_angleData.Atot){
		int4 angle = c_angleData.d_angles[d_i];
		float4 r1 = tex1Dfetch(t_coord, angle.x);
		float4 r2 = tex1Dfetch(t_coord, angle.y);
		float4 r3 = tex1Dfetch(t_coord, angle.z);
		float2 par = tex1Dfetch(t_angleTypes, angle.w);
		float3 dr12, dr32;
		dr12.x = r1.x - r2.x;
		dr12.y = r1.y - r2.y;
		dr12.z = r1.z - r2.z;
		DO_PBC(dr12);
		dr32.x = r3.x - r2.x;
		dr32.y = r3.y - r2.y;
		dr32.z = r3.z - r2.z;
		DO_PBC(dr32);
		float r12inv = 1.0f/sqrtf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z);
		float r32inv = 1.0f/sqrtf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z);
		float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
		if(costheta > 1.0f){
			costheta = 1.0f;
		} else
		if(costheta < -1.0f){
			costheta = -1.0f;
		}
		float theta = acos(costheta);
		float diff = theta - par.y;
		c_angleData.d_angleEnergies[d_i] = par.x*diff*diff;
	}
}

inline void computeEnergy(){
	harmonicAnglePotentialEnergy_kernel<<<angleBlockCount, angleBlockSize>>>();
	cudaMemcpy(angleData.h_angleEnergies, angleData.d_angleEnergies,
					angleData.Atot*sizeof(float), cudaMemcpyDeviceToHost);
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		float pot = 0.0f;
		for(i = 0; i < angleData.A; i++){
			pot += angleData.h_angleEnergies[i + traj*angleData.A];
		}
		energyOutput.values[traj] = pot;
	}
	checkCUDAError("angle energy");
}

void destroy(){

}

#undef LOG

} // namespace angle_potential
