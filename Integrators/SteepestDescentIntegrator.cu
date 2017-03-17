/*
 * SteepestDescentIntegrator.cu
 *
 *  Created on: Nov 13, 2010
 *      Author: zhmurov
 */


#include "SteepestDescentIntegrator.cuh"
#include "../Util/mdunitsconverter.h"
#include "../Util/Log.h"
#include "../Core/global.h"

namespace sd_integrator {

class Log: public ILog {
    virtual void Write(const char* message) const {
        std::cout << makeTimePrefix() << "<sd_integrator> " << message << std::endl;
    }
} log;

#define LOG LogStream(log)

void create(){
	steepestDescentIntegrator.h = getFloatParameter(PARAMETER_TIMESTEP);
	steepestDescentIntegrator.integrate = &computeSteepestDescentIntegrator;
	steepestDescentIntegrator.finalize = &finalizeSteepestDescentIntegrator;
	sprintf(steepestDescentIntegrator.name, "Steepest Descent integrator");
	integrator = &steepestDescentIntegrator;
	init();

	if(getYesNoParameter(PARAMETER_MIN_CHANGE_TIMESTEP, 0)){
		steepestDescentTimestepUpdater.update = updateTimestep;
		steepestDescentTimestepUpdater.destroy = destroyUpdater;
		steepestDescentTimestepUpdater.frequency = getIntegerParameter(PARAMETER_MIN_UPDATE_TIMESTEP_FREQ, 1000);
		sprintf(steepestDescentTimestepUpdater.name, "SD timestep updater");
		updaters[updatersCount] = &steepestDescentTimestepUpdater;
		updatersCount ++;
	}
}

void init(){
	blockCount = gsystem.Ntot/BLOCK_SIZE + 1;
	blockSize = BLOCK_SIZE;
	steepestDescent.h = steepestDescentIntegrator.h;
	steepestDescent.maxForce = getFloatParameter(PARAMETER_MIN_MAXFORCE, DEFAULT_MIN_MAXFORCE);
	if(getYesNoParameter(PARAMETER_MIN_CHANGE_TIMESTEP, 0)){
    	steepestDescent.timestepChangeVel = getFloatParameter(PARAMETER_MIN_CHANGE_TIMESTEP_AT);
    	steepestDescent.timestepFactor = getFloatParameter(PARAMETER_MIN_TIMESTEP_FACTOR);
    	steepestDescent.finalTimestep = getFloatParameter(PARAMETER_MIN_FINAL_TIMESTEP);
    }
	cudaMemcpyToSymbol(c_steepestDescent, &steepestDescent, sizeof(SteepestDescent), 0, cudaMemcpyHostToDevice);
	LOG << "Initialized";
}



__global__ void ComputeSteepestDescent_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 midcoord = c_gsystem.d_midcoord[d_i];
		float4 vel = c_gsystem.d_vel[d_i];
		float4 f = c_gsystem.d_forces[d_i];
		float mult = c_steepestDescent.maxForce;
		f.w = sqrt(f.x*f.x + f.y*f.y + f.z*f.z);
		if(f.w > mult){
			mult /= f.w;
			f.x *= mult;
			f.y *= mult;
			f.z *= mult;
		}


		f.w = tex1Dfetch(t_m, (int)midcoord.w);
		mult = c_steepestDescent.h/f.w;

		vel.x = f.x*mult;
		vel.y = f.y*mult;
		vel.z = f.z*mult;

		vel.w += vel.x*vel.x + vel.y*vel.y + vel.z*vel.z;

		midcoord.x = vel.x*c_steepestDescent.h;
		midcoord.y = vel.y*c_steepestDescent.h;
		midcoord.z = vel.z*c_steepestDescent.h;
        //DO_PBC(coord);

		c_gsystem.d_midcoord[d_i] = midcoord;
		c_gsystem.d_vel[d_i] = vel;
		f.x = 0.0f;
		f.y = 0.0f;
		f.z = 0.0f;
		c_gsystem.d_forces[d_i] = f;
	}
}

__global__ void FinalizeSteepestDescent_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 coord    = c_gsystem.d_coord[d_i];
		float4 midcoord = c_gsystem.d_midcoord[d_i];
		float4 vel      = c_gsystem.d_vel[d_i];

		vel.x = midcoord.x/c_steepestDescent.h;
		vel.y = midcoord.y/c_steepestDescent.h;
		vel.z = midcoord.z/c_steepestDescent.h;

		coord.x += midcoord.x;
		coord.y += midcoord.y;
		coord.z += midcoord.z;
        DO_PBC(coord);

		c_gsystem.d_coord[d_i] = coord;
		c_gsystem.d_vel[d_i] = vel;

	}
}

void inline computeSteepestDescentIntegrator(){
	ComputeSteepestDescent_kernel<<<blockCount, blockSize>>>();
}

void inline finalizeSteepestDescentIntegrator(){
	FinalizeSteepestDescent_kernel<<<blockCount, blockSize>>>();
}

void inline updateTimestep(){
	copyVelocitiesFromGPU();
	if(steepestDescentIntegrator.h < steepestDescent.finalTimestep){
		int i;
		int doChange = 1;
		for(i = 0; i < gsystem.Ntot; i++){
			if(gsystem.h_vel[i].w > steepestDescent.timestepChangeVel){
				doChange = 0;
				break;
			}
		}
		if(doChange){
			steepestDescent.h *= steepestDescent.timestepFactor;
			cudaMemcpyToSymbol(c_steepestDescent, &steepestDescent, sizeof(SteepestDescent), 0, cudaMemcpyHostToDevice);
		}
	}
	if(steepestDescentIntegrator.h > steepestDescent.finalTimestep){
		steepestDescentIntegrator.h = steepestDescent.finalTimestep;
		cudaMemcpyToSymbol(c_steepestDescent, &steepestDescent, sizeof(SteepestDescent), 0, cudaMemcpyHostToDevice);
	}
}


void destroyUpdater(){

}

#undef LOG

} // namespace sd_integrator
