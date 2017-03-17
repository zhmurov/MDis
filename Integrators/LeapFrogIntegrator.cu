/*
 * LeapFrogIntegrator.cu
 *
 *  Created on: Jul 29, 2010
 *      Author: zhmurov
 */
#include "LeapFrogIntegrator.cuh"
#include "../Util/mdunitsconverter.h"
#include "../Util/Log.h"
#include "../Core/global.h"

namespace leapfrog_integrator {

class Log: public ILog {
    virtual void Write(const char* message) const {
        std::cout << makeTimePrefix() << "<leapfrog_integrator> " << message << std::endl;
    }
} log;

#define LOG LogStream(log)

void create(){
	leapFrogIntegrator.h = getFloatParameter(PARAMETER_TIMESTEP);
	leapFrogIntegrator.integrate = &computeLeapFrogIntegrator;
	leapFrogIntegrator.finalize = &finalizeLeapFrogIntegrator;
	sprintf(leapFrogIntegrator.name, "Leap Frog integrator");
	integrator = &leapFrogIntegrator;
	init();
}

void init(){
	blockCount = gsystem.Ntot/BLOCK_SIZE + 1;
	blockSize = BLOCK_SIZE;
	leapFrog.h = leapFrogIntegrator.h;
	leapFrog.gamma = getFloatParameter(PARAMETER_DAMPING, DEFAULT_DAMPING);
	//initRand(getLongIntegerParameter(PARAMETER_RSEED), gsystem.Ntot);
	cudaMemcpyToSymbol(c_leapFrog, &leapFrog, sizeof(LeapFrog), 0, cudaMemcpyHostToDevice);
	//LOG << "Initialized with timestep " << leapFrog.h << " ps and damping " << leapFrog.gamma << "1/ps";
	LOG << "Initialized with timestep " << leapFrog.h << " ps";
}


__global__ void computeLeapFrog_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 midcoord = c_gsystem.d_midcoord[d_i];
		float4 vel 		= c_gsystem.d_vel[d_i];
		float4 f 		= c_gsystem.d_forces[d_i];
		f.w = tex1Dfetch(t_m, (int)midcoord.w);
		float mult = 1.0f*c_leapFrog.h/f.w;

		vel.x += mult*f.x;
		vel.y += mult*f.y;
		vel.z += mult*f.z;

/*		vel.w += vel.x*vel.x + vel.y*vel.y + vel.z*vel.z;

		vel.x += mult*f.x;
		vel.y += mult*f.y;
		vel.z += mult*f.z;
*/
		midcoord.x = vel.x*c_leapFrog.h;
		midcoord.y = vel.y*c_leapFrog.h;
		midcoord.z = vel.z*c_leapFrog.h;
        //DO_PBC(coord);

		c_gsystem.d_midcoord[d_i] = midcoord;
//		c_gsystem.d_vel[d_i] = vel;
		f.x = 0.0f;
		f.y = 0.0f;
		f.z = 0.0f;
		c_gsystem.d_forces[d_i] = f;
	}
}

__global__ void finalizeLeapFrog_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 coord    = c_gsystem.d_coord[d_i];
		float4 midcoord = c_gsystem.d_midcoord[d_i];
		float4 vel_old  = c_gsystem.d_vel[d_i];
		float4 vel 		= vel_old;

		vel.x = midcoord.x/c_leapFrog.h;
		vel.y = midcoord.y/c_leapFrog.h;
		vel.z = midcoord.z/c_leapFrog.h;
//!!
		vel_old.x = 0.5* (vel_old.x + vel.x);
		vel_old.y = 0.5* (vel_old.y + vel.y);
		vel_old.z = 0.5* (vel_old.z + vel.z);
		vel.w    += vel_old.x*vel_old.x + vel_old.y*vel_old.y + vel_old.z*vel_old.z;
//!!
		coord.x += midcoord.x;
		coord.y += midcoord.y;
		coord.z += midcoord.z;
        DO_PBC(coord);

		c_gsystem.d_coord[d_i] = coord;
		c_gsystem.d_vel[d_i] = vel;

	}
}

void inline computeLeapFrogIntegrator(){
	computeLeapFrog_kernel<<<blockCount, blockSize>>>();
}

void inline finalizeLeapFrogIntegrator(){
	finalizeLeapFrog_kernel<<<blockCount, blockSize>>>();
}

#undef LOG

} // namespace leapfrog_integrator
