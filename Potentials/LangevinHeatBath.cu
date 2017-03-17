/*
 * LangevinHeatBath.cu
 *
 *  Created on: Nov 14, 2010
 *      Author: zhmurov
 */

#include "../Core/global.h"
#include "../Util/Log.h"
#include "LangevinHeatBath.cuh"
#include "HybridTaus.cu"

namespace langevin_heat_bath {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<langevin_heat_bath> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void create(){
	if(getYesNoParameter(PARAMETER_REMD_ENABLED, 0)) {
		LOG << "REMD will be used instead of Langevin heat bath";
		return;
	}
	LOG << "Initializing Langevin heat bath";
	if(getYesNoParameter(PARAMETER_HEATING, DEFAULT_HEATING)){
		LOG << "Heating simulations requested.";
		updater.update = updateTemperature;
		updater.destroy = destroyUpdater;
		updater.frequency = getIntegerParameter(PARAMETER_TEMPERATURE_UPDATE_FREQ);
		sprintf(updater.name, "Temperature control");
		updaters[updatersCount] = &updater;
		updatersCount ++;
		initTemperatureControl();
	} else {
		initConstantTemperature();
	}
	potential.compute = &compute;
	potential.destroy = &destroy;
	sprintf(potential.name, "Langevin heatbath");
	potentials[potentialsCount] = &potential;
	potentialsCount ++;
	init();
	LOG << "Done initializing Langevin heat bath";
}

void init(){
	langevinHeatBathBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;
	langevinHeatBathBlockSize = BLOCK_SIZE;
	hybrid_taus::initRand(parameters.rseed, gsystem.Ntot);
}

__global__ void langevinHeatBath_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 f = c_gsystem.d_forces[d_i];
		float4 v = c_gsystem.d_vel[d_i];
		int at = c_gsystem.d_atomTypes[d_i];
		f.w = tex1Dfetch(t_m, at);
		float mult = c_var*sqrtf(f.w);
		float4 rf = hybrid_taus::rforce(d_i);
		float mgamma = f.w*c_gamma;
		f.x += mult*rf.x - mgamma*v.x;
		f.y += mult*rf.y - mgamma*v.y;
		f.z += mult*rf.z - mgamma*v.z;
		c_gsystem.d_forces[d_i] = f;
	}
}

void inline compute(){
	langevinHeatBath_kernel<<<langevinHeatBathBlockCount, langevinHeatBathBlockSize>>>();
}

void destroy(){

}

void initConstantTemperature(){
	float T = getFloatParameter(PARAMETER_TEMPERATURE, DEFAULT_TEMPERATURE);
	float h = getFloatParameter(PARAMETER_TIMESTEP);
	temperatureControl.gamma = getFloatParameter(PARAMETER_DAMPING, DEFAULT_DAMPING);
	LOG << "Simulations will be held at constant temperature of " << T << "K";
	temperatureControl.var = sqrtf(2.0f*temperatureControl.gamma*Kb_MD*T/h);
	cudaMemcpyToSymbol(c_var, &(temperatureControl.var), sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_gamma, &(temperatureControl.gamma), sizeof(float), 0, cudaMemcpyHostToDevice);
}

void initTemperatureControl(){
	temperatureControl.initialT = getFloatParameter(PARAMETER_INITIAL_TEMPERATURE, DEFAULT_INITIAL_TEMPERATURE);
	temperatureControl.finalT = getFloatParameter(PARAMETER_FINAL_TEMPERATURE, DEFAULT_FINAL_TEMPERATURE);
	temperatureControl.deltaT = getFloatParameter(PARAMETER_TEMPERATURE_INCREMENT, 0.0);
	if(temperatureControl.deltaT == 0.0){
			temperatureControl.deltaT = (temperatureControl.finalT - temperatureControl.initialT)/
					((parameters.totalTime-parameters.initialTime)/(updater.frequency*getFloatParameter(PARAMETER_TIMESTEP)));
		}
	LOG << "Temperature will be increased by " << temperatureControl.deltaT << " every " << updater.frequency << " steps starting from " << temperatureControl.initialT << "K untill value of " << temperatureControl.finalT << "K is reached.";
	temperatureControl.T = temperatureControl.initialT;
	float h = getFloatParameter(PARAMETER_TIMESTEP);
	temperatureControl.gamma = getFloatParameter(PARAMETER_DAMPING, DEFAULT_DAMPING);
	temperatureControl.mult = sqrtf(2.0f*temperatureControl.gamma*Kb_MD/h);
	cudaMemcpyToSymbol(c_gamma, &(temperatureControl.gamma), sizeof(float), 0, cudaMemcpyHostToDevice);
}

void updateTemperature(){
	if( (temperatureControl.T < temperatureControl.finalT && temperatureControl.deltaT > 0) || 
	    (temperatureControl.T > temperatureControl.finalT && temperatureControl.deltaT < 0)){
		temperatureControl.T = temperatureControl.initialT + temperatureControl.deltaT*(step/updater.frequency);
		temperatureControl.var = sqrtf(temperatureControl.T)*temperatureControl.mult;
		LOG << "Updating temperature to " << temperatureControl.T << " at step "<< step;
		cudaMemcpyToSymbol(c_var, &(temperatureControl.var), sizeof(float), 0, cudaMemcpyHostToDevice);
	}
}

void destroyUpdater(){

}

#undef LOG

} // namespace langevin_heat_bath
