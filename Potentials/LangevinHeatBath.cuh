/*
 * LangevinHeatBath.cuh
 *
 *  Created on: Nov 14, 2010
 *      Author: zhmurov
 */

#pragma once

namespace langevin_heat_bath {

int langevinHeatBathBlockSize;
int langevinHeatBathBlockCount;

typedef struct {
	float initialT;
	float finalT;
	float deltaT;
	float T;
	float mult;
	float var;
	float gamma;
} TemperatureControl;

TemperatureControl temperatureControl;
__device__ __constant__ float c_var;
__device__ __constant__ float c_gamma;

Potential potential;
Updater updater;

void init();
void inline compute();
void destroy();
void initConstantTemperature();
void initTemperatureControl();
void inline updateTemperature();
void destroyUpdater();

} // namespace langevin_heat_bath
