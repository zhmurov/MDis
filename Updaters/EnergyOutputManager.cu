/*
 * EnergyOutputManager.cu
 *
 *  Created on: Aug 2, 2010
 *      Author: zhmurov
 */

#include "../Core/global.h"
#include "../Util/Log.h"
#include "EnergyOutputManager.cuh"

energy_output::EnergyOutputData energyOutputData;
int constrCount = 0;//needed for proper temperature calculation

namespace energy_output {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<energy_output> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

char** energyOutputFilename;

void create(){
	updater.update = save;
	updater.destroy = destroy;
	updater.frequency = getIntegerParameter(PARAMETER_ENERGYOUTPUT_FREQ);
	sprintf(updater.name, "Energy output manager");
	updaters[updatersCount] = &updater;
	updatersCount ++;
	init();
}


void init(){
	energyOutputFilename = (char**)calloc(parameters.Ntr, sizeof(char*));
	int traj;
	char trajnum[10];
	for(traj = 0; traj < parameters.Ntr; traj++){
		sprintf(trajnum, "%d", traj + parameters.firstrun);
		energyOutputFilename[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
		getMaskedParameterWithReplacement(energyOutputFilename[traj], PARAMETER_ENERGYOUTPUT_FILE,
				trajnum, "<run>");
		energyOutputFile = safe_fopen(energyOutputFilename[traj], "w");
		fclose(energyOutputFile);
	}
	energyOutputData.T = (float*)calloc(parameters.Ntr, sizeof(float));
	energyOutputData.E = (float*)calloc(parameters.Ntr, sizeof(float));
	char legendFilename[PARAMETER_LENGTH];
	getMaskedParameter(legendFilename, PARAMETER_ENERGYOUTPUT_LEGEND_FILE, "none");
	int i;
	if(strcmp(legendFilename, "none") != 0){
		FILE* legendFile = safe_fopen(legendFilename, "w");
		fprintf(legendFile, "%s\n%s\nMD/LD\n%s\n",
				ENERGY_OUTPUT_NAME_TIMESTEP,
				ENERGY_OUTPUT_NAME_TIME,
				ENERGY_OUTPUT_NAME_TEMPERATURE);
		for(i = 0; i < energyOutputsCount; i++){
			fprintf(legendFile, "%s\n", energyOutputs[i]->name);
		}
		fclose(legendFile);
	}
}

void saveRigidBody(int traj){
	energyOutputFile = safe_fopen(energyOutputFilename[traj], "a");
	fprintf(energyOutputFile, "%lld\t%f\tLD\t", step, trajectoryTime[traj]);
	fprintf(energyOutputFile, "0\t");
	for(int i = 0; i < energyOutputsCount; i++){
		fprintf(energyOutputFile, "0\t");
	}
	fprintf(energyOutputFile, "\n");
	fclose(energyOutputFile);
}

void save(){
	copyVelocitiesFromGPU();
	int i, itot, traj;
	int naned = 0;
	int lookForNaNs = 1;
	if(lookForNaNs){
		copyCoordinatesFromGPU(0);
	}
	for(traj = 0; traj < parameters.Ntr; traj++){
		energyOutputData.T[traj] = 0.0f;
		for(i = 0; i < gsystem.N; i++){
			itot = traj*gsystem.N + i;
			energyOutputData.T[traj] += gsystem.h_vel[itot].w*topology.atoms[i].mass;
			float velocity = sqrtf(gsystem.h_vel[itot].w)/((float)updater.frequency);
			if(velocity > 100.0f){
				printf("WARNING! Atom velocity is > 100 on atom %d (%s%d) in trajectory #%d\n",
						topology.atoms[i].id,
						topology.atoms[i].resName,
						topology.atoms[i].resid,
						parameters.firstrun + traj);
			}
			if(isnan(velocity)){
				printf("ERROR: Atom velocity is NaN on atom %d (%s%d) in trajectory #%d\n",
						topology.atoms[i].id,
						topology.atoms[i].resName,
						topology.atoms[i].resid,
						parameters.firstrun + traj);
				naned = 1;
			}
			if(lookForNaNs){
					if(isnan(gsystem.h_coord[itot].x) || isnan(gsystem.h_coord[itot].y) || isnan(gsystem.h_coord[itot].z)){
						printf("ERROR: Atom coordinate is NaN on atom %d (%s%d) in trajectory #%d\n",
									topology.atoms[i].id,
									topology.atoms[i].resName,
									topology.atoms[i].resid,
									parameters.firstrun + traj);
						naned = 1;
					}
			}
			gsystem.h_vel[itot].w = 0.0f;
		}
	}
	if(naned){
		exit(0);
	}
	printf("%*s%-*s%*s%*s%*s",
			ENERGY_OUTPUT_WIDTH, "",
			ENERGY_OUTPUT_WIDTH, ENERGY_OUTPUT_NAME_TRAJECTORY,
			ENERGY_OUTPUT_WIDTH, ENERGY_OUTPUT_NAME_TIMESTEP,
			ENERGY_OUTPUT_WIDTH, ENERGY_OUTPUT_NAME_TIME,
			ENERGY_OUTPUT_WIDTH, ENERGY_OUTPUT_NAME_TEMPERATURE);
	for(i = 0; i < energyOutputsCount; i++){
		printf("%*s", ENERGY_OUTPUT_WIDTH, energyOutputs[i]->name);
		energyOutputs[i]->computeValues();
	}
	printf("%*s", ENERGY_OUTPUT_WIDTH, ENERGY_OUTPUT_NAME_POTENTIAL);
	printf("\n");

	for(traj = 0; traj < parameters.Ntr; traj++){
		energyOutputFile = safe_fopen(energyOutputFilename[traj], "a");
		printf("%-*s#%-*d%*lld%*f",
				ENERGY_OUTPUT_WIDTH, ENERGY_OUTPUT_TITLE,
				ENERGY_OUTPUT_WIDTH-1, traj + parameters.firstrun,
				ENERGY_OUTPUT_WIDTH, step,
				ENERGY_OUTPUT_WIDTH, trajectoryTime[traj]);
		fprintf(energyOutputFile, "%lld\t%f\tMD\t", step, trajectoryTime[traj]);
		energyOutputData.T[traj] /= ((float)gsystem.N*3.0f-constrCount)*((float)updater.frequency)*Kb_MD;
		printf("%*f", ENERGY_OUTPUT_WIDTH, energyOutputData.T[traj]);
		fprintf(energyOutputFile, "%f\t", energyOutputData.T[traj]);
		energyOutputData.E[traj] = 0.0f;
		for(i = 0; i < energyOutputsCount; i++){
			energyOutputData.E[traj] += energyOutputs[i]->values[traj];
			printf("%*f", ENERGY_OUTPUT_WIDTH, energyOutputs[i]->values[traj]*KCALL_PER_KJ);
			fprintf(energyOutputFile, "%f\t", energyOutputs[i]->values[traj]*KCALL_PER_KJ);
		}
		printf("%*f", ENERGY_OUTPUT_WIDTH, energyOutputData.E[traj]*KCALL_PER_KJ);
		printf("\n");
		fprintf(energyOutputFile, "%f", energyOutputData.E[traj]*KCALL_PER_KJ);
		fprintf(energyOutputFile, "\n");

		fclose(energyOutputFile);
	}
	printTime(step-firststep);
	printEstimatedTimeleft(((float)(step-firststep))/((float)(parameters.numsteps-firststep)));
	cudaMemcpy(gsystem.d_vel, gsystem.h_vel, gsystem.Ntot*sizeof(float4), cudaMemcpyHostToDevice);
}

void destroy(){
}

#undef LOG

} // namespace energy_output
