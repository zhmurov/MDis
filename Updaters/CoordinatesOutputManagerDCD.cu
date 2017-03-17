/*
 * CoordinatesOutputManagerDCD.cu
 *
 *  Created on: Aug 2, 2010
 *      Author: zhmurov
 */
#include "../Core/global.h"
#include "../IO/dcdio.h"
#include "../Util/Log.h"
#include "CoordinatesOutputManagerDCD.cuh"

namespace coordinates_output_dcd {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<coordinates_output_dcd> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

char** dcdOutputFilename;

void create(){
	updater.update = save;
	updater.destroy = destroy;
	updater.frequency = getIntegerParameter(PARAMETER_DCDOUTPUT_FREQ);
	sprintf(updater.name, "DCD output manager");
	updaters[updatersCount] = &updater;
	updatersCount ++;
	init();
}

void init(){
	LOG << "Initializing...";
	allocateCPU((void**)&dcdOutput.X, gsystem.N*sizeof(float));
	allocateCPU((void**)&dcdOutput.Y, gsystem.N*sizeof(float));
	allocateCPU((void**)&dcdOutput.Z, gsystem.N*sizeof(float));
	dcdOutputFilename = (char**)calloc(parameters.Ntr, sizeof(char*));
	int traj;
	char trajnum[10];
	for(traj = 0; traj < parameters.Ntr; traj++){
		sprintf(trajnum, "%d", traj + parameters.firstrun);
		dcdOutputFilename[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
		getMaskedParameterWithReplacement(dcdOutputFilename[traj], PARAMETER_DCDOUTPUT_FILE,
				trajnum, "<run>");
		dcdOutput.file = dcd_open_write(dcdOutput.file, dcdOutputFilename[traj]);
		dcd_write_header(dcdOutput.file, dcdOutputFilename[traj], gsystem.N,
				parameters.numsteps/updater.frequency + 1,
				1, 1, integrator->h*updater.frequency);
		dcd_close(dcdOutput.file);
	}
	LOG << "Done intializing DCD output";
}

void save(){
	checkCUDAError("DCD output manager: before memcpy");
	copyCoordinatesFromGPU();
	checkCUDAError("DCD output manager: memcpy");
	for(int traj = 0; traj < parameters.Ntr; traj ++){
		saveTrajectory(traj);
	}
}

void saveTrajectory(int traj) {
	for(int i = 0; i < gsystem.N; i++){
		dcdOutput.X[i] = gsystem.h_coord[i + traj*gsystem.N].x*10.0f;
		dcdOutput.Y[i] = gsystem.h_coord[i + traj*gsystem.N].y*10.0f;
		dcdOutput.Z[i] = gsystem.h_coord[i + traj*gsystem.N].z*10.0f;
	}
	/*printf("%f\t%f\t%f\n",
			gsystem.h_coord[traj*gsystem.N].x,
			gsystem.h_coord[traj*gsystem.N].y,
			gsystem.h_coord[traj*gsystem.N].z);*/
	dcdOutput.file = dcd_open_append(dcdOutput.file, dcdOutputFilename[traj]);
	dcd_write_frame(dcdOutput.file, gsystem.N, dcdOutput.X, dcdOutput.Y, dcdOutput.Z);
	dcd_close(dcdOutput.file);
}

void destroy(){
}

#undef LOG

} // namespace coordinates_output_dcd

