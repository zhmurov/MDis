/*
 * PDBOutputManager.cu
 *
 *  Created on: Nov 13, 2010
 *      Author: zhmurov
 */

#include "../Core/global.h"
#include "../IO/pdbio.h"
#include "../IO/xyzio.h"
#include "../Util/Log.h"
#include "RestartOutputManager.cuh"
#include "../Util/wrapper.h"

namespace restart_output {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<restart_output> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void (*saveCoordVel)();
int restartFileFormat;

char** coordOutputFilename;
char** velOutputFilename;

void create(){
	char tempFilename[PARAMETER_LENGTH];
	getMaskedParameter(tempFilename, PARAMETER_RESTART_COORD_OUTPUT_FILE);
	restartFileFormat = getFileType(tempFilename);
	if(restartFileFormat == FILETYPE_PDB){
		LOG << "Restart/final coordinates and velocities will be saved in PDB format";
		saveCoordVel = savePDBs;
		updater.update = savePDBs;
	} else
	if(restartFileFormat == FILETYPE_XYZ){
		LOG << "Restart/final coordinates and velocities will be saved in XYZ format";
		saveCoordVel = saveXYZs;
		updater.update = saveXYZs;
	} else {
		DIE("Unknown format for restart coordinates/velocities output.");
	}
	updater.destroy = destroy;
	updater.frequency = getIntegerParameter(PARAMETER_RESTARTOUTPUT_FREQ);
	sprintf(updater.name, "Restart output manager");
	updaters[updatersCount] = &updater;
	updatersCount ++;
	init();
}

void writeRestartKey() {
    char keyFilename[PARAMETER_LENGTH];
    getMaskedParameter(keyFilename, PARAMETER_RESTART_KEY, PARAMETER_STRING_UNDEFINED);
    if (! strncmp(keyFilename,PARAMETER_STRING_UNDEFINED,PARAMETER_LENGTH)) return;
    LOG << "Saving restart key to " << keyFilename;
    FILE *f = safe_fopen(keyFilename, "w");
	fprintf(f, "%d %lld\n", getIntegerParameter(PARAMETER_RESTART_COUNT), step);
	for (int traj = 0; traj < parameters.Ntr; traj++) {
		fprintf(f, "%.6f ", trajectoryTime[traj]);
	}
	fprintf(f, "\n");
	for (int i = 0; i < restartersCount; i++) {
		fprintf(f, "\n%s\n", restarters[i]->name);
		restarters[i]->save(f);
	}
	fclose(f);
}


void init(){
	LOG << "Intializing restart/final coordinates and velocities output...";
	coordOutputFilename = (char**)calloc(parameters.Ntr, sizeof(char*));
	velOutputFilename = (char**)calloc(parameters.Ntr, sizeof(char*));
	int traj;
	char trajnum[10];
	for(traj = 0; traj < parameters.Ntr; traj++){
		sprintf(trajnum, "%d", traj + parameters.firstrun);
		coordOutputFilename[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
		velOutputFilename[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
		getMaskedParameterWithReplacement(coordOutputFilename[traj], PARAMETER_RESTART_COORD_OUTPUT_FILE,
				trajnum, "<run>");
		getMaskedParameterWithReplacement(velOutputFilename[traj], PARAMETER_RESTART_VEL_OUTPUT_FILE,
						trajnum, "<run>");
	}

	int i;
	//if(restartFileFormat == FILETYPE_PDB){ // Does not work if filetypes are different for restart and final coordinates/velocities
		pdbOutputData.atomCount = gsystem.N;
		pdbOutputData.ssCount = 0;
		allocateCPU((void**)&pdbOutputData.atoms, gsystem.N*sizeof(PDBAtom));
		for(i = 0; i < gsystem.N; i++){
			pdbOutputData.atoms[i].id = i+1;
			sprintf(pdbOutputData.atoms[i].name, "%s", topology.atoms[i].name);
			pdbOutputData.atoms[i].resid = topology.atoms[i].resid;
			sprintf(pdbOutputData.atoms[i].resName, "%s", topology.atoms[i].resName);
			pdbOutputData.atoms[i].chain = ' ';
			pdbOutputData.atoms[i].altLoc = ' ';
			pdbOutputData.atoms[i].beta = 0.0;
			pdbOutputData.atoms[i].occupancy = 0.0;
			pdbOutputData.atoms[i].x = 0.0;
			pdbOutputData.atoms[i].y = 0.0;
			pdbOutputData.atoms[i].z = 0.0;
		}
	//} else
	//if(restartFileFormat == FILETYPE_XYZ){
		xyzOutputData.atomCount = gsystem.N;
		allocateCPU((void**)&xyzOutputData.atoms, gsystem.N*sizeof(XYZAtom));
		for(i = 0; i < gsystem.N; i++){
			xyzOutputData.atoms[i].name = topology.atoms[i].name[0];
		}
	//}
	LOG << "Done intializing restart/final coordinates and velocities output.";

}

void savePDBs(){
	LOG << "Saving coordinates/velocities pdbs at step " << step;
	checkCUDAError("PDB output manager: before memcpy");
	copyCoordinatesFromGPU();
	copyVelocitiesFromGPU();
	checkCUDAError("PDB output manager: memcpy");
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj ++){
		for(i = 0; i < gsystem.N; i++){
			pdbOutputData.atoms[i].x = gsystem.h_coord[i + traj*gsystem.N].x*10.0f;
			pdbOutputData.atoms[i].y = gsystem.h_coord[i + traj*gsystem.N].y*10.0f;
			pdbOutputData.atoms[i].z = gsystem.h_coord[i + traj*gsystem.N].z*10.0f;
		}
		/*printf("%f\t%f\t%f\n",
				gsystem.h_coord[traj*gsystem.N].x,
				gsystem.h_coord[traj*gsystem.N].y,
				gsystem.h_coord[traj*gsystem.N].z);*/
		writePDB(coordOutputFilename[traj], &pdbOutputData);
		for(i = 0; i < gsystem.N; i++){
			pdbOutputData.atoms[i].x = gsystem.h_vel[i + traj*gsystem.N].x*10.0f;
			pdbOutputData.atoms[i].y = gsystem.h_vel[i + traj*gsystem.N].y*10.0f;
			pdbOutputData.atoms[i].z = gsystem.h_vel[i + traj*gsystem.N].z*10.0f;
		}
		writePDB(velOutputFilename[traj], &pdbOutputData);
	}
    writeRestartKey();
}

void saveXYZs(){
	LOG << "Saving coordinates/velocities xyzs at step " << step;
	checkCUDAError("XYZ output manager: before memcpy");
	copyCoordinatesFromGPU();
	copyVelocitiesFromGPU();
	checkCUDAError("XYZ output manager: memcpy");
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj ++){
		for(i = 0; i < gsystem.N; i++){
			xyzOutputData.atoms[i].x = gsystem.h_coord[i + traj*gsystem.N].x*10.0f;
			xyzOutputData.atoms[i].y = gsystem.h_coord[i + traj*gsystem.N].y*10.0f;
			xyzOutputData.atoms[i].z = gsystem.h_coord[i + traj*gsystem.N].z*10.0f;
		}
		/*printf("%f\t%f\t%f\n",
				gsystem.h_coord[traj*gsystem.N].x,
				gsystem.h_coord[traj*gsystem.N].y,
				gsystem.h_coord[traj*gsystem.N].z);*/
		writeXYZ(coordOutputFilename[traj], &xyzOutputData);
		for(i = 0; i < gsystem.N; i++){
			xyzOutputData.atoms[i].x = gsystem.h_vel[i + traj*gsystem.N].x*10.0f;
			xyzOutputData.atoms[i].y = gsystem.h_vel[i + traj*gsystem.N].y*10.0f;
			xyzOutputData.atoms[i].z = gsystem.h_vel[i + traj*gsystem.N].z*10.0f;
		}
		writeXYZ(velOutputFilename[traj], &xyzOutputData);
	}
    writeRestartKey();
}

void destroy(){
	LOG << "Saving final velocities/coordinates.";
	int traj;
	char trajnum[10];
	for(traj = 0; traj < parameters.Ntr; traj++){
		sprintf(trajnum, "%d", traj + parameters.firstrun);
		coordOutputFilename[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
		velOutputFilename[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
		getMaskedParameterWithReplacement(coordOutputFilename[traj], PARAMETER_FINAL_COORD_OUTPUT_FILE,
				trajnum, "<run>");
		getMaskedParameterWithReplacement(velOutputFilename[traj], PARAMETER_FINAL_VEL_OUTPUT_FILE,
						trajnum, "<run>");
	}
    int coord_file_format = getFileType(coordOutputFilename[0]); // FIX: I guess, MDis would fail long before if we have < 1 trajectory. Though, I'm still a bad person
    if (coord_file_format == FILETYPE_PDB)
        savePDBs(); 
    else if (coord_file_format == FILETYPE_XYZ) 
        saveXYZs();
    else
    	saveCoordVel(); // Use the same format as used for restarts
	LOG << "Done saving final velocities/coordinates.";
}

#undef LOG

} // namespace restart_output
