/*
 * parameters.c
 *
 *  Created on: Jul 5, 2009
 *      Author: zhmurov
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sstream>
#include "parameters.h"
#include "../IO/configreader.h"
#include "../Util/wrapper.h"

extern double* trajectoryTime;

void parseParametersFileNAMD(char* filename, Parameters* parameters){
	parameters->device = getIntegerParameter(PARAMETER_MPI_DEV_CUR, getIntegerParameter(PARAMETER_DEVICE, -1));

	parameters->firststep = getLongIntegerParameter(PARAMETER_FIRSTSTEP, 0);
	parameters->numsteps = getLongIntegerParameter(PARAMETER_NUMSTEPS, -1);
	parameters->totalTime = getFloatParameter(PARAMETER_TOTALTIME, -1.0f);
	float timeStep = getFloatParameter(PARAMETER_TIMESTEP);
	parameters->initialTime = ((float)parameters->firststep)*timeStep;
	if(parameters->numsteps == -1 && parameters->totalTime == -1.0f){
		DIE("Either '%s' or '%s' should be specified.",
				PARAMETER_NUMSTEPS,
				PARAMETER_TOTALTIME);
	} else
	if(parameters->totalTime == -1.0f){
		parameters->totalTime = ((float)parameters->numsteps)*timeStep;
	} else
	if(parameters->numsteps == -1){
		parameters->numsteps = (int)(parameters->totalTime/timeStep);
	} else {
		printf("Both '%s' and '%s' are specified. '%s' will be used to redefine %s'.\n",
				PARAMETER_NUMSTEPS,	PARAMETER_TOTALTIME,
				PARAMETER_TOTALTIME, PARAMETER_NUMSTEPS);
		parameters->numsteps = (int)(parameters->totalTime/timeStep);
		printf("New value for '%s' is %lld'.\n", PARAMETER_NUMSTEPS, parameters->numsteps);
	}
	printf("Total simulation time is %fps.\n", parameters->totalTime);
	int run = getIntegerParameter(PARAMETER_RUN, -1);
	if(run == -1){
		parameters->firstrun = getIntegerParameter(PARAMETER_FIRSTRUN, -1);
		if(parameters->firstrun == -1){
			DIE("Either '%s' or '%s' and '%s' should be specified in configuration file.",
					PARAMETER_RUN,
					PARAMETER_FIRSTRUN,
					PARAMETER_RUNNUM);
		}
		parameters->Ntr = getIntegerParameter(PARAMETER_RUNNUM);
		printf("Will run %d trajectories, starting from trajectory %d\n.",
				parameters->Ntr,
				parameters->firstrun);
	} else {
		parameters->firstrun = run;
		parameters->Ntr = 1;
		printf("Will run one trajectory (#%d).", parameters->firstrun);
	}

	trajectoryTime = (double*)calloc(parameters->Ntr, sizeof(double));

	parameters->rseed = getIntegerParameter(PARAMETER_RSEED, (int)time(NULL))*parameters->firstrun;
#ifdef USE_MPI
    parameters->rseed *= getIntegerParameter(PARAMETER_MPI_RANK) + 1;
#endif

	getMaskedParameter(parameters->topologyFilename, PARAMETER_TOPOLOGY_FILE);
	getMaskedParameter(parameters->coordFilename, PARAMETER_COORD_FILE);
	getMaskedParameter(parameters->velFilename, PARAMETER_VEL_FILE, PARAMETER_STRING_UNDEFINED);
	getMaskedParameter(parameters->forceFieldFilename, PARAMETER_FF_FILE);
	getMaskedParameter(parameters->topologiesFilename, PARAMETER_TOP_FILE, PARAMETER_STRING_UNDEFINED);
	getMaskedParameter(parameters->parametersType, PARAMETER_FF_TYPE);
	getMaskedParameter(parameters->restartKey, PARAMETER_RESTART_KEY, PARAMETER_STRING_UNDEFINED);

}

FILE* parametersPrepareRestart(Parameters* parameters) {
    if (strlen(parameters->restartKey) != 0) {
        FILE *keyf = safe_fopen(parameters->restartKey, "r");
        int restartcount = 0;
		if (fscanf(keyf,"%d %lld",&restartcount,&parameters->firststep) != 2)
            DIE("Reading restartkey %s: unable to get restartcout and firststep", parameters->restartKey);
		for (int traj = 0; traj < parameters->Ntr; traj++) {
			if (fscanf(keyf, "%lf", &trajectoryTime[traj]) != 1)
                DIE("Reading restartkey %s: unable to read trajectory's #%d time", parameters->restartKey, traj);
		}
        printf("Initiating %d-th restart, from step %lld\n",restartcount+1, parameters->firststep);
    	getMaskedParameterWithReplacementT(parameters->coordFilename, PARAMETER_RESTART_COORD_OUTPUT_FILE, restartcount, "<restartcount>");
    	getMaskedParameterWithReplacementT(parameters->velFilename  , PARAMETER_RESTART_VEL_OUTPUT_FILE  , restartcount, "<restartcount>");
//        parameters->restartcount ++;
        addParameterT(PARAMETER_RESTART_COUNT, restartcount+1);
        setParameterT(PARAMETER_FIRSTSTEP,  parameters->firststep, 1);
        return keyf;
    } else {
         DIE("Trying to restart simulation, but restartkey is not specified!");
    }
}

char *parameter_mpi_device(int i) {
    std::stringstream s;
    s << "mpi_device_" << i;
    char *a;
    a = new char[s.str().size()+1];
    strcpy(a,s.str().c_str());
    return a;
}


