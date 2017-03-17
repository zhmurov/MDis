#include "parameters.h"
#include "../Core/global.h"
#include "../Util/wrapper.h"
#include <unistd.h>
#ifdef USE_MPI
#include <mpi.h>
#endif

char* configFile;
bool is_restart = false;

Parameters parameters;
Topology topology;
ForceField ff;

extern void initGPU();
extern void compute();
extern void parametrizeSOP();

int main(int argc, char *argv[]){

	printf("===============================\n");
	printf("MDis version trunk\n");
	printf("(built %s %s)\n", __DATE__, __TIME__);
#ifdef USE_MPI
	printf("[With MPI support]\n");
#endif
#ifdef USE_PBC
	printf("[With PBC support]\n");
#endif
	printf("===============================\n");
    char hostname[1024];
    gethostname(hostname, 1024);
    time_t t_now = time(NULL);
    printf("Running on host '%s', %s\n", hostname, ctime(&t_now));
    printf("Args:");
    for (int i = 0; i < argc; ++i) printf(" %s", argv[i]);
    printf("\n");
	printf("===============================\n");

#ifdef USE_MPI
    MPI::Init(argc, argv);
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();
    const int mpi_size = MPI::COMM_WORLD.Get_size();
#endif

    if(argc < 2){
		die("Parameters file should be specified.\n");
	}
	configFile = argv[1];
	parseParametersFile(configFile, argc, argv);

    // Now we initialize some fake variables...
    if (argc > 2) {
        for (int i = 2; i < argc; ++i)
            is_restart = is_restart || (strncmp(argv[i], "--restart", 9) == 0);
    }

#ifdef USE_MPI
    addParameterT<int>(PARAMETER_MPI_RANK, mpi_rank);
    int mpi_dev_cur;
    if (getYesNoParameter(PARAMETER_MPI_DEV_AUTO, 0)) {
        int mpi_dpn = getIntegerParameter(PARAMETER_MPI_DEVPERNODE);
        mpi_dev_cur = mpi_rank % mpi_dpn;
    } else if (mpi_size > 1) {
        mpi_dev_cur = getIntegerParameter(PARAMETER_MPI_DEVICE(mpi_rank));
    } else {
        mpi_dev_cur = getIntegerParameter(PARAMETER_MPI_DEVICE(mpi_rank),getIntegerParameter(PARAMETER_DEVICE));
    }
    addParameterT<int>(PARAMETER_MPI_DEV_CUR, mpi_dev_cur);
#endif

	parseParametersFileNAMD(configFile, &parameters);

#ifdef USE_MPI
    if (getYesNoParameter(PARAMETER_MPI_FIRSTRUN, 1)) // We adjust the 'firstrun' parameter so the trajectory numbers are different for each thread
        parameters.firstrun = mpi_rank * parameters.Ntr + parameters.firstrun;
#endif

	FILE* restartkey = 0;
    if (is_restart) {
        printf("Trying to restart computation from checkpoint...\n");
		restartkey = parametersPrepareRestart(&parameters);
    } else {
        addParameter(PARAMETER_RESTART_COUNT, "0");
    }
    if(strncmp("CHARMM", parameters.parametersType, 6) == 0){
		readCHARMMFF(parameters.forceFieldFilename, &ff, parameters.parametersType);
		//if(strcmp(FF_TYPE_CHARMM19, parameters.parametersType) == 0){
		//	convertAtomTypesCHARMM19(parameters.topologiesFilename, atoms, atomCount);
		//}
	} else {
		die("FF parameters type other then CHARMM is not yet supported.\n");
	}

    // Redirect STDERR/STDOUT to files if needed
    char r_stdout[512], r_stderr[512];
    getMaskedParameter(r_stdout,PARAMETER_REDIRECT_STDOUT,PARAMETER_STRING_UNDEFINED);
    getMaskedParameter(r_stderr,PARAMETER_REDIRECT_STDERR,PARAMETER_STRING_UNDEFINED);
    bool r_stdout_flag = strcmp(r_stdout,PARAMETER_STRING_UNDEFINED);
    bool r_stderr_flag = strcmp(r_stderr,PARAMETER_STRING_UNDEFINED);
    if (r_stderr_flag) {
        printf ("Redirecting STDERR to %s\n", r_stderr);
        check(freopen(r_stderr, "w", stderr), "redirecting stderr");
    }
    if (r_stdout_flag) {
        printf ("Redirecting STDOUT to %s\n", r_stdout);
        check(freopen(r_stdout, "w", stdout), "redirecting stdout");
    }

    initTimer();

	initTopology();

	initAtomTypes();
	initAngleTypes();

	initGPU();

	if (is_restart) {
		void launchRestarters(FILE*);
		launchRestarters(restartkey);
	}
	splitTimer();
	compute();

	/*if(getYesNoParameter(PARAMETER_PARAMETRIZE_SOP, 0)){
		parametrizeSOP();
	} else {
		compute();
	}*/
    if (r_stdout_flag) fclose(stdout);
    if (r_stderr_flag) fclose(stderr);
#ifdef USE_MPI
    MPI::COMM_WORLD.Barrier();
    MPI::Finalize();
#endif
    return 0;
}

