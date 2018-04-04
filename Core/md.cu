/*
 * md.cu
 *
 *  Created on: Jul 28, 2010
 *      Author: zhmurov
 */
#include "md.cuh"
#include "../Util/Cuda.h"
#include "../Util/memory.h"
#include "../Util/mystl.h"
#include "global.h"
#include "../Potentials/PeriodicBoundary.cu"
#include "../Integrators/LeapFrogIntegrator.cu"
#include "../Integrators/SteepestDescentIntegrator.cu"
#include "../Potentials/HarmonicPotential.cu"
#include "../Potentials/AnglePotential.cu"
#include "../Potentials/DihedralPotential.cu"
#include "../Potentials/ImproperPotential.cu"
#include "../Potentials/NonBondedPotential.cu"
#include "../Potentials/SASAPotential.cu"
#include "../Potentials/GenBornPotential.cu"
#include "../Potentials/LangevinHeatBath.cu"
#include "../Potentials/ReplicaExchange.cu"
#include "../Potentials/FixForcePotential.cu"
#include "../Potentials/UmbrellaSampling.cu"
#include "../Potentials/RepulsiveBoundaryPotential.cu"
#include "../Potentials/HarmonicConstraints.cu"
//#include "../Potentials/GBSWPotential.cu"
#include "../Potentials/GBPotential.cu"
#include "../Potentials/PullingPlanePotential.cu"
#include "../Potentials/PushingPlanePotential.cu"
#include "../Potentials/DrumPotential.cu"
#include "../Potentials/FragmemPotential.cu"
#include "../Updaters/CoordinatesOutputManagerDCD.cu"
#include "../Updaters/EnergyOutputManager.cu"
#include "../Updaters/PairsListsUpdater.cu"
#include "../Updaters/RestartOutputManager.cu"
#include "../Updaters/RigidBody.cu"
#include "../Updaters/PairsMeshUpdater.cu"
#include "../ConstrAlgorithms/shake.cu"
#include "../ConstrAlgorithms/ccma.cu"
//#include "SOPGPU/SOPGPUParameterizer.cu"

class LogMD: public ILog {
    virtual void Write(const char* message) const {
        std::cout << makeTimePrefix() << "<md_core> " << message << std::endl;
    }
} log_md;

#define LOG LogStream(log_md)

Potential** potentials;
Updater** updaters;
Integrator* integrator;
ConstrAlg** constrAlgs;
EnergyOutput** energyOutputs;
Restarter** restarters;

int potentialsCount;
int updatersCount;
int constrAlgsCount;
int energyOutputsCount;
int restartersCount;

long int lastStepCoordCopied = -1;
long int lastStepVelCopied = -1;

void initAtomTypesOnGPU();
void copyCoordinatesToGPU(int traj);
void copyVelocitiesToGPU(int traj);
void copyMassesToGPU();
void copyForcesToGPU();
void copyAtomTypesToGPU();
void copyCoordinatesFromGPU(bool force);
void copyVelocitiesFromGPU();

void addRestarter(const char* name, void (*save)(FILE*), void (*load)(FILE*)){
	for (const char* i = name; *i; i++) {
		if (isspace(*i)) {
			DIE("internal error: restarter name contain whitespace");
		}
	}
	Restarter* r = (Restarter*)malloc(sizeof(Restarter));
	if (!r) {
		DIE("Out of memory");
	}
	r->save = save;
	r->load = load;
	if (snprintf(r->name, sizeof(r->name), "%s", name) >= (int) sizeof(r->name)) {
		DIE("500 Program internal error");
	}
	restarters[restartersCount++] = r;
}

void launchRestarters(FILE* keyf) {
	for (int i = 0; i < restartersCount; i++) {
		char name[4096];
		if (fscanf(keyf, "\n%s\n", name) != 1)
            DIE("Cannot read restarter #%d", i);
		if (strcmp(name, restarters[i]->name) != 0) {
			DIE("Got wrong restarter name '%s', expected '%s'", name, restarters[i]->name);
		}
		restarters[i]->load(keyf);
	}
}

void dumpTOP(){
	TOPData top;
	top.atomCount = topology.atomCount;
	top.atoms = (TOPAtom*)calloc(top.atomCount, sizeof(TOPAtom));
	int i;
	for(i = 0; i < topology.atomCount; i++){
		top.atoms[i].id = topology.atoms[i].id;
		sprintf(top.atoms[i].type, "%d", topology.atoms[i].typeId);
		top.atoms[i].resid = topology.atoms[i].resid;
		sprintf(top.atoms[i].resName, "%s", topology.atoms[i].resName);
		sprintf(top.atoms[i].name, "%s", topology.atoms[i].name);
		top.atoms[i].chain = topology.atoms[i].segment[0];
		top.atoms[i].charge = topology.atoms[i].charge;
		top.atoms[i].mass = topology.atoms[i].mass;
	}
	top.bondCount = topology.bondCount;
	top.bonds = (TOPPair*)calloc(top.bondCount, sizeof(TOPPair));
	int b;
	for(b = 0; b < topology.bondCount; b++){
		top.bonds[b].i = topology.bonds[b].i;
		top.bonds[b].j = topology.bonds[b].j;
		top.bonds[b].func = 1;
		top.bonds[b].c0 = topology.bonds[b].b0;
		top.bonds[b].c1 = topology.bonds[b].kb;
	}
	top.angleCount = topology.angleCount;
	top.angles = (TOPAngle*)calloc(top.angleCount, sizeof(TOPAngle));
	int a;
	for(a = 0; a < top.angleCount; a++){
		top.angles[a].i = topology.angles[a].i;
		top.angles[a].j = topology.angles[a].j;
		top.angles[a].k = topology.angles[a].k;
		top.angles[a].func = 1;
		top.angles[a].c0 = topology.angles[a].theta0*180.0/M_PI;
		top.angles[a].c1 = topology.angles[a].ktheta;
	}
	int dihedralCount = 0;
	int d;
	for(d = 0; d < topology.dihedralCount; d++){
		dihedralCount += topology.dihedrals[d].multiplicity;
	}
	for(d = 0; d < topology.improperCount; d++){
		dihedralCount += topology.impropers[d].multiplicity;
	}
	top.dihedralCount = dihedralCount;
	top.dihedrals = (TOPDihedral*)calloc(top.dihedralCount, sizeof(TOPDihedral));

	dihedralCount = 0;
	for(d = 0; d < topology.dihedralCount; d++){
		for(i = 0; i < topology.dihedrals[d].multiplicity; i++){
			top.dihedrals[dihedralCount].i = topology.dihedrals[d].i;
			top.dihedrals[dihedralCount].j = topology.dihedrals[d].j;
			top.dihedrals[dihedralCount].k = topology.dihedrals[d].k;
			top.dihedrals[dihedralCount].l = topology.dihedrals[d].l;
			top.dihedrals[dihedralCount].func = 1;
			top.dihedrals[dihedralCount].parCount = 3;
			top.dihedrals[dihedralCount].c0 = topology.dihedrals[d].delta[i]*180.0/M_PI;
			top.dihedrals[dihedralCount].c1 = topology.dihedrals[d].kchi[i];
			top.dihedrals[dihedralCount].c2 = topology.dihedrals[d].n[i];
			dihedralCount++;
		}
	}
	for(d = 0; d < topology.improperCount; d++){
		for(i = 0; i < topology.impropers[d].multiplicity; i++){
			top.dihedrals[dihedralCount].i = topology.impropers[d].i;
			top.dihedrals[dihedralCount].j = topology.impropers[d].j;
			top.dihedrals[dihedralCount].k = topology.impropers[d].k;
			top.dihedrals[dihedralCount].l = topology.impropers[d].l;
			top.dihedrals[dihedralCount].func = 2;
			top.dihedrals[dihedralCount].parCount = 2;
			top.dihedrals[dihedralCount].c0 = topology.impropers[d].psi0[i]*180.0/M_PI;
			top.dihedrals[dihedralCount].c1 = topology.impropers[d].kpsi[i];
			//top.dihedrals[dihedralCount].c2 = topology.impropers[d].n[i];
			dihedralCount++;
		}
	}
	top.exclusionCount = topology.exclusionsCount;
	top.exclusions = (TOPExclusion*)calloc(top.exclusionCount, sizeof(TOPExclusion));
	int e;
	for(e = 0; e < top.exclusionCount; e++){
		top.exclusions[e].i = topology.exclusions[e].i;
		top.exclusions[e].j = topology.exclusions[e].j;
		top.exclusions[e].func = 1;
	}
	top.pairsCount = 0;
	writeTOP("dump.top", &top);
}


void initGPU(){
	dumpTOP();
	LOG << "Preparing system on a GPU...";
	if (parameters.device >= 0) {
		cudaSetDevice(parameters.device);
	} else {
		cudaGetDevice(&parameters.device); // Allow automagic to do all the work
	}
	cudaGetDeviceProperties(&deviceProps, parameters.device);
	gsystem.N = topology.atomCount;
	gsystem.Nsim = getIntegerParameter(PARAMETER_NSIM, gsystem.N);
	if(gsystem.Nsim % BLOCK_SIZE != 0){
		gsystem.widthSim = (gsystem.Nsim/16 + 1)*16;
	} else {
		gsystem.widthSim = gsystem.Nsim;
	}
	gsystem.Ntot = gsystem.N*parameters.Ntr;
	if(gsystem.Ntot % BLOCK_SIZE != 0){
		gsystem.widthTot = (gsystem.Ntot/BLOCK_SIZE + 1)*BLOCK_SIZE;
	} else {
		gsystem.widthTot = gsystem.Ntot;
	}
	LOG << "Preparing data for " << parameters.Ntr << " trajectories (" << gsystem.N << " atoms in a system, " << gsystem.Ntot << " total atoms).";
	LOG << "Arrays will be aligned to the width of " << gsystem.widthSim << " and total width of " << gsystem.widthTot;

	allocateCPU((void**)&gsystem.h_coord, gsystem.Ntot*sizeof(float4));
	allocateCPU((void**)&gsystem.h_vel, gsystem.Ntot*sizeof(float4));
	allocateCPU((void**)&gsystem.h_m, atomTypesCount*sizeof(float));
	allocateCPU((void**)&gsystem.h_forces, gsystem.Ntot*sizeof(float4));
	allocateCPU((void**)&gsystem.h_atomTypes, gsystem.Ntot*sizeof(int));
	allocateGPU((void**)&gsystem.d_coord, gsystem.Ntot*sizeof(float4));
	allocateGPU((void**)&gsystem.d_midcoord, gsystem.Ntot*sizeof(float4));
	allocateGPU((void**)&gsystem.d_vel, gsystem.Ntot*sizeof(float4));
	allocateGPU((void**)&gsystem.d_m, atomTypesCount*sizeof(float));
	allocateGPU((void**)&gsystem.d_forces, gsystem.Ntot*sizeof(float4));
	allocateGPU((void**)&gsystem.d_atomTypes, gsystem.Ntot*sizeof(int));

	cudaMemset(gsystem.d_forces, 0, gsystem.Ntot*sizeof(float4));
	cudaMemset(gsystem.d_vel, 0, gsystem.Ntot*sizeof(float4));

	potentialsCount = 0;
	updatersCount = 0;
	energyOutputsCount = 0;
	restartersCount = 0;
	potentials = (Potential**)calloc(MAX_POTENTIALS_COUNT, sizeof(Potential*));
	int p, u, i;
	for(p = 0; p < MAX_POTENTIALS_COUNT; p++){
		potentials[p] = (Potential*)malloc(sizeof(Potential*));
	}
	updaters = (Updater**)calloc(MAX_UPDATERS_COUNT, sizeof(Updater*));
	for(u = 0; u < MAX_UPDATERS_COUNT; u++){
		updaters[u] = (Updater*)malloc(sizeof(Updater*));
	}
	constrAlgs = (ConstrAlg**)calloc(MAX_CONSTR_ALGS_COUNT, sizeof(ConstrAlg*));
	for(i = 0; i < MAX_CONSTR_ALGS_COUNT; i++){
		constrAlgs[i] = (ConstrAlg*)malloc(sizeof(ConstrAlg*));
	}
	energyOutputs = (EnergyOutput**)calloc(MAX_ENERGY_OUTPUT_COUNT, sizeof(EnergyOutput*));
	for(i = 0; i < MAX_ENERGY_OUTPUT_COUNT; i++){
		energyOutputs[i] = (EnergyOutput*)malloc(sizeof(EnergyOutput*));
	}
	restarters = (Restarter**)calloc(MAX_RESTARTERS_COUNT, sizeof(Restarter*));
	for(i = 0; i < MAX_RESTARTERS_COUNT; i++){
		restarters[i] = (Restarter*)malloc(sizeof(Restarter*));
	}

	int traj;
	char trajnum[10];
	char trajFilename[100];
	for(traj = 0; traj < parameters.Ntr; traj++){
		sprintf(trajnum, "%d", traj + parameters.firstrun);
		replaceString(trajFilename, parameters.coordFilename, trajnum, "<run>");
		readCoordinates(trajFilename);
		if(strcmp(parameters.velFilename, PARAMETER_STRING_UNDEFINED) != 0){
			replaceString(trajFilename, parameters.velFilename, trajnum, "<run>");
			readVelocities(trajFilename);
		} else {
			float T = getFloatParameter(PARAMETER_INITIAL_TEMPERATURE, -1.0f);
			if(T == -1.0f){
				T = getFloatParameter(PARAMETER_TEMPERATURE, 0.0f);
			}
			LOG << "generating velocities with T=" << T << " K";
			generateVelocities(T);
		}
		copyCoordinatesToGPU(traj);
		copyVelocitiesToGPU(traj);
	}
	checkCUDAError("copying coordinate/velocities to GPU");

	copyMassesToGPU();
	checkCUDAError("copying masses to GPU");
	copyForcesToGPU();
	checkCUDAError("copying forces to GPU");
	copyAtomTypesToGPU();
	checkCUDAError("copying atoms to GPU");
	cudaMemcpyToSymbol(c_gsystem, &gsystem, sizeof(GSystem), 0, cudaMemcpyHostToDevice);
	checkCUDAError("init c_gsystem");

	char integratorName[100];
	getMaskedParameter(integratorName, PARAMETER_INTEGRATOR, PARAMETER_VALUE_INTEGRATOR_LEAPFROG);
	if(strcmp(integratorName, PARAMETER_VALUE_INTEGRATOR_STEEPEST_DESCENT) == 0){
		LOG << "Steepest Descent integrator requested. Energy minimization will be performed.";
		sd_integrator::create();
	} else
	if(strcmp(integratorName, PARAMETER_VALUE_INTEGRATOR_LEAPFROG) == 0){
		LOG << "Leap-Frog integrator will be used.";
		leapfrog_integrator::create();
	}
	checkCUDAError("init integrator");

	char rigidbonds[100];
	getMaskedParameter(rigidbonds, PARAMETER_RIGIDBONDS, PARAMETER_VALUE_RIGIDBONDS_NONE);
	if(strcmp(rigidbonds, PARAMETER_VALUE_RIGIDBONDS_HYDROGEN) == 0){
			LOG << "Hbond lengths will be constrained with SHAKE algorithm\n";
			shake_constrAlg::create();
	}
	else if(strcmp(rigidbonds, PARAMETER_VALUE_RIGIDBONDS_ALL) == 0){
			LOG << "All bond lengths will be constrained with CCMA algorithm\n";
			ccma_constrAlg::create();
	}
	else if(strcmp(integratorName, PARAMETER_VALUE_RIGIDBONDS_NONE) == 0){
			LOG << "No constraints will be appied.\n";
	}
	checkCUDAError("init constraint algorithms");





    umbrella_sampling::create();
    checkCUDAError("Init umbrella sampling");

	harmonic_potential::create();
	checkCUDAError("init harmonic potential");
	angle_potential::create();
	checkCUDAError("init angle potential");
	dihedral_potential::create();
	checkCUDAError("init dihedral potential");
	improper_potential::create();
	checkCUDAError("init improper potential");
	non_bonded_potential::create();
	checkCUDAError("init nonbonded potential");
	sasa_potential::create();
	checkCUDAError("init SASA potential");
	if(getYesNoParameter(PARAMETER_LANGEVIN_ON, DEFAULT_LANGEVIN_ON)){
		langevin_heat_bath::create();
		checkCUDAError("init Langevin Heat-Bath potential");
	}
	repulsive_boundary::create();
	checkCUDAError("init repulsive boundary");
	periodic_boundary::create();
	checkCUDAError("init periodic boundary");

	harmonic_constraints::create();
	checkCUDAError("init harmonic constraints boundary");
	

	//gbsw_potential::create();
	//checkCUDAError("init GBSW potential");
	gb_potential::create();
	checkCUDAError("init GB potential");

	// Pair list must be updated after we finish LD
	rigid_body::create();
	checkCUDAError("init Rigid Body");

	pairslist_updater::create();
	checkCUDAError("init pairlist updater");

	// Since GenBorn Potential uses pairlists in its initialization
	genborn_potential::create();
	checkCUDAError("init GenBorn potential");

	fixforce_potential::create();
	checkCUDAError("init Fix-Force potential");

	pulling_plane_potential::create();
	checkCUDAError("init Pulling-Plane potential");

	pushing_plane_potential::create();
	checkCUDAError("init Pushing-Plane potential");

	drum_potential::create();
	checkCUDAError("init Drum potential");

	fragmem_potential::create();
	checkCUDAError("init FragMem potential");

	coordinates_output_dcd::create();
	checkCUDAError("init DCD output manager");
	energy_output::create();
	checkCUDAError("init energy output manager");
	// Replica exchange updater must be called after output manager, bacause latter computes energies
	replica_exchange::create();
	checkCUDAError("init replica exchange");
	restart_output::create();
	checkCUDAError("init output manager");
	pairsmesh::create();
	checkCUDAError("init mesh updater");


	cudaBindTexture(0, t_coord, gsystem.d_coord, gsystem.Ntot*sizeof(float4));
	cudaBindTexture(0, t_m, gsystem.d_m, atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_atomTypes, gsystem.d_atomTypes, gsystem.Ntot*sizeof(float));

	checkCUDAError("init");
	LOG << "Done preparing system on a GPU.";
	printMemoryUsed();
}

void compute(){
	int p, u, nav;
	nav = updaters[0]->frequency;
	for(u = 0; u < updatersCount; u++){
		nav = GCD(nav, updaters[u]->frequency);
	}
	int i, traj;
	firststep = getLongIntegerParameter(PARAMETER_FIRSTSTEP, 0);
	step = firststep;
	// This is now implemented in parameters.cpp and read from restartkey file
	//for(traj = 0; traj < parameters.Ntr; traj++){
	//	trajectoryTime[traj] = firststep*integrator->h;
	//}
	printTime(step - firststep);

	while(step < parameters.numsteps){
		for(u = 0; u < updatersCount; u++){
			if((step - firststep) % updaters[u]->frequency == 0){
				cudaThreadSynchronize();
				updaters[u]->update();
				checkCUDAError(updaters[u]->name);
			}
		}
		for(i = 0; i < nav; i++){
			for(p = 0; p < potentialsCount; p++){
				//cudaThreadSynchronize();
				potentials[p]->compute();
				checkCUDAError(potentials[p]->name);
			}
			/*cudaMemcpy(gsystem.h_forces, gsystem.d_forces, gsystem.Ntot*sizeof(float4), cudaMemcpyDeviceToHost);
			FILE* out = fopen("forces.dat", "w");
			for(i = 0; i < gsystem.N; i++){
				fprintf(out, "%d\t%f\t%f\t%f\n", i, gsystem.h_forces[i].x, gsystem.h_forces[i].y, gsystem.h_forces[i].z);
			}
			fclose(out);
			exit(0);*/
            step++;
			cudaThreadSynchronize();
			integrator->integrate();
			cudaThreadSynchronize();
			for(p = 0; p < constrAlgsCount; p++){
				constrAlgs[p]->compute();
				checkCUDAError(constrAlgs[p]->name);
				cudaThreadSynchronize();
			}
			integrator->finalize();
			//checkCUDAError(integrator->name);
			cudaThreadSynchronize();
		}
		for(traj = 0; traj < parameters.Ntr; traj++){
			trajectoryTime[traj] += nav*integrator->h;
		}
	}
	for(u = 0; u < updatersCount; u++){
		if(step % updaters[u]->frequency == 0){
			cudaThreadSynchronize();
			updaters[u]->update();
			checkCUDAError(updaters[u]->name);
		}
	}
	checkCUDAError("finalizing");

	for(p = 0; p < potentialsCount; p++){
		potentials[p]->destroy();
	}
	for(u = 0; u < updatersCount; u++){
		updaters[u]->destroy();
	}
}

void copyCoordinatesToGPU(int traj, int reset = 1){
	int i;
	for(i = 0; i < gsystem.N && reset; i++){
		gsystem.h_coord[i + traj*gsystem.N].x = topology.atoms[i].x;
		gsystem.h_coord[i + traj*gsystem.N].y = topology.atoms[i].y;
		gsystem.h_coord[i + traj*gsystem.N].z = topology.atoms[i].z;
		gsystem.h_coord[i + traj*gsystem.N].w = (float)topology.atoms[i].typeId;
	}
	cudaMemcpy(&gsystem.d_coord[traj*gsystem.N], &gsystem.h_coord[traj*gsystem.N],
			gsystem.N*sizeof(float4), cudaMemcpyHostToDevice);
	cudaMemcpy(&gsystem.d_midcoord[traj*gsystem.N], &gsystem.h_coord[traj*gsystem.N],
				gsystem.N*sizeof(float4), cudaMemcpyHostToDevice);//the coordinates will be overwritten on the first step, we do it to get atomid's in midcoord
}
void copyCoordinatesToGPU(int traj) { copyCoordinatesToGPU(traj, 1); }

void copyVelocitiesToGPU(int traj, int reset = 1){
	int i;
	for(i = 0; i < gsystem.N && reset; i++){
		gsystem.h_vel[i + traj*gsystem.N].x = topology.atoms[i].vx;
		gsystem.h_vel[i + traj*gsystem.N].y = topology.atoms[i].vy;
		gsystem.h_vel[i + traj*gsystem.N].z = topology.atoms[i].vz;
		gsystem.h_vel[i + traj*gsystem.N].w = 0.0f;
	}
	cudaMemcpy(&gsystem.d_vel[traj*gsystem.N], &gsystem.h_vel[traj*gsystem.N],
			gsystem.N*sizeof(float4), cudaMemcpyHostToDevice);
}
void copyVelocitiesToGPU(int traj) { copyVelocitiesToGPU(traj, 1); }

void copyMassesToGPU(){
	int i;
	for(i = 0; i < atomTypesCount; i++){
		gsystem.h_m[i] = atomTypes[i].mass;
	}
	cudaMemcpy(gsystem.d_m, gsystem.h_m, atomTypesCount*sizeof(float), cudaMemcpyHostToDevice);
}

void copyForcesToGPU(){
	int i, itot, traj;
	for(i = 0; i < gsystem.N; i++){
		gsystem.h_forces[i].x = 0.0f;
		gsystem.h_forces[i].y = 0.0f;
		gsystem.h_forces[i].z = 0.0f;
		gsystem.h_forces[i].w = topology.atoms[i].mass;
	}
	for(traj = 1; traj < parameters.Ntr; traj ++){
		for(i = 0; i < gsystem.N; i++){
			itot = gsystem.N*traj + i;
			gsystem.h_forces[itot] = gsystem.h_forces[i];
		}
	}
	cudaMemcpy(gsystem.d_forces, gsystem.h_forces, gsystem.Ntot*sizeof(float4), cudaMemcpyHostToDevice);
}

void copyAtomTypesToGPU(){
	int i, itot, traj;
	for(i = 0; i < gsystem.N; i++){
		gsystem.h_atomTypes[i] = topology.atoms[i].typeId;
	}
	for(traj = 1; traj < parameters.Ntr; traj ++){
		for(i = 0; i < gsystem.N; i++){
			itot = gsystem.N*traj + i;
			gsystem.h_atomTypes[itot] = gsystem.h_atomTypes[i];
		}
	}
	cudaMemcpy(gsystem.d_atomTypes, gsystem.h_atomTypes, gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);
}

void copyCoordinatesFromGPU(bool force){
	if(force || lastStepCoordCopied != step){
		cudaMemcpy(gsystem.h_coord, gsystem.d_coord, gsystem.Ntot*sizeof(float4), cudaMemcpyDeviceToHost);
		lastStepCoordCopied = step;
	}
}

void copyVelocitiesFromGPU(){
	if(lastStepVelCopied != step){
		cudaMemcpy(gsystem.h_vel, gsystem.d_vel, gsystem.Ntot*sizeof(float4), cudaMemcpyDeviceToHost);
		lastStepVelCopied = step;
	}
}
