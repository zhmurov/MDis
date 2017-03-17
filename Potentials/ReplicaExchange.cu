/*
 * ReplicaExchange.cu
 *
 *  Created on: Aug 10, 2011
 *      Author: serxa
 */

#include <iostream>
#include <algorithm>
#include <fstream>
#include "../Util/Log.h"
#include <map>
#include "../Core/global.h"
#include "ReplicaExchange.cuh"
#include "../Core/parameters.h"
#include "../Util/ran2.h"
#include "../Updaters/EnergyOutputManager.cuh"

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace replica_exchange {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<replica_exchange> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

int numberOfTrajectories() {
#ifdef USE_MPI
	static const int result = parameters.Ntr * MPI::COMM_WORLD.Get_size();
#else
	static const int result = parameters.Ntr;
#endif
	return result;
}

int mpiRank() {
#ifdef USE_MPI
	static const int result = MPI::COMM_WORLD.Get_rank();
#else
	static const int result = 0;
#endif
	return result;
}

int trajId(int traj, int rank) {
	return traj + rank * parameters.Ntr;
}

void create() {
	LOG << "create";

	debug = (bool)getYesNoParameter(PARAMETER_REMD_DEBUG, 0);
	if (!getYesNoParameter(PARAMETER_REMD_ENABLED, 0))
		return;

	// Initialize globals
	isExchangesDisabled = (bool)getYesNoParameter(PARAMETER_REMD_DISABLE_EXCHANGES, 0);
	heatUpdatesLimit = getIntegerParameter(PARAMETER_REMD_HEATSTEPS, 0) / getIntegerParameter(PARAMETER_REMD_FREQ);
	blockCount = (int)ceil((float)gsystem.Ntot/BLOCK_SIZE);
	blockSize = BLOCK_SIZE;
	hybrid_taus::initRand(getLongIntegerParameter(PARAMETER_RSEED)*parameters.firstrun*(mpiRank()+1), gsystem.Ntot);
	rseed = getIntegerParameter(PARAMETER_RSEED);

	// Initialize replica-related stuff
	float Tmin = getFloatParameter(PARAMETER_REMD_MIN_TEMPERATURE);
	float Tmax = getFloatParameter(PARAMETER_REMD_MAX_TEMPERATURE);
	float dT = (Tmax - Tmin) / numberOfTrajectories();
	replica.resize(numberOfTrajectories());
	for (int traj = 0; traj < numberOfTrajectories(); traj++) {
		replica[traj].id = traj;
		replica[traj].Tmax = Tmin + dT * traj;
		if (isHeatModeEnabled())
			replica[traj].T = 0;
		else
			replica[traj].T = replica[traj].Tmax;
	}

	localReplica = &replica[mpiRank() * parameters.Ntr];

	refreshVars();

	// Initialize potential
	potential.compute = compute;
	potential.destroy = destroy;
	sprintf(potential.name, "Replica exchange");
	potentials[potentialsCount++] = &potential;

	// Initialize updater
	updater.update = update;
	updater.destroy = destroy;
	updater.frequency = getIntegerParameter(PARAMETER_REMD_FREQ);
	sprintf(updater.name, "Replica exchange");
	updaters[updatersCount++] = &updater;

	// Initialize restarter
	addRestarter("ReplicaExchange", save, load);

	LOG << "create done";
}

__global__ void compute_kernel(float* var) {
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < c_gsystem.Ntot) {
		float4 f = c_gsystem.d_forces[i];
		int at = c_gsystem.d_atomTypes[i];
		f.w = c_gsystem.d_m[at]; // Mass is now here.
		// TODO: This should be optimize with use of constant/texture memory. Artem, what do you think will be best here?
		// The best would be constant memory, since most of the times all threads in a warp will access the same 'var'.
		// However, I'm not sure that it is possible to allocate constant memory dynamically - we don't know in advance how many trajectories we will have
		// Texture memory is cached spatially, so some cache will be wasted if texture is used here.
		// I think texture is the best choice here though.
		float mult = var[i / c_gsystem.N] * sqrtf(f.w);
		float4 rf = hybrid_taus::rforce(i);
		f.x += mult * rf.x;
		f.y += mult * rf.y;
		f.z += mult * rf.z;
		c_gsystem.d_forces[i] = f;
	}
}

void inline compute() {
	compute_kernel<<<blockCount, blockSize>>>(d_var);
}

void destroy() {}

float computeReplicaEnergy(int traj) {
	// NOTE: this works only if called after energyOutputManager.update()
	return energyOutputData.E[traj];
}

double exchangeProbability(const Replica& a, const Replica& b) {
	return std::min(1.0, exp((a.E - b.E) * (1.0/(Kb_MD * a.T) - 1.0/(Kb_MD * b.T))));
}

void refreshVars() {
	// Allocate memory for the first time
	if (h_var.empty()) {
		h_var.resize(parameters.Ntr);
		cudaMalloc((void**)&d_var, sizeof(float) * parameters.Ntr);
	}

	// Compute coefficients
	static float h = getFloatParameter(PARAMETER_TIMESTEP);
	static float gamma = getFloatParameter(PARAMETER_DAMPING, DEFAULT_DAMPING);
	for (int traj = 0; traj < parameters.Ntr; traj++)
		h_var[traj] = sqrtf(2.0f*gamma*Kb_MD*localReplica[traj].T/h);

	// Write coefficients to device
	cudaMemcpy(d_var, &h_var[0], sizeof(float) * parameters.Ntr, cudaMemcpyHostToDevice);
}

__global__ void normalizeVelocities_kernel(int base, int size, float coef) {
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < size)
		c_gsystem.d_vel[base + i] *= coef;
}

void logExchangeReplicas(int i, int j) {
	static bool first = true;
	static std::string fileName;
	if (first) {
		first = false;
		// Get and save filename for exchanges
		fileName = getParameterAs<std::string>(PARAMETER_REMD_EXCHANGESFILE, "");
		// Clear file (we are going to append text to the end of it)
		if (!fileName.empty())
			std::ofstream(fileName.c_str());
	}

	if (!fileName.empty()) {
		std::ofstream ofs(fileName.c_str(), std::ios::app);
		ofs << step << '\t' << replica[i].id << '\t' << replica[j].id << '\n';
	}
}

int localTrajId(int traj) {
	return traj - mpiRank() * parameters.Ntr;
}

int rankForTraj(int traj) {
	return traj / parameters.Ntr;
}

void normalizeVelocities(int traj, float coef) {
	if (mpiRank() == rankForTraj(traj)) {
		unsigned gridSize = (unsigned)ceil((float)gsystem.N/blockSize);
		normalizeVelocities_kernel<<<gridSize, blockSize>>>(gsystem.N * localTrajId(traj), gsystem.N, coef);
	}
}

void exchangeReplicas(int i, int j) {
	// Actually exchange replicas/trajectories on CPU and GPU
	float coef = sqrt(replica[i].T/replica[j].T);
	normalizeVelocities(i, 1.0f/coef);
	normalizeVelocities(j, coef);
	std::swap(replica[i], replica[j]);
}

void broadcastExchangeInfo(int i, int j) {
#ifdef USE_MPI
	int size = MPI::COMM_WORLD.Get_size();
	for (int dst_rank = 1; dst_rank < size; dst_rank++) {
		if (MPI_Send(&i, 1, MPI_INT, dst_rank, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
			DIE("Unable to send traj num to exchange: traj=%d dst_rank=%d", i, dst_rank);
		if (MPI_Send(&j, 1, MPI_INT, dst_rank, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
			DIE("Unable to send traj num to exchange: traj=%d dst_rank=%d", j, dst_rank);
	}
#endif
}

void broadcastEndOfExchangesMarker() {
#ifdef USE_MPI
	int size = MPI::COMM_WORLD.Get_size();
	int marker = -1;
	for (int dst_rank = 1; dst_rank < size; dst_rank++)
		if (MPI_Send(&marker, 1, MPI_INT, dst_rank, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
			DIE("Unable to send marker: dst_rank=%d", dst_rank);
#endif
}

void update() {
	updatesCount++;
	if (isHeatMode())
		heatMode();
	else
		exchangeMode();
}

void heatMode() {
	LOG << "update (heat mode)";

	// Update temperatures of all replicas
	for (int i = 0; i < numberOfTrajectories(); i++)
		replica[i].T = double(updatesCount) / double(heatUpdatesLimit) * replica[i].Tmax;

	// Update environment
	refreshVars();
}

struct ReplicaTemperatureComp {
	bool operator()(int l, int r) const {
		return replica[l].T < replica[r].T;
	}
};

void sendReceiveEnergy() {
#ifdef USE_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	MPI_Status status;
	if (rank == 0) {
		// Receive energies
		for (int src_rank = 1; src_rank < size; src_rank++)
			for (int traj = 0; traj < parameters.Ntr; traj++)
				if (MPI_Recv(&replica[trajId(traj, src_rank)].E, 1, MPI_FLOAT, src_rank, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
					DIE("Unable to receive energy: src_rank=%d traj=%d error_code=%d", src_rank, traj, status.MPI_ERROR);
	} else {
		// Send energies
		for (int traj = 0; traj < parameters.Ntr; traj++)
			if (MPI_Send(&localReplica[traj].E, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
				DIE("Unable to send energy: traj=%d", traj);
	}
#endif
}

void exchangeMode() {
	LOG << "update (exchange mode)";
	static bool checked = false;
	if (!checked) {
		checked = true;
		if (updater.frequency % energy_output::updater.frequency != 0) {
			DIE("REMD updater frequency (%d) is not a multiple of energy output manager frequency (%d)", updater.frequency, energy_output::updater.frequency);
		}
	}

	// Compute local replicas energies
	for (int i = 0; i < parameters.Ntr; i++)
		localReplica[i].energySum += (localReplica[i].E = computeReplicaEnergy(i));

	// Send/receive replicas energies
	sendReceiveEnergy();

	if (mpiRank() == 0) {
		// Create sorted by replica temperature list of replica indices
		std::vector<int> idx;
		for (int i = 0, e = numberOfTrajectories(); i < e; i++)
			idx.push_back(i);
		std::sort(idx.begin(), idx.end(), ReplicaTemperatureComp());

		// Attempt replica exchanges
		int totalExchanges = 0, successfulExchanges = 0;
		std::stringstream sxchg;
		for (int n = int(ran2::ran2(&rseed) * 2); n < numberOfTrajectories() - 1; n += 2) {
			int i = idx[n];
			int j = idx[n + 1];
			replica[i].exchangeAttempts++;
			replica[j].exchangeAttempts++;
			totalExchanges++;
			double p = exchangeProbability(replica[i], replica[j]);
			int success = (ran2::ran2(&rseed) < p);
			if (success && !isExchangesDisabled) {
				replica[i].successfulExchanges++;
				replica[j].successfulExchanges++;

				logExchangeReplicas(i, j);
				broadcastExchangeInfo(i, j);
				exchangeReplicas(i, j);
				successfulExchanges++;
			}
			sxchg << " (" << i << ", " << j << ", " << p << ", " << success << ")";
		}

		broadcastEndOfExchangesMarker();

		// Create replica_id -> replica map
		std::map<int, Replica*> replicaIdx;
		std::stringstream smap;
		for (int i = 0; i < numberOfTrajectories(); i++) {
			smap << " " << replica[i].id;
			replicaIdx[replica[i].id] = &replica[i];
		}

		std::stringstream srates;
		std::stringstream senergy;
		for (int i = 0; i < numberOfTrajectories(); i++) {
			srates << " " << float(replicaIdx[i]->successfulExchanges) / replicaIdx[i]->exchangeAttempts;
			senergy << " " << (replicaIdx[i]->energySum / updatesCount) * KCALL_PER_KJ;
		}

		LOG << "Xchg:" << sxchg.str();
		LOG << successfulExchanges << "/" << totalExchanges << " exchanges have been done";
		LOG << "Map:" << smap.str();
		LOG << "Rates:" << srates.str();
		LOG << "AvgEnergy:" << senergy.str();
	} else {
#ifdef USE_MPI
		// Receive list of traj numers to exchange
		MPI_Status status;
		while (true) {
			int i, j;
			if (MPI_Recv(&i, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
				DIE("Unable to receive the first traj num to exchange");
			if (i == -1) // End-Of-Exchanges-List marker
				break;
			if (MPI_Recv(&j, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
				DIE("Unable to receive the second traj num to exchange");
			exchangeReplicas(i, j);
		}
#endif
	}

	// Update environment
	refreshVars();
}

void save(FILE* f) {
	LOG << "restart:save";
	fprintf(f, "%lu %d\n", updatesCount, rseed);
	for (size_t i = 0; i < (size_t) numberOfTrajectories(); ++i) {
		Replica& r = replica[i];
		fprintf(f, "%.10e %.10e %.10e %d %d %d %.10e\n", r.Tmax, r.T, r.E, r.id, r.exchangeAttempts, r.successfulExchanges, r.energySum);
	}
}

void load(FILE* f) {
	LOG << "restart:load";
	if (fscanf(f, " %lu %d ", &updatesCount, &rseed) != 2)
		DIE("Loading restart from file: unable to get updatesCount and rseed");
	printf("%lu %d\n", updatesCount, rseed);
	for (size_t i = 0; i < (size_t) numberOfTrajectories(); ++i) {
		Replica& r = replica[i];
		int ret;
		ret = fscanf(f, "%e %e %e %d %d %d %le ", &r.Tmax, &r.T, &r.E, &r.id, &r.exchangeAttempts, &r.successfulExchanges, &r.energySum);
		if (ret != 7) 
			DIE("Loading restart from file: unable to load data for traj #%d", i);
		printf("%.10e %.10e %.10e %d %d %d %.10e\n", r.Tmax, r.T, r.E, r.id, r.exchangeAttempts, r.successfulExchanges, r.energySum);
	}
	refreshVars();
}

#undef LOG

} // namespace replica_exchange
