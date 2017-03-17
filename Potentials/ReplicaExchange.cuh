/*
 * ReplicaExchange.cuh
 *
 *  Created on: Aug 10, 2011
 *      Author: serxa
 */

#pragma once

#include <vector>

namespace replica_exchange {

int blockSize;
int blockCount;

struct Replica {
	Replica(): T(0.0f), E(0.0f), id(-1), exchangeAttempts(0), successfulExchanges(0), energySum(0.0) {}
	float Tmax; // Max temperature (is equal to T after heating)
	float T; // Current temperature [K]
	float E; // Energy [kJ]
	int id; // Replica id at start (corresponds to its temperature)
	int exchangeAttempts;
	int successfulExchanges;
	double energySum; // Need to evaluate averageEnergy (i.e. devided by number of update() calls)
};

// array of all replica distributed by MPI
// NOTE: indexed by trajectory number, that is NOT a Replica::id
std::vector<Replica> replica;

// This is a pointer to the local part of replica array
Replica* localReplica;

std::vector<float> h_var;
float* d_var;
int rseed;
bool debug;
bool isExchangesDisabled; // [for debug] Do not exchange any replicas if set
size_t updatesCount; // number of calls to update()
size_t heatUpdatesLimit; // number of calls to update() for heating

Potential potential;
Updater updater;

inline bool isHeatModeEnabled() { return heatUpdatesLimit; }
inline bool isHeatMode() { return isHeatModeEnabled() && updatesCount <= heatUpdatesLimit; }
void create();
void inline compute();
void update();
void exchangeMode();
void heatMode();
void destroy();
void refreshVars();

void save(FILE*);
void load(FILE*);

}
