/*
 * PairsListsUpdater.cuh
 *
 *  Created on: Aug 5, 2010
 *      Author: zhmurov
 */

#pragma once

namespace pairslist_updater {

int pairsListsBlockSize;
int pairsListsBlockCount;
int possiblePairsListBlockSize;
int possiblePairsListBlockCount;

typedef struct {

	int maxExcludedPerAtom;
	int* h_pairsExclusionListCount;
	int* d_pairsExclusionListCount;
	int* h_pairsExclusionList;
	int* d_pairsExclusionList;

	int maxPairs12PerAtom;
	int* h_pairs12ListCount;
	int* d_pairs12ListCount;
	int* h_pairs12List;
	int* d_pairs12List;

	int maxPairs13PerAtom;
	int* h_pairs13ListCount;
	int* d_pairs13ListCount;
	int* h_pairs13List;
	int* d_pairs13List;

	int maxPairs14PerAtom;
	int* h_pairs14ListCount;
	int* d_pairs14ListCount;
	int* h_pairs14List;
	int* d_pairs14List;

	int maxPairsLJPerAtom;
	int* h_pairsLJListCount;
	int* d_pairsLJListCount;
	int* h_pairsLJList;
	int* d_pairsLJList;
	float ljPairsCutoff;

	int maxPairsCoulombPerAtom;
	int* h_pairsCoulombListCount;
	int* d_pairsCoulombListCount;
	int* h_pairsCoulombList;
	int* d_pairsCoulombList;
	float coulombPairsCutoff;

	int maxPossiblePairsPerAtom;
	int* h_possiblePairsListCount;
	int* d_possiblePairsListCount;
	int* h_possiblePairsList;
	int* d_possiblePairsList;
	float possiblePairsCutoff;

} PairsListsData;

void create();
void init();
void updatePairsLists();
void updatePossiblePairsList();
void destroyPairsListsUpdater();
void destroyPossiblePairsListUpdater();

} // namespace pairslist_updater

pairslist_updater::PairsListsData pairsListsData;
int pairlistsExtension;
__device__ __constant__ pairslist_updater::PairsListsData c_pairsListsData;

Updater pairsListsUpdater;
Updater possiblePairsListUpdater;