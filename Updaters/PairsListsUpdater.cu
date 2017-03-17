/*
 * PairsListsUpdater.cu
 *
 *  Created on: Aug 5, 2010
 *      Author: zhmurov
 */
#include "../Core/global.h"
#include "../Core/md.cuh"
#include "PairsListsUpdater.cuh"

namespace pairslist_updater {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<pairslist_updater> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

int checkAtomForAromaticExclusions(Atom atom);
int checkAtomsForExplicitExclusion(int a1, int a2);
int addAtomToListIfNotYetThere(int i, int j, int* counts, int* list, int maxCount, int width);
void sortList(int* counts, int* list, int N, int width);
int checkAtomsBonded(int i, int j);
int checkAtoms13Bonded(int i, int j);
int checkAtoms14Bonded(int i, int j);
int checkAtomsInExclusionList(int i, int j);

void create(){

	possiblePairsListUpdater.update = updatePossiblePairsList;
	possiblePairsListUpdater.destroy = destroyPossiblePairsListUpdater;
	possiblePairsListUpdater.frequency = getIntegerParameter(PARAMETER_POSSIBLEPAIRS_FREQ);
	sprintf(possiblePairsListUpdater.name, "Possible Pairs List updater");
	updaters[updatersCount] = &possiblePairsListUpdater;
	updatersCount ++;

	pairsListsUpdater.update = updatePairsLists;
	pairsListsUpdater.destroy = destroyPairsListsUpdater;
	pairsListsUpdater.frequency = getIntegerParameter(PARAMETER_PAIRS_FREQ);
	sprintf(pairsListsUpdater.name, "Verlet list updater");
	updaters[updatersCount] = &pairsListsUpdater;
	updatersCount ++;

	init();
}

void init(){
	LOG << "Initializing pairs lists...";

	pairsListsBlockSize = BLOCK_SIZE;
	pairsListsBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;
	possiblePairsListBlockSize = BLOCK_SIZE;
	possiblePairsListBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;

	pairlistsExtension = getIntegerParameter(PARAMETER_PAIRSLISTS_EXTENSION, DEFAULT_PAIRSLISTS_EXTENSION);

	pairsListsData.ljPairsCutoff = getFloatParameter(PARAMETER_PAIRS_CUTOFF_LJ, -1.0f);
	if(pairsListsData.ljPairsCutoff == -1.0f){
		pairsListsData.ljPairsCutoff = getFloatParameter(PARAMETER_PAIRS_CUTOFF_NB, DEFAULT_PAIRS_CUTOFF_LJ);
	}
	pairsListsData.coulombPairsCutoff = getFloatParameter(PARAMETER_PAIRS_CUTOFF_NB, -1.0f);
	if(pairsListsData.coulombPairsCutoff == -1.0f){
		pairsListsData.coulombPairsCutoff = getFloatParameter(PARAMETER_PAIRS_CUTOFF_NB, DEFAULT_PAIRS_CUTOFF_COULOMB);
	}
	pairsListsData.possiblePairsCutoff = getFloatParameter(PARAMETER_POSSIBLEPAIRS_CUTOFF, DEFAULT_POSSIBLEPAIRS_CUTOFF);


	LOG << "Searching for 1-2, 1-3 and 1-4 pairs for exclusion list and 1-4 list...";

	LOG << "Counting pairs...";
	int count12 = 0;
	int count13 = 0;
	int count14 = 0;
	int countLJ = 0;
	int countExplicitExclusions = 0;
	int countExcluded = 0;
	int countCoulomb = 0;
	int countPP = 0;
	int i, j, k, b, b1, b2, a, d;
	int traj, itot, jtot, p;


	allocateCPU((void**)&pairsListsData.h_pairs12ListCount, gsystem.Ntot*sizeof(int));
	pairsListsData.maxPairs12PerAtom = getIntegerParameter(PARAMETER_MAX_PAIRS_12, 0);
	if(pairsListsData.maxPairs12PerAtom == 0){
		LOG << "Counting covalent bonds...";
		for(i = 0; i < gsystem.N; i++){
			pairsListsData.h_pairs12ListCount[i] = 0;
		}
		for(b = 0; b < topology.bondCount; b++){
			Bond bond = topology.bonds[b];
			i = bond.i;
			j = bond.j;
			pairsListsData.h_pairs12ListCount[i] ++;
			pairsListsData.h_pairs12ListCount[j] ++;
		}
		pairsListsData.maxPairs12PerAtom = 0;
		for(i = 0; i < gsystem.N; i++){
			if(pairsListsData.maxPairs12PerAtom < pairsListsData.h_pairs12ListCount[i]){
				pairsListsData.maxPairs12PerAtom = pairsListsData.h_pairs12ListCount[i];
			}
			count12 += pairsListsData.h_pairs12ListCount[i];
		}
		LOG << "Maximum covalent bonds per atom is " << pairsListsData.maxPairs12PerAtom << ". Total is " << count12/2;
	} else {
		LOG << "Using pre-defined 1-2 listsize of " <<  pairsListsData.maxPairs12PerAtom;
	}

	LOG << "Constructing list of covalent bonds...";
	allocateCPU((void**)&pairsListsData.h_pairs12List,
			gsystem.widthTot*pairsListsData.maxPairs12PerAtom*sizeof(int));
	for(i = 0; i < gsystem.N; i++){
		pairsListsData.h_pairs12ListCount[i] = 0;
	}
	for(b = 0; b < topology.bondCount; b++){
		Bond bond = topology.bonds[b];
		i = bond.i;
		j = bond.j;
		pairsListsData.h_pairs12List[pairsListsData.h_pairs12ListCount[i]*gsystem.widthTot + i] = j;
		pairsListsData.h_pairs12List[pairsListsData.h_pairs12ListCount[j]*gsystem.widthTot + j] = i;
		pairsListsData.h_pairs12ListCount[i] ++;
		pairsListsData.h_pairs12ListCount[j] ++;
	}
	count12 = 0;
	for(i = 0; i < gsystem.N; i++){
		count12 += pairsListsData.h_pairs12ListCount[i];
	}
	LOG << "Done constructing list of covalent bonds.";




	allocateCPU((void**)&pairsListsData.h_pairs13ListCount, gsystem.Ntot*sizeof(int));
	pairsListsData.maxPairs13PerAtom = getIntegerParameter(PARAMETER_MAX_PAIRS_13, 0);
	if(pairsListsData.maxPairs13PerAtom == 0){
		LOG << "Counting 1-3 pairs...";
		for(i = 0; i < gsystem.N; i++){
			pairsListsData.h_pairs13ListCount[i] = 0;
		}
		for(i = 0; i < gsystem.N; i++){
			for(b1 = 0; b1 < pairsListsData.h_pairs12ListCount[i]; b1++){
				j = pairsListsData.h_pairs12List[b1*gsystem.widthTot + i];
				for(b2 = 0; b2 < pairsListsData.h_pairs12ListCount[j]; b2++){
					k = pairsListsData.h_pairs12List[b2*gsystem.widthTot + j];
					if(i != k){
						pairsListsData.h_pairs13ListCount[i] ++;
					}
				}
			}
		}
		pairsListsData.maxPairs13PerAtom = 0;
		for(i = 0; i < gsystem.N; i++){
			if(pairsListsData.maxPairs13PerAtom < pairsListsData.h_pairs13ListCount[i]){
				pairsListsData.maxPairs13PerAtom = pairsListsData.h_pairs13ListCount[i];
			}
			count13 += pairsListsData.h_pairs13ListCount[i];
		}
		LOG << "Maximum 1-3 pairs per atom is " << pairsListsData.maxPairs13PerAtom << ". Total is " << count13/2;
	} else {
		LOG << "Using pre-defined 1-3 listsize of " << pairsListsData.maxPairs13PerAtom;
	}

	LOG << "Constructing a 1-3 pairs list...";
	allocateCPU((void**)&pairsListsData.h_pairs13List,
				gsystem.widthTot*pairsListsData.maxPairs13PerAtom*sizeof(int));
	for(i = 0; i < gsystem.N; i++){
		pairsListsData.h_pairs13ListCount[i] = 0;
	}
	for(i = 0; i < gsystem.N; i++){
		for(b1 = 0; b1 < pairsListsData.h_pairs12ListCount[i]; b1++){
			j = pairsListsData.h_pairs12List[b1*gsystem.widthTot + i];
			for(b2 = 0; b2 < pairsListsData.h_pairs12ListCount[j]; b2++){
				k = pairsListsData.h_pairs12List[b2*gsystem.widthTot + j];
				if(i != k && !checkAtomsBonded(i, k)){
					addAtomToListIfNotYetThere(i, k, pairsListsData.h_pairs13ListCount,
							pairsListsData.h_pairs13List, pairsListsData.maxPairs13PerAtom, gsystem.widthTot);
					//pairsListsData.h_pairs13List[pairsListsData.h_pairs13ListCount[i]*gsystem.widthTot + i] = k;
					//pairsListsData.h_pairs13ListCount[i] ++;
				}
			}
		}
	}
	count13 = 0;
	for(i = 0; i < gsystem.N; i++){
		count13 += pairsListsData.h_pairs13ListCount[i];
	}
	LOG << "Done constructing list of 1-3 pairs.";




	allocateCPU((void**)&pairsListsData.h_pairs14ListCount, gsystem.Ntot*sizeof(int));
	pairsListsData.maxPairs14PerAtom = getIntegerParameter(PARAMETER_MAX_PAIRS_14, 0);
	if(pairsListsData.maxPairs14PerAtom == 0){
		LOG << "Counting list of 1-4 pairs...";
		for(i = 0; i < gsystem.N; i++){
			pairsListsData.h_pairs14ListCount[i] = 0;
		}
		for(i = 0; i < gsystem.N; i++){
			for(b = 0; b < pairsListsData.h_pairs12ListCount[i]; b++){
				j = pairsListsData.h_pairs12List[b*gsystem.widthTot + i];
				for(a = 0; a < pairsListsData.h_pairs13ListCount[j]; a++){
					k = pairsListsData.h_pairs13List[a*gsystem.widthTot + j];
					if(i != k && !checkAtomsBonded(i, k) && !checkAtoms13Bonded(i, k)){
						//if(!(checkAtomForAromaticExclusions(atoms[i]) && checkAtomForAromaticExclusions(atoms[k]))){
						if(!checkAtomsForExplicitExclusion(i, k)){
							pairsListsData.h_pairs14ListCount[i] ++;
						}
					}
				}
			}
		}
		pairsListsData.maxPairs14PerAtom = 0;
		for(i = 0; i < gsystem.N; i++){
			if(pairsListsData.maxPairs14PerAtom < pairsListsData.h_pairs14ListCount[i]){
				pairsListsData.maxPairs14PerAtom = pairsListsData.h_pairs14ListCount[i];
			}
			count14 += pairsListsData.h_pairs14ListCount[i];
		}
		LOG << "Maximum 1-4 pairs per atom is " << pairsListsData.maxPairs14PerAtom << ". Total is " << count14/2 << " (Note that explicit exclusions are taken into account).";
	} else {
		LOG << "Using pre-defined 1-4 listsize of " << pairsListsData.maxPairs14PerAtom;
	}


	LOG << "Constructing a 1-4 pairs lists...";
	allocateCPU((void**)&pairsListsData.h_pairs14List,
				gsystem.widthTot*pairsListsData.maxPairs14PerAtom*sizeof(int));

	for(i = 0; i < gsystem.N; i++){
		pairsListsData.h_pairs14ListCount[i] = 0;
	}

	for(i = 0; i < gsystem.N; i++){
		for(b = 0; b < pairsListsData.h_pairs12ListCount[i]; b++){
			j = pairsListsData.h_pairs12List[b*gsystem.widthTot + i];
			for(a = 0; a < pairsListsData.h_pairs13ListCount[j]; a++){
				k = pairsListsData.h_pairs13List[a*gsystem.widthTot + j];
				if(!checkAtomsBonded(i, k) && !checkAtoms13Bonded(i, k)){
					//if(!(checkAtomForAromaticExclusions(atoms[i]) && checkAtomForAromaticExclusions(atoms[k]))){
					/*if((checkAtomForAromaticExclusions(atoms[i]) && checkAtomForAromaticExclusions(atoms[k]))){
						if(!checkAtomsForExplicitExclusion(i,k)){
							printf("%s%d %s - %s%d %s\n",
									atoms[i].resName, atoms[i].resid, atoms[i].name,
									atoms[k].resName, atoms[k].resid, atoms[k].name);
						}
					}*/
					if(!checkAtomsForExplicitExclusion(i, k)){
						addAtomToListIfNotYetThere(i, k, pairsListsData.h_pairs14ListCount,
									pairsListsData.h_pairs14List, pairsListsData.maxPairs14PerAtom, gsystem.widthTot);
						/*pairsListsData.h_pairs14List[pairsListsData.h_pairs14ListCount[i]*gsystem.widthTot + i] = k;
						pairsListsData.h_pairs14ListCount[i] ++;*/
					}
				}
			}
		}
	}

	count14 = 0;
	for(i = 0; i < gsystem.N; i++){
		count14 += pairsListsData.h_pairs14ListCount[i];
	}

	LOG << "Done constructing list of 1-4 pairs.";




	allocateCPU((void**)&pairsListsData.h_pairsExclusionListCount, gsystem.Ntot*sizeof(int));
	pairsListsData.maxExcludedPerAtom = getIntegerParameter(PARAMETER_MAX_PAIRS_EXCLUDED, 0);
	if(pairsListsData.maxExcludedPerAtom == 0){

		for(i = 0; i < gsystem.N; i++){
			pairsListsData.h_pairsExclusionListCount[i] = 0;
		}

		/*for(i = 0; i < gsystem.N; i++){
			for(b = 0; b < pairsListsData.h_pairs12ListCount[i]; b++){
				j = pairsListsData.h_pairs12List[b*gsystem.widthTot + i];
				for(a = 0; a < pairsListsData.h_pairs13ListCount[j]; a++){
					k = pairsListsData.h_pairs13List[a*gsystem.widthTot + j];
					if(!checkAtomsBonded(i, k) && !checkAtoms13Bonded(i, k) && !checkAtoms14Bonded(i, k)){
						//if(checkAtomForAromaticExclusions(atoms[i]) && checkAtomForAromaticExclusions(atoms[k])){
						if(checkAtomsForExplicitExclusion(i, k)){
							pairsListsData.h_pairsExclusionListCount[i] ++;
							countExplicitExclusions ++;
						}
					}
				}
			}
		}*/
		for(b = 0; b < topology.exclusionsCount; b++){
			i = topology.exclusions[b].i;
			j = topology.exclusions[b].j;
			pairsListsData.h_pairsExclusionListCount[i] ++;
			pairsListsData.h_pairsExclusionListCount[j] ++;
			countExplicitExclusions ++;
		}

		for(i = 0; i < gsystem.N; i++){
			pairsListsData.h_pairsExclusionListCount[i] ++;
			pairsListsData.h_pairsExclusionListCount[i] += pairsListsData.h_pairs12ListCount[i];
			pairsListsData.h_pairsExclusionListCount[i] += pairsListsData.h_pairs13ListCount[i];
			pairsListsData.h_pairsExclusionListCount[i] += pairsListsData.h_pairs14ListCount[i];
		}

		pairsListsData.maxExcludedPerAtom = 0;
		for(i = 0; i < gsystem.N; i++){
			if(pairsListsData.maxExcludedPerAtom < pairsListsData.h_pairsExclusionListCount[i]){
				pairsListsData.maxExcludedPerAtom = pairsListsData.h_pairsExclusionListCount[i];
			}
			countExcluded += pairsListsData.h_pairsExclusionListCount[i];
		}
		LOG << "Maximum excluded pairs per atom is estimated as " << pairsListsData.maxExcludedPerAtom;
	} else {
		LOG << "Using pre-defined exclusions list size of " << pairsListsData.maxExcludedPerAtom;
	}

	allocateCPU((void**)&pairsListsData.h_pairsExclusionList,
						gsystem.widthTot*pairsListsData.maxExcludedPerAtom*sizeof(int));

	for(i = 0; i < gsystem.N; i++){
		pairsListsData.h_pairsExclusionListCount[i] = 0;
	}

	/*for(i = 0; i < gsystem.N; i++){
		for(b = 0; b < pairsListsData.h_pairs12ListCount[i]; b++){
			j = pairsListsData.h_pairs12List[b*gsystem.widthTot + i];
			for(a = 0; a < pairsListsData.h_pairs13ListCount[j]; a++){
				k = pairsListsData.h_pairs13List[a*gsystem.widthTot + j];
				if(!checkAtomsBonded(i, k) && !checkAtoms13Bonded(i, k) && !checkAtoms14Bonded(i, k)){
					//if(checkAtomForAromaticExclusions(atoms[i]) && checkAtomForAromaticExclusions(atoms[k])){
					if(checkAtomsForExplicitExclusion(i, k)){
						addAtomToListIfNotYetThere(i, k, pairsListsData.h_pairsExclusionListCount,
								pairsListsData.h_pairsExclusionList, pairsListsData.maxExcludedPerAtom);
					}
				}
			}
		}
	}*/
	countExplicitExclusions = 0;
	for(b = 0; b < topology.exclusionsCount; b++){
		i = topology.exclusions[b].i;
		j = topology.exclusions[b].j;
		if(addAtomToListIfNotYetThere(i, j,
				pairsListsData.h_pairsExclusionListCount,
				pairsListsData.h_pairsExclusionList, pairsListsData.maxExcludedPerAtom, gsystem.widthTot)){
			countExplicitExclusions ++;
		}
	}

	for(i = 0; i < gsystem.N; i++){
		addAtomToListIfNotYetThere(i, i, pairsListsData.h_pairsExclusionListCount,
					pairsListsData.h_pairsExclusionList, pairsListsData.maxExcludedPerAtom, gsystem.widthTot);
		for(b = 0; b < pairsListsData.h_pairs12ListCount[i]; b++){
			j = pairsListsData.h_pairs12List[b*gsystem.widthTot + i];
			addAtomToListIfNotYetThere(i, j, pairsListsData.h_pairsExclusionListCount,
					pairsListsData.h_pairsExclusionList, pairsListsData.maxExcludedPerAtom, gsystem.widthTot);
		}
		for(a = 0; a < pairsListsData.h_pairs13ListCount[i]; a++){
			j = pairsListsData.h_pairs13List[a*gsystem.widthTot + i];
			addAtomToListIfNotYetThere(i, j, pairsListsData.h_pairsExclusionListCount,
					pairsListsData.h_pairsExclusionList, pairsListsData.maxExcludedPerAtom, gsystem.widthTot);
		}
		for(d = 0; d < pairsListsData.h_pairs14ListCount[i]; d++){
			j = pairsListsData.h_pairs14List[d*gsystem.widthTot + i];
			addAtomToListIfNotYetThere(i, j, pairsListsData.h_pairsExclusionListCount,
					pairsListsData.h_pairsExclusionList, pairsListsData.maxExcludedPerAtom, gsystem.widthTot);
		}
	}

	countExcluded = 0;
	for(i = 0; i < gsystem.N; i++){
		countExcluded += pairsListsData.h_pairsExclusionListCount[i];
	}

	sortList(pairsListsData.h_pairsExclusionListCount, pairsListsData.h_pairsExclusionList,
			gsystem.N, gsystem.widthTot);

	LOG << "Done constructing list of excluded pairs.";



	allocateCPU((void**)&pairsListsData.h_pairsLJListCount, gsystem.Ntot*sizeof(int));
	allocateCPU((void**)&pairsListsData.h_pairsCoulombListCount, gsystem.Ntot*sizeof(int));
	allocateCPU((void**)&pairsListsData.h_possiblePairsListCount, gsystem.Ntot*sizeof(int));

	pairsListsData.maxPairsLJPerAtom = getIntegerParameter(PARAMETER_MAX_PAIRS_LJ, 0);
	pairsListsData.maxPairsCoulombPerAtom = getIntegerParameter(PARAMETER_MAX_PAIRS_COULOMB, 0);
	pairsListsData.maxPossiblePairsPerAtom = getIntegerParameter(PARAMETER_MAX_POSSIBLEPAIRS, 0);
	if((pairsListsData.maxPairsLJPerAtom != 0 ||
		pairsListsData.maxPairsCoulombPerAtom != 0 ||
		pairsListsData.maxPossiblePairsPerAtom != 0) &&
		(pairsListsData.maxPairsLJPerAtom *
				pairsListsData.maxPairsCoulombPerAtom *
				pairsListsData.maxPossiblePairsPerAtom == 0)){
		LOG << "Sizes for all changeable pairlists (LJ, Coulomb and PossiblePairs) should be defined. "
				"Re-counting will be forced.";
	}
	if(pairsListsData.maxPairsLJPerAtom == 0 ||
			pairsListsData.maxPairsCoulombPerAtom == 0 ||
			pairsListsData.maxPossiblePairsPerAtom == 0){

		LOG << "Counting members of current LJ/Coulomb Verlet and pairs lists...";

		for(i = 0; i < gsystem.Ntot; i++){
			pairsListsData.h_pairsLJListCount[i] = 0;
			pairsListsData.h_pairsCoulombListCount[i] = 0;
			pairsListsData.h_possiblePairsListCount[i] = 0;
		}

		float4 r;
		for(traj = 0; traj < parameters.Ntr; traj++){
			for(i = 0; i < gsystem.N; i++){
				for(j = 0; j < i; j++){
					itot = i + traj*gsystem.N;
					jtot = j + traj*gsystem.N;
					if(!checkAtomsInExclusionList(i, j)){
						r.x = gsystem.h_coord[itot].x - gsystem.h_coord[jtot].x;
						r.y = gsystem.h_coord[itot].y - gsystem.h_coord[jtot].y;
						r.z = gsystem.h_coord[itot].z - gsystem.h_coord[jtot].z;
						DO_PBC(r);
						r.w = sqrtf(r.x*r.x + r.y*r.y + r.z*r.z);
						if(r.w < pairsListsData.possiblePairsCutoff){
							pairsListsData.h_possiblePairsListCount[itot]++;
							pairsListsData.h_possiblePairsListCount[jtot]++;
							countPP++;
						}
						if(r.w < pairsListsData.coulombPairsCutoff){
							if(r.w < pairsListsData.ljPairsCutoff){
								pairsListsData.h_pairsLJListCount[itot]++;
								pairsListsData.h_pairsLJListCount[jtot]++;
								countLJ++;
							} else {
								pairsListsData.h_pairsCoulombListCount[itot]++;
								pairsListsData.h_pairsCoulombListCount[jtot]++;
								countCoulomb++;
							}
						}
					}
				}
			}
		}

		pairsListsData.maxPairsLJPerAtom = 0;
		pairsListsData.maxPairsCoulombPerAtom = 0;
		pairsListsData.maxPossiblePairsPerAtom = 0;
		for(i = 0; i < gsystem.Ntot; i++){
			if(pairsListsData.maxPairsLJPerAtom < pairsListsData.h_pairsLJListCount[i]){
				pairsListsData.maxPairsLJPerAtom = pairsListsData.h_pairsLJListCount[i];
			}
			if(pairsListsData.maxPairsCoulombPerAtom < pairsListsData.h_pairsCoulombListCount[i]){
				pairsListsData.maxPairsCoulombPerAtom = pairsListsData.h_pairsCoulombListCount[i];
			}
			if(pairsListsData.maxPossiblePairsPerAtom < pairsListsData.h_possiblePairsListCount[i]){
				pairsListsData.maxPossiblePairsPerAtom = pairsListsData.h_possiblePairsListCount[i];
			}
		}

		LOG << "Done counting members of current LJ/Coulomb Verlet and pairs lists.";
	} else {
		LOG << "Using pre-defined LJ, Coulomb and Possible Pairs list sizes of " << pairsListsData.maxPairsLJPerAtom << ", " << pairsListsData.maxPairsCoulombPerAtom << " and " << pairsListsData.maxPossiblePairsPerAtom;
		LOG << "Actual numbers of pairs in changeable lists were not yet computed and will be zeros in following table.";
		pairsListsData.maxPairsCoulombPerAtom -= pairsListsData.maxPairsLJPerAtom;
	}
	
	const int tw = 10; // Width of table field
	LOG << "";
	LOG.table(tw,3) << "Pairs" << "Total" << "Max/Atom";
	LOG.table(tw,3) << "1-2" << count12/2 << pairsListsData.maxPairs12PerAtom;
	LOG.table(tw,3) << "1-3" << count13/2 << pairsListsData.maxPairs13PerAtom;
	LOG.table(tw,3) << "1-4" << count14/2 << pairsListsData.maxPairs14PerAtom;
	LOG.table(tw,3) << "Excluded" << (countExcluded-gsystem.Nsim)/2  << pairsListsData.maxExcludedPerAtom << " (including " << countExplicitExclusions << " explicit exclusions)";
	LOG.table(tw,3) << "LJ" << countLJ << pairsListsData.maxPairsLJPerAtom;
	LOG.table(tw,3) << "Coulomb" << countCoulomb+countLJ << pairsListsData.maxPairsLJPerAtom+pairsListsData.maxPairsCoulombPerAtom;
	LOG.table(tw,3) << "Pairslist" << countPP << pairsListsData.maxPossiblePairsPerAtom;
	LOG << "";
	LOG << "Total number of changeable lists has been counted for all " << parameters.Ntr << " trajectories";
	LOG << "Adding " << pairlistsExtension << " to the length of all changeable lists to avoid memory conflicts";
	pairsListsData.maxPairsLJPerAtom += pairlistsExtension;
	pairsListsData.maxPairsCoulombPerAtom += pairlistsExtension;
	pairsListsData.maxPossiblePairsPerAtom += pairlistsExtension;

	allocateGPU((void**)&pairsListsData.d_pairs12ListCount, gsystem.Ntot*sizeof(int));
	allocateGPU((void**)&pairsListsData.d_pairs12List,
			pairsListsData.maxPairs12PerAtom*gsystem.widthTot*sizeof(int));

	allocateGPU((void**)&pairsListsData.d_pairs13ListCount, gsystem.Ntot*sizeof(int));
	allocateGPU((void**)&pairsListsData.d_pairs13List,
			pairsListsData.maxPairs13PerAtom*gsystem.widthTot*sizeof(int));

	allocateGPU((void**)&pairsListsData.d_pairs14ListCount, gsystem.Ntot*sizeof(int));
	allocateGPU((void**)&pairsListsData.d_pairs14List,
			pairsListsData.maxPairs14PerAtom*gsystem.widthTot*sizeof(int));

	allocateGPU((void**)&pairsListsData.d_pairsExclusionListCount, gsystem.Ntot*sizeof(int));
	allocateGPU((void**)&pairsListsData.d_pairsExclusionList,
			pairsListsData.maxExcludedPerAtom*gsystem.widthTot*sizeof(int));

	allocateGPU((void**)&pairsListsData.d_pairsLJListCount, gsystem.Ntot*sizeof(int));
	allocateCPU((void**)&pairsListsData.h_pairsLJList,
			pairsListsData.maxPairsLJPerAtom*gsystem.widthTot*sizeof(int));
	allocateGPU((void**)&pairsListsData.d_pairsLJList,
			pairsListsData.maxPairsLJPerAtom*gsystem.widthTot*sizeof(int));

	allocateGPU((void**)&pairsListsData.d_pairsCoulombListCount, gsystem.Ntot*sizeof(int));
	allocateCPU((void**)&pairsListsData.h_pairsCoulombList,
			pairsListsData.maxPairsCoulombPerAtom*gsystem.widthTot*sizeof(int));
	allocateGPU((void**)&pairsListsData.d_pairsCoulombList,
			pairsListsData.maxPairsCoulombPerAtom*gsystem.widthTot*sizeof(int));

	allocateGPU((void**)&pairsListsData.d_possiblePairsListCount, gsystem.Ntot*sizeof(int));
	allocateCPU((void**)&pairsListsData.h_possiblePairsList,
			pairsListsData.maxPossiblePairsPerAtom*gsystem.widthTot*sizeof(int));
	allocateGPU((void**)&pairsListsData.d_possiblePairsList,
			pairsListsData.maxPossiblePairsPerAtom*gsystem.widthTot*sizeof(int));

	for(traj = 1; traj < parameters.Ntr; traj++){
		for(i = 0; i < gsystem.N; i++){
			itot = traj*gsystem.N + i;
			pairsListsData.h_pairsExclusionListCount[itot] = pairsListsData.h_pairsExclusionListCount[i];
			pairsListsData.h_pairs12ListCount[itot] = pairsListsData.h_pairs12ListCount[i];
			pairsListsData.h_pairs13ListCount[itot] = pairsListsData.h_pairs13ListCount[i];
			pairsListsData.h_pairs14ListCount[itot] = pairsListsData.h_pairs14ListCount[i];
			pairsListsData.h_pairsLJListCount[itot] = pairsListsData.h_pairsLJListCount[i];
			pairsListsData.h_pairsCoulombListCount[itot] = pairsListsData.h_pairsCoulombListCount[i];
			pairsListsData.h_possiblePairsListCount[itot] = pairsListsData.h_possiblePairsListCount[i];
			for(p = 0; p < pairsListsData.maxExcludedPerAtom; p++){
				pairsListsData.h_pairsExclusionList[p*gsystem.widthTot + itot] =
						pairsListsData.h_pairsExclusionList[p*gsystem.widthTot + i] + traj*gsystem.N;
			}
			for(p = 0; p < pairsListsData.maxPairs12PerAtom; p++){
				pairsListsData.h_pairs12List[p*gsystem.widthTot + itot] =
						pairsListsData.h_pairs12List[p*gsystem.widthTot + i] + traj*gsystem.N;
			}
			for(p = 0; p < pairsListsData.maxPairs13PerAtom; p++){
				pairsListsData.h_pairs13List[p*gsystem.widthTot + itot] =
						pairsListsData.h_pairs13List[p*gsystem.widthTot + i] + traj*gsystem.N;
			}
			for(p = 0; p < pairsListsData.maxPairs14PerAtom; p++){
				pairsListsData.h_pairs14List[p*gsystem.widthTot + itot] =
						pairsListsData.h_pairs14List[p*gsystem.widthTot + i] + traj*gsystem.N;
			}
			for(p = 0; p < pairsListsData.maxPairsLJPerAtom; p++){
				pairsListsData.h_pairsLJList[p*gsystem.widthTot + itot] =
						pairsListsData.h_pairsLJList[p*gsystem.widthTot + i] + traj*gsystem.N;
			}
			for(p = 0; p < pairsListsData.maxPairsCoulombPerAtom; p++){
				pairsListsData.h_pairsCoulombList[p*gsystem.widthTot + itot] =
						pairsListsData.h_pairsCoulombList[p*gsystem.widthTot + i] + traj*gsystem.N;
			}
			for(p = 0; p < pairsListsData.maxPossiblePairsPerAtom; p++){
				pairsListsData.h_possiblePairsList[p*gsystem.widthTot + itot] =
						pairsListsData.h_possiblePairsList[p*gsystem.widthTot + i] + traj*gsystem.N;
			}
		}
	}


	//int j;
	/*printf("\n\nList of excluded pairs:\n");
	for(i = 0; i < gsystem.N; i++){
		printf("%d: ", i);
		for(j = 0; j < pairsListsData.h_pairsExclusionListCount[i]; j++){
			printf("%d ",
					pairsListsData.h_pairsExclusionList[j*gsystem.widthTot + i]);
		}
		printf("\n");
	}
	exit(0);*/
	/*printf("\n\nList of 1-4 pairs:\n");
	for(i = 0; i < gsystem.N; i++){
		printf("%d: ", i);
		for(j = 0; j < pairsListsData.h_pairs14ListCount[i]; j++){
			printf("%d ",
					pairsListsData.h_pairs14List[j*gsystem.widthTot + i]);
		}
		printf("\n");
	}*/

	cudaMemcpy(pairsListsData.d_pairs12ListCount, pairsListsData.h_pairs12ListCount,
				gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(pairsListsData.d_pairs12List, pairsListsData.h_pairs12List,
				gsystem.widthTot*pairsListsData.maxPairs12PerAtom*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(pairsListsData.d_pairs13ListCount, pairsListsData.h_pairs13ListCount,
				gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(pairsListsData.d_pairs13List, pairsListsData.h_pairs13List,
				gsystem.widthTot*pairsListsData.maxPairs13PerAtom*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(pairsListsData.d_pairs14ListCount, pairsListsData.h_pairs14ListCount,
				gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(pairsListsData.d_pairs14List, pairsListsData.h_pairs14List,
				gsystem.widthTot*pairsListsData.maxPairs14PerAtom*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(pairsListsData.d_pairsExclusionListCount, pairsListsData.h_pairsExclusionListCount,
				gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(pairsListsData.d_pairsExclusionList, pairsListsData.h_pairsExclusionList,
				gsystem.widthTot*pairsListsData.maxExcludedPerAtom*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpyToSymbol(c_pairsListsData, &pairsListsData,
				sizeof(PairsListsData), 0, cudaMemcpyHostToDevice);

	LOG << "Done initializing pairs lists";
}

__global__ void updatePairsLists_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 r1 = tex1Dfetch(t_coord, d_i);
		float4 r2;
		int i, j;
		int ljCount = 0;
		int coulombCount = 0;
		for(i = 0; i < c_pairsListsData.d_possiblePairsListCount[d_i]; i++){
			j = c_pairsListsData.d_possiblePairsList[i*c_gsystem.widthTot + d_i];
			r2 = tex1Dfetch(t_coord, j);
			r2 -= r1;
			DO_PBC(r2);
			r2.w = sqrtf(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);
			if(r2.w < c_pairsListsData.ljPairsCutoff){
				c_pairsListsData.d_pairsLJList[ljCount*c_gsystem.widthTot + d_i] = j;
				ljCount++;
			} else if(r2.w < c_pairsListsData.coulombPairsCutoff){
				c_pairsListsData.d_pairsCoulombList[coulombCount*c_gsystem.widthTot + d_i] = j;
				coulombCount++;
			}
		}
		c_pairsListsData.d_pairsLJListCount[d_i] = ljCount;
		c_pairsListsData.d_pairsCoulombListCount[d_i] = coulombCount;
	}
}

void updatePairsLists(){
	checkCUDAError("before Verlet list");
	updatePairsLists_kernel<<<pairsListsBlockCount, pairsListsBlockSize>>>();
	checkCUDAError("Verlet list, before copying data");

	/*cudaMemcpy(pairsListsData.h_pairsLJListCount, pairsListsData.d_pairsLJListCount,
				gsystem.Ntot*sizeof(int), cudaMemcpyDeviceToHost);

	cudaMemcpy(pairsListsData.h_pairsCoulombListCount, pairsListsData.d_pairsCoulombListCount,
						gsystem.Ntot*sizeof(int), cudaMemcpyDeviceToHost);

	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		for(i = 0; i < gsystem.N; i++){
			if(pairsListsData.h_pairsLJListCount[i] > pairsListsData.maxPairsLJPerAtom){
				DIE("ERROR: Number of LJ pairs on atom %d (trajectory %d) is %d, which exceeds the limit of %d.",
						i, traj+parameters.firstrun,
						pairsListsData.h_pairsLJListCount[i], pairsListsData.maxPairsLJPerAtom);
			}
			if(pairsListsData.h_pairsCoulombListCount[i] > pairsListsData.maxPairsCoulombPerAtom){
				DIE("ERROR: Number of Coulomb pairs on atom %d (trajectory %d) is %d, which exceeds the limit of %d.",
						i, traj+parameters.firstrun,
						pairsListsData.h_pairsCoulombListCount[i], pairsListsData.maxPairsCoulombPerAtom);
			}
		}
	}*/
	/*int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		int totalPairsLJ = 0;
		int totalPairsCoulomb = 0;
		for(i = 0; i < gsystem.N; i++){
			totalPairsLJ += pairsListsData.h_pairsLJListCount[gsystem.N*traj + i];
			totalPairsCoulomb += pairsListsData.h_pairsCoulombListCount[gsystem.N*traj + i];
		}
		printf("Total pairs for trajectory #%d: %d(LJ), %d(Coulomb)\n",
				traj, totalPairsLJ, totalPairsCoulomb);
	}*/


	/*cudaMemcpy(pairsListsData.h_pairsLJList, pairsListsData.d_pairsLJList,
				gsystem.widthTot*pairsListsData.maxPairsLJPerAtom*sizeof(int), cudaMemcpyDeviceToHost);


	cudaMemcpy(pairsListsData.h_pairsCoulombList, pairsListsData.d_pairsCoulombList,
					gsystem.widthTot*pairsListsData.maxPairsCoulombPerAtom*sizeof(int), cudaMemcpyDeviceToHost);

	int j;
	int totalPairs = 0;
	printf("\n\nList of LJ pairs:\n");
	for(i = 0; i < gsystem.Ntot; i++){
		printf("%d (%d): ", i, pairsListsData.h_pairsLJListCount[i]);
		for(j = 0; j < pairsListsData.h_pairsLJListCount[i]; j++){
			printf("%d ",
					pairsListsData.h_pairsLJList[j*gsystem.widthTot + i]);
		}
		totalPairs += pairsListsData.h_pairsLJListCount[i];
		printf("\n");
	}
	printf("Total: %d\n", totalPairs);
	printf("\n\nList of Coulomb pairs:\n");
	for(i = 0; i < gsystem.Ntot; i++){
		printf("%d (%d): ", i, pairsListsData.h_pairsCoulombListCount[i]);
		for(j = 0; j < pairsListsData.h_pairsCoulombListCount[i]; j++){
			printf("%d ",
					pairsListsData.h_pairsCoulombList[j*gsystem.widthTot + i]);
		}
		printf("\n");
	}*/
	//exit(0);
}

__global__ void updatePossiblePairsLists_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		int i;
		float4 r1 = tex1Dfetch(t_coord, d_i);
		float4 r2;
		int possiblePairsCount = 0;
		int traj = d_i/c_gsystem.N;
		int currentExclusion = c_pairsListsData.d_pairsExclusionList[d_i];;
		int currentExclusionID = 0;
		int totalExcluded = c_pairsListsData.d_pairsExclusionListCount[d_i];
		for(i = traj*c_gsystem.N; i < (traj + 1)*c_gsystem.N; i++){
			if(i != currentExclusion){
				r2 = tex1Dfetch(t_coord, i);
				r2 -= r1;
				DO_PBC(r2);
				r2.w = sqrtf(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);
				if(r2.w < c_pairsListsData.possiblePairsCutoff){
					c_pairsListsData.d_possiblePairsList[possiblePairsCount*c_gsystem.widthTot + d_i] = i;
					possiblePairsCount++;
				}
			} else {
				currentExclusionID ++;
				if(currentExclusionID == totalExcluded){
					currentExclusionID --;
				}
				currentExclusion = c_pairsListsData.d_pairsExclusionList[currentExclusionID*c_gsystem.widthTot + d_i];
			}
		}
		c_pairsListsData.d_possiblePairsListCount[d_i] = possiblePairsCount;
	}
}


void updatePossiblePairsList(){
	checkCUDAError("before possible pairs");
#ifdef CUDA_USE_L1
	cudaFuncSetCacheConfig(updatePossiblePairsLists_kernel,cudaFuncCachePreferL1);
#endif
	updatePossiblePairsLists_kernel<<<possiblePairsListBlockCount, possiblePairsListBlockSize>>>();
	checkCUDAError("possible pairs before copying data");

/*	cudaMemcpy(pairsListsData.h_possiblePairsListCount, pairsListsData.d_possiblePairsListCount,
					gsystem.Ntot*sizeof(int), cudaMemcpyDeviceToHost);*/

	/*int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		for(i = 0; i < gsystem.N; i++){
			if(pairsListsData.h_possiblePairsListCount[i] > pairsListsData.maxPossiblePairsPerAtom){
				DIE("ERROR: Number of possible pairs on atom %d (trajectory %d) is %d, which exceeds the limit of %d.",
						i, traj+parameters.firstrun,
						pairsListsData.h_possiblePairsListCount[i], pairsListsData.maxPossiblePairsPerAtom);
			}
		}
	}*/

/*	cudaMemcpy(pairsListsData.h_possiblePairsList, pairsListsData.d_possiblePairsList,
					gsystem.Ntot*pairsListsData.maxPossiblePairsPerAtom*sizeof(int), cudaMemcpyDeviceToHost);*/

/*	int i, j;
	//printf("\n\nList of possible pairs:\n");
	FILE* ppfile = fopen("pp.txt", "w");
	for(i = 0; i < gsystem.Ntot; i++){
		fprintf(ppfile, "%d (%d):", i, pairsListsData.h_possiblePairsListCount[i]);
		for(j = 0; j < pairsListsData.h_possiblePairsListCount[i]; j++){
			fprintf(ppfile, "%d ",
					pairsListsData.h_possiblePairsList[j*gsystem.widthTot + i]);
		}
		fprintf(ppfile, "\n");
	}
	fclose(ppfile);
	exit(0);*/
}

void destroyPairsListsUpdater(){

}

void destroyPossiblePairsListUpdater(){

}

int isHeavySidechainAtom(Atom atom){
	if(atom.name[0] == 'H'){
		return 0;
	}
	if(strcmp(atom.name, "H") == 0 ||
			strcmp(atom.name, "C") == 0 ||
			strcmp(atom.name, "CA") == 0 ||
			strcmp(atom.name, "N") == 0 ||
			strcmp(atom.name, "O") == 0){
		return 0;
	}
	return 1;
}

int checkAtomForAromaticExclusions(Atom atom){
	//if(pairsListsData.checkAromaticExclusions){
		if(strcmp(atom.resName, "PHE") == 0){
			/*if(isHeavySidechainAtom(atom) && strcmp(atom.name, "CB") != 0){
				return 1;
			}*/
			if(strcmp(atom.name, "CD1") == 0 ||
					strcmp(atom.name, "CD2") == 0 ||
					strcmp(atom.name, "CE1") == 0 ||
					strcmp(atom.name, "CE2") == 0 ||
					strcmp(atom.name, "CG") == 0 ||
					strcmp(atom.name, "CZ") == 0){
				/*printf("Excluding atom %d as %s side-chain atom: %s (%s)\n",
						atom.id, atom.resName, atom.name, atom.type);*/
				return 1;
			} else {
				return 0;
			}

		}
		if(strcmp(atom.resName, "TYR") == 0){
			/*if(isHeavySidechainAtom(atom) && strcmp(atom.name, "CB") != 0 && strcmp(atom.name, "OH") != 0){
				return 1;
			}*/
			if(strcmp(atom.name, "CD1") == 0 ||
					strcmp(atom.name, "CD2") == 0 ||
					strcmp(atom.name, "CE1") == 0 ||
					strcmp(atom.name, "CE2") == 0 ||
					strcmp(atom.name, "CG") == 0 ||
					strcmp(atom.name, "CZ") == 0){
				/*printf("Excluding atom %d as %s side-chain atom: %s (%s)\n",
						atom.id, atom.resName, atom.name, atom.type);*/
				return 1;
			} else {
				return 0;
			}
		}
		if(strcmp(atom.resName, "TRP") == 0){
			/*if(isHeavySidechainAtom(atom) && strcmp(atom.name, "CB") != 0){
				return 1;
			}*/
			if(strcmp(atom.name, "CD2") == 0 ||
					strcmp(atom.name, "CE2") == 0 ||
					strcmp(atom.name, "CE3") == 0 ||
					strcmp(atom.name, "CH2") == 0 ||
					strcmp(atom.name, "CZ2") == 0 ||
					strcmp(atom.name, "CZ3") == 0){
				/*printf("Excluding atom %d as %s side-chain atom: %s (%s)\n",
						atom.id, atom.resName, atom.name, atom.type);*/
				return 1;
			} else {
				return 0;
			}
		}
		/*if(strcmp(atom.resName, "HSD") == 0 ||
				strcmp(atom.resName, "HSE") == 0 ||
				strcmp(atom.resName, "HIS") == 0 ||
				strcmp(atom.resName, "HSC") == 0){
			if(isHeavySidechainAtom(atom) && strcmp(atom.name, "CB") != 0){*/
				/*printf("Excluding atom %d as %s side-chain atom: %s (%s)\n",
						atom.id, atom.resName, atom.name, atom.type);*/
		/*		return 1;
			} else {
				return 0;
			}
		}
		if(strcmp(atom.resName, "PRO") == 0){
			if(strcmp(atom.name, "CA") == 0 ||
					strcmp(atom.name, "CB") == 0 ||
					strcmp(atom.name, "CG") == 0 ||
					strcmp(atom.name, "CD") == 0 ||
					strcmp(atom.name, "N") == 0){*/
			/*	printf("Excluding atom %d as %s ring member: %s (%s)\n",
						atom.id, atom.resName, atom.name, atom.type);*/
		/*		return 1;
			} else {
				return 0;
			}
		}*/
		return 0;
	/*} else {
		return 0;
	}*/
}

int checkAtomsForExplicitExclusion(int a1, int a2){
	int i;
	for(i = 0; i < topology.exclusionsCount; i++){
		if((a1 == topology.exclusions[i].i && a2 == topology.exclusions[i].j) ||
				(a1 == topology.exclusions[i].j && a2 == topology.exclusions[i].i)){
					return 1;
				}
	}
	return 0;
}


int checkAtomsBonded(int i, int j){
	int b;
	for(b = 0; b < pairsListsData.h_pairs12ListCount[i]; b++){
		if(pairsListsData.h_pairs12List[b*gsystem.widthTot + i] == j){
			return 1;
		}
	}
	return 0;
}

int checkAtoms13Bonded(int i, int j){
	int a;
	for(a = 0; a < pairsListsData.h_pairs13ListCount[i]; a++){
		if(pairsListsData.h_pairs13List[a*gsystem.widthTot + i] == j){
			return 1;
		}
	}
	return 0;
}

int checkAtoms14Bonded(int i, int j){
	int d;
	for(d = 0; d < pairsListsData.h_pairs14ListCount[i]; d++){
		if(pairsListsData.h_pairs14List[d*gsystem.widthTot + i] == j){
			return 1;
		}
	}
	return 0;
}

int checkAtomsInExclusionList(int i, int j){
	int d;
	for(d = 0; d < pairsListsData.h_pairsExclusionListCount[i]; d++){
		if(pairsListsData.h_pairsExclusionList[d*gsystem.widthTot + i] == j){
			return 1;
		}
	}
	return 0;
}

int addAtomToListIfNotYetThere(int i, int j, int* counts, int* list, int maxCount, int width){
	int k = 0;
	int found = 0;
	while(k < counts[i] && found == 0){
		if(list[k*width + i] == j){
			found = 1;
		}
		k++;
	}
	if(found == 0){
		list[counts[i]*width + i] = j;
		counts[i]++;
		if(i != j){
			list[counts[j]*width+ j] = i;
			counts[j]++;
		}
		if(counts[i] > maxCount || counts[j] > maxCount){
			DIE("ERROR: Number of pairs in a list is larger then a maximum of %d.", maxCount);
		}
		return 1;
	} else {
		return 0;
	}
}

void sortList(int* counts, int* list, int N, int width){
	int i, j, k;
	int maxim = 0;
	int pos;
	//Find maximum count
	for(i = 0; i < N; i++){
		if(counts[i] > maxim){
			maxim = counts[i];
		}
	}
	//Allocate memory for temporary array
	int* temp = (int*)calloc(maxim, sizeof(int));
	// Sorting
	for(i = 0; i < N; i++){
		//Copy data
		for(j = 0; j < counts[i]; j++){
			temp[j] = list[j*width + i];
		}
		for(j = 0; j < counts[i]; j++){
			pos = 0;
			for(k = 0; k < counts[i]; k++){
				if(temp[k] < temp[j]){
					pos++;
				}
			}
			list[pos*width + i] = temp[j];
		}
	}
	free(temp);
}

#undef LOG

} // namespace pairslist_updater
