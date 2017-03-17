/*
 * SASAPotential.cu
 *
 *  Created on: Aug 6, 2010
 *      Author: zhmurov
 */
#include "../Core/global.h"
#include "../Core/md.cuh"
#include "sasa.h"
#include "SASAPotential.cuh"

namespace sasa_potential {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<sasa_potential> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

//__global__ void computeSASAPotentialB_kernel();

void create(){
	if (!getYesNoParameter(PARAMETER_SASA_ON, 0))
		return;
	potential.compute = &compute;
	potential.destroy = &destroy;
	sprintf(potential.name, "SASA potential");
	potentials[potentialsCount] = &potential;
	potentialsCount ++;

	sasaPairsListUpdater.update = updateSASAPairsList;
	sasaPairsListUpdater.destroy = destroySASAPairsListUpdater;
	sasaPairsListUpdater.frequency = getIntegerParameter(PARAMETER_SASA_PAIRS_FREQ);
	sprintf(sasaPairsListUpdater.name, "SASA updater");
	updaters[updatersCount] = &sasaPairsListUpdater;
	updatersCount ++;

	energyOutput.computeValues = &computeEnergy;
	allocateCPU((void**)&energyOutput.values, parameters.Ntr*sizeof(float));
	strcpy(energyOutput.name, ENERGY_OUTPUT_NAME_SASA);
	energyOutputs[energyOutputsCount] = &energyOutput;
	energyOutputsCount ++;

	//if(deviceProps.major == 2){
		//cudaFuncSetCacheConfig(computeSASAPotentialB_kernel, cudaFuncCachePreferL1);
	//}

	init();
}


void init(){
	LOG << "Initializing SASA implicit solvent potential...";

	int i, j;

	pairlistsExtension = getIntegerParameter(PARAMETER_PAIRSLISTS_EXTENSION, DEFAULT_PAIRSLISTS_EXTENSION);

	sasaData.Rprobe = getFloatParameter(PARAMETER_SASA_RPROBE, DEFAULT_SASA_RPROBE);
	sasaData.pij_cov = getFloatParameter(PARAMETER_SASA_PIJ_COV, DEFAULT_SASA_PIJ_COV);
	sasaData.pij_nb = getFloatParameter(PARAMETER_SASA_PIJ_NB, DEFAULT_SASA_PIJ_NB);

	sasaData.sasaPairsListCutoff = getFloatParameter(PARAMETER_SASA_PAIRSLIST_CUTOFF);

	allocateCPU((void**)&sasaData.h_threadAtom, gsystem.Ntot*sizeof(int));
	allocateGPU((void**)&sasaData.d_threadAtom, gsystem.Ntot*sizeof(int));
	allocateCPU((void**)&sasaData.h_atomThread, gsystem.Ntot*sizeof(int));
	allocateGPU((void**)&sasaData.d_atomThread, gsystem.Ntot*sizeof(int));

	int hydrogenCount = 0;
	sasaData.threadsCount = 0;
	for(i = 0; i < gsystem.N; i++){
		if(topology.atoms[i].type[0] == 'H'){
			hydrogenCount ++;
			sasaData.h_atomThread[i] = -1;
		} else {
			sasaData.h_atomThread[i] = sasaData.threadsCount;
			sasaData.h_threadAtom[sasaData.threadsCount] = i;
			sasaData.threadsCount ++;
		}
	}
/*	for(i = 0; i < sasaData.threadsCount; i++){
		printf("%d: %s - %d\n", i,
				topology.atoms[sasaData.h_threadAtom[i]].name,
				sasaData.h_threadAtom[i]);
	}*/
	sasaData.threadsCountTot = sasaData.threadsCount*parameters.Ntr;

	sasaBlockSize = BLOCK_SIZE;
	sasaBlockCount = sasaData.threadsCountTot/BLOCK_SIZE + 1;
	sasaPairsListBlockSize = BLOCK_SIZE;
	sasaPairsListBlockCount = sasaData.threadsCountTot/BLOCK_SIZE + 1;

	if(sasaData.threadsCountTot % 16 == 0){
		sasaData.widthTot = sasaData.threadsCountTot;
	} else {
		sasaData.widthTot = (sasaData.threadsCountTot/16)*16 + 16;
	}

	LOG << "Found " << hydrogenCount << " hydrogens, " << gsystem.N - hydrogenCount << " non-hydrogens";
	LOG << "Arrays will be aligned to the width of " <<  sasaData.widthTot;

	//allocateCPU((void**)&sasaData.h_sasaParameters, atomTypesCount*sizeof(GSASAParameters));
	//allocateGPU((void**)&sasaData.d_sasaParameters, atomTypesCount*sizeof(GSASAParameters));


	allocateCPU((void**)&sasaData.h_sasaRipr, atomTypesCount*sizeof(float));
	allocateCPU((void**)&sasaData.h_sasaRipr2, atomTypesCount*sizeof(float));
	allocateCPU((void**)&sasaData.h_sasaPiOverSi, atomTypesCount*sizeof(float));
	allocateCPU((void**)&sasaData.h_sasaSigmaSi, atomTypesCount*sizeof(float));
	allocateCPU((void**)&sasaData.h_sasaSi, atomTypesCount*sizeof(float));

	allocateGPU((void**)&sasaData.d_sasaRipr, atomTypesCount*sizeof(float));
	allocateGPU((void**)&sasaData.d_sasaRipr2, atomTypesCount*sizeof(float));
	allocateGPU((void**)&sasaData.d_sasaPiOverSi, atomTypesCount*sizeof(float));
	allocateGPU((void**)&sasaData.d_sasaSigmaSi, atomTypesCount*sizeof(float));
	allocateGPU((void**)&sasaData.d_sasaSi, atomTypesCount*sizeof(float));

	char sasaFilename[100];
	getMaskedParameter(sasaFilename, PARAMETER_SASA_PARAMETERS_FILE);
	readSASAParameters(sasaFilename);
	for(i = 0; i < atomTypesCount; i++){
		SASAParameters param = getSASAParametersCHARMM(atomTypes[i].name);
		float Ripr = (param.R + sasaData.Rprobe);
		float Si = 4.0f*M_PI*Ripr*Ripr;
		//sasaData.h_sasaParameters[i].Ripr = Ripr;
		//sasaData.h_sasaParameters[i].Ripr2 = Ripr*Ripr;
		//sasaData.h_sasaParameters[i].piOverSi = param.p/Si;
		//sasaData.h_sasaParameters[i].sigmaSi = param.sigma*Si;
		sasaData.h_sasaRipr[i] = Ripr;//sasaData.h_sasaParameters[i].Ripr;
		sasaData.h_sasaRipr2[i] = Ripr*Ripr;//sasaData.h_sasaParameters[i].Ripr2;
		sasaData.h_sasaPiOverSi[i] = param.p/Si;//sasaData.h_sasaParameters[i].piOverSi;
		sasaData.h_sasaSigmaSi[i] = param.sigma*Si;//sasaData.h_sasaParameters[i].sigmaSi;
		sasaData.h_sasaSi[i] = Si;
	}

	allocateCPU((void**)&sasaData.h_pairs12Counts, sasaData.threadsCountTot*sizeof(int));
	allocateGPU((void**)&sasaData.d_pairs12Counts, sasaData.threadsCountTot*sizeof(int));

	allocateCPU((void**)&sasaData.h_pairsListCount, sasaData.threadsCountTot*sizeof(int));
	allocateGPU((void**)&sasaData.d_pairsListCount, sasaData.threadsCountTot*sizeof(int));

	allocateCPU((void**)&sasaData.h_sasaListCount, sasaData.threadsCountTot*sizeof(int));
	allocateGPU((void**)&sasaData.d_sasaListCount, sasaData.threadsCountTot*sizeof(int));

	LOG << "Counting members of current SASA pairlist...";
	int totalSASAPairlistCount = 0;
	int totalSASAListCount = 0;
	for(i = 0; i < sasaData.threadsCountTot; i++){
		sasaData.h_pairsListCount[i] = 0;
		sasaData.h_sasaListCount[i] = 0;
	}

	sasaData.maxSASAPairsPerAtom = getIntegerParameter(PARAMETER_MAX_PAIRS_SASA, 0);
	sasaData.maxPairsListItemsPerAtom = getIntegerParameter(PARAMETER_MAX_PAIRLIST_SASA, 0);
	if((sasaData.maxSASAPairsPerAtom != 0 ||
		sasaData.maxPairsListItemsPerAtom != 0) &&
		(sasaData.maxSASAPairsPerAtom *	sasaData.maxPairsListItemsPerAtom == 0)){
		LOG << "Sizes for both SASA lists should be defined. "
				"Re-counting will be forced.";
	}
	if(sasaData.maxSASAPairsPerAtom == 0 ||
			sasaData.maxPairsListItemsPerAtom == 0){
		float4 r1, r2;
		int a1, a2;
		for(i = 0; i < sasaData.threadsCount; i++){
			a1 = sasaData.h_threadAtom[i];
			if(topology.atoms[a1].type[0] != 'H'){
				r1 = gsystem.h_coord[a1];
				for(j = 0; j < i; j++){
					a2 = sasaData.h_threadAtom[j];
					if(topology.atoms[a2].type[0] != 'H'){
						r2 = gsystem.h_coord[a2];
						r2.x = r2.x - r1.x;
						r2.y = r2.y - r1.y;
						r2.z = r2.z - r1.z;
						r2.w = sqrtf(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);
						if(r2.w < sasaData.sasaPairsListCutoff){
							sasaData.h_pairsListCount[i] ++;
							sasaData.h_pairsListCount[j] ++;
							totalSASAPairlistCount ++;
							if(r2.w < sasaData.h_sasaRipr[gsystem.h_atomTypes[a1]] + sasaData.h_sasaRipr[gsystem.h_atomTypes[a2]]){
								sasaData.h_sasaListCount[i] ++;
								sasaData.h_sasaListCount[j] ++;
								totalSASAListCount ++;
							}
						}
					}
				}
			}
		}

		sasaData.maxPairsListItemsPerAtom = 0;
		sasaData.maxSASAPairsPerAtom = 0;
		for(i = 0; i < sasaData.threadsCount; i++){
			if(sasaData.maxPairsListItemsPerAtom < sasaData.h_pairsListCount[i]){
				sasaData.maxPairsListItemsPerAtom = sasaData.h_pairsListCount[i];
			}
			if(sasaData.maxSASAPairsPerAtom < sasaData.h_sasaListCount[i]){
				sasaData.maxSASAPairsPerAtom = sasaData.h_sasaListCount[i];
			}
		}
	    const int tw = 10; // Width of table field
	    LOG << "";
    	LOG.table(tw,3) << "Pairs" << "Total" << "Max/Atom";
    	LOG.table(tw,3) << "SASA" << totalSASAListCount << sasaData.maxSASAPairsPerAtom;
    	LOG.table(tw,3) << "Pairlist" << totalSASAPairlistCount << sasaData.maxPairsListItemsPerAtom;
        LOG << "";
	} else {
		LOG << "Using pre-defined list sizes of " << sasaData.maxSASAPairsPerAtom << " for SASA list and " << sasaData.maxPairsListItemsPerAtom << " for SASA pairlist";
		LOG << "Actual numbers of pairs in changeable lists were not yet computed and will be zeros in following table.";
	}

	DPRINTF("Pairs:        Total:      Max/Atom:\n");
	DPRINTF("SASA List %10d     %10d\n", totalSASAListCount, sasaData.maxSASAPairsPerAtom);
	DPRINTF("Pairlist  %10d     %10d\n", totalSASAPairlistCount, sasaData.maxPairsListItemsPerAtom);
	LOG << "Adding " << pairlistsExtension << " to the length of all changeable lists to avoid memory conflicts";
	sasaData.maxSASAPairsPerAtom += 4; //For covalent bonds
	sasaData.maxSASAPairsPerAtom += pairlistsExtension;
	sasaData.maxPairsListItemsPerAtom += pairlistsExtension;

	allocateCPU((void**)&sasaData.h_pairsList,
			sasaData.widthTot*sasaData.maxPairsListItemsPerAtom*sizeof(int));
	allocateGPU((void**)&sasaData.d_pairsList,
			sasaData.widthTot*sasaData.maxPairsListItemsPerAtom*sizeof(int));

	allocateCPU((void**)&sasaData.h_BGrad,
			sasaData.widthTot*sasaData.maxSASAPairsPerAtom*sizeof(float4));
	allocateGPU((void**)&sasaData.d_BGrad,
			sasaData.widthTot*sasaData.maxSASAPairsPerAtom*sizeof(float4));
	allocateCPU((void**)&sasaData.h_B,
			sasaData.widthTot*sasaData.maxSASAPairsPerAtom*sizeof(float));
	allocateGPU((void**)&sasaData.d_B,
			sasaData.widthTot*sasaData.maxSASAPairsPerAtom*sizeof(float));

	allocateCPU((void**)&sasaData.h_BGradT,
			sasaData.widthTot*sasaData.maxSASAPairsPerAtom*sizeof(float4));
	allocateGPU((void**)&sasaData.d_BGradT,
			sasaData.widthTot*sasaData.maxSASAPairsPerAtom*sizeof(float4));
	allocateCPU((void**)&sasaData.h_BT,
			sasaData.widthTot*sasaData.maxSASAPairsPerAtom*sizeof(float));
	allocateGPU((void**)&sasaData.d_BT,
			sasaData.widthTot*sasaData.maxSASAPairsPerAtom*sizeof(float));

	allocateCPU((void**)&sasaData.h_sasaList,
			sasaData.widthTot*sasaData.maxSASAPairsPerAtom*sizeof(int));
	allocateGPU((void**)&sasaData.d_sasaList,
			sasaData.widthTot*sasaData.maxSASAPairsPerAtom*sizeof(int));

	for(i = 0; i < sasaData.threadsCount; i++){
		sasaData.h_pairs12Counts[i] = 0;
	}



	int b;
	for(b = 0; b < topology.bondCount; b++){
		Bond bond = topology.bonds[b];
		i = sasaData.h_atomThread[bond.i];
		j = sasaData.h_atomThread[bond.j];
		if(topology.atoms[bond.i].type[0] != 'H' && topology.atoms[bond.j].type[0] != 'H'){
			sasaData.h_sasaList[sasaData.h_pairs12Counts[i]*sasaData.widthTot + i] = bond.j;
			sasaData.h_sasaList[sasaData.h_pairs12Counts[j]*sasaData.widthTot + j] = bond.i;
			sasaData.h_pairs12Counts[i] ++;
			sasaData.h_pairs12Counts[j] ++;
		}
	}

	/*for(i = 0; i < sasaData.threadsCount; i++){
		a1 = sasaData.h_threadAtom[i];
		sasaData.h_pairs12Counts[i] = pairsListsData.h_pairs12ListCount[a1];
		//printf("(%s)%d\n", topology.atoms[a1].name, sasaData.h_pairs12Counts[i]);
		for(j = 0; j < pairsListsData.h_pairs12ListCount[a1]; j++){
			sasaData.h_sasaList[j*sasaData.widthTot + i] =
					pairsListsData.h_pairs12List[j*gsystem.widthTot + a1];
		}
	}*/

	allocateCPU((void**)&sasaData.h_sasa, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&sasaData.d_sasa, gsystem.Ntot*sizeof(float));

	allocateCPU((void**)&sasaData.h_sasaEnergies, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&sasaData.d_sasaEnergies, gsystem.Ntot*sizeof(float));

	for(i = 0; i < gsystem.Ntot; i++){
		sasaData.h_sasaEnergies[i] = 0.0f;
	}


	/*printf("\n\nList of Covalent Bonds (in SASA):\n");
	for(i = 0; i < sasaData.threadsCount; i++){
		printf("%d (%d): ", i, sasaData.h_pairs12Counts[i]);
		for(j = 0; j < sasaData.h_pairs12Counts[i]; j++){
			printf("%d ",
					sasaData.h_sasaList[j*sasaData.widthTot + i]);
		}
		printf("\n");
	}*/

	int traj, itot;
	for(traj = 1; traj < parameters.Ntr; traj++){
		for(i = 0; i < sasaData.threadsCount; i++){
			itot = traj*sasaData.threadsCount + i;
			sasaData.h_threadAtom[itot] = sasaData.h_threadAtom[i] + traj*gsystem.N;
			sasaData.h_pairs12Counts[itot] = sasaData.h_pairs12Counts[i];
			for(j = 0; j < sasaData.h_pairs12Counts[i]; j++){
				sasaData.h_sasaList[j*sasaData.widthTot + itot] =
						sasaData.h_sasaList[j*sasaData.widthTot + i] + traj*gsystem.N;
			}
		}
	}


	DPRINTF("SASA Parameters:\n");
	DPRINTF("Atom type:\tR_ipr: \t\t(R_i+R_pr)^2: \t\tp_i/S_i: \t\tsigma*S_i\n");
	for(i = 0; i < atomTypesCount; i++){
		DPRINTF("%d.%s\t\t%6.4f\t\t%6.4f\t\t%6.4f\t\t%6.4f\n", i, atomTypes[i].name,
				sasaData.h_sasaRipr[i],
				sasaData.h_sasaRipr2[i],
				sasaData.h_sasaPiOverSi[i],
				sasaData.h_sasaSigmaSi[i]);
	}


	//cudaMemcpy(sasaData.d_sasaParameters, sasaData.h_sasaParameters,
	//				atomTypesCount*sizeof(GSASAParameters), cudaMemcpyHostToDevice);

	cudaMemcpy(sasaData.d_sasaRipr, sasaData.h_sasaRipr,
						atomTypesCount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(sasaData.d_sasaRipr2, sasaData.h_sasaRipr2,
						atomTypesCount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(sasaData.d_sasaPiOverSi, sasaData.h_sasaPiOverSi,
						atomTypesCount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(sasaData.d_sasaSigmaSi, sasaData.h_sasaSigmaSi,
						atomTypesCount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(sasaData.d_sasaSi, sasaData.h_sasaSi,
							atomTypesCount*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(sasaData.d_threadAtom, sasaData.h_threadAtom,
						sasaData.threadsCountTot*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(sasaData.d_pairs12Counts, sasaData.h_pairs12Counts,
						sasaData.threadsCountTot*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(sasaData.d_sasaList, sasaData.h_sasaList,
							sasaData.widthTot*sasaData.maxSASAPairsPerAtom*sizeof(int),
							cudaMemcpyHostToDevice);

	cudaMemcpy(sasaData.d_sasaEnergies, sasaData.h_sasaEnergies,
					gsystem.Ntot*sizeof(float), cudaMemcpyHostToDevice);

	/*cudaMemcpyToSymbol(c_sasaParameters, sasaData.h_sasaParameters,
			atomTypesCount*sizeof(GSASAParameters), 0);*/
	cudaMemcpyToSymbol(c_sasaData, &sasaData, sizeof(GSASAData), 0);

	cudaBindTexture(0, t_sasaRipr, sasaData.d_sasaRipr,
			atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_sasaRipr2, sasaData.d_sasaRipr2,
			atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_sasaPiOverSi, sasaData.d_sasaPiOverSi,
			atomTypesCount*sizeof(float));
	cudaBindTexture(0, t_sasaSigmaSi, sasaData.d_sasaSigmaSi,
			atomTypesCount*sizeof(float));

	/*cudaBindTexture(0, t_sasaBGrad, sasaData.d_BGrad,
			gsystem.widthTot*sasaData.maxSASAPairsPerAtom*sizeof(float4));
	cudaBindTexture(0, t_sasaB, sasaData.d_B,
			gsystem.widthTot*sasaData.maxSASAPairsPerAtom*sizeof(float));*/

	LOG << "Done initializing SASA implicit solvent potential.";
}

__global__ void generateSASAList_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_sasaData.threadsCountTot){
		int i = c_sasaData.d_threadAtom[d_i];
		float4 r1 = tex1Dfetch(t_coord, i);
		r1.w = tex1Dfetch(t_sasaRipr, (int)r1.w);//c_sasaParameters[(int)r1.w].Ripr;

		int sasaCount = c_sasaData.d_pairs12Counts[d_i];
		for(i = 0; i < c_sasaData.d_pairsListCount[d_i]; i++){
			float mult;
			float4 r2;

			int j = c_sasaData.d_pairsList[i*c_sasaData.widthTot + d_i];
			r2 = tex1Dfetch(t_coord, j);//abs(j));
			mult = tex1Dfetch(t_sasaRipr, (int)r2.w);//c_sasaParameters[(int)r2.w].Ripr;
			mult += r1.w;

			r2 -= r1;
			DO_PBC(r2);
			r2.w = sqrtf(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);

			if(r2.w < mult){
				c_sasaData.d_sasaList[sasaCount*c_sasaData.widthTot + d_i] = j;
				/*c_sasaData.d_sasaListPij[sasaCount*c_sasaData.widthTot + d_i] =
						c_sasaData.d_pairsListPij[i*c_sasaData.widthTot + d_i];*/
				//c_sasaData.d_sasadrs[sasaCount*c_sasaData.widthTot + d_i] = r2;
				sasaCount ++;
			}
		}
		c_sasaData.d_sasaListCount[d_i] = sasaCount;
	}
}

__global__ void computeSASAPotentialB_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_sasaData.threadsCountTot){
		int i = c_sasaData.d_threadAtom[d_i];
		int ati = c_gsystem.d_atomTypes[i];
		float4 B;
		float4 dr;
		s_Ri[threadIdx.x] = tex1Dfetch(t_coord, i);
		int atj;
		float pij;
		int covalentCount = c_sasaData.d_pairs12Counts[d_i];
		for(i = 0; i < c_sasaData.d_sasaListCount[d_i]; i++){

			atj = c_sasaData.d_sasaList[i*c_sasaData.widthTot + d_i];
			if(i < covalentCount){
				pij = c_sasaData.pij_cov;
			} else {
				pij = c_sasaData.pij_nb;
			}
			//atj = abs(atj);
			//float4 dr = c_sasaData.d_sasadrs[i*c_sasaData.widthTot + d_i];
			B = s_Ri[threadIdx.x];
			dr = tex1Dfetch(t_coord, atj);
			dr -= B;
			DO_PBC(dr);
			dr.w = sqrtf(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);


			atj = c_gsystem.d_atomTypes[atj];
			B.w = dr.w*dr.w;
			B.x = tex1Dfetch(t_sasaRipr2, ati);
			B.y = tex1Dfetch(t_sasaRipr2, atj);
			B.w = (B.y - B.x)/B.w;
			//B.w = (c_sasaParameters[atj].Ripr2 - c_sasaParameters[ati].Ripr2)/B.w;
			float mult = (1.0f + B.w)/dr.w;
			B.x = tex1Dfetch(t_sasaPiOverSi, ati);
			B.y = tex1Dfetch(t_sasaRipr, ati);
			mult *= B.x;//c_sasaParameters[ati].piOverSi;
			mult *= pij;
			mult *= M_PI;
			mult *= B.y;//c_sasaParameters[ati].Ripr;

			B.x = mult*dr.x;
			B.y = mult*dr.y;
			B.z = mult*dr.z;

			c_sasaData.d_BGrad[i*c_sasaData.widthTot + d_i] = B;

			B.x = tex1Dfetch(t_sasaPiOverSi, atj);
			B.y = tex1Dfetch(t_sasaRipr, atj);
			mult = (1.0f - B.w)/dr.w;
			mult *= B.x;//c_sasaParameters[atj].piOverSi;
			mult *= pij;
			mult *= M_PI;
			mult *= B.y;//c_sasaParameters[atj].Ripr;

			B.x = mult*dr.x;
			B.y = mult*dr.y;
			B.z = mult*dr.z;

			c_sasaData.d_BGradT[i*c_sasaData.widthTot + d_i] = B;

			dr.x = tex1Dfetch(t_sasaRipr, ati);
			dr.y = tex1Dfetch(t_sasaRipr, atj);
			dr.z = tex1Dfetch(t_sasaPiOverSi, ati);

			//B.x = c_sasaParameters[ati].Ripr;
			//B.x += c_sasaParameters[atj].Ripr;
			B.x = dr.x + dr.y;
			B.x -= dr.w;
			B.x *= M_PI;
			B.x *= pij;
			//B.y = c_sasaParameters[atj].Ripr;
			//B.y -= c_sasaParameters[ati].Ripr;
			B.y = dr.y - dr.x;
			B.y /= dr.w;

			dr.w = tex1Dfetch(t_sasaPiOverSi, atj);

			B.w = 1.0f + B.y;
			B.w *= B.x;
			B.w *= dr.x;//c_sasaParameters[ati].Ripr;
			B.w *= dr.z;//c_sasaParameters[ati].piOverSi;
			B.w = 1.0f - B.w;

			c_sasaData.d_BGrad[i*c_sasaData.widthTot + d_i].w = B.w;
			c_sasaData.d_B[i*c_sasaData.widthTot + d_i] = B.w;

			//c_sasaData.d_BGradT[i*c_sasaData.widthTot + d_i].w =
			//		1.0f - c_sasaParameters[atj].piOverSi*c_sasaParameters[atj].Ripr*B.x*(1.0f - B.y);
			c_sasaData.d_BGradT[i*c_sasaData.widthTot + d_i].w =
					1.0f - dr.w*dr.y*B.x*(1.0f - B.y);

			//c_sasaData.d_BT[sasaCount*c_gsystem.N + d_i] = B.w;
		}
	}
}

__global__ void computeSASAPotentialEnergy_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_sasaData.threadsCountTot){
		int a1 = c_sasaData.d_threadAtom[d_i];
		int j;
		int ati = c_gsystem.d_atomTypes[a1];
		int sasaCount = c_sasaData.d_sasaListCount[d_i];
		float pot = tex1Dfetch(t_sasaSigmaSi, ati);//c_sasaParameters[ati].sigmaSi;
		for(j = 0; j < sasaCount; j++){
			pot *= c_sasaData.d_B[j*c_sasaData.widthTot + d_i];//tex1Dfetch(t_sasaB, j*c_gsystem.N + d_i);
		}
		c_sasaData.d_sasaEnergies[a1] = pot;
	}
}

__global__ void computeSASA_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_sasaData.threadsCountTot){
		int a1 = c_sasaData.d_threadAtom[d_i];
		int j;
		int ati = c_gsystem.d_atomTypes[a1];
		int sasaCount = c_sasaData.d_sasaListCount[d_i];
		float pot = c_sasaData.d_sasaSi[ati];//c_sasaParameters[ati].sigmaSi;
		for(j = 0; j < sasaCount; j++){
			pot *= c_sasaData.d_B[j*c_sasaData.widthTot + d_i];//tex1Dfetch(t_sasaB, j*c_gsystem.N + d_i);
		}
		c_sasaData.d_sasa[a1] = pot;
	}
}

__global__ void computeSASAPotential_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_sasaData.threadsCountTot){
		int j, k;
		int a1 =  c_sasaData.d_threadAtom[d_i];
		float4 f = c_gsystem.d_forces[a1];
		float mult;
		float4 B;
		float Vi = c_sasaData.d_sasaEnergies[a1];
		for(k = 0; k < c_sasaData.d_sasaListCount[d_i]; k++){
			j = c_sasaData.d_sasaList[k*c_sasaData.widthTot + d_i];
			B = c_sasaData.d_BGrad[k*c_sasaData.widthTot + d_i];
			mult = Vi/B.w;
			f.x += mult*B.x;
			f.y += mult*B.y;
			f.z += mult*B.z;
			B = c_sasaData.d_BGradT[k*c_sasaData.widthTot + d_i];
			mult = c_sasaData.d_sasaEnergies[j];
			mult /= B.w;
			f.x += mult*B.x;
			f.y += mult*B.y;
			f.z += mult*B.z;
		}
		c_gsystem.d_forces[a1] = f;
	}
}

inline void compute(){
	generateSASAList_kernel<<<sasaBlockCount, sasaBlockSize>>>();
	computeSASAPotentialB_kernel<<<sasaBlockCount, sasaBlockSize>>>();
	/*cudaMemcpy(sasaData.h_sasaListCount, sasaData.d_sasaListCount,
			gsystem.Ntot*sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(sasaData.h_BGrad, sasaData.d_BGrad,
			gsystem.widthTot*sasaData.maxSASAPairsPerAtom*sizeof(float4), cudaMemcpyDeviceToHost);
	int i, j;
	printf("\n\nB's:\n");
		for(i = 0; i < gsystem.Ntot; i++){
			printf("%d (%d): ", i, sasaData.h_sasaListCount[i]);
			for(j = 0; j < sasaData.h_sasaListCount[i]; j++){
				printf("%5.2f ",
						sasaData.h_BGrad[j*gsystem.widthTot + i].w);
			}
			printf("\n");
		}*/

	computeSASAPotentialEnergy_kernel<<<sasaBlockCount, sasaBlockSize>>>();
	computeSASAPotential_kernel<<<sasaBlockCount, sasaBlockSize>>>();
	/*cudaMemcpy(gsystem.h_forces, gsystem.d_forces, gsystem.Ntot*sizeof(float4), cudaMemcpyDeviceToHost);
	//int i;
	float3 force = make_float3(0.0f, 0.0f, 0.0f);
	for(i = 0; i < gsystem.Ntot; i++){
		force.x += gsystem.h_forces[i].x;
		force.y += gsystem.h_forces[i].y;
		force.z += gsystem.h_forces[i].z;
		printf("%d: (%f, %f, %f) %f\n", i, gsystem.h_forces[i].x, gsystem.h_forces[i].y, gsystem.h_forces[i].z,
					sqrtf(gsystem.h_forces[i].x*gsystem.h_forces[i].x + gsystem.h_forces[i].y*gsystem.h_forces[i].y + gsystem.h_forces[i].z*gsystem.h_forces[i].z));
	}
	printf("Net force (sasa): (%f, %f, %f) %f\n", force.x, force.y, force.z,
			sqrtf(force.x*force.x + force.y*force.y + force.z*force.z));*/

	/*computeSASA_kernel<<<sasaBlockCount, sasaBlockSize>>>();
	cudaMemcpy(sasaData.h_sasa, sasaData.d_sasa, gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	float totalSASA = 0.0;
	FILE* file = fopen("sasa.dat", "w");
	int i;
	for(i = 0; i < gsystem.Ntot; i++){
		if(sasaData.h_atomThread[i] != -1){
			fprintf(file, "%d\t%f\n", i, sasaData.h_sasa[i]);
			totalSASA += sasaData.h_sasa[i];
		}
	}
	fclose(file);
	printf("Total SASA: %f\n", totalSASA);*/
	//exit(0);
}

inline void computeEnergy(){
	if(step == 0){
		updateSASAPairsList();
		generateSASAList_kernel<<<sasaBlockCount, sasaBlockSize>>>();
		computeSASAPotentialB_kernel<<<sasaBlockCount, sasaBlockSize>>>();
		computeSASAPotentialEnergy_kernel<<<sasaBlockCount, sasaBlockSize>>>();
	}
	cudaMemcpy(sasaData.h_sasaEnergies, sasaData.d_sasaEnergies,
				gsystem.Ntot*sizeof(float), cudaMemcpyDeviceToHost);
	checkCUDAError("sasa - Energy copy");
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		float pot = 0.0f;
		for(i = 0; i < gsystem.N; i++){
			pot += sasaData.h_sasaEnergies[i + gsystem.N*traj];
		}
		energyOutput.values[traj] = pot;
	}
	/*cudaMemcpy(sasaData.h_pairs12Counts, sasaData.d_pairs12Counts,
						sasaData.threadsCountTot*sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(sasaData.h_sasaListCount, sasaData.d_sasaListCount,
							sasaData.threadsCountTot*sizeof(int), cudaMemcpyDeviceToHost);

	cudaMemcpy(sasaData.h_sasaList, sasaData.d_sasaList,
							sasaData.widthTot*sasaData.maxSASAPairsPerAtom*sizeof(int),
							cudaMemcpyDeviceToHost);

	printf("\n");
	for(i = 0; i < sasaData.threadsCount; i++){
		int a1 = sasaData.h_threadAtom[i];
		int j;
		printf("%d (%s, %d in list, %d cov):  ", a1, topology.atoms[a1].name,
				sasaData.h_sasaListCount[i], sasaData.h_pairs12Counts[i]);
		for(j = 0; j < sasaData.h_sasaListCount[i]; j++){
			printf("%d  ", sasaData.h_sasaList[j*sasaData.widthTot + i]);
		}
		printf("\n");
	}*/

}

void destroy(){

}

__global__ void updateSASAPairsList_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_sasaData.threadsCountTot){
		int i, j, found;
		int a1, a2;
		a1 = c_sasaData.d_threadAtom[d_i];
		float4 r1 = tex1Dfetch(t_coord, a1);
		float4 r2;
		int sasaPairsCount = 0;
		int traj = d_i/c_sasaData.threadsCount;
		for(i = traj*c_sasaData.threadsCount; i < (traj + 1)*c_sasaData.threadsCount; i++){
			a2 = c_sasaData.d_threadAtom[i];
			if(a2 != a1){
				r2 = tex1Dfetch(t_coord, a2);
				r2 -= r1;
				DO_PBC(r2);
				r2.w = sqrtf(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);
				if(r2.w < c_sasaData.sasaPairsListCutoff){
					int covalentCount = c_sasaData.d_pairs12Counts[d_i];
					j = 0;
					found = 0;
					while(found == 0 && j < covalentCount){
						if(a2 == c_sasaData.d_sasaList[j*c_sasaData.widthTot + d_i]){
							found = 1;
						}
						j++;
					}


					//float pij;
					if(found == 0){
						//pij = c_sasaData.pij_nb;
					/*	j = -a2;
					} else {*/
					//	j = a2;
						//pij = c_sasaData.pij_cov;
						c_sasaData.d_pairsList[sasaPairsCount*c_sasaData.widthTot + d_i] = a2;
						sasaPairsCount ++;
					}
					//c_sasaData.d_pairsListPij[sasaPairsCount*c_sasaData.widthTot + d_i] = pij;



					/*if(found == 0){
						c_sasaData.d_pairsList[sasaPairsCount*c_sasaData.widthTot + d_i] = i;
						sasaPairsCount ++;
					}*/
				}
			}
		}
		c_sasaData.d_pairsListCount[d_i] = sasaPairsCount;
	}
}

inline void updateSASAPairsList(){
	updateSASAPairsList_kernel<<<sasaPairsListBlockCount, sasaPairsListBlockSize>>>();
	cudaMemcpy(sasaData.h_pairsListCount, sasaData.d_pairsListCount,
					sasaData.threadsCountTot*sizeof(int), cudaMemcpyDeviceToHost);
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		for(i = 0; i < sasaData.threadsCount; i++){
			if(sasaData.h_pairsListCount[i] > MAX_SASA_PAIRSLIST_ITEMS_PER_ATOM){
				DIE("Number of SASA pairs on atom %d (trajectory %d) is %d, which exceeds the limit of %d.",
						sasaData.h_threadAtom[i], traj+parameters.firstrun,
						sasaData.h_pairsListCount[i], MAX_SASA_PAIRSLIST_ITEMS_PER_ATOM);
			}
		}
	}
	cudaThreadSynchronize();
	checkCUDAError("update sasa pairs");
	cudaMemcpy(sasaData.h_pairsListCount, sasaData.d_pairsListCount,
			sasaData.threadsCountTot*sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(sasaData.h_pairsList, sasaData.d_pairsList,
			sasaData.widthTot*sasaData.maxPairsListItemsPerAtom*sizeof(int), cudaMemcpyDeviceToHost);
	/*int j;
	printf("\n\nList of SASA pairs:\n");
	for(i = 0; i < sasaData.threadsCountTot; i++){
		if(i % sasaData.threadsCount == 0){
		printf("%d. %d (%d): ", i, sasaData.h_threadAtom[i],
				sasaData.h_pairsListCount[i]);
		for(j = 0; j < sasaData.h_pairsListCount[i]; j++){
			printf("%d ",
					sasaData.h_pairsList[j*sasaData.widthTot + i]);
		}
		printf("\n");
		}
	}
	for(j = 0; j < sasaData.threadsCount; j++){
		if(sasaData.h_pairsListCount[j] != sasaData.h_pairsListCount[j+traj*sasaData.threadsCount]){
			DIE("Sukanah!");
		}
	}*/
}

void destroySASAPairsListUpdater(){

}

#undef LOG

} // namespace sasa_potential

