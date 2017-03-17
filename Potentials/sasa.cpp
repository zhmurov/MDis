/*
 * sasa.cpp
 *
 *  Created on: May 22, 2009
 *      Author: zhmurov
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sasa.h"
#include "../Util/mdunitsconverter.h"
#include "../Util/Log.h"
#include "../Util/wrapper.h"

namespace sasa_potential {

SASAParameters* sasaParameters;
int* sasaRef;
int sasaParamCount;


void convertSASAUnitsToMD();

SASAParameters getSASAParameters(char* atomType){
	return sasaParameters[getRef(atomType)];
}

SASAParameters getSASAParametersCHARMM(char* atomType){
	return sasaParameters[getRefCHARMM(atomType)];
}

void readSASAParameters(char* filename){
	sasaParamCount = 0;
	FILE* sasaFile = safe_fopen(filename, "r");
	char buffer[buf_size];
	while(fgets(buffer, buf_size, sasaFile) != NULL){
		if(buffer[0] != '#' && buffer[0] != ' ' && buffer[0] != '\n' && buffer[0] != '\t'){
			//printf("%s\n", buffer);
			sasaParamCount++;
		}
	}
	sasaParameters = (SASAParameters*)malloc(sasaParamCount*sizeof(SASAParameters));
	int i = 0;
	rewind(sasaFile);
	while(fgets(buffer, buf_size, sasaFile) != NULL){
		if(buffer[0] != '#' && buffer[0] != ' ' && buffer[0] != '\n' && buffer[0] != '\t'){
			char* pch = strtok(buffer, " \t");
			if (strlen(pch) > SASA_ATOMTYPE_SIZE - 1) // mind the \0
				DIE("Too long atom type name: %s", pch);
			strcpy(sasaParameters[i].type, pch);
			pch = strtok(NULL, " \t\n");
			sasaParameters[i].Rmin = atof(pch);
			pch = strtok(NULL, " \t\n");
			sasaParameters[i].R = atof(pch);
			pch = strtok(NULL, " \t\n");
			sasaParameters[i].p = atof(pch);
			pch = strtok(NULL, " \t\n");
			sasaParameters[i].sigma = atof(pch);
			DPRINTF("%s\t%f\t%f\t%f\t%f\n", sasaParameters[i].type, sasaParameters[i].Rmin, sasaParameters[i].R, sasaParameters[i].p, sasaParameters[i].sigma);
			i++;
		}
	}
	fclose(sasaFile);
	convertSASAUnitsToMD();
}

int getRef(char* atomType){
	FILE* refFile = safe_fopen("charmm-charmm.map", "r");
	char buffer[buf_size];
	int i = 0;
	char charmmType[SASA_ATOMTYPE_SIZE];
	while(fgets(buffer, buf_size, refFile) != NULL){
		if(buffer[0] != '#' && buffer[0] != ' ' && buffer[0] != '\n' && buffer[0] != '\t'){
			//printf("%s", buffer);
			char* pch = strtok(buffer, " \t");
			if(strcmp(atomType, pch) == 0){
				pch = strtok(NULL, " \t\n");
				strcpy(charmmType, pch);
				fclose(refFile);
				for(i = 0; i < sasaParamCount; i++){
					if(strcmp(charmmType, sasaParameters[i].type) == 0){
						return i;
					}
				}
			}
		}
	}
	fclose(refFile);
	DIE("Cannot find reference for type '%s'", atomType);
}

int getRefCHARMM(char* atomType){
	int i;
	for(i = 0; i < sasaParamCount; i++){
		if(strcmp(atomType, sasaParameters[i].type) == 0){
			return i;
		}
	}
	DIE("Cannot find reference for type '%s'", atomType);
}


void convertSASAUnitsToMD(){
	int i;
	for(i = 0; i < sasaParamCount; i++){
		sasaParameters[i].R = convertDistanceCHARMMtoMD(sasaParameters[i].R);
		sasaParameters[i].Rmin = convertDistanceCHARMMtoMD(sasaParameters[i].Rmin);
		sasaParameters[i].sigma = convertKbCHARMMtoMD(sasaParameters[i].sigma);
	}
}

} // namespace sasa_potential
