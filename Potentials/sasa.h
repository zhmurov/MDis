/*
 * sasa.h
 *
 *  Created on: May 22, 2009
 *      Author: zhmurov
 */

#pragma once

namespace sasa_potential {

#define SASA_ATOMTYPE_SIZE 8

#define buf_size 256

#define max_types 100

typedef struct {
	char type[SASA_ATOMTYPE_SIZE];
	float Rmin;
	float R;
	float p;
	float sigma;
} SASAParameters;

void readSASAParameters(char* filename);
SASAParameters getSASAParameters(char* atomType);
SASAParameters getSASAParametersCHARMM(char* atomType);
int getRef(char* atomType);
int getRefCHARMM(char* atomType);

} // namespace sasa_potential
