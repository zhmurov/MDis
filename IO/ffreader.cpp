/*
 * ffreader.c
 *
 *  Created on: Jul 4, 2009
 *      Author: zhmurov
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "../Util/wrapper.h"
#include <cmath>
#include "../Core/forcefield.h"
#include "ffreader.h"
#include "../Util/mdunitsconverter.h"
#include "../Core/topology.h"

#define BUF_SIZE 1024
//#define DEBUGFF

#define FF_SECTION_COUNT		6
#define FF_ALL_SECTION_COUNT	14

#define FF_TYPENAME_SIZE	5

#define FF_TYPEID_CHARMM22			0
#define FF_TYPEID_CHARMM19			1

#define FF_ATOM_NAMES_COUNT		7

const char ffAtomNames[FF_ATOM_NAMES_COUNT] = {'C', 'H', 'O', 'N', 'S', 'F', 'Z'};

const char ffSections[FF_SECTION_COUNT][2][16] = {
		{"BONDS", "BOND"},
		{"ANGLES", "THETAS"},
		{"DIHEDRALS", "PHI"},
		{"IMPROPER", "IMPHI"},
		{"CMAP", "CMAP"},
		{"NONBONDED", "NONBONDED"}
};

const char sectionsAll[FF_ALL_SECTION_COUNT][16] = {
		"BONDS",
		"ANGLES",
		"DIHEDRALS",
		"IMPROPER",
		"NONBONDED"
		"BOND",
		"THETAS",
		"PHI",
		"IMPHI",
		"CMAP",
		"HBOND",
		"END",
		"NBFIX",
		"cutnb",

};

int ffTypeID;

/*
 * Private methods
 */
void countFFEntries(FILE* file, ForceField* ff);
void readBonds(FILE* file, ForceField* ff);
void readAngles(FILE* file, ForceField* ff);
void readDihedrals(FILE* file, ForceField* ff);
void readImpropers(FILE* file, ForceField* ff);
void readCMAPs(FILE* file, ForceField* ff);
void readNonBonded(FILE* file, ForceField* ff);
void convertUnits(ForceField* ff);

void findFFSection(int sectionID, FILE* file);
int checkIfFFSectionEnded(char* buffer);

char buffer[BUF_SIZE];

/*
 * Parses NAMD force-field file (usually called par_all*.inp).
 * Parameters:
 * 		filename: name of the file to parse
 * 		ff: pointer to the object to save data into
 */
void readCHARMMFF(const char* filename, ForceField* ff, char* fftype){
	printf("Reading force-field data from '%s'.\n", filename);
	if(strcmp(fftype, FF_TYPE_CHARMM19) == 0){
		ffTypeID = FF_TYPEID_CHARMM19;
	} else {
		ffTypeID = FF_TYPEID_CHARMM22;
	}
	FILE* file = safe_fopen(filename, "r");
	if(file != NULL){

		//Get the numbers of entities in the file
		countFFEntries(file, ff);

		//Output the numbers of entities
		printf("Found:\n");
		printf("%d bond types.\n", ff->bondTypesCount);
		printf("%d angle types.\n", ff->angleTypesCount);
		printf("%d dihedral types.\n", ff->dihedralTypesCount);
		printf("%d improper types.\n", ff->improperTypesCount);
		printf("%d CMAP corrections types.\n", ff->cmapTypesCount);
		printf("%d non-bonded types.\n", ff->nonbondedTypesCount);

		//Allocate memory for entities
		ff->bondTypes = (FFBondType*)calloc(ff->bondTypesCount, sizeof(FFBondType));
		ff->angleTypes = (FFAngleType*)calloc(ff->angleTypesCount, sizeof(FFAngleType));
		ff->dihedralTypes = (FFDihedralType*)calloc(ff->dihedralTypesCount, sizeof(FFDihedralType));
		ff->improperTypes = (FFImproperType*)calloc(ff->improperTypesCount, sizeof(FFImproperType));
		ff->cmapTypes = (FFCMAPType*)calloc(ff->cmapTypesCount, sizeof(FFCMAPType));
		ff->nonbondedTypes = (FFNonbondedType*)calloc(ff->nonbondedTypesCount, sizeof(FFNonbondedType));

		//Start from the beginning of the file
		rewind(file);

		//Read data
		readBonds(file, ff);
		readAngles(file, ff);
		readDihedrals(file, ff);
		readImpropers(file, ff);
		readCMAPs(file, ff);
		readNonBonded(file, ff);

		//Convert units
		convertUnits(ff);

		printf("Done reading force-field data.\n");

	} else {
		DIE("ERROR: Force-field file '%s' can not be found.", filename);
	}
	fclose(file);
}

/*
 * Counts the number of bonds, angles, dihedral, impropers and non-bonded types in
 * a force-field file
 */
void countFFEntries(FILE* file, ForceField* ff){

	ff->bondTypesCount = 0;
	ff->angleTypesCount = 0;
	ff->dihedralTypesCount = 0;
	ff->improperTypesCount = 0;
	ff->cmapTypesCount = 0;

	char * pch;
	int i;
	char* eof;
	findFFSection(FF_SECTION_BONDS, file);
	if(feof(file) == 0){
		eof = safe_fgets(buffer, BUF_SIZE, file);

		while (!checkIfFFSectionEnded(buffer) && eof != NULL) {
			pch = strtok(buffer, " \t");
			if(pch[0] != '!' && pch[0] != ' ' && pch[0] != '\t' && pch[0] != '\n'){
				ff->bondTypesCount++;
			}
			eof = safe_fgets(buffer, BUF_SIZE, file);
		}
	}

	rewind(file);
	findFFSection(FF_SECTION_ANGLES, file);
	if(feof(file) == 0){
		eof = safe_fgets(buffer, BUF_SIZE, file);
		while (!checkIfFFSectionEnded(buffer) && eof != NULL){
			pch = strtok(buffer, " \t");
			if(pch[0] != '!' && pch[0] != ' ' && pch[0] != '\t' && pch[0] != '\n'){
				ff->angleTypesCount++;
			}
			eof = safe_fgets(buffer, BUF_SIZE, file);
		}
	}

	rewind(file);
	findFFSection(FF_SECTION_DIHEDRALS, file);
	if(feof(file) == 0){
		eof = safe_fgets(buffer, BUF_SIZE, file);
		while (!checkIfFFSectionEnded(buffer) && eof != NULL){
			pch = strtok(buffer, " \t");
			if(pch[0] != '!' && pch[0] != ' ' && pch[0] != '\t' && pch[0] != '\n'){
				ff->dihedralTypesCount++;
			}
			eof = safe_fgets(buffer, BUF_SIZE, file);
		}
	}

	rewind(file);
	findFFSection(FF_SECTION_IMPROPERS, file);
	if(feof(file) == 0){
		eof = safe_fgets(buffer, BUF_SIZE, file);
		while (!checkIfFFSectionEnded(buffer) && eof != NULL){
			pch = strtok(buffer, " \t");
			if(pch[0] != '!' && pch[0] != ' ' && pch[0] != '\t' && pch[0] != '\n'){
				ff->improperTypesCount++;
			}
			eof = safe_fgets(buffer, BUF_SIZE, file);
		};
	}

	rewind(file);
	findFFSection(FF_SECTION_CMAP, file);
	if(feof(file) == 0){
		eof = safe_fgets(buffer, BUF_SIZE, file);
		while (!checkIfFFSectionEnded(buffer) && eof != NULL){
			pch = strtok(buffer, " \t");
			for(i = 0; i < FF_ATOM_NAMES_COUNT; i++){
				if(pch[0] == ffAtomNames[i]){
					ff->cmapTypesCount++;
				}
			}
			eof = safe_fgets(buffer, BUF_SIZE, file);
		};
	}

	rewind(file);
	findFFSection(FF_SECTION_NONBONDED, file);
	if(feof(file) == 0){
		eof = safe_fgets(buffer, BUF_SIZE, file);
		while (!checkIfFFSectionEnded(buffer) && eof != NULL){
			eof = safe_fgets(buffer, BUF_SIZE, file);
			pch = strtok(buffer, " \t");
			if(pch[0] != '!' && pch[0] != ' ' && pch[0] != '\t' && pch[0] != '\n'){
				ff->nonbondedTypesCount++;
			}
		}
	}
}

/*
 * Reads the bond parameters (kb and b0) from all entries in 'bonds' section
 *
 * V(bond) = kb*(b - b0)^2
 * kb: kcal/mole/A^2
 * b0: A
 *
 * Parameters:
 * 		file: pointer to the FILE object to read from
 * 		ff: ForceField object to save data into
 */
void readBonds(FILE* file, ForceField* ff){
	int index = 0;
	char * pch;
	findFFSection(FF_SECTION_BONDS, file);
	if(feof(file) == 0){
		char* eof = safe_fgets(buffer, BUF_SIZE, file);
		if(eof == NULL){
			return;
		}
		do{
			pch = strtok(buffer, " \t");
			if(pch[0] != '!' && pch[0] != ' ' && pch[0] != '\t' && pch[0] != '\n'){
				strcpy(ff->bondTypes[index].atomType1.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->bondTypes[index].atomType2.name, pch);
				pch = strtok(NULL, " \t");
				ff->bondTypes[index].kb = atof(pch);
				pch = strtok(NULL, " \t");
				ff->bondTypes[index].b0 = atof(pch);
	#ifdef DEBUGFF
				printf("Bond: %s %s %f %f\n",
						ff->bondTypes[index].atomType1.name,
						ff->bondTypes[index].atomType1.name,
						ff->bondTypes[index].kb,
						ff->bondTypes[index].b0);
	#endif
				index ++;
			}
			safe_fgets(buffer, BUF_SIZE, file);
		} while (!checkIfFFSectionEnded(buffer));
	}
	rewind(file);
}

/*
 * Reads the angles parameters (ktheta, theta0, kub and s0) from all entries in 'angles' section
 *
 * V(angle) = ktheta*(theta - theta0)^2
 * V(Urey-Bradley) = kub*(s - s0)^2
 * ktheta: kcal/mole/rad^2
 * theta0: degrees
 * kub: kcal/mole/A^2 (Urey-Bradley)
 * s0: A
 *
 * Parameters:
 * 		file: pointer to the FILE object to read from
 * 		ff: ForceField object to save data into
 */
void readAngles(FILE* file, ForceField* ff){
	int index = 0;
	char * pch;
	findFFSection(FF_SECTION_ANGLES, file);
	if(feof(file) == 0){
		char* eof = safe_fgets(buffer, BUF_SIZE, file);
		if(eof == NULL){
			return;
		}
		do{
			pch = strtok(buffer, " \t");
			if(pch[0] != '!' && pch[0] != ' ' && pch[0] != '\t' && pch[0] != '\n'){
				strcpy(ff->angleTypes[index].atomType1.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->angleTypes[index].atomType2.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->angleTypes[index].atomType3.name, pch);
				pch = strtok(NULL, " \t");
				ff->angleTypes[index].ktheta = atof(pch);
				pch = strtok(NULL, " \t");
				ff->angleTypes[index].theta0 = atof(pch);
				pch = strtok(NULL, " \t");
				if(pch != NULL){
					ff->angleTypes[index].kub = atof(pch);
					pch = strtok(NULL, " \t");
					if(pch != NULL){
						ff->angleTypes[index].s0 = atof(pch);
					} else {
						ff->angleTypes[index].s0 = 0.0f;
					}
				} else {
					ff->angleTypes[index].kub = 0.0f;
					ff->angleTypes[index].s0 = 0.0f;
				}
	#ifdef DEBUGFF
				printf("Angle: %s %s %s %f %f %f %f\n",
						ff->angleTypes[index].atomType1.name,
						ff->angleTypes[index].atomType2.name,
						ff->angleTypes[index].atomType3.name,
						ff->angleTypes[index].theta0,
						ff->angleTypes[index].ktheta,
						ff->angleTypes[index].kub,
						ff->angleTypes[index].s0);
	#endif

				index++;
			}
			safe_fgets(buffer, BUF_SIZE, file);
		} while (!checkIfFFSectionEnded(buffer));
	}
	rewind(file);
}

/*
 * Reads the dihedral parameters (kchi, n and delta) from all entries in 'dihedrals' section
 *
 * V(dihedral) = kchi*(1 + cos(n(chi) - delta))
 * Kchi: kcal/mole
 * n: multiplicity
 * delta: degrees
 *
 * Parameters:
 * 		file: pointer to the FILE object to read from
 * 		ff: ForceField object to save data into
 */
void readDihedrals(FILE* file, ForceField* ff){
	int index = 0;
	char * pch;
	findFFSection(FF_SECTION_DIHEDRALS, file);
	if(feof(file) == 0){
		char* eof = safe_fgets(buffer, BUF_SIZE, file);
		if(eof == NULL){
			return;
		}
		do{
			pch = strtok(buffer, " \t");
			if(pch[0] != '!' && pch[0] != ' ' && pch[0] != '\t' && pch[0] != '\n'){
				strcpy(ff->dihedralTypes[index].atomType1.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->dihedralTypes[index].atomType2.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->dihedralTypes[index].atomType3.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->dihedralTypes[index].atomType4.name, pch);
				pch = strtok(NULL, " \t");
				ff->dihedralTypes[index].kchi = atof(pch);
				pch = strtok(NULL, " \t");
				ff->dihedralTypes[index].n = atoi(pch);
				pch = strtok(NULL, " \t");
				ff->dihedralTypes[index].delta = atof(pch);
	#ifdef DEBUGFF
				printf("Dihedral: %s %s %s %s %f %d %f\n",
						ff->dihedralTypes[index].atomType1.name,
						ff->dihedralTypes[index].atomType2.name,
						ff->dihedralTypes[index].atomType3.name,
						ff->dihedralTypes[index].atomType4.name,
						ff->dihedralTypes[index].kchi,
						ff->dihedralTypes[index].n,
						ff->dihedralTypes[index].delta);
	#endif
				index++;
			}
			safe_fgets(buffer, BUF_SIZE, file);
		} while (!checkIfFFSectionEnded(buffer));
	}
	rewind(file);
	/*int i;
	for(i = 0; i < ff->dihedralTypesCount-1; i++){
		if(strncmp(ff->dihedralTypes[i].atomType1.name, ff->dihedralTypes[i+1].atomType1.name, 5) == 0 &&
				strncmp(ff->dihedralTypes[i].atomType2.name, ff->dihedralTypes[i+1].atomType2.name, 5) == 0 &&
				strncmp(ff->dihedralTypes[i].atomType3.name, ff->dihedralTypes[i+1].atomType3.name, 5) == 0 &&
				strncmp(ff->dihedralTypes[i].atomType4.name, ff->dihedralTypes[i+1].atomType4.name, 5) == 0){
			printf("Dihedral (%s - %s - %s - %s) has multiple parameters.\n",
					ff->dihedralTypes[i].atomType1.name,
					ff->dihedralTypes[i].atomType2.name,
					ff->dihedralTypes[i].atomType3.name,
					ff->dihedralTypes[i].atomType4.name);
		}
	}*/

}

/*
 * Reads the impropers parameters (kpsi and psi0) from all entries in 'impropers' section
 *
 * V(improper) = kpsi*(psi - psi0)^2
 * kpsi: kcal/mole/rad^2
 * psi0: degrees
 *
 * Parameters:
 * 		file: pointer to the FILE object to read from
 * 		ff: ForceField object to save data into
 */
void readImpropers(FILE* file, ForceField* ff){
	int index = 0;
	char * pch;
	findFFSection(FF_SECTION_IMPROPERS, file);
	if(feof(file) == 0){
		char* eof = safe_fgets(buffer, BUF_SIZE, file);
		if(eof == NULL){
			return;
		}
		do{
			pch = strtok(buffer, " \t");
			if(pch[0] != '!' && pch[0] != ' ' && pch[0] != '\t' && pch[0] != '\n'){
				strcpy(ff->improperTypes[index].atomType1.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->improperTypes[index].atomType2.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->improperTypes[index].atomType3.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->improperTypes[index].atomType4.name, pch);
				pch = strtok(NULL, " \t");
				ff->improperTypes[index].kpsi = atof(pch);
				pch = strtok(NULL, " \t");
				ff->improperTypes[index].n = atoi(pch);
				pch = strtok(NULL, " \t");
				ff->improperTypes[index].psi0 = atof(pch);
	#ifdef DEBUGFF
				printf("Improper: %s %s %s %s %f %d %f\n",
						ff->improperTypes[index].atomType1.name,
						ff->improperTypes[index].atomType2.name,
						ff->improperTypes[index].atomType3.name,
						ff->improperTypes[index].atomType4.name,
						ff->improperTypes[index].kpsi,
						ff->improperTypes[index].n,
						ff->improperTypes[index].psi0);
	#endif
				index++;
			}
			safe_fgets(buffer, BUF_SIZE, file);
		} while (!checkIfFFSectionEnded(buffer));
	}
	rewind(file);
}


/*
 * Reads the CMAP tables
 *
 * Parameters:
 * 		file: pointer to the FILE object to read from
 * 		ff: ForceField object to save data into
 */
void readCMAPs(FILE* file, ForceField* ff){
	int index = 0;
	char * pch;
	int i, j;
	findFFSection(FF_SECTION_CMAP, file);
	if(feof(file) == 0){
		char* eof = safe_fgets(buffer, BUF_SIZE, file);
		if(eof == NULL){
			return;
		}
		do{
			pch = strtok(buffer, " \t");
			if(pch[0] != '!' && pch[0] != ' ' && pch[0] != '\t' && pch[0] != '\n'){
				strcpy(ff->cmapTypes[index].atomType1.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->cmapTypes[index].atomType2.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->cmapTypes[index].atomType3.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->cmapTypes[index].atomType4.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->cmapTypes[index].atomType5.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->cmapTypes[index].atomType6.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->cmapTypes[index].atomType7.name, pch);
				pch = strtok(NULL, " \t");
				strcpy(ff->cmapTypes[index].atomType8.name, pch);
				pch = strtok(NULL, " \t");
				int gridSize = atoi(pch);
				ff->cmapTypes[index].gridSize = gridSize;
				ff->cmapTypes[index].data = (float*)calloc(gridSize*gridSize, sizeof(float));
				for(i = 0; i < gridSize; i++){
					do{
						safe_fgets(buffer, BUF_SIZE, file);
						//printf("%s", buffer);
						pch = strtok(buffer, " \t");
					} while(pch[0] == '!' || pch[0] == ' ' || pch[0] == '\t' || pch[0] == '\n');
					for(j = 0; j < gridSize; j++){
						int ij = i*gridSize + j;
						ff->cmapTypes[index].data[ij] = atof(pch);
						if((j+1) % 5 == 0){
							safe_fgets(buffer, BUF_SIZE, file);
							pch = strtok(buffer, " \t");
						} else {
							pch = strtok(NULL, " \t");
						}
					}
				}
	#ifdef DEBUGFF
				printf("CMAP: %s %s %s %s %s %s %s %s %d\n",
						ff->cmapTypes[index].atomType1.name,
						ff->cmapTypes[index].atomType2.name,
						ff->cmapTypes[index].atomType3.name,
						ff->cmapTypes[index].atomType4.name,
						ff->cmapTypes[index].atomType5.name,
						ff->cmapTypes[index].atomType6.name,
						ff->cmapTypes[index].atomType7.name,
						ff->cmapTypes[index].atomType8.name,
						ff->cmapTypes[index].gridSize);
				for(i = 0; i < gridSize; i++){
					for(j = 0; j < gridSize; j++){
						int ij = i*gridSize + j;
						printf("%5.3f  ", ff->cmapTypes[index].data[ij]);
					}
					printf("\n");
				}
	#endif
				index++;
			}
			safe_fgets(buffer, BUF_SIZE, file);
		} while (!checkIfFFSectionEnded(buffer));
	}
	rewind(file);
}

/*
 * Reads the non-bonded parameters (epsilon (1-4) and Rmin/2 (1-4)) from all entries in 'nonbonded' section
 *
 * V(Lennard-Jones) = Eps,i,j*[(Rmin,i,j/ri,j)^12 - 2*(Rmin,i,j/ri,j)^6]
 * epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
 * Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
 *
 * Parameters:
 * 		file: pointer to the FILE object to read from
 * 		ff: ForceField object to save data into
 */
void readNonBonded(FILE* file, ForceField* ff){
	int index = 0;
	char * pch;
	findFFSection(FF_SECTION_NONBONDED, file);
	if(feof(file) == 0){
		safe_fgets(buffer, BUF_SIZE, file);
		do{
			pch = strtok(buffer, " \t");
			if(pch[0] != '!' && pch[0] != ' ' && pch[0] != '\t' && pch[0] != '\n'){
				strcpy(ff->nonbondedTypes[index].atomType.name, pch);
				pch = strtok(NULL, " \t");
				ff->nonbondedTypes[index].ignored = atof(pch);
				pch = strtok(NULL, " \t");
				ff->nonbondedTypes[index].epsilon = atof(pch);
				pch = strtok(NULL, " \t");
				ff->nonbondedTypes[index].RminOver2 = atof(pch);
				pch = strtok(NULL, " \t");
				if(pch != NULL){
					ff->nonbondedTypes[index].ignored_14 = atof(pch);
					pch = strtok(NULL, " \t");
					if(pch != NULL){
						ff->nonbondedTypes[index].epsilon_14 = atof(pch);
						pch = strtok(NULL, " \t");
						if(pch != NULL){
							ff->nonbondedTypes[index].RminOver2_14 = atof(pch);
						} else {
							ff->nonbondedTypes[index].RminOver2_14 = ff->nonbondedTypes[index].RminOver2;
						}
					} else {
						ff->nonbondedTypes[index].epsilon_14 = ff->nonbondedTypes[index].epsilon;
						ff->nonbondedTypes[index].RminOver2_14 = ff->nonbondedTypes[index].RminOver2;
					}
				} else {
					ff->nonbondedTypes[index].ignored_14 = ff->nonbondedTypes[index].ignored;
					ff->nonbondedTypes[index].epsilon_14 = ff->nonbondedTypes[index].epsilon;
					ff->nonbondedTypes[index].RminOver2_14 = ff->nonbondedTypes[index].RminOver2;
				}
				if(ff->nonbondedTypes[index].epsilon_14 == 0.0f){
					ff->nonbondedTypes[index].ignored_14 = ff->nonbondedTypes[index].ignored;
					ff->nonbondedTypes[index].epsilon_14 = ff->nonbondedTypes[index].epsilon;
					ff->nonbondedTypes[index].RminOver2_14 = ff->nonbondedTypes[index].RminOver2;
				}
				/*if(ffTypeID == FF_TYPEID_CHARMM19){
					ff->nonbondedTypes[index].RminOver2 /= 2.0f;
					ff->nonbondedTypes[index].RminOver2_14 /= 2.0f;
				}*/
	#ifdef DEBUGFF
				printf("Non-bonded: %s: %f %f %f - %f %f %f\n",
						ff->nonbondedTypes[index].atomType.name,
						ff->nonbondedTypes[index].ignored,
						ff->nonbondedTypes[index].epsilon,
						ff->nonbondedTypes[index].RminOver2,
						ff->nonbondedTypes[index].ignored_14,
						ff->nonbondedTypes[index].epsilon_14,
						ff->nonbondedTypes[index].RminOver2_14);
	#endif
				index++;
			}
			safe_fgets(buffer, BUF_SIZE, file);
		} while (!checkIfFFSectionEnded(buffer));
	}
	rewind(file);
}

/*
 * Converts all the parameters readed into MD units. See Gromacs Users Manual, Chapter I for details.
 * Parameters:
 * 		ff: ForceField object to convert
 */
void convertUnits(ForceField* ff){
	int i;
	for(i = 0; i < ff->bondTypesCount; i++){
		ff->bondTypes[i].kb = convertKbCHARMMtoMD(ff->bondTypes[i].kb);
		ff->bondTypes[i].b0 = convertDistanceCHARMMtoMD(ff->bondTypes[i].b0);
	}
	for(i = 0; i < ff->angleTypesCount; i++){
		ff->angleTypes[i].ktheta = convertKthetaCHARMMtoMD(ff->angleTypes[i].ktheta);
		ff->angleTypes[i].theta0 = convertAngleCHARMMtoMD(ff->angleTypes[i].theta0);
		ff->angleTypes[i].kub = convertKubCHARMMtoMD(ff->angleTypes[i].kub);
		ff->angleTypes[i].s0 = convertDistanceCHARMMtoMD(ff->angleTypes[i].s0);
	}
	for(i = 0; i < ff->dihedralTypesCount; i++){
		ff->dihedralTypes[i].kchi = convertKchiCHARMMtoMD(ff->dihedralTypes[i].kchi);
		ff->dihedralTypes[i].delta = convertAngleCHARMMtoMD(ff->dihedralTypes[i].delta);
	}
	for(i = 0; i < ff->improperTypesCount; i++){
		ff->improperTypes[i].kpsi = convertKpsiCHARMMtoMD(ff->improperTypes[i].kpsi);
		ff->improperTypes[i].psi0 = convertAngleCHARMMtoMD(ff->improperTypes[i].psi0);
	}
	for(i = 0; i < ff->nonbondedTypesCount; i++){
		ff->nonbondedTypes[i].epsilon = convertEpsilonCHARMMtoMD(ff->nonbondedTypes[i].epsilon);
		ff->nonbondedTypes[i].epsilon_14 = convertEpsilonCHARMMtoMD(ff->nonbondedTypes[i].epsilon_14);
		ff->nonbondedTypes[i].RminOver2 = convertDistanceCHARMMtoMD(ff->nonbondedTypes[i].RminOver2);
		ff->nonbondedTypes[i].RminOver2_14 = convertDistanceCHARMMtoMD(ff->nonbondedTypes[i].RminOver2_14);
	}

}

void findFFSection(int sectionID, FILE* file){

	char* eof = safe_fgets(buffer, BUF_SIZE, file);
	while(strncmp(buffer, ffSections[sectionID][ffTypeID], 3) != 0 && eof != NULL){
		//printf("Searching for %s or %s in: %s\n", sectionsNAMD[sectionID], sectionsCHARMM[sectionID], buffer);
		eof = safe_fgets(buffer, BUF_SIZE, file);
	}
	//printf("Found %s or %s in: %s\n", sectionsNAMD[sectionID], sectionsCHARMM[sectionID], buffer);
}

int checkIfFFSectionEnded(char* buffer){
	//printf("Checking: %s\n", buffer);
	int ended = 0;
	int i;
	for(i = 0; i < FF_ALL_SECTION_COUNT; i++){
		if(strncmp(buffer, sectionsAll[i], 3) == 0){
			ended = 1;
		}
	}
	return ended;
}

void convertAtomTypesCHARMM19(const char* topFilename, Atom* atoms, int atomCount){
	printf("Searching CHARMM atom types in '%s'...\n", topFilename);
	FILE* file = safe_fopen(topFilename, "r");
	int typesCount = 0;
	char* pch;
	char* result;
	if(file != NULL){
		result = safe_fgets(buffer, BUF_SIZE, file);
		while(result != NULL){
			pch = strtok(buffer, " \t");
			if(strncmp(pch, "MASS", 4) == 0){
				typesCount ++;
			}
			result = safe_fgets(buffer, BUF_SIZE, file);
		}
		printf("Found %d CHARMM atom types:\n", typesCount);
		int* typesIDs = (int*)calloc(typesCount, sizeof(int));
		char** typeNames = (char**)calloc(typesCount, sizeof(char*));
		int i, j;
		for(i = 0; i < typesCount; i++){
			typeNames[i] = (char*)calloc(FF_TYPENAME_SIZE, sizeof(char*));
		}
		rewind(file);
		typesCount = 0;
		result = safe_fgets(buffer, BUF_SIZE, file);
		while(result != NULL){
			pch = strtok(buffer, " \t");
			if(strncmp(pch, "MASS", 4) == 0){
				pch = strtok(NULL, " \t");
				typesIDs[typesCount] = atoi(pch);
				pch = strtok(NULL, " \t");
				strcpy(typeNames[typesCount], pch);
				printf("%d - %s\n", typesIDs[typesCount], typeNames[typesCount]);
				typesCount ++;
			}
			result = safe_fgets(buffer, BUF_SIZE, file);
		}
		printf("Converting type IDs to type names for all atoms.\n");
		for(i = 0; i < atomCount; i++){
			for(j = 0; j < typesCount; j++){
				int typeID = atoi(atoms[i].type);
				if(typeID == typesIDs[j]){
					strcpy(atoms[i].type, typeNames[j]);
				}
			}
		}
	} else {
		DIE("ERROR: File '%s' not found.", topFilename);
	}
	fclose(file);
	printf("Done converting force-field atom types.\n");
}
