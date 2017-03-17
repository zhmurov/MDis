/*
 * topology.cpp
 *
 *  Created on: Jul 29, 2010
 *      Author: zhmurov
 */
#include "topology.h"
#include "../Core/global.h"
#include "../Util/ran2.h"
#include "../Util/wrapper.h"

int* atomIDs;

void initTopology(){
	PSF psf;
	readPSF(parameters.topologyFilename, &psf);
	printf("Initializing topology...\n");

	topology.atomCount = psf.natom;
	topology.bondCount = psf.nbond;
	topology.angleCount = psf.ntheta;
	topology.ureyBradleyCount = 0;
	topology.dihedralCount = psf.nphi;
	topology.improperCount = psf.nimphi;
	topology.cmapCount = psf.ncmap;
	topology.exclusionsCount = psf.nnb;

	topology.atoms = (Atom*)calloc(topology.atomCount, sizeof(Atom));
	topology.bonds = (Bond*)calloc(topology.bondCount, sizeof(Bond));
	topology.angles = (Angle*)calloc(topology.angleCount, sizeof(Angle));
	topology.dihedrals = (Dihedral*)calloc(topology.dihedralCount, sizeof(Dihedral));
	topology.impropers = (Improper*)calloc(topology.improperCount, sizeof(Improper));
	if(topology.cmapCount > 0){
		topology.cmaps = (CMAPCorrection*)calloc(topology.cmapCount, sizeof(CMAPCorrection));
	}
	if(topology.exclusionsCount > 0){
		topology.exclusions = (ExplicitExclusion*)calloc(topology.exclusionsCount, sizeof(ExplicitExclusion));
	}

	int i, j;
	int largestID = 0;
	for(i = 0; i < topology.atomCount; i++){
		topology.atoms[i].id = psf.atoms[i].id;
		strcpy(topology.atoms[i].name, psf.atoms[i].name);
		strcpy(topology.atoms[i].type, psf.atoms[i].type);
		topology.atoms[i].resid = psf.atoms[i].resid;
		strcpy(topology.atoms[i].resName, psf.atoms[i].resName);
		strcpy(topology.atoms[i].segment, psf.atoms[i].segment);
		topology.atoms[i].charge = psf.atoms[i].q;
		topology.atoms[i].mass = psf.atoms[i].m;
		if(largestID < psf.atoms[i].id){
			largestID = psf.atoms[i].id;
		}
	}
	atomIDs = (int*)calloc(largestID+1, sizeof(int));
	for(i = 0; i < largestID+1; i++){
		atomIDs[i] = -1;
	}
	for(i = 0; i < topology.atomCount; i++){
		atomIDs[topology.atoms[i].id] = i;
	}

	if(strncmp("CHARMM", parameters.parametersType, 6) == 0){
		convertAtomTypesCHARMM19(parameters.topologiesFilename, topology.atoms, topology.atomCount);
	}

	for(i = 0; i < topology.atomCount; i++){
		FFNonbondedType ffpar = findNonbondedType(topology.atoms[i].type, &ff);
		topology.atoms[i].ignored = ffpar.ignored;
		topology.atoms[i].epsilon = ffpar.epsilon;
		topology.atoms[i].RminOver2 = ffpar.RminOver2;
		topology.atoms[i].ignored_14 = ffpar.ignored_14;
		topology.atoms[i].epsilon_14 = ffpar.epsilon_14;
		topology.atoms[i].RminOver2_14 = ffpar.RminOver2_14;
	}

	for(i = 0; i < topology.bondCount; i++){
		topology.bonds[i].i = findAtom(psf.bonds[i].i);
		topology.bonds[i].j = findAtom(psf.bonds[i].j);
		FFBondType bt = findBondType(
				topology.atoms[topology.bonds[i].i].type,
				topology.atoms[topology.bonds[i].j].type,
				&ff);
		topology.bonds[i].kb = bt.kb;
		topology.bonds[i].b0 = bt.b0;
	}

	//printf("Angles in topology:\n");
	for(i = 0; i < topology.angleCount; i++){
		topology.angles[i].i = findAtom(psf.angles[i].i);
		topology.angles[i].j = findAtom(psf.angles[i].j);
		topology.angles[i].k = findAtom(psf.angles[i].k);
		//printf("top: %d - %d - %d\n", topology.angles[i].i, topology.angles[i].j, topology.angles[i].k);
		FFAngleType at = findAngleType(
				topology.atoms[topology.angles[i].i].type,
				topology.atoms[topology.angles[i].j].type,
				topology.atoms[topology.angles[i].k].type,
				&ff);
		topology.angles[i].ktheta = at.ktheta;
		topology.angles[i].theta0 = at.theta0;
		if(at.kub != 0){
			topology.ureyBradleyCount ++;
		}
	}
	topology.ureyBradley = (UreyBradley*)calloc(topology.ureyBradleyCount, sizeof(UreyBradley));
	topology.ureyBradleyCount = 0;
	for(i = 0; i < topology.angleCount; i++){
		topology.angles[i].i = findAtom(psf.angles[i].i);
		topology.angles[i].j = findAtom(psf.angles[i].j);
		topology.angles[i].k = findAtom(psf.angles[i].k);
		//printf("top: %d - %d - %d\n", topology.angles[i].i, topology.angles[i].j, topology.angles[i].k);
		FFAngleType at = findAngleType(
				topology.atoms[topology.angles[i].i].type,
				topology.atoms[topology.angles[i].j].type,
				topology.atoms[topology.angles[i].k].type,
				&ff);
		if(at.kub != 0){
			topology.ureyBradley[topology.ureyBradleyCount].i = findAtom(psf.angles[i].i);
			topology.ureyBradley[topology.ureyBradleyCount].j = findAtom(psf.angles[i].k);
			topology.ureyBradley[topology.ureyBradleyCount].kub = at.kub;
			topology.ureyBradley[topology.ureyBradleyCount].s0 = at.s0;

			/*printf("U-B: %d(%s)-%d(%s):\tkub=%f\ts0=%f\n",
					topology.ureyBradley[topology.ureyBradleyCount].i,
					topology.atoms[topology.angles[i].i].type,
					topology.ureyBradley[topology.ureyBradleyCount].j,
					topology.atoms[topology.angles[i].k].type,
					topology.ureyBradley[topology.ureyBradleyCount].kub,
					topology.ureyBradley[topology.ureyBradleyCount].s0);
*/
			topology.ureyBradleyCount ++;
		}
	}

	for(i = 0; i < topology.dihedralCount; i++){
		topology.dihedrals[i].i = findAtom(psf.dihedrals[i].i);
		topology.dihedrals[i].j = findAtom(psf.dihedrals[i].j);
		topology.dihedrals[i].k = findAtom(psf.dihedrals[i].k);
		topology.dihedrals[i].l = findAtom(psf.dihedrals[i].l);
		FFDihedralType dt[MAX_DIHEDRAL_MULTIPLICITY];
		int multiplicity = findDihedralType(
				topology.atoms[topology.dihedrals[i].i].type,
				topology.atoms[topology.dihedrals[i].j].type,
				topology.atoms[topology.dihedrals[i].k].type,
				topology.atoms[topology.dihedrals[i].l].type,
				&ff, dt);
		topology.dihedrals[i].multiplicity = multiplicity;
		if(multiplicity > 1){
			/*printf("Multiplicity is > 1 (%d) for dihedral (%s - %s - %s - %s) with parameters:\n",
					multiplicity,
					atoms[topology.dihedrals[i].i].type,
					atoms[topology.dihedrals[i].j].type,
					atoms[topology.dihedrals[i].k].type,
					atoms[topology.dihedrals[i].l].type);*/
		}
		for(j = 0; j < multiplicity; j++){
			topology.dihedrals[i].kchi[j] = dt[j].kchi;
			topology.dihedrals[i].n[j] = dt[j].n;
			topology.dihedrals[i].delta[j] = dt[j].delta;
			/*printf("%d: K = %f; n = %d; delta = %f; multiplicity = %d\n",
				i, dt[j].kchi, dt[j].n, dt[j].delta,
				topology.dihedrals[i].multiplicity);*/
		}

	}

	for(i = 0; i < topology.improperCount; i++){
		topology.impropers[i].i = findAtom(psf.impropers[i].i);
		topology.impropers[i].j = findAtom(psf.impropers[i].j);
		topology.impropers[i].k = findAtom(psf.impropers[i].k);
		topology.impropers[i].l = findAtom(psf.impropers[i].l);
		FFImproperType it[MAX_IMPROPER_MULTIPLICITY];
		int multiplicity = findImproperType(
				topology.atoms[topology.impropers[i].i].type,
				topology.atoms[topology.impropers[i].j].type,
				topology.atoms[topology.impropers[i].k].type,
				topology.atoms[topology.impropers[i].l].type,
				&ff, it);
		topology.impropers[i].multiplicity = multiplicity;
		if(multiplicity > 1){
		/*printf("Multiplicity is > 1 (%d) for improper (%s - %s - %s - %s) with parameters:\n",
				multiplicity,
				atoms[topology.impropers[i].i].type,
				atoms[topology.impropers[i].j].type,
				atoms[topology.impropers[i].k].type,
				atoms[topology.impropers[i].l].type);*/
		}
		for(j = 0; j < multiplicity; j++){
			topology.impropers[i].kpsi[j] = it[j].kpsi;
			topology.impropers[i].n[j] = it[j].n;
			topology.impropers[i].psi0[j] = it[j].psi0;
			/*printf("%d: K = %f; n = %d; delta = %f; multiplicity = %d\n",
				i, it[j].kpsi, it[j].n, it[j].psi0,
				topology.impropers[i].multiplicity);*/
		}
	}

	for(i = 0; i < topology.cmapCount; i++){
		topology.cmaps[i].i1 = findAtom(psf.cmaps[i].i1);
		topology.cmaps[i].j1 = findAtom(psf.cmaps[i].j1);
		topology.cmaps[i].k1 = findAtom(psf.cmaps[i].k1);
		topology.cmaps[i].l1 = findAtom(psf.cmaps[i].l1);
		topology.cmaps[i].i2 = findAtom(psf.cmaps[i].i2);
		topology.cmaps[i].j2 = findAtom(psf.cmaps[i].j2);
		topology.cmaps[i].k2 = findAtom(psf.cmaps[i].k2);
		topology.cmaps[i].l2 = findAtom(psf.cmaps[i].l2);
		FFCMAPType cmap = findCMAPType(
				topology.atoms[topology.cmaps[i].i1].type,
				topology.atoms[topology.cmaps[i].j1].type,
				topology.atoms[topology.cmaps[i].k1].type,
				topology.atoms[topology.cmaps[i].l1].type,
				topology.atoms[topology.cmaps[i].i2].type,
				topology.atoms[topology.cmaps[i].j2].type,
				topology.atoms[topology.cmaps[i].k2].type,
				topology.atoms[topology.cmaps[i].l2].type,
				&ff);
		//!TODO Should be CMAP type or something
	}

	j = 0;
	for(i = 0; i < topology.atomCount; i++){
		while(j < psf.nbExclusionsCounts[i]){
			topology.exclusions[j].i = i;
			topology.exclusions[j].j = findAtom(psf.nbExclusions[j]);
			j++;
		}
	}

	printf("Done initializing topology...\n");
}

void readCoordinates(char* filename){
	int fileType = getFileType(filename);
	double* x = (double*)calloc(topology.atomCount, sizeof(double));
	double* y = (double*)calloc(topology.atomCount, sizeof(double));
	double* z = (double*)calloc(topology.atomCount, sizeof(double));
	switch (fileType) {
		case FILETYPE_PDB:
			readCoordinatesFromPDB(filename, x, y, z, topology.atomCount);
			break;
		case FILETYPE_XYZ:
			readCoordinatesFromXYZ(filename, x, y, z, topology.atomCount);
			break;
		default:
			DIE("Unsupported file extension.\n");
	}
	int i;
	for(i = 0; i < topology.atomCount; i++){
		topology.atoms[i].x = x[i]/10.0; //Coordinates are in A in pdb file
		topology.atoms[i].y = y[i]/10.0;
		topology.atoms[i].z = z[i]/10.0;
	}

	free(x);
	free(y);
	free(z);
}

void readVelocities(char* filename){
	int fileType = getFileType(filename);
	double* x = (double*)calloc(topology.atomCount, sizeof(double));
	double* y = (double*)calloc(topology.atomCount, sizeof(double));
	double* z = (double*)calloc(topology.atomCount, sizeof(double));
	switch (fileType) {
		case FILETYPE_PDB:
			readCoordinatesFromPDB(filename, x, y, z, topology.atomCount);
			break;
		case FILETYPE_XYZ:
			readCoordinatesFromXYZ(filename, x, y, z, topology.atomCount);
			break;
		default:
			DIE("Unsupported file extension.\n");
	}
	int i;
	for(i = 0; i < topology.atomCount; i++){
		topology.atoms[i].vx = x[i]/10.0; //Coordinates are in A in pdb file
		topology.atoms[i].vy = y[i]/10.0;
		topology.atoms[i].vz = z[i]/10.0;
	}

	free(x);
	free(y);
	free(z);
}

void generateVelocities(float T){
	//printf("Generating velocities at temperature T=%fK.\n", T);
	int i;
	if(T < 0){
		DIE("Negative value for temperature is set (T = %fK).", T);
	} else
	if(T == 0){
		for(i = 0; i < topology.atomCount; i++){
			topology.atoms[i].vx = 0.0;
			topology.atoms[i].vy = 0.0;
			topology.atoms[i].vz = 0.0;
		}
	} else {
		for(i = 0; i < topology.atomCount; i++){
			double var = sqrt(Kb_MD*T/topology.atoms[i].mass);
			topology.atoms[i].vx = var*ran2::gasdev(&parameters.rseed);
			topology.atoms[i].vy = var*ran2::gasdev(&parameters.rseed);
			topology.atoms[i].vz = var*ran2::gasdev(&parameters.rseed);
		}
	}
}

int findAtom(int atomId){
	int id = atomIDs[atomId];
	if(id != -1){
		return id;
	} else {
		DIE("ERROR: Atom with id %d can't be found.", atomId);
	}
}

void printAtom(Atom atom){
	printf("Atom %d: %s (%s) from %s%d\n",
			atom.id,
			atom.name,
			atom.type,
			atom.resName,
			atom.resid);
}
