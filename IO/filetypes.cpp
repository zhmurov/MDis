/*
 * filetypes.cpp
 *
 *  Created on: Apr 7, 2011
 *      Author: zhmurov
 */
#include "filetypes.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../Util/wrapper.h"

struct filetypes {
	const char* extension;
	int id;
} filetypestab[] = {
        {".pdb", FILETYPE_PDB},
        {".xyz", FILETYPE_XYZ},
        {".gro", FILETYPE_GRO},
        {".crd", FILETYPE_CRD},
        {".dcd", FILETYPE_DCD},
        {".trr", FILETYPE_TRR}
};


int getFileTypeByExtension(char* extension);

int getFileType(char* filename){
	char* extension = strrchr(filename, '.');
	//printf("%s extension is '%s'\n", filename, extension);
	return getFileTypeByExtension(extension);

}

int getFileTypeByExtension(char* extension){
	int i;
	for(i = 0; i < FILETYPES_COUNT; i++){
		if(strcmp(extension, filetypestab[i].extension) == 0){
			return filetypestab[i].id;
		}
	}
	DIE("ERROR: File type '%s' is not supported.", extension);
}
