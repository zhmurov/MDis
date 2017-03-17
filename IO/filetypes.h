/*
 * filetypes.h
 *
 *  Created on: Apr 7, 2011
 *      Author: zhmurov
 */

#pragma once

#define FILETYPES_COUNT		6

#define FILETYPE_PDB		1
#define FILETYPE_XYZ		2
#define FILETYPE_GRO		3
#define FILETYPE_CRD		4

#define FILETYPE_DCD		10
#define FILETYPE_TRR		11

int getFileType(char* filename);
