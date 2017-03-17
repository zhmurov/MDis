/*
 * gbswradiiio.h
 *
 *  Created on: Oct 7, 2011
 *      Author: zhmurov
 */

#ifndef GBSWRADIIIO_H_
#define GBSWRADIIIO_H_

typedef struct {
	char name[10];
	char resName[10];
	float value;
} PBRFileEntry;

typedef struct {
	int rCount;
	PBRFileEntry* r;
} PBRFileData;

void readPBRadiiFile(char* filename, PBRFileData* data);
int comparePBRadiiWithWC(char* str, char* mask);

#endif /* GBSWRADIIIO_H_ */
