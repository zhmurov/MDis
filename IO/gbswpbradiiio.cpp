/*
 * gbswradiiio.cpp
 *
 *  Created on: Oct 7, 2011
 *      Author: zhmurov
 */
#include "gbswpbradiiio.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>

#define buf_size 256

void readPBRadiiFile(char* filename, PBRFileData* data){
	FILE* file = fopen(filename, "r");

	data->rCount = 0;
	char buffer[buf_size];
	while(fgets(buffer, buf_size, file) != NULL){
		if(buffer[0] != '!' && buffer[0] != ' ' && buffer[0] != '\n' && buffer[0] != '\t'){
			//printf("%s\n", buffer);
			data->rCount++;
		}
	}
	data->r = (PBRFileEntry*)calloc(data->rCount, sizeof(PBRFileEntry));
	int i = 0;
	rewind(file);
	while(fgets(buffer, buf_size, file) != NULL){
		if(buffer[0] != '!' && buffer[0] != ' ' && buffer[0] != '\n' && buffer[0] != '\t'){
			char* pch = strtok(buffer, " \t");
			strcpy(data->r[i].name, pch);
			pch = strtok(NULL, " \t\n");
			strcpy(data->r[i].resName, pch);
			pch = strtok(NULL, " \t\n");
			data->r[i].value = atof(pch);
			printf("%s\t%s\t%f\n", data->r[i].name, data->r[i].resName, data->r[i].value);
			i++;
		}
	}

	fclose(file);
}

int comparePBRadiiWithWC(char* str, char* mask){
	unsigned int i;
	unsigned int len = strlen(str);
	if(len == strlen(mask)){
		for(i = 0; i < len-1; i++){
			if(str[i] != mask[i]){
				return 0;
			}
		}
		if(str[len-1] == mask[len-1] || mask[len-1] == '*'){
			return 1;
		} else {
			return 0;
		}
	}
	if(len > strlen(mask)){
		for(i = 0; i < strlen(mask)-1; i++){
			if(str[i] != mask[i]){
				return 0;
			}
		}
		if(mask[strlen(mask)-1] == '*'){
			return 1;
		}
	}
	if(len == strlen(mask)-1){
		for(i = 0; i < len; i++){
			if(str[i] != mask[i]){
				return 0;
			}
		}
		if(mask[strlen(mask)-1] == '*'){
			return 1;
		}
	}
	return 0;
}
