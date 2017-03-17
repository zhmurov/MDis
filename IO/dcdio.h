/*
 * dcdio.h
 *
 *  Created on: Feb 27, 2009
 *      Author: zhmurov
 */

#pragma once

#include <stdio.h>

FILE* dcd_open_write(FILE* dcd_file, char *dcd_filename);
FILE* dcd_open_append(FILE* dcd_file, char *dcd_filename);
FILE* dcd_open_read(FILE* dcd_file, char *dcd_filename);

void dcd_write_header(FILE* dcd_file, char *dcd_filename, int N, int NFILE, int NPRIV, int NSAVC, double DELTA);
void dcd_write_frame(FILE* dcd_file, int N, float *X, float *Y, float *Z);
void dcd_read_header(FILE* dcd_file, int* N, int* NFILE, int* NPRIV, int* NSAVC, float* DELTA);
int dcd_read_frame(FILE* dcd_file, int N, float *X, float *Y, float *Z);

void dcd_close(FILE* dcd_file);
