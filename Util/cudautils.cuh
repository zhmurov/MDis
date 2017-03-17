/*
 * cudautils.cuh
 *
 *  Created on: Apr 14, 2011
 *      Author: zhmurov
 */

#pragma once

#define REDUCE_BLOCKSIZE	256

float* d_sums;
float* h_sums;

float reduce(float* d_data, int N);
