/*
 * global.h
 *
 *  Created on: Jul 28, 2010
 *      Author: zhmurov
 */

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>

#include "../Util/physconstants.h"
#include "../IO/pdbio.h"
#include "../IO/xyzio.h"
#include "../IO/psfio.h"
#include "../IO/topio.h"
#include "../IO/dcdio.h"
#include "../IO/configreader.h"
#include "../Core/topology.h"
#include "../Core/parameters.h"
#include "../Core/forcefield.h"
#include "../Core/ffutils.h"
#include "../Core/atomTypes.h"
#include "../Core/angleTypes.h"
#include "../Util/timer.h"
#include "../IO/filetypes.h"


extern char* configFile;

extern Parameters parameters;
extern ForceField ff;

/*
// This optimization turned out to be useless, so never use it
#if __CUDACC__
# if CUDA_VERSION >= 4000 && __CUDA_ARCH__ >= 200
#  define CUDA_USE_L1 
# else
#  if __CUDA_ARCH__ >= 200
#   warning For devices with CC >= 2.0 it is recommended to use CUDA 4.0 or later
#  endif
# endif
#endif
*/
