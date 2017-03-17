/*
 * PeriodicBoundary.cu
 *
 *  Created on: Jul 12, 2011
 *      Author: alekseenko
 */

#include "../Util/Log.h"
#include "PeriodicBoundary.cuh"

namespace periodic_boundary {
	
class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<periodic_boundary> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

#ifdef USE_PBC
#define DO_PBC(x) d_applyPBC(x)
template<typename T>
__inline__ __device__ __host__ void d_applyPBC(T &a) {
	if (c_gpbc.w == 0.0f) return;
#define NORMALIZE_COMPONENT(Q) \
	if (a.Q >= c_gpbc.Q) \
		a.Q -= 2*c_gpbc.Q; \
	else if (a.Q <= -c_gpbc.Q) \
		a.Q += 2*c_gpbc.Q;
	NORMALIZE_COMPONENT(x);
	NORMALIZE_COMPONENT(y);
	NORMALIZE_COMPONENT(z);
#undef NORMALIZE_COMPONENT
}
#else
#define DO_PBC(x)
#endif

void create() {
	init();
}

void init() {
#ifdef USE_PBC
	LOG << "Initializing Periodic boundary Conditions (PBC)";
	float4 h_pbc;
	if(getYesNoParameter(PARAMETER_PBC_ON, 0)){
		h_pbc.w = 1.0f;
		h_pbc.x = getFloatParameter(PARAMETER_PBC_LX);
		h_pbc.y = getFloatParameter(PARAMETER_PBC_LX, h_pbc.x);
		h_pbc.z = getFloatParameter(PARAMETER_PBC_LX, h_pbc.x);
		LOG << "Using PBC with box " << make_float3(h_pbc);
		float pbc_min = min(min(h_pbc.x, h_pbc.y), h_pbc.z);
		if (pbc_min < getFloatParameter(PARAMETER_POSSIBLEPAIRS_CUTOFF)) {
			DIE("Periodic Box size is less than %s! Since it certainly won't work, terminating...", PARAMETER_POSSIBLEPAIRS_CUTOFF);
		}
		if (getYesNoParameter(PARAMETER_REPULSIVE_BOUNDARY_ON, 0)){
			DIE("Periodic Boundary Conditions are used simultaneously with Repulsive Boundary. Most likely, you are doing something wrong. Terminating. But if you really need such ridiculous thing, patch the source code yourself, it's somewhere around %s:%d. I think the best way is to add some parameter to allow such setup, but I think it's highly unlikely someone will ever need it.\n", __FILE__, __LINE__);
		}
	} else {
		h_pbc.x = h_pbc.y = h_pbc.z = h_pbc.w = 0.0f;
	}
	cudaMemcpyToSymbol(c_gpbc, &h_pbc, sizeof(float4), 0, cudaMemcpyHostToDevice);
	checkCUDAError("init PBC");
#else
	if(getYesNoParameter(PARAMETER_PBC_ON, 0)){
		DIE("You have specified to use PBC in your config file. However, this binary is built without pbc support. Please, use appropriately built binary file, or set 'pbc off' in your config.\n");
	}
#endif
}

#undef LOG

} // namespace periodic_boundary 
