#pragma once

namespace periodic_boundary {

#ifdef USE_PBC
__device__ __constant__ float4 c_gpbc; // Size of periodic boundary
#endif

void create();
void init();

} // namespace periodic_boundary 
