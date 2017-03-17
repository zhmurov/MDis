/*
 * UmbrellaSampling.cu
 *
 *  Created on: Jul 22, 2011
 *      Author: alekseenko
 *  
 */
#include <Core/global.h>
#include <Util/Log.h>
#include <Util/wrapper.h>
#include <Util/mystl.h>
#include <Util/atomfilter.h>

#include "UmbrellaSampling.cuh"
#include <set>

namespace umbrella_sampling {


class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<umbrella_sampling> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

#include "UmbrellaSampling_com.cu"
//#include "UmbrellaSampling_rot.cu"

void create() {
	LOG << "create";
	if(getYesNoParameter(PARAMETER_UMBRELLA, 0)) {
		// Initialize all necessary structures...
		potential.destroy = &destroy;
		sprintf(potential.name, "Umbrella potential");
		debug = getYesNoParameter(PARAMETER_UMBRELLA_DEBUG, 0);
		init();
		// read parameters...
		std::string stage, method;
		stage = getParameterAs<std::string> (PARAMETER_UMBRELLA_STAGE);
		method = getParameterAs<std::string> (PARAMETER_UMBRELLA_METHOD);
		sampling_params.win_step = getFloatParameter (PARAMETER_UMBRELLA_WINSTEP);
		sampling_params.energyfile = getMaskedParameterAs<std::string> (PARAMETER_UMBRELLA_OUTFILE);
		for (int i = 0; i < parameters.Ntr; ++i)
			fclose(safe_fopen(string_replace(sampling_params.energyfile,"<run>",any2str(i + parameters.firstrun)).c_str(),"w")); // Clear file. Ugly way, but I'm too lazy to implement smth. better

		// Allocate some memory...
		allocateCPU((void**)&h_data.atomGroups, gsystem.Ntot*sizeof(int));
		allocateGPU((void**)&d_data.atomGroups, gsystem.Ntot*sizeof(int));
		allocateCPU((void**)&sampling_params.rcs, parameters.Ntr*sizeof(float4));
		allocateCPU((void**)&sampling_params.energy, parameters.Ntr*sizeof(float2));

		// Get K-spring for movement along reaction coordinate
		h_data.ks = d_data.ks = getFloatParameter(PARAMETER_UMBRELLA_KS);
		// Get K-spring for orthogonal movement. Default is 0 (orthogonal movement is unrestricted)
		h_data.ks_ortho = d_data.ks_ortho = getFloatParameter(PARAMETER_UMBRELLA_ORTHOGONAL_KS, 0.0f);
		// Set up updater freqency: how often we move cantilever (if it's moved) and how often we output energies.
		updater.frequency = getIntegerParameter(PARAMETER_UMBRELLA_FREQ);
		// Specify stage of sampling..
		if (stage == PARAMETER_VALUE_UMBRELLA_STAGE_PREPARE) {
			LOG << "Preparation (pulling) will be performed";
			sampling_params.stage = UMBRELLA_STAGE_PREPARATION;
			// Set initial RCs to zero
			for (int i = 0; i < parameters.Ntr; ++i)
				sampling_params.rcs[i] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
		} else {
			if (stage == PARAMETER_VALUE_UMBRELLA_STAGE_PRODUCTION) {
				LOG << "Production run will be performed";
				sampling_params.stage = UMBRELLA_STAGE_PRODUCTION;
				// Load trajectories (we just overwrite previous data...). Not very elegant, but works...
				for (int i = 0; i < parameters.Ntr; ++i) {
					int traj = i + parameters.firstrun;
					char trajFilename[1024];
					LOG << "Loading window # " << traj;
					getMaskedParameterWithReplacementT(trajFilename, PARAMETER_UMBRELLA_FILE_FRAME, PARAMETER_STRING_UNDEFINED, traj, "<run>");
					// Read coordinates to host buffer
					readCoordinates(trajFilename);
					// Screw velocities. Random ones are ok.
					float T = getFloatParameter(PARAMETER_INITIAL_TEMPERATURE, -1.0f);
					if(T == -1.0f){
						T = getFloatParameter(PARAMETER_TEMPERATURE, 0.0f);
					}
					generateVelocities ( T );
					// Copy data to appropriate place on GPU
					copyCoordinatesToGPU(i, 1);
					copyVelocitiesToGPU(i, 1);
					sampling_params.rcs[i] = make_float4(0.0f, 0.0f, 0.0f, traj); // .w-th component is window number
				}
			} else {
				LOG << "Error: Stage parameter should be either '" << PARAMETER_VALUE_UMBRELLA_STAGE_PREPARE << "' for generating initial positions" << \
					"or '" << PARAMETER_VALUE_UMBRELLA_STAGE_PRODUCTION << "' for production run. It is '" << stage << "'";
				DIE("Wrong parameter!");
			}
		}
		// Prepare run for specific method
		if (method == PARAMETER_VALUE_UMBRELLA_METHOD_COM) {
			initCoM();
		/* else if (method == PARAMETER_VALUE_UMBRELLA_METHOD_ROT) {
			initRot(); */
		} else {
			LOG << "Error: Umbrella sampling supports only '" << PARAMETER_VALUE_UMBRELLA_METHOD_COM << "' method right now. We're very sorry.";
			DIE("Wrong parameter!");
		}
		cudaMemcpy(d_data.atomGroups, h_data.atomGroups, gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);
		checkCUDAError("Copying atom groups to device");
		// Copy common data...
//		cudaMemcpyToSymbol("umbrella_sampling_c_data", &d_data, sizeof(Data), 0, cudaMemcpyHostToDevice);
//		checkCUDAError("Initializing c_data");

		// Announce that we have updater and potential
		updaters[updatersCount++] = &updater;
		potentials[potentialsCount++] = &potential;
	} else {
		LOG << "Umbrella sampling is disabled";
	}
	LOG << "Done initializing";
}

void init() {
	blockSize = BLOCK_SIZE;
	blockCount = gsystem.Ntot/BLOCK_SIZE + 1;
	blockCountTr = gsystem.N / BLOCK_SIZE + 1;
}

/*
//   Unused now, can be used later... Somewhere.. May be...
template <typename T>
__device__ __host__ void findPerpendiculars(const <T> &in, <T> &out1, <T> &out2) {
#define X in.x
#define Y in.y
#define Z in.z
    out1 = T(-Y, X, 0.0f);
    if (X == Y && Y == 0.0f) {
        out1 = T(1.0f, 1.0f, 0.0f);
    }
    out2 = T(-X*Z, -Y*Z, X*X + Y*Y);
#undef X
#undef Y
#undef Z
}
*/


void destroy() {
	// Freeing allocated memory is too mainstream
}

void updaterDestroy() {}

#undef LOG

} // namespace umbrella_sampling
