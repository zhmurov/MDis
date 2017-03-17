/*
 * UmbrellaSampling_com.cu
 *
 *  Created on: Jul 22, 2011
 *      Author: alekseenko
 */

void getCoMFiltered(int group_id, float4 *res, bool oneres = false, bool force_cpu = false);

void initCoM() {
	sampling_params.method = UMBRELLA_METHOD_COM;
	allocateCPU((void**)&com_params.positions, parameters.Ntr*sizeof(float4), 0);
	com_params.use_gpu = getYesNoParameter(PARAMETER_UMBRELLA_COM_USEGPU, 0);
	if (com_params.use_gpu) {
		LOG << "Using GPU for CoM calculations...";
		allocateGPU((void**)&com_params.d_com, blockCountTr * sizeof(float4));
		allocateCPU((void**)&com_params.h_com, blockCountTr * sizeof(float4));
	}
	// Alternative behaviour of fixed atoms
	// Instead of fixing each atom, we fix their center of mass
	com_params.fix_com = getYesNoParameter(PARAMETER_UMBRELLA_COM_FIX_COM, 0);

	// Set atom groups...
	// Atom group specify whether atom in question will be pulled, fixed or ignored
	std::ostringstream oss_f, oss_p;
	// Configure fixed atoms
	if (hasParameter(PARAMETER_UMBRELLA_COM_FIXED_MASK)) {
		oss_f << getParameterAs<std::string>(PARAMETER_UMBRELLA_COM_FIXED_MASK);
	} else {
		oss_f << "segid " << getParameterAs<std::string>(PARAMETER_UMBRELLA_COM_CHAIN_FIXED)[0] << " and name CA ";
		if (!getYesNoParameter(PARAMETER_UMBRELLA_COM_FIXED_ALL, 0)) {
			oss_f << " and resid " << getParameterAs<std::string>(PARAMETER_UMBRELLA_COM_FIXED_LIST);
		}
	}
	// Configure pulled atoms
	if (hasParameter(PARAMETER_UMBRELLA_COM_PULLED_MASK)) {
		oss_p << getParameterAs<std::string>(PARAMETER_UMBRELLA_COM_PULLED_MASK);
	} else {
		oss_p << "segid " << getParameterAs<std::string>(PARAMETER_UMBRELLA_COM_CHAIN_PULLED)[0];
	}
	// Masquerade
	AtomFilter af_f(oss_f.str()), af_p(oss_p.str());
	memset(h_data.atomGroups, UMBRELLA_ATOM_GROUP_OTHER, gsystem.N * sizeof(int));
	af_f.masquerade(topology.atoms, gsystem.N, h_data.atomGroups, UMBRELLA_ATOM_GROUP_FIXED);
	af_p.masquerade(topology.atoms, gsystem.N, h_data.atomGroups, UMBRELLA_ATOM_GROUP_PULLED);
#ifdef DEBUG
	LOG << "Atom states:";
	for (int i = 0; i < gsystem.N; ++i)
		LOG << "  " << topology.atoms[i].id << "\t" << topology.atoms[i].name << "\t" << topology.atoms[i].segment << topology.atoms[i].resid << "\t" << topology.atoms[i].resName << "\t" << (h_data.atomGroups[i] == UMBRELLA_ATOM_GROUP_PULLED ? "PUL" : (h_data.atomGroups[i] == UMBRELLA_ATOM_GROUP_FIXED ? "FIX" : "==="));
#endif
	// Repeat for all trajectories
	for (int i = 1; i < parameters.Ntr; ++i) memcpy(h_data.atomGroups + gsystem.N * i, h_data.atomGroups, gsystem.N * sizeof(int));

	// Make some statistics
	int count_f = 0, count_p = 0;
	for (int i = 0; i < gsystem.N; ++i) {
		if (h_data.atomGroups[i] == UMBRELLA_ATOM_GROUP_FIXED)  count_f ++;
		if (h_data.atomGroups[i] == UMBRELLA_ATOM_GROUP_PULLED) count_p ++;
	}
	LOG << "Fixed: " << count_f << " atoms, pulled: " << count_p << " atoms";

	h_data.ks_rest = getFloatParameter(PARAMETER_UMBRELLA_RESTRAIN_KS); // kJ/mol/nm^2

	FILE *fileTemp; // file specified by PARAMETER_UMBRELLA_FILE_COM_TEMP and used for preserving some essential info between two stages of umbrella sampling
#define XYZ_CS(a) a.x, a.y, a.z
	// Set up the pulling direction...
	if (sampling_params.stage == UMBRELLA_STAGE_PREPARATION) {
		// Get initial positions...
		if (getYesNoParameter(PARAMETER_UMBRELLA_COM_LOAD_VECTOR, 0)) { // Load from file
			fileTemp = safe_fopen(getMaskedParameterAs<std::string>(PARAMETER_UMBRELLA_FILE_COM_TEMP).c_str(), "r");
			if (fscanf(fileTemp, "%f %f %f %f %f %f", XYZ_CS(&com_params.initial_fixpos), XYZ_CS(&com_params.initial_pos)) != 6)
				DIE("Loading data from <##PARAMETER_UMBRELLA_FILE_COM_TEMP##> file");
			fclose(fileTemp);
		} else { // Get from first trajectory
			float4 com_fixed, com_pulled;
			getCoMFiltered(UMBRELLA_ATOM_GROUP_PULLED, &com_pulled, true, true);
			getCoMFiltered(UMBRELLA_ATOM_GROUP_FIXED, &com_fixed, true, true);
			com_params.initial_pos = com_pulled - com_fixed;
			com_params.initial_fixpos = com_fixed;
			// save initial CoM-CoM vector to temp-file...
			fileTemp = safe_fopen(getMaskedParameterAs<std::string>(PARAMETER_UMBRELLA_FILE_COM_TEMP).c_str(), "w");
			fprintf(fileTemp, "%f %f %f %f %f %f\n", XYZ_CS(com_params.initial_fixpos), XYZ_CS(com_params.initial_pos));
			fclose(fileTemp);
		}
		LOG << "Initial vector : " << com_params.initial_pos;
		// Set up initial positions
		for (int i = 0; i < parameters.Ntr; ++i)
			com_params.positions[i] = com_params.initial_pos;
	} else { // UMBRELLA_STAGE_PRODUCTION
		if (sampling_params.stage != UMBRELLA_STAGE_PRODUCTION)
			DIE("What? I don't even..."); // It will print line number anyway, so why bother with verbose error messages?
		// Loading initial positions from file...
		fileTemp = safe_fopen(getMaskedParameterAs<std::string>(PARAMETER_UMBRELLA_FILE_COM_TEMP).c_str(), "r");
		if (fscanf(fileTemp, "%f %f %f %f %f %f", XYZ_CS(&com_params.initial_fixpos), XYZ_CS(&com_params.initial_pos)) != 6)
			DIE("Loading data from <##PARAMETER_UMBRELLA_FILE_COM_TEMP##> file");
		fclose(fileTemp);
	}
#undef XYZ_CS
	// Initialize pulling direction and velocity
	if (hasParameter(PARAMETER_UMBRELLA_DIRECTION)) {
		h_data.d = getFloat3Parameter(PARAMETER_UMBRELLA_DIRECTION);
	} else {
		LOG << "Using initial CoM-CoM vector as pulling direction";
		h_data.d = com_params.initial_pos;
	}
	h_data.d.w = 0.0f;
	h_data.v = normalize(h_data.d);
	if (hasParameter(PARAMETER_UMBRELLA_VELOCITY)) {
		h_data.v = getFloatParameter(PARAMETER_UMBRELLA_VELOCITY);
	}
	LOG << "Pulling with velocity v = " << h_data.d * h_data.v << " (|v| = " << h_data.v << " um/s = " << h_data.v * 1e-9 << " nm/ps)";
	d_data.d = h_data.d;
	d_data.v = h_data.v;

	com_params.rc_zero = getFloatParameter(PARAMETER_UMBRELLA_COM_COORD_ZERO, 0.0f); // Use this to artificially displace initial "position"
	for (int i = 0; i < parameters.Ntr; ++i)
		com_params.positions[i] = com_params.initial_pos + h_data.d * com_params.rc_zero;
	if (sampling_params.stage == UMBRELLA_STAGE_PRODUCTION) {
		// Set up initial positions...
		for (int i = 0; i < parameters.Ntr; ++i) {
			com_params.positions[i] = com_params.initial_pos + h_data.d * (i + parameters.firstrun - 1) * sampling_params.win_step;
			LOG << "Position for traj #" << i << " : " << com_params.positions[i];
		}
	}

	//allocate some memory
	allocateCPU((void**)&com_params.h_dfp, parameters.Ntr * sizeof(float4));
	allocateGPU((void**)&com_params.d_dfp, parameters.Ntr * sizeof(float4));
	// Create copy of initial coordinates for applying restraints on receptor
	allocateGPU((void**)&com_params.d_initial_coords, gsystem.Ntot * sizeof(float4));
	cudaMemcpy(com_params.d_initial_coords, gsystem.d_coord, gsystem.Ntot * sizeof(float4), cudaMemcpyDeviceToDevice);
	checkCUDAError("Making copy of initial positions");
	potential.compute = computeUmbrellaCoM;
	updater.update = updateUmbrellaCoM;
	updater.destroy = updaterDestroy;
	LOG << "Umbrella-CoM initialized";
}


/* 
 * This kernel adds additional force, which results in same additional acceleration for all atoms in the group
 *
 * @param a Array float4[Ntr], which specifys the acceleration to be given to chosen atom group for each trajectory. .w-component is unused
 * @param flags Array int[Ntot] filled with flags. If flag matches flag_ok, then the atom is considered chosen
 * @param flag_ok Flag value, which marks atoms chosen for procedure
 *
 */
__global__ void computeAcceleration_kernel(float4 *a, const int *flags, int flag_ok) {
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot) {
		if (flags[d_i] == flag_ok) {
			const int trid = d_i / c_gsystem.N;
			const int type = tex1Dfetch(t_coord, d_i).w;
			c_gsystem.d_forces[d_i] += a[trid] * tex1Dfetch(t_m, type);
		}
	}
}

/* 
 * This kernel applies harmonic restraint to all atoms marked with flag to keep them in desired position
 *
 * @param coord0 Array float4[Ntot], with desired coordinates for chosen atoms.
 * @param flags Array int[Ntot] filled with flags. If flag matches flag_ok, then the atom is considered chosen
 * @param flag_ok Flag value, which marks atoms chosen for procedure
 * @param ks Spring constant, kJ/nm^2
 *
 */
__global__ void computeRestraints_kernel(const float4 *coords0, const int *flags, int flag_ok, float ks) {
	int d_i = blockIdx.x * blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot) {
		if (flags[d_i] == flag_ok) {
			float4 dr = coords0[d_i] - tex1Dfetch(t_coord, d_i);
			c_gsystem.d_forces[d_i] += dr * ks;
		}
	}
}

void computeUmbrellaCoM(int dry) {
	if (!com_params.use_gpu)
		copyCoordinatesFromGPU();
	std::vector<float4> com_f(parameters.Ntr), com_p(parameters.Ntr);
	getCoMFiltered(UMBRELLA_ATOM_GROUP_FIXED, &com_f[0]);
	getCoMFiltered(UMBRELLA_ATOM_GROUP_PULLED, &com_p[0]);
	// Calculate accelerations..
	for (int i = 0; i < parameters.Ntr; ++i) {
		com_params.h_dfp[i] = (com_params.positions[i] - (com_p[i] - com_f[i]));
		com_params.h_dfp[i].w = 0.0f;
		if (step % 100 == 0 && debug)
			LOG << "[" << step << "] CoM distance = " << (com_p[i] - com_f[i]) << " || " << com_params.positions[i];
		// And now some linear algebra
		float4 dfp_d = h_data.d * dot(h_data.d, com_params.h_dfp[i]); // Projection to h_data.d
		float4 dfpXd = com_params.h_dfp[i] - dfp_d; // Projection on the plane perpendicular to h_data.d
		sampling_params.energy[i].x = h_data.ks       * abs2(dfp_d) / 2;
		sampling_params.energy[i].y = h_data.ks_ortho * abs2(dfpXd) / 2;
		com_params.h_dfp[i] = dfp_d * h_data.ks + dfpXd * h_data.ks_ortho;
		com_params.h_dfp[i] /= com_p[i].w; // Actually, we need to pass acceleration to the kernel, not force...
		// sampling_params.rcs[i].x = abs(dfp_d);
		// sampling_params.rcs[i].y = abs(dfpXd);
		sampling_params.rcs[i].x = abs(com_p[i] - com_f[i]);  // dCoM
		sampling_params.rcs[i].y = dot(com_p[i] - com_f[i], h_data.d); 
		sampling_params.rcs[i].z = dot(com_p[i] - com_f[i] - com_params.initial_pos, h_data.d);
	}
	if (!dry) {
		// Copy acceleration to device
		cudaMemcpy(com_params.d_dfp, com_params.h_dfp, sizeof(float4) * parameters.Ntr, cudaMemcpyHostToDevice);
		checkCUDAError("Copying dfp to GPU");
		//Run kernels
		computeAcceleration_kernel<<<blockCount, blockSize>>>(com_params.d_dfp, d_data.atomGroups, UMBRELLA_ATOM_GROUP_PULLED);
		checkCUDAError("computeAcceleration_kernel (for pulled)");
		if (com_params.fix_com) { // We fix fixed atoms' CoM
			for (int i = 0; i < parameters.Ntr; ++i) {
				float4 dx = com_params.initial_fixpos - com_f[i];
				com_params.h_dfp[i] = dx * h_data.ks_rest;
			}
			cudaMemcpy(com_params.d_dfp, com_params.h_dfp, sizeof(float4) * parameters.Ntr, cudaMemcpyHostToDevice);
			computeAcceleration_kernel<<<blockCount, blockSize>>>(com_params.d_dfp, d_data.atomGroups, UMBRELLA_ATOM_GROUP_FIXED);
			checkCUDAError("computeAcceleration_kernel (for fixed)");
		} else { // We simply fix each fixed atom individually
			computeRestraints_kernel<<<blockCount, blockSize>>>(com_params.d_initial_coords, d_data.atomGroups, UMBRELLA_ATOM_GROUP_FIXED, h_data.ks_rest);
			checkCUDAError("computeRestraints_kernel (for fixed)");
		}
	}
}

void updateUmbrellaCoM() {
	if (debug) LOG << "Updater called!";
	computeUmbrellaCoM(1);
	std::vector<float> rc0(parameters.Ntr);
	if (sampling_params.stage == UMBRELLA_STAGE_PREPARATION) {
        for (int i = 0; i < parameters.Ntr; ++i) {
		rc0[i] = h_data.v * step * 1e-9 * integrator->h;
		com_params.positions[i] = com_params.initial_pos + h_data.d * rc0[i];
	}
	} else {
		for (int i = 0; i < parameters.Ntr; ++i)
			rc0[i] = (i + parameters.firstrun - 1) * sampling_params.win_step + com_params.rc_zero;
	}
	printf("%*s%-*s%*s%*s%*s%*s%*s%*s%*s%*s\n",
			PULLING_OUTPUT_WIDTH, "",
			PULLING_OUTPUT_WIDTH, PULLING_OUTPUT_NAME_TRAJECTORY,
			PULLING_OUTPUT_WIDTH, PULLING_OUTPUT_NAME_STEP,
			PULLING_OUTPUT_WIDTH, "E_p",
			PULLING_OUTPUT_WIDTH, "E_o",
			PULLING_OUTPUT_WIDTH, "CoM_RC",
			PULLING_OUTPUT_WIDTH, "CoM_RC_0",
			PULLING_OUTPUT_WIDTH, "CoM",
			PULLING_OUTPUT_WIDTH, "CoM_p",
			PULLING_OUTPUT_WIDTH, "CoM_p_0");
	for (int i = 0; i < parameters.Ntr; ++i) {
		printf("%-*s%-*d%*lld%*f%*f%*f%*f%*f%*f%*f\n",
			PULLING_OUTPUT_WIDTH, "",
			PULLING_OUTPUT_WIDTH, i + parameters.firstrun,
			PULLING_OUTPUT_WIDTH, step,
			PULLING_OUTPUT_WIDTH, sampling_params.energy[i].x * KCALL_PER_KJ,
			PULLING_OUTPUT_WIDTH, sampling_params.energy[i].y * KCALL_PER_KJ,
			PULLING_OUTPUT_WIDTH, sampling_params.rcs[i].z,
			PULLING_OUTPUT_WIDTH, rc0[i],
			PULLING_OUTPUT_WIDTH, sampling_params.rcs[i].x,
			PULLING_OUTPUT_WIDTH, sampling_params.rcs[i].y,
			PULLING_OUTPUT_WIDTH, dot(com_params.positions[i],h_data.d));
		FILE *fout;
		fout = safe_fopen(string_replace(sampling_params.energyfile,"<run>",any2str(i + parameters.firstrun)).c_str(), "a");
		fprintf(fout, "%-*s%-*d%*lld%*f%*f%*f%*f%*f%*f%*f\n",
			PULLING_OUTPUT_WIDTH, "",
			PULLING_OUTPUT_WIDTH, i + parameters.firstrun,
			PULLING_OUTPUT_WIDTH, step,
			PULLING_OUTPUT_WIDTH, sampling_params.energy[i].x,
			PULLING_OUTPUT_WIDTH, sampling_params.energy[i].y,
			PULLING_OUTPUT_WIDTH, sampling_params.rcs[i].z,
			PULLING_OUTPUT_WIDTH, rc0[i],
			PULLING_OUTPUT_WIDTH, sampling_params.rcs[i].x,
			PULLING_OUTPUT_WIDTH, sampling_params.rcs[i].y,
			PULLING_OUTPUT_WIDTH, dot(com_params.positions[i],h_data.d));
		fclose(fout);
	}
}

// This code has been taken from nVidia CUDA SDK 3.1 and rewritten to suit our needs
/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and
 * proprietary rights in and to this software and related documentation.
 * Any use, reproduction, disclosure, or distribution of this software
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA)
 * associated with this source code for terms and conditions that govern
 * your use of this NVIDIA software.
 *
 */
// I've used reduce2 kernel, because it's rather fast and is easy to incorporate filtering into

__global__ void getCoMFiltered_kernel(int group_id, int first_index, int end_index, float4 *ret, const int *groups) {
	__shared__ float4 mr[BLOCK_SIZE];

	const unsigned int tid = threadIdx.x;
	const unsigned int i = blockIdx.x*BLOCK_SIZE + threadIdx.x + first_index;
	mr[tid] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);

	if (i < end_index) {
		if (groups[i] == group_id) {
			float4 coords = tex1Dfetch(t_coord, i);
			float m = tex1Dfetch(t_m, coords.w);
			coords *= m; coords.w = m;
			mr[tid] = coords;
		}
	}
	__syncthreads();

	// do reduction in shared mem
	for(unsigned int s = BLOCK_SIZE/2; s > 0; s >>= 1) {
		if (tid < s) {
			mr[tid] += mr[tid + s];
		}
		__syncthreads();
	}

	// write result for this block to global mem
	if (tid == 0) ret[blockIdx.x] = mr[0];
}

void getCoMFiltered(int group_id, float4 *res, bool oneres, bool force_cpu) {
	bool gpu = com_params.use_gpu && !force_cpu;
	for (int tr = 0; tr < (oneres ? 1 : parameters.Ntr); tr++) res[tr] = float4();
	if (gpu) {
		for (int tr = 0; tr < (oneres ? 1 : parameters.Ntr); tr++) {
			int idx = tr * gsystem.N;
			getCoMFiltered_kernel<<<blockCountTr, BLOCK_SIZE>>> (group_id, idx, idx + gsystem.N, com_params.d_com, d_data.atomGroups);
			checkCUDAError("getCoMFiltered_kernel");
			cudaMemcpy(com_params.h_com, com_params.d_com, blockCountTr * sizeof(float4), cudaMemcpyDeviceToHost);
			checkCUDAError("cudaMemcpy after getCoMFiltered_kernel");
			for (int i = 0; i < blockCountTr; ++i)
				res[tr] += com_params.h_com[i];
		}
	} else {
		for (int i = 0; i < (oneres ? gsystem.N : gsystem.Ntot); ++i) {
			int ntr = i / gsystem.N;
			int idx = i - ntr * gsystem.N;
			if (h_data.atomGroups[idx] == group_id) {
				float m = gsystem.h_m[gsystem.h_atomTypes[i]];
				float4 t = gsystem.h_coord[i] * m;
				t.w = m;
				res[ntr] += t;
			}
		}
	}
	for(int tr = 0; tr < (oneres ? 1 : parameters.Ntr); tr++) { 
		float m = res[tr].w;
		if (m == 0) DIE("For trajectory %d and group %d we got zero mass, which is not good", tr, group_id);
		res[tr] /= m;
		res[tr].w = m;
	}
}
