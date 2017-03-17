/*
 * FixForcePotential.cu
 *
 *  Created on: Nov 14, 2010
 *      Author: zhmurov
 */
#include "../Core/global.h"
#include "../Util/Log.h"
#include "FixForcePotential.cuh"

namespace fixforce_potential {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<fixforce_potential> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void create(){
	LOG << "Initializing Fix-Force potential...";
	if(getYesNoParameter(PARAMETER_PULLING, DEFAULT_PULLING)){
		fixForcePotential.destroy = &destroy;
		sprintf(fixForcePotential.name, "External force");
		potentials[potentialsCount] = &fixForcePotential;
		potentialsCount ++;
		init();
		char pullProtocol[PARAMETER_LENGTH];
		getMaskedParameter(pullProtocol, PARAMETER_PULLING_PROTOCOL);
		if(strcmp(pullProtocol, PARAMETER_VALUE_PULLING_PROTOCOL_FCONST) == 0){
			LOG << "Constant force pulling will be performed.";
			fixForcePotential.compute = &computeFConst;
			initFConst();
		} else
		if(strcmp(pullProtocol, PARAMETER_VALUE_PULLING_PROTOCOL_FRAMP) == 0){
			LOG << "Constant speed (force-ramp) pulling will be performed.";
			fixForcePotential.compute = &computeSMD;
			smdForceUpdater.update = updatePullingForceSMD;
			smdForceUpdater.destroy = destroySMDForceUpdater;
			smdForceUpdater.frequency = getIntegerParameter(PARAMETER_PULLING_SMD_UPDATE_FREQ);
			sprintf(smdForceUpdater.name, "SMD Force updater");
			updaters[updatersCount] = &smdForceUpdater;
			updatersCount ++;
			initSMD();
		} else {
			DIE("Pulling protocol parameter should be either %s for constant force pulling"
					"or %s for force-ramp (constant speed) pulling.",
					PARAMETER_VALUE_PULLING_PROTOCOL_FCONST,
					PARAMETER_VALUE_PULLING_PROTOCOL_FRAMP);
		}
	} else {
		LOG << "No atoms will be fixed/pulled.";
	}
	LOG << "Done initializing Fix-Force potential.";
}

void init(){
	fixForceBlockSize = BLOCK_SIZE;
	fixForceBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;
	allocateCPU((void**)&fixForceData.h_atomMasks, gsystem.Ntot*sizeof(int));
	allocateGPU((void**)&fixForceData.d_atomMasks, gsystem.Ntot*sizeof(int));
	allocateCPU((void**)&fixForceData.h_extForces, gsystem.Ntot*sizeof(float4));
	allocateGPU((void**)&fixForceData.d_extForces, gsystem.Ntot*sizeof(float4));
	PDB fixForceRefPDBData;
	char fixForcePDBFilename[100];
	int traj, i, itot;
	char trajnum[10];
	for(traj = 0; traj < parameters.Ntr; traj++){
		sprintf(trajnum, "%d", traj + parameters.firstrun);
		getMaskedParameterWithReplacement(fixForcePDBFilename, PARAMETER_PULLING_REFFILE,
				trajnum, "<run>");
		readPDB(fixForcePDBFilename, &fixForceRefPDBData);
		for(i = 0; i < gsystem.N; i++){
			itot = traj*gsystem.N + i;
			if(fixForceRefPDBData.atoms[i].beta > 0){
				DPRINTF("Atom #%d (%s%d-%s) in trajectory #%d will be fixed.\n",
					i, topology.atoms[i].resName, topology.atoms[i].resid, topology.atoms[i].name, traj+parameters.firstrun);
				fixForceData.h_atomMasks[itot] = FIXFORCE_ATOM_FIXED;
			} else
			if(fixForceRefPDBData.atoms[i].occupancy > 0){
				fixForceData.h_extForces[itot].x = fixForceRefPDBData.atoms[i].x;
				fixForceData.h_extForces[itot].y = fixForceRefPDBData.atoms[i].y;
				fixForceData.h_extForces[itot].z = fixForceRefPDBData.atoms[i].z;
				fixForceData.h_extForces[itot].w = fixForceRefPDBData.atoms[i].occupancy;
				fixForceData.h_atomMasks[itot] = FIXFORCE_ATOM_PULLED;
			} else {
				fixForceData.h_atomMasks[itot] = FIXFORCE_ATOM_FREE;
			}
		}
	}
	cudaMemcpyToSymbol(c_fixForceData, &fixForceData,
				sizeof(FixForceData), 0, cudaMemcpyHostToDevice);
}

void initFConst(){
	int traj, i, itot;
	for(traj = 0; traj < parameters.Ntr; traj++){
		for(i = 0; i < gsystem.N; i++){
			itot = traj*gsystem.N + i;
			fixForceData.h_extForces[itot].x *= fixForceData.h_extForces[itot].w;
			fixForceData.h_extForces[itot].y *= fixForceData.h_extForces[itot].w;
			fixForceData.h_extForces[itot].z *= fixForceData.h_extForces[itot].w;
			fixForceData.h_extForces[itot].w = sqrtf(
					fixForceData.h_extForces[itot].x*fixForceData.h_extForces[itot].x +
					fixForceData.h_extForces[itot].y*fixForceData.h_extForces[itot].y +
					fixForceData.h_extForces[itot].z*fixForceData.h_extForces[itot].z
					);
			if(fixForceData.h_atomMasks[itot] == FIXFORCE_ATOM_PULLED){
				DPRINTF("Atom #%d (%s%d-%s) in trajectory #%d will be pulled with "
					"constant external force of %5.3f pN (%5.3f, %5.3f, %5.3f).\n",
					i, topology.atoms[i].resName, topology.atoms[i].resid, topology.atoms[i].name, traj+parameters.firstrun,
					fixForceData.h_extForces[itot].w,
					fixForceData.h_extForces[itot].x,
					fixForceData.h_extForces[itot].y,
					fixForceData.h_extForces[itot].z);
			}
		}
	}
	cudaMemcpy(fixForceData.d_atomMasks, fixForceData.h_atomMasks,
				gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(fixForceData.d_extForces, fixForceData.h_extForces,
				gsystem.Ntot*sizeof(float4), cudaMemcpyHostToDevice);
}

void initSMD(){

	allocateCPU((void**)&fixForceOutputFilename, parameters.Ntr*sizeof(char*));
	int traj, i, itot;
	char trajnum[10];
	for(traj = 0; traj < parameters.Ntr; traj++){
		sprintf(trajnum, "%d", traj + parameters.firstrun);
		fixForceOutputFilename[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
		getMaskedParameterWithReplacement(fixForceOutputFilename[traj], PARAMETER_PULLING_SMD_FILE,
				trajnum, "<run>");
		fixForceOutputFile = safe_fopen(fixForceOutputFilename[traj], "w");
		fclose(fixForceOutputFile);

	}

	allocateCPU((void**)&smdData.h_initialTipPositions, gsystem.Ntot*sizeof(float4));
	allocateGPU((void**)&smdData.d_initialTipPositions, gsystem.Ntot*sizeof(float4));
	allocateCPU((void**)&smdData.h_tipPositions, gsystem.Ntot*sizeof(float4));
	allocateGPU((void**)&smdData.d_tipPositions, gsystem.Ntot*sizeof(float4));
	allocateCPU((void**)&smdData.h_ks, gsystem.Ntot*sizeof(float));
	allocateGPU((void**)&smdData.d_ks, gsystem.Ntot*sizeof(float));
	allocateCPU((void**)&smdData.h_pullVelocities, gsystem.Ntot*sizeof(float4));
	allocateGPU((void**)&smdData.d_pullVelocities, gsystem.Ntot*sizeof(float4));
	float mult = 1.0e-9f*integrator->h*smdForceUpdater.frequency;
	char initialTipPDBFilename[256];
	PDB initialTipPDB;
	int restart = getYesNoParameter(PARAMETER_PULLING_SMD_RESTART, 0);
	float4 initialTipPosition;
	getVectorParameter(PARAMETER_PULLING_SMD_INITIAL_TIP_VECTOR,
			&initialTipPosition.x, &initialTipPosition.y, &initialTipPosition.z,
			-666.0f, -666.0f, -666.0f);
	if(initialTipPosition.x == -666.0f && initialTipPosition.y == -666.0f && initialTipPosition.z == -666.0f){
		initialTipPosition.w = 0.0f;
	} else {
		initialTipPosition.w = 1.0f;
	}
	float rescaleVelocity = getFloatParameter(PARAMETER_PULLING_SMD_RVEL, 1.0f);
	float rescaleKS       = getFloatParameter(PARAMETER_PULLING_SMD_RKS , 1.0f);
	for(traj = 0; traj < parameters.Ntr; traj++){
		sprintf(trajnum, "%d", traj + parameters.firstrun);
		if(restart){
			getMaskedParameterWithReplacement(initialTipPDBFilename, PARAMETER_PULLING_SMD_INITIAL_TIP_FILE,
						trajnum, "<run>");
			readPDB(initialTipPDBFilename, &initialTipPDB);
		}

		for(i = 0; i < gsystem.N; i++){
			itot = traj*gsystem.N + i;
			smdData.h_pullVelocities[itot].x = fixForceData.h_extForces[itot].x * rescaleVelocity;
			smdData.h_pullVelocities[itot].y = fixForceData.h_extForces[itot].y * rescaleVelocity;
			smdData.h_pullVelocities[itot].z = fixForceData.h_extForces[itot].z * rescaleVelocity;
			smdData.h_ks[itot] = fixForceData.h_extForces[itot].w * rescaleKS;
			fixForceData.h_extForces[itot].w = sqrtf(
					fixForceData.h_extForces[itot].x*fixForceData.h_extForces[itot].x +
					fixForceData.h_extForces[itot].y*fixForceData.h_extForces[itot].y +
					fixForceData.h_extForces[itot].z*fixForceData.h_extForces[itot].z
					) * rescaleVelocity;

			if(fixForceData.h_atomMasks[itot] == FIXFORCE_ATOM_PULLED){
				DPRINTF("Atom #%d (%s%d-%s) in trajectory #%d will be pulled with "
					"constant velocity of %5.3f um/s (%5.3f, %5.3f, %5.3f) and Ks=%5.3f pN/nm\n",
					i, topology.atoms[i].resName, topology.atoms[i].resid, topology.atoms[i].name, traj+parameters.firstrun,
					fixForceData.h_extForces[itot].w,
					smdData.h_pullVelocities[itot].x,
					smdData.h_pullVelocities[itot].y,
					smdData.h_pullVelocities[itot].z,
					smdData.h_ks[itot]);
				if(restart == 0){
					if(initialTipPosition.w == 0.0f){
						smdData.h_initialTipPositions[itot].x = gsystem.h_coord[itot].x;
						smdData.h_initialTipPositions[itot].y = gsystem.h_coord[itot].y;
						smdData.h_initialTipPositions[itot].z = gsystem.h_coord[itot].z;
					} else {
						smdData.h_initialTipPositions[itot].x = initialTipPosition.x;
						smdData.h_initialTipPositions[itot].y = initialTipPosition.y;
						smdData.h_initialTipPositions[itot].z = initialTipPosition.z;
					}
				} else {
					smdData.h_initialTipPositions[itot].x = initialTipPDB.atoms[i].x/10.0f;
					smdData.h_initialTipPositions[itot].y = initialTipPDB.atoms[i].y/10.0f;
					smdData.h_initialTipPositions[itot].z = initialTipPDB.atoms[i].z/10.0f;
				}
				smdData.h_initialTipPositions[itot].w = 0.0f;
			} else {
				smdData.h_initialTipPositions[itot].x = 0.0f;
				smdData.h_initialTipPositions[itot].y = 0.0f;
				smdData.h_initialTipPositions[itot].z = 0.0f;
				smdData.h_initialTipPositions[itot].w = 0.0f;
			}

			smdData.h_pullVelocities[itot].x *= mult;
			smdData.h_pullVelocities[itot].y *= mult;
			smdData.h_pullVelocities[itot].z *= mult;
			float deltax = sqrtf(
					smdData.h_pullVelocities[itot].x*smdData.h_pullVelocities[itot].x +
					smdData.h_pullVelocities[itot].y*smdData.h_pullVelocities[itot].y +
					smdData.h_pullVelocities[itot].z*smdData.h_pullVelocities[itot].z
					);

			if(fixForceData.h_atomMasks[itot] == FIXFORCE_ATOM_PULLED){
				DPRINTF("   (i.e. it will be shifted by %f nm every %d steps)", deltax, smdForceUpdater.frequency);
				DPRINTF("   Initial tip position for this atom: %5.3f %5.3f %5.3f\n",
						smdData.h_initialTipPositions[itot].x, 
						smdData.h_initialTipPositions[itot].y, 
						smdData.h_initialTipPositions[itot].z, 
				       );
			}
			if(fixForceData.h_atomMasks[itot] == FIXFORCE_ATOM_PULLED){
				smdData.h_tipPositions[itot].x = smdData.h_initialTipPositions[itot].x;
				smdData.h_tipPositions[itot].y = smdData.h_initialTipPositions[itot].y;
				smdData.h_tipPositions[itot].z = smdData.h_initialTipPositions[itot].z;
			} else {
				smdData.h_tipPositions[itot].x = 0.0f;
				smdData.h_tipPositions[itot].y = 0.0f;
				smdData.h_tipPositions[itot].z = 0.0f;
			}
			fixForceData.h_extForces[itot].x = 0.0f;
			fixForceData.h_extForces[itot].y = 0.0f;
			fixForceData.h_extForces[itot].z = 0.0f;
			fixForceData.h_extForces[itot].w = 0.0f;
		}
	}

	cudaMemcpy(fixForceData.d_atomMasks, fixForceData.h_atomMasks,
					gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(fixForceData.d_extForces, fixForceData.h_extForces,
				gsystem.Ntot*sizeof(float4), cudaMemcpyHostToDevice);
	cudaMemcpy(smdData.d_initialTipPositions, smdData.h_initialTipPositions,
					gsystem.Ntot*sizeof(float4), cudaMemcpyHostToDevice);
	cudaMemcpy(smdData.d_tipPositions, smdData.h_tipPositions,
					gsystem.Ntot*sizeof(float4), cudaMemcpyHostToDevice);
	cudaMemcpy(smdData.d_ks, smdData.h_ks,
						gsystem.Ntot*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(smdData.d_pullVelocities, smdData.h_pullVelocities,
					gsystem.Ntot*sizeof(float4), cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_smdData, &smdData,
					sizeof(SMDData), 0, cudaMemcpyHostToDevice);
}

__global__ void computeFixForcePotentialFConst_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		int mask = c_fixForceData.d_atomMasks[d_i];
		if(mask == FIXFORCE_ATOM_FIXED){
			float4 vel = c_gsystem.d_vel[d_i];
			float4 f = c_gsystem.d_forces[d_i];
			vel.x = 0.0f;
			vel.y = 0.0f;
			vel.z = 0.0f;
			c_gsystem.d_vel[d_i] = vel;
			f.x = 0.0f;
			f.y = 0.0f;
			f.z = 0.0f;
			c_gsystem.d_forces[d_i] = f;
		} else
		if(mask == FIXFORCE_ATOM_PULLED){
			float4 f = c_gsystem.d_forces[d_i];
			float4 df = c_fixForceData.d_extForces[d_i];
			f.x += df.x;
			f.y += df.y;
			f.z += df.z;
			c_gsystem.d_forces[d_i] = f;
		}
	}
}

inline void computeFConst(){
	computeFixForcePotentialFConst_kernel<<<fixForceBlockCount, fixForceBlockSize>>>();
}

__global__ void computeFixForcePotentialSMD_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		int mask = c_fixForceData.d_atomMasks[d_i];
		if(mask == FIXFORCE_ATOM_FIXED){
			float4 vel = c_gsystem.d_vel[d_i];
			float4 f = c_gsystem.d_forces[d_i];
			vel.x = 0.0f;
			vel.y = 0.0f;
			vel.z = 0.0f;
			c_gsystem.d_vel[d_i] = vel;
			f.x = 0.0f;
			f.y = 0.0f;
			f.z = 0.0f;
			c_gsystem.d_forces[d_i] = f;
		} else
		if(mask == FIXFORCE_ATOM_PULLED){
			float4 f = c_gsystem.d_forces[d_i];
			float4 coord = c_gsystem.d_coord[d_i];
			float4 tip = c_smdData.d_tipPositions[d_i];
			float ks = c_smdData.d_ks[d_i];
			float4 extF = c_fixForceData.d_extForces[d_i];
			float4 df;
			df.x = tip.x - coord.x;
			df.y = tip.y - coord.y;
			df.z = tip.z - coord.z;
			DO_PBC(df);
			df.x *= ks;
			df.y *= ks;
			df.z *= ks;
			extF.x += df.x;
			extF.y += df.y;
			extF.z += df.z;
			c_fixForceData.d_extForces[d_i] = extF;
			f.x += df.x;
			f.y += df.y;
			f.z += df.z;
			c_gsystem.d_forces[d_i] = f;
		}
	}
}

inline void computeSMD(){
	computeFixForcePotentialSMD_kernel<<<fixForceBlockCount, fixForceBlockSize>>>();
}



void destroy(){

}

__global__ void updatePullingForceSMD_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		if(c_fixForceData.d_atomMasks[d_i] == FIXFORCE_ATOM_PULLED){
			float4 r0 = c_smdData.d_initialTipPositions[d_i];
			float4 vf = c_smdData.d_pullVelocities[d_i];
			vf.x *= c_fixForceTimesMoved;
			vf.y *= c_fixForceTimesMoved;
			vf.z *= c_fixForceTimesMoved;
			r0.x += vf.x;
			r0.y += vf.y;
			r0.z += vf.z;
			DO_PBC(r0);
			c_smdData.d_tipPositions[d_i] = r0;
			c_fixForceData.d_extForces[d_i] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
		}
	}
}

void updatePullingForceSMD(){
	fixForceTimesMoved = step/smdForceUpdater.frequency;
	//printf("%d\n", fixForceTimesMoved);
	cudaMemcpyToSymbol(c_fixForceTimesMoved, &fixForceTimesMoved,
						sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpy(fixForceData.h_extForces, fixForceData.d_extForces,
					gsystem.Ntot*sizeof(float4), cudaMemcpyDeviceToHost);
	cudaMemcpy(smdData.h_tipPositions, smdData.d_tipPositions,
						gsystem.Ntot*sizeof(float4), cudaMemcpyDeviceToHost);
	int traj, i, itot;
	printf("%*s%-*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s\n",
			PULLING_OUTPUT_WIDTH, "",
			PULLING_OUTPUT_WIDTH, PULLING_OUTPUT_NAME_TRAJECTORY,
			PULLING_OUTPUT_WIDTH, PULLING_OUTPUT_NAME_STEP,
			PULLING_OUTPUT_WIDTH, PULLING_OUTPUT_NAME_ATOM,
			PULLING_OUTPUT_WIDTH, PULLING_OUTPUT_NAME_FEXT,
			PULLING_OUTPUT_WIDTH, PULLING_OUTPUT_NAME_FX,
			PULLING_OUTPUT_WIDTH, PULLING_OUTPUT_NAME_FY,
			PULLING_OUTPUT_WIDTH, PULLING_OUTPUT_NAME_FZ,
			PULLING_OUTPUT_WIDTH, PULLING_OUTPUT_NAME_R,
			PULLING_OUTPUT_WIDTH, PULLING_OUTPUT_NAME_RX,
			PULLING_OUTPUT_WIDTH, PULLING_OUTPUT_NAME_RY,
			PULLING_OUTPUT_WIDTH, PULLING_OUTPUT_NAME_RZ);

	for(traj = 0; traj < parameters.Ntr; traj++){
		for(i = 0; i < gsystem.N; i++){
			itot = traj*gsystem.N + i;
			if(fixForceData.h_atomMasks[itot] == FIXFORCE_ATOM_PULLED){
				smdData.h_tipPositions[itot].x -= smdData.h_initialTipPositions[itot].x;
				smdData.h_tipPositions[itot].y -= smdData.h_initialTipPositions[itot].y;
				smdData.h_tipPositions[itot].z -= smdData.h_initialTipPositions[itot].z;
				smdData.h_tipPositions[itot].w = sqrtf(
						smdData.h_tipPositions[itot].x*smdData.h_tipPositions[itot].x +
						smdData.h_tipPositions[itot].y*smdData.h_tipPositions[itot].y +
						smdData.h_tipPositions[itot].z*smdData.h_tipPositions[itot].z
						);
				fixForceData.h_extForces[itot].x /= ((float)smdForceUpdater.frequency);
				fixForceData.h_extForces[itot].y /= ((float)smdForceUpdater.frequency);
				fixForceData.h_extForces[itot].z /= ((float)smdForceUpdater.frequency);
				fixForceData.h_extForces[itot].w = sqrtf(
						fixForceData.h_extForces[itot].x*fixForceData.h_extForces[itot].x +
						fixForceData.h_extForces[itot].y*fixForceData.h_extForces[itot].y +
						fixForceData.h_extForces[itot].z*fixForceData.h_extForces[itot].z
						);
				printf("%-*s#%-*d%*lld%*d%*f%*f%*f%*f%*f%*f%*f%*f\n",
						PULLING_OUTPUT_WIDTH, PULLING_OUTPUT_TITLE,
						PULLING_OUTPUT_WIDTH-1, traj+parameters.firstrun,
						PULLING_OUTPUT_WIDTH, step,
						PULLING_OUTPUT_WIDTH, topology.atoms[i].id,
						PULLING_OUTPUT_WIDTH, fixForceData.h_extForces[itot].w,
						PULLING_OUTPUT_WIDTH, fixForceData.h_extForces[itot].x,
						PULLING_OUTPUT_WIDTH, fixForceData.h_extForces[itot].y,
						PULLING_OUTPUT_WIDTH, fixForceData.h_extForces[itot].z,
						PULLING_OUTPUT_WIDTH, smdData.h_tipPositions[itot].w,
						PULLING_OUTPUT_WIDTH, smdData.h_tipPositions[itot].x,
						PULLING_OUTPUT_WIDTH, smdData.h_tipPositions[itot].y,
						PULLING_OUTPUT_WIDTH, smdData.h_tipPositions[itot].z);
				fixForceOutputFile = safe_fopen(fixForceOutputFilename[traj], "a");
				fprintf(fixForceOutputFile, "Atom#%d\t%lld\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%f\t%f\t%f\t%f\n",
										topology.atoms[i].id, step,
										fixForceData.h_extForces[itot].w,
										fixForceData.h_extForces[itot].x,
										fixForceData.h_extForces[itot].y,
										fixForceData.h_extForces[itot].z,
										smdData.h_tipPositions[itot].w,
										smdData.h_tipPositions[itot].x,
										smdData.h_tipPositions[itot].y,
										smdData.h_tipPositions[itot].z);
				fclose(fixForceOutputFile);
			}
		}
	}
	updatePullingForceSMD_kernel<<<fixForceBlockCount, fixForceBlockSize>>>();

}
void destroySMDForceUpdater(){

}

#undef LOG

} // namespace fixforce_potential

