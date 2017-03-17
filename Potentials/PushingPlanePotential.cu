/*
 * PushingPlanePotential.cu
 *
 *  Created on: Apr 16, 2016
 *      Author: kir_min
 */
#include "../Core/global.h"
#include "../Util/Log.h"
#include "PushingPlanePotential.cuh"

namespace pushing_plane_potential
{

class Log: public ILog
{
	virtual void Write(const char* message) const
	{
		std::cout << makeTimePrefix() << "<pushing_plane_potential> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void create()
{
	if(getYesNoParameter(PARAMETER_PUSHING_PLANE, DEFAULT_PUSHING_PLANE))
	{
		potential.compute = &computePlanePushing;
		potential.destroy = &destroy;
		sprintf(potential.name, "Pushing-Plane potential");
		potentials[potentialsCount] = &potential;
		potentialsCount++;
		sprintf(planeLocationUpdater.name, "Pushing-Plane updater");
		planeLocationUpdater.update = updatePlaneLocation;
		planeLocationUpdater.destroy = destroyPlaneLocationUpdater;
		planeLocationUpdater.frequency = getIntegerParameter(PARAMETER_PUSHING_PLANE_UPDATE_FREQ);
		updaters[updatersCount] = &planeLocationUpdater;
		updatersCount ++;
		init();
		LOG << "Done initializing Pushing-Plane potential.";
	}
	else
	{
		LOG << "No Plane pushing will be performed.";
	}
}
void init()
{
	LOG << "Initializing Pulling-Plane potential...";
	
	blockSize = BLOCK_SIZE;
	blockCount = (gsystem.Ntot-1)/BLOCK_SIZE + 1;

	potentialData.planeCount = getIntegerParameter(PARAMETER_PUSHING_PLANE_COUNT, 2, 1);
	potentialData.sigma = getFloatParameter(PARAMETER_PUSHING_PLANE_SIGMA, 1.0, 1);
	potentialData.epsilon = getFloatParameter(PARAMETER_PUSHING_PLANE_EPSILON, 1.0, 1);

	allocateCPU((void**)&potentialData.h_planeNormal, potentialData.planeCount*sizeof(float4));
	allocateGPU((void**)&potentialData.d_planeNormal, potentialData.planeCount*sizeof(float4));

	allocateCPU((void**)&potentialData.h_planePosition, parameters.Ntr*potentialData.planeCount*sizeof(float4));
	allocateGPU((void**)&potentialData.d_planePosition, parameters.Ntr*potentialData.planeCount*sizeof(float4));

	allocateCPU((void**)&potentialData.h_planePosition0, parameters.Ntr*potentialData.planeCount*sizeof(float4));
	allocateGPU((void**)&potentialData.d_planePosition0, parameters.Ntr*potentialData.planeCount*sizeof(float4));

	allocateCPU((void**)&potentialData.r, parameters.Ntr*sizeof(float));

	allocateCPU((void**)&potentialData.h_forces, potentialData.planeCount*gsystem.Ntot*sizeof(float4));
	allocateGPU((void**)&potentialData.d_forces, potentialData.planeCount*gsystem.Ntot*sizeof(float4));

	int i;
	for(i = 0; i < potentialData.planeCount*gsystem.Ntot; i++)
	{
		potentialData.h_forces[i].x = 0.0;
		potentialData.h_forces[i].y = 0.0;
		potentialData.h_forces[i].z = 0.0;
	}

	getVectorParameter(PARAMETER_PUSHING_PLANE_NORM_FIX, &potentialData.h_planeNormal[FIX_PLANE].x, &potentialData.h_planeNormal[FIX_PLANE].y, &potentialData.h_planeNormal[FIX_PLANE].z);
	getVectorParameter(PARAMETER_PUSHING_PLANE_POSITION_FIX, &potentialData.h_planePosition0[FIX_PLANE].x, &potentialData.h_planePosition0[FIX_PLANE].y, &potentialData.h_planePosition0[FIX_PLANE].z);
	getVectorParameter(PARAMETER_PUSHING_PLANE_NORM_PUSH, &potentialData.h_planeNormal[PUSHING_PLANE].x, &potentialData.h_planeNormal[PUSHING_PLANE].y, &potentialData.h_planeNormal[PUSHING_PLANE].z);
	getVectorParameter(PARAMETER_PUSHING_PLANE_POSITION_PUSH, &potentialData.h_planePosition0[PUSHING_PLANE].x, &potentialData.h_planePosition0[PUSHING_PLANE].y, &potentialData.h_planePosition0[PUSHING_PLANE].z);

	potentialData.pushSpeed = getFloatParameter(PARAMETER_PUSHING_PLANE_PUSHSPEED);
	potentialData.springConstant = getFloatParameter(PARAMETER_PUSHING_PLANE_PUSHSPRING, 100.0, 1);
	potentialData.dampingConstant = getFloatParameter(PARAMETER_PUSHING_PLANE_PUSHDAMPING, 100.0, 1);
	potentialData.updateFreq = getLongIntegerParameter(PARAMETER_PUSHING_PLANE_UPDATE_FREQ);
	potentialData.outputFreq = getLongIntegerParameter(PARAMETER_PUSHING_PLANE_OUTPUT_FREQ);
	
	int p;
	for(p = 0; p < potentialData.planeCount; p++){
		float4 h_n = potentialData.h_planeNormal[p];
		float h_n2 = h_n.x*h_n.x + h_n.y*h_n.y + h_n.z*h_n.z;
		float h_n_mod = sqrt(h_n2);
		h_n.x /= h_n_mod;
		h_n.y /= h_n_mod;
		h_n.z /= h_n_mod;
		potentialData.h_planeNormal[p] = h_n;
	}
	
	//potentialData.x0 = 0.0;

	/*getMaskedParameter(outputfilename, PARAMETER_PUSHING_PLANE_OUTPUT_FILE);
	FILE* datout = fopen(outputfilename, "w");
	fclose(datout);*/

	allocateCPU((void**)&outputfilename, parameters.Ntr*sizeof(char*));
	char trajnum[10];
	int traj;
	for(traj = 0; traj < parameters.Ntr; traj++)
	{
		potentialData.r[traj] = 0.0f;

		potentialData.h_planePosition[FIX_PLANE + traj*potentialData.planeCount].x = potentialData.h_planePosition0[FIX_PLANE].x;
		potentialData.h_planePosition[FIX_PLANE + traj*potentialData.planeCount].y = potentialData.h_planePosition0[FIX_PLANE].y;
		potentialData.h_planePosition[FIX_PLANE + traj*potentialData.planeCount].z = potentialData.h_planePosition0[FIX_PLANE].z;

		potentialData.h_planePosition[PUSHING_PLANE + traj*potentialData.planeCount].x = potentialData.h_planePosition0[PUSHING_PLANE].x;
		potentialData.h_planePosition[PUSHING_PLANE + traj*potentialData.planeCount].y = potentialData.h_planePosition0[PUSHING_PLANE].y;
		potentialData.h_planePosition[PUSHING_PLANE + traj*potentialData.planeCount].z = potentialData.h_planePosition0[PUSHING_PLANE].z;

		if (traj != 0)
		{
			potentialData.h_planePosition0[FIX_PLANE + traj*potentialData.planeCount].x = potentialData.h_planePosition0[FIX_PLANE].x;
			potentialData.h_planePosition0[FIX_PLANE + traj*potentialData.planeCount].y = potentialData.h_planePosition0[FIX_PLANE].y;
			potentialData.h_planePosition0[FIX_PLANE + traj*potentialData.planeCount].z = potentialData.h_planePosition0[FIX_PLANE].z;

			potentialData.h_planePosition0[PUSHING_PLANE + traj*potentialData.planeCount].x = potentialData.h_planePosition0[PUSHING_PLANE].x;
			potentialData.h_planePosition0[PUSHING_PLANE + traj*potentialData.planeCount].y = potentialData.h_planePosition0[PUSHING_PLANE].y;
			potentialData.h_planePosition0[PUSHING_PLANE + traj*potentialData.planeCount].z = potentialData.h_planePosition0[PUSHING_PLANE].z;
		}
		
		sprintf(trajnum, "%d", traj + parameters.firstrun);
		outputfilename[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
		getMaskedParameterWithReplacement(outputfilename[traj], PARAMETER_PULLING_PLANE_OUTPUT_FILE, trajnum, "<run>");
		FILE* datout = fopen(outputfilename[traj], "w");
		fclose(datout);
	}

	getMaskedParameter(reffilename, PARAMETER_PUSHING_PLANE_REF_FILE);
	readPDB(reffilename, &refpdb);

	cudaMemcpy(potentialData.d_planeNormal, potentialData.h_planeNormal, potentialData.planeCount*sizeof(float4), cudaMemcpyHostToDevice);
	cudaMemcpy(potentialData.d_planePosition, potentialData.h_planePosition, parameters.Ntr*potentialData.planeCount*sizeof(float4), cudaMemcpyHostToDevice);
	cudaMemcpy(potentialData.d_forces, potentialData.h_forces, potentialData.planeCount*gsystem.Ntot*sizeof(float4), cudaMemcpyHostToDevice);

	
	cudaMemcpyToSymbol(c_potentialData, &potentialData, sizeof(PotentialData), 0, cudaMemcpyHostToDevice);
}
__global__ void computePushingPlane_kernel()
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < c_gsystem.Ntot){
		int traj = i/c_gsystem.N;
		float4 dr = c_gsystem.d_coord[i];
		float4 f = c_gsystem.d_forces[i];
		int p;
		for(p = 0; p < c_potentialData.planeCount; p++){
			float4 fpl = c_potentialData.d_forces[i + p*c_gsystem.Ntot];
			float4 dr0 = c_potentialData.d_planePosition[p + c_potentialData.planeCount*traj];
			float4 dn = c_potentialData.d_planeNormal[p];
			float mod_dr = (dr.x-dr0.x)*dn.x + (dr.y-dr0.y)*dn.y + (dr.z-dr0.z)*dn.z;
			float sor = c_potentialData.sigma/mod_dr;
			float sor2 = sor*sor;
			float sor6 = sor2*sor2*sor2;
			float mul = 6.0*c_potentialData.epsilon*sor6/mod_dr;
			float4 df;
			df.x = mul*dn.x;
			df.y = mul*dn.y;
			df.z = mul*dn.z;
			fpl.x += df.x;
			fpl.y += df.y;
			fpl.z += df.z;
			f.x += df.x;
			f.y += df.y;
			f.z += df.z;
			c_potentialData.d_forces[i + p*c_gsystem.Ntot] = fpl;
		}
		c_gsystem.d_forces[i] = f;
	}
}
inline void computePlanePushing()
{	
	computePushingPlane_kernel<<<blockSize, blockCount>>>();
}
void updatePlaneLocation()
{
	if(step % potentialData.updateFreq == 0 || step == 0){
		cudaMemcpy(potentialData.h_forces, potentialData.d_forces, potentialData.planeCount*gsystem.Ntot*sizeof(float4), cudaMemcpyDeviceToHost);
		int traj;
		for(traj = 0; traj < parameters.Ntr; traj++)
		{
			float4 n = potentialData.h_planeNormal[PUSHING_PLANE];
			float k = 0.6*potentialData.springConstant;
			float v = potentialData.pushSpeed/1000000000.0;
			float t = step*integrator->h;
			float pushforce = k*(v*t-potentialData.r[traj]);
			float atomforce = 0.0;
			int i;
			for(i = 0; i < gsystem.N; i++){
				int itot = i + traj*gsystem.N;
				float4 f = potentialData.h_forces[itot + PUSHING_PLANE*gsystem.Ntot];
				//atomforce += f.x*n.x + f.y*n.y + f.z*n.z;
				float f2 = f.x*f.x + f.y*f.y + f.z*f.z;
				float mod_f = sqrt(f2);
				atomforce += mod_f;
				potentialData.h_forces[itot + PUSHING_PLANE*gsystem.Ntot].x = 0.0;
				potentialData.h_forces[itot + PUSHING_PLANE*gsystem.Ntot].y = 0.0;
				potentialData.h_forces[itot + PUSHING_PLANE*gsystem.Ntot].z = 0.0;
			}
			atomforce /= potentialData.updateFreq;
			float force = pushforce - atomforce;
			float tau = integrator->h;
			potentialData.r[traj] += tau*potentialData.updateFreq*force/potentialData.dampingConstant;
			potentialData.h_planePosition[PUSHING_PLANE + traj*potentialData.planeCount].x = potentialData.h_planePosition0[PUSHING_PLANE + traj*potentialData.planeCount].x + potentialData.r[traj]*n.x;
			potentialData.h_planePosition[PUSHING_PLANE + traj*potentialData.planeCount].y = potentialData.h_planePosition0[PUSHING_PLANE + traj*potentialData.planeCount].y + potentialData.r[traj]*n.y;
			potentialData.h_planePosition[PUSHING_PLANE + traj*potentialData.planeCount].z = potentialData.h_planePosition0[PUSHING_PLANE + traj*potentialData.planeCount].z + potentialData.r[traj]*n.z;
		
			//End-to-end
			copyCoordinatesFromGPU(false);
			double TrueExtension;
			double x_fixed = 0.0;
			double y_fixed = 0.0;
			double z_fixed = 0.0;
			double x_pulled = 0.0;
			double y_pulled = 0.0;
			double z_pulled = 0.0;
			int count_fixed = 0;
			int count_pulled = 0;
			int itot;
			for(i = 0; i < gsystem.N; i++){
				itot = i + traj*gsystem.N;
				if(refpdb.atoms[i].occupancy == 1.0){
					x_fixed	+= gsystem.h_coord[itot].x;
					y_fixed	+= gsystem.h_coord[itot].y;
					z_fixed	+= gsystem.h_coord[itot].z;
					count_fixed ++;
				}
				if(refpdb.atoms[i].beta == 1.0){
					x_pulled += gsystem.h_coord[itot].x;
					y_pulled += gsystem.h_coord[itot].y;
					z_pulled += gsystem.h_coord[itot].z;
					count_pulled ++;
				}
			}
			x_fixed /= count_fixed;
			y_fixed /= count_fixed;
			z_fixed /= count_fixed;
			x_pulled /= count_pulled;
			y_pulled /= count_pulled;
			z_pulled /= count_pulled;
			double deltaX = x_pulled - x_fixed;
			double deltaY = y_pulled - y_fixed;
			double deltaZ = z_pulled - z_fixed;
			double X = deltaX * deltaX;
			double Y = deltaY * deltaY;
			double Z = deltaZ * deltaZ;
			TrueExtension = sqrt(X + Y + Z);

			if(step % potentialData.outputFreq == 0 || step == 0){
				FILE* datout = fopen(outputfilename[traj], "a");
				fprintf(datout, "%lld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", step,
					potentialData.h_planePosition[PUSHING_PLANE + traj*potentialData.planeCount].x,
					potentialData.h_planePosition[PUSHING_PLANE + traj*potentialData.planeCount].y,
					potentialData.h_planePosition[PUSHING_PLANE + traj*potentialData.planeCount].z,
					pushforce, atomforce, force, TrueExtension);
				fclose(datout);
			}
			cudaMemcpy(potentialData.d_planePosition, potentialData.h_planePosition, parameters.Ntr*potentialData.planeCount*sizeof(float4), cudaMemcpyHostToDevice);
			cudaMemcpy(potentialData.d_forces, potentialData.h_forces, potentialData.planeCount*gsystem.Ntot*sizeof(float4), cudaMemcpyHostToDevice);
		}//close for
	}
			
}
void destroyPlaneLocationUpdater()
{

}

void destroy()
{

}

}
