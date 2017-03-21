/*
 * PullingPlanePotential.cu
 *
 *  Created on: Sep 19, 2013
 *      Author: zip
 */
#include "../Core/global.h"
#include "../Util/Log.h"
#include "PullingPlanePotential.cuh"

namespace pulling_plane_potential
{

class Log: public ILog
{
	virtual void Write(const char* message) const
	{
		std::cout << makeTimePrefix() << "<pulling_plane_potential> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void create()
{
	if(getYesNoParameter(PARAMETER_PULLING_PLANE, DEFAULT_PULLING_PLANE))
	{
		potential.destroy = &destroy;
		sprintf(potential.name, "Pulling-Plane potential");
		potentials[potentialsCount] = &potential;
		potentialsCount++;
		init();
		char pullProtocol[PARAMETER_LENGTH];
		getMaskedParameter(pullProtocol, PARAMETER_PULLING_PLANE_PROTOCOL);
		if(strcmp(pullProtocol, PARAMETER_VALUE_PULLING_PLANE_PROTOCOL_FCONST) == 0)
		{
			potential.compute = &computePlanePulling;
			potentialData.pullForce = getFloatParameter(PARAMETER_PULLING_PLANE_PULLFORCE, 0, 0);
			LOG << "Constant Force Plane pulling will be performed. Pulling force value is set to " << potentialData.pullForce << "pN.";
			//convert pN to kJ/(mol*nm)
			potentialData.pullForce = 0.6*potentialData.pullForce;
			potentialData.pullSpring = 0;
			potentialData.pullSpeed = 0;
			sprintf(planeLocationUpdater.name, "Pulling-Plane Constant Force updater");
			planeLocationUpdater.update = updatePlaneLocation;
			planeLocationUpdater.destroy = destroyPlaneLocationUpdater;
			planeLocationUpdater.frequency = getIntegerParameter(PARAMETER_PULLING_PLANE_UPDATE_FREQ);
			updaters[updatersCount] = &planeLocationUpdater;
			updatersCount ++;
		}
		else
			if(strcmp(pullProtocol, PARAMETER_VALUE_PULLING_PLANE_PROTOCOL_FRAMP) == 0)
			{
				potential.compute = &computePlanePulling;
				potentialData.pullForce = 0;
				potentialData.pullSpring = getFloatParameter(PARAMETER_PULLING_PLANE_PULLSPRING, 0, 0);
				potentialData.pullSpeed = getFloatParameter(PARAMETER_PULLING_PLANE_PULLSPEED, 0, 0);
				LOG << "Constant Velocity Plane pulling will be performed. Pulling speed is set to " << potentialData.pullSpeed << " um/s, spring constant is set to " << potentialData.pullSpring << " pN/nm.";
				//convert um/s to nm/ps
				potentialData.pullSpeed = potentialData.pullSpeed / 1000000000.0;
				//convert pN/nm to kJ/(mol*nm^2)
				potentialData.pullSpring = 0.6*potentialData.pullSpring;
				sprintf(planeLocationUpdater.name, "Pulling-Plane Constant Velocity updater");
				planeLocationUpdater.update = updatePlaneLocation;
				planeLocationUpdater.destroy = destroyPlaneLocationUpdater;
				planeLocationUpdater.frequency = getIntegerParameter(PARAMETER_PULLING_PLANE_UPDATE_FREQ);
				updaters[updatersCount] = &planeLocationUpdater;
				updatersCount ++;
			}
			else
			{
				DIE("Pulling-Plane protocol parameter should be either %s for constant force pulling or %s for force-ramp (constant speed) pulling.",
						PARAMETER_VALUE_PULLING_PLANE_PROTOCOL_FCONST,
						PARAMETER_VALUE_PULLING_PLANE_PROTOCOL_FRAMP);
			}
		cudaMemcpy(potentialData.d_workList, potentialData.h_workList, potentialData.Ntotal*sizeof(WorkList), cudaMemcpyHostToDevice);
		cudaMemcpy(potentialData.d_force, potentialData.h_force, potentialData.Ntotal*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(potentialData.d_planeDisplacement, potentialData.h_planeDisplacement, potentialData.Nplane*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(c_potentialData, &potentialData, sizeof(PotentialData), 0, cudaMemcpyHostToDevice);
		//init output
		allocateCPU((void**)&outputFilename, parameters.Ntr*sizeof(char*));
		char trajnum[10];
		int traj;
		for(traj = 0; traj < parameters.Ntr; traj++)
		{
			sprintf(trajnum, "%d", traj + parameters.firstrun);
			outputFilename[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
			getMaskedParameterWithReplacement(outputFilename[traj], PARAMETER_PULLING_PLANE_OUTPUT_FILE, trajnum, "<run>");
			outputFile = safe_fopen(outputFilename[traj], "w");
			fclose(outputFile);
		}
		outputFreq = getIntegerParameter(PARAMETER_PULLING_PLANE_OUTPUT_FREQ);
		logOutputFreq = getIntegerParameter(PARAMETER_ENERGYOUTPUT_FREQ);
		planePullingGridSize = potentialData.Ntotal/BLOCK_SIZE + 1;
		planePullingBlockSize = BLOCK_SIZE;
		LOG << "Done initializing Pulling-Plane potential.";
	}
	else
	{
		LOG << "No Plane pulling will be performed.";
	}
}

void init()
{
	LOG << "Initializing Pulling-Plane potential...";
	//1tr
	//allocateCPU((void**)&potentialData.prc, gsystem.N*sizeof(int));
	//ntr
	allocateCPU((void**)&potentialData.prc, gsystem.Ntot*sizeof(int));
	potentialData.Npulled = 0;
	potentialData.Nfixed = 0;
	int traj, i, itot;
	char refFilename[100];
	char trajnum[10];
	PDB refFile;
	//1tr
	//traj = 0;
	//if(parameters.Ntr != 1)
	//	DIE("Pulling-Plane protocol can not yet perform pulling for more than 1 trajectory.");
	//ntr
	for(traj = 0; traj < parameters.Ntr; traj ++)
	{
		sprintf(trajnum, "%d", traj + parameters.firstrun);
		getMaskedParameterWithReplacement(refFilename, PARAMETER_PULLING_PLANE_REFFILE, trajnum, "<run>");
		readPDB(refFilename, &refFile);
		for(i = 0; i < gsystem.N; i++)
		{
			itot = i + traj*gsystem.N;
			if(refFile.atoms[i].beta > 0)
			{
				potentialData.prc[itot] = PLANE_ATOM_FIXED;
				potentialData.Nfixed++;
			}
			else
				if(refFile.atoms[i].occupancy > 0)
				{
					potentialData.prc[itot] = PLANE_ATOM_PULLED;
					potentialData.Npulled++;
				}
				else
				{
					potentialData.prc[itot] = PLANE_ATOM_FREE;
				}
		}
	}

	if(potentialData.Npulled + potentialData.Nfixed != 0)
	{
		potentialData.Ntotal = potentialData.Npulled + potentialData.Nfixed;
		LOG << "Total amount of atoms, included in Plane-pulling worklist is " << potentialData.Ntotal << ".";
		potentialData.Nplane = 2*parameters.Ntr;
		LOG << "Total amount of planes is " << potentialData.Nplane << ".";
		getVectorParameter(PARAMETER_PULLING_PLANE_NORM, &potentialData.planeNorm.x, &potentialData.planeNorm.y, &potentialData.planeNorm.z, 0, 0, 0, 0);
		potentialData.planeNorm.x = potentialData.planeNorm.x / sqrt(potentialData.planeNorm.x*potentialData.planeNorm.x + potentialData.planeNorm.y*potentialData.planeNorm.y + potentialData.planeNorm.z*potentialData.planeNorm.z);
		potentialData.planeNorm.y = potentialData.planeNorm.y / sqrt(potentialData.planeNorm.x*potentialData.planeNorm.x + potentialData.planeNorm.y*potentialData.planeNorm.y + potentialData.planeNorm.z*potentialData.planeNorm.z);
		potentialData.planeNorm.z = potentialData.planeNorm.z / sqrt(potentialData.planeNorm.x*potentialData.planeNorm.x + potentialData.planeNorm.y*potentialData.planeNorm.y + potentialData.planeNorm.z*potentialData.planeNorm.z);
		LOG << "Plane normal is set to (" << potentialData.planeNorm.x << " " << potentialData.planeNorm.y << " " << potentialData.planeNorm.z << ").";
		getVectorParameter(PARAMETER_PULLING_PLANE_POINT, &potentialData.planePoint.x, &potentialData.planePoint.y, &potentialData.planePoint.z, 0, 0, 0, 0);
		LOG << "Plane point is set to (" << potentialData.planePoint.x << " " << potentialData.planePoint.y << " " << potentialData.planePoint.z << ").";
		potentialData.planeMass = getFloatParameter(PARAMETER_PULLING_PLANE_MASS, 0, 0);
		LOG << "Plane mass is set to " << potentialData.planeMass << ".";
		potentialData.fixSpring = getFloatParameter(PARAMETER_PULLING_PLANE_FIXSPRING, 0, 0);
		LOG << "Plane fixing spring is set to " << potentialData.fixSpring << " pN/nm.";
		potentialData.fixSpring = 0.6*potentialData.fixSpring;
		allocateCPU((void**)&potentialData.h_force, potentialData.Ntotal*sizeof(float));
		allocateGPU((void**)&potentialData.d_force, potentialData.Ntotal*sizeof(float));
		for(i = 0; i < potentialData.Ntotal; i++)
		{
			potentialData.h_force[i] = 0;
		}
		allocateCPU((void**)&potentialData.h_planeDisplacement, potentialData.Nplane*sizeof(float));
		allocateGPU((void**)&potentialData.d_planeDisplacement, potentialData.Nplane*sizeof(float));
		for(i = 0; i < potentialData.Nplane; i++)
		{
			potentialData.h_planeDisplacement[i] = 0;
		}
		allocateCPU((void**)&potentialData.h_workList, potentialData.Ntotal*sizeof(WorkList));
		allocateGPU((void**)&potentialData.d_workList, potentialData.Ntotal*sizeof(WorkList));
		initWorkList();
	}
	else
	{
		DIE("No atoms were selected for Pulling-Plane protocol, while it was switched on.");
	}
}

void initWorkList()
{
	int i, j, traj, itot;
	float4 d;
	j = 0;
	//1tr
	//traj = 0;
	//ntr

	for(traj = 0; traj < parameters.Ntr; traj++)
	{
		for(i = 0; i < gsystem.N; i++)
		{
			itot = i + traj*gsystem.N;
			if(potentialData.prc[i] == PLANE_ATOM_FIXED)
			{
				DPRINTF("Atom %d-%d\t(%s%d-%s) in trajectory %d will be connected to the fixed plane.\n",
						itot, i, topology.atoms[i].resName, topology.atoms[i].resid, topology.atoms[i].name, traj+parameters.firstrun);
				LOG << "Atom " << itot << "-" << i << "\t(" << topology.atoms[i].resName << topology.atoms[i].resid << " - " << topology.atoms[i].name << ") " << " will be connected to the fixed plane.";
				potentialData.h_workList[j].atomID = itot;
				potentialData.h_workList[j].planeID = traj*2;
				j++;
			}
			if(potentialData.prc[i] == PLANE_ATOM_PULLED)
			{
				DPRINTF("Atom %d-%d\t(%s%d-%s) in trajectory %d will be connected to the pulled plane.\n",
						itot, i, topology.atoms[i].resName, topology.atoms[i].resid, topology.atoms[i].name, traj+parameters.firstrun);
				LOG << "Atom " << itot << "-" << i << "\t(" << topology.atoms[i].resName << topology.atoms[i].resid << " - " << topology.atoms[i].name << ") " << " will be connected to the pulled plane.";
				potentialData.h_workList[j].atomID = itot;
				potentialData.h_workList[j].planeID = traj*2 + 1;
				j++;
			}
		}
	}

	if(getYesNoParameter(PARAMETER_PULLING_PLANE_USE_PDB, 0))
	{
		float3* refCoords = (float3*)calloc(gsystem.Ntot, sizeof(float3));
		PDB refPDB;
		char trajnum[10];
		char* refPDBFilename = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
		for(traj = 0; traj < parameters.Ntr; traj++)
		{
			sprintf(trajnum, "%d", traj + parameters.firstrun);
			getMaskedParameterWithReplacement(refPDBFilename, PARAMETER_PULLING_PLANE_COORDS_PDB, trajnum, "<run>");
			readPDB(refPDBFilename, &refPDB);
			for(i = 0; i < gsystem.N; i++)
			{
				itot = i + traj*gsystem.N;
				refCoords[itot].x = refPDB.atoms[i].x;
				refCoords[itot].y = refPDB.atoms[i].y;
				refCoords[itot].z = refPDB.atoms[i].z;
			}
		}
		for(i = 0; i < potentialData.Ntotal; i++)
		{
			j = potentialData.h_workList[i].atomID;
			d.x = potentialData.planePoint.x - refCoords[j].x;
			d.y = potentialData.planePoint.y - refCoords[j].y;
			d.z = potentialData.planePoint.z - refCoords[j].z;
			potentialData.h_workList[i].bDistance = d.x*potentialData.planeNorm.x + d.y*potentialData.planeNorm.y + d.z*potentialData.planeNorm.z;
		}
	} else {
		for(i = 0; i < potentialData.Ntotal; i++)
		{
			j = potentialData.h_workList[i].atomID;
			d.x = potentialData.planePoint.x - gsystem.h_coord[j].x;
			d.y = potentialData.planePoint.y - gsystem.h_coord[j].y;
			d.z = potentialData.planePoint.z - gsystem.h_coord[j].z;
			potentialData.h_workList[i].bDistance = d.x*potentialData.planeNorm.x + d.y*potentialData.planeNorm.y + d.z*potentialData.planeNorm.z;
		}
	}

}

__global__ void computePlanePulling_kernel()
{
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_potentialData.Ntotal)
	{
		int atomID = c_potentialData.d_workList[d_i].atomID;
		int planeID = c_potentialData.d_workList[d_i].planeID;
		float4 coord = tex1Dfetch(t_coord, atomID);
		float4 currentpoint = c_potentialData.planePoint;
		float4 norm = c_potentialData.planeNorm;
		//here dis is a plane displacement
		float dis = c_potentialData.d_planeDisplacement[planeID];
		//calculating current plane position
		currentpoint.x += dis*norm.x;
		currentpoint.y += dis*norm.y;
		currentpoint.z += dis*norm.z;
		//here dis is the distance between atom and plane
		dis = (currentpoint.x - coord.x)*norm.x + (currentpoint.y - coord.y)*norm.y + (currentpoint.z - coord.z)*norm.z;
		//here dis is atom's displacement from its balanced state
		dis = dis - c_potentialData.d_workList[d_i].bDistance;
		//here dis becomes a force value
		dis = dis*c_potentialData.fixSpring;
		//here we write the force value to use it for integrating plane movement
		c_potentialData.d_force[d_i] += dis;
		//here we add the forces to all the other forces
		//we use currentpoint to store them
		currentpoint = c_gsystem.d_forces[atomID];
		currentpoint.x += dis*norm.x;
		currentpoint.y += dis*norm.y;
		currentpoint.z += dis*norm.z;
		c_gsystem.d_forces[atomID] = currentpoint;
	}
}

inline void computePlanePulling()
{
	computePlanePulling_kernel<<<planePullingGridSize, planePullingBlockSize>>>();
}

void updatePlaneLocation()
{
	//float4 force;
	cudaMemcpy(potentialData.h_force, potentialData.d_force, potentialData.Ntotal*sizeof(float),cudaMemcpyDeviceToHost);
	float force, pullforce, atomforce;
	int i, j, traj;
	traj = 0;
	for(i = 1; i < potentialData.Nplane; i = i + 2)
	{
		force = 0;
		atomforce = 0;
		for(j = 0; j < potentialData.Ntotal; j++)
		{
			if(potentialData.h_workList[j].planeID == i)
			{
				atomforce += potentialData.h_force[j];
				potentialData.h_force[j] = 0;
			}
		}
		//here take - instead of +, because this is a force applied to the plane
		//LOG << -force*potentialData.planeNorm.x/planeLocationUpdater.frequency <<" "<< -force*potentialData.planeNorm.y/planeLocationUpdater.frequency <<" "<< -force*potentialData.planeNorm.z/planeLocationUpdater.frequency <<" ";
		//LOG << "Plane " << i << ": atomic force increased to\t"
		//		<< -force*potentialData.planeNorm.x/planeLocationUpdater.frequency <<" "<< -force*potentialData.planeNorm.y/planeLocationUpdater.frequency <<" "<< -force*potentialData.planeNorm.z/planeLocationUpdater.frequency <<" ";
		//compute force for CF
		if(potentialData.pullForce != 0)
		{
			atomforce = -atomforce/planeLocationUpdater.frequency;
			force = atomforce + potentialData.pullForce;
		}
		//compute force for CV
		else
		{
			pullforce = potentialData.pullSpring*(potentialData.pullSpeed*step*integrator->h - potentialData.h_planeDisplacement[i]);
			atomforce = -atomforce/planeLocationUpdater.frequency;
			force = atomforce + pullforce;
		}
		//here force somehow becomes speed
		force = force/potentialData.planeMass;
		//compute current displacement
		potentialData.h_planeDisplacement[i] += force*planeLocationUpdater.frequency*integrator->h;
		//LOG << "Plane " << i << ": relocating point position to\t"
		//		<< potentialData.planePoint.x + potentialData.planeNorm.x*potentialData.h_planeDisplacement[i] << "(" << potentialData.planeNorm.x*potentialData.h_planeDisplacement[i] << ") "
		//		<< potentialData.planePoint.y + potentialData.planeNorm.y*potentialData.h_planeDisplacement[i] << "(" << potentialData.planeNorm.y*potentialData.h_planeDisplacement[i] << ") "
		//		<< potentialData.planePoint.z + potentialData.planeNorm.z*potentialData.h_planeDisplacement[i] << "(" << potentialData.planeNorm.z*potentialData.h_planeDisplacement[i] << ").";
		outputFile = safe_fopen(outputFilename[traj], "a");
		if(potentialData.pullForce != 0)
		{
			if(step%outputFreq == 0 && step != 0)
			{
				fprintf(outputFile, "%lld\t%f\t%f\t%f\t%f\t%f\n",
					step,
					potentialData.planePoint.x + potentialData.planeNorm.x*potentialData.h_planeDisplacement[i],
					potentialData.planePoint.y + potentialData.planeNorm.y*potentialData.h_planeDisplacement[i],
					potentialData.planePoint.z + potentialData.planeNorm.z*potentialData.h_planeDisplacement[i],
					potentialData.h_planeDisplacement[i],
					atomforce/0.6);
				if((step%logOutputFreq - outputFreq == 0) && traj == 0)
				{
					printf("Pulling-plane CF output:\n");
					printf("%*s%*s%*s%*s%*s%*s%*s\n",
							ENERGY_OUTPUT_WIDTH, "Step",
							ENERGY_OUTPUT_WIDTH, "Plane",
							ENERGY_OUTPUT_WIDTH, "X",
							ENERGY_OUTPUT_WIDTH, "Y",
							ENERGY_OUTPUT_WIDTH, "Z",
							ENERGY_OUTPUT_WIDTH, "Displacement",
							ENERGY_OUTPUT_WIDTH, "Atom Force");
				}
				printf("%*lld%*d%*f%*f%*f%*f%*f\n",
						ENERGY_OUTPUT_WIDTH, step,
						ENERGY_OUTPUT_WIDTH, i,
						ENERGY_OUTPUT_WIDTH, potentialData.planePoint.x + potentialData.planeNorm.x*potentialData.h_planeDisplacement[i],
						ENERGY_OUTPUT_WIDTH, potentialData.planePoint.y + potentialData.planeNorm.y*potentialData.h_planeDisplacement[i],
						ENERGY_OUTPUT_WIDTH, potentialData.planePoint.z + potentialData.planeNorm.z*potentialData.h_planeDisplacement[i],
						ENERGY_OUTPUT_WIDTH, potentialData.h_planeDisplacement[i],
						ENERGY_OUTPUT_WIDTH, atomforce/0.6);
			}
		}
		else
		{
			if(step%outputFreq == 0 && step != 0)
			{
				fprintf(outputFile, "%lld\t%f\t%f\t%f\t%f\t%f\t%f\n",
					step,
					potentialData.planePoint.x + potentialData.planeNorm.x*potentialData.h_planeDisplacement[i],
					potentialData.planePoint.y + potentialData.planeNorm.y*potentialData.h_planeDisplacement[i],
					potentialData.planePoint.z + potentialData.planeNorm.z*potentialData.h_planeDisplacement[i],
					potentialData.h_planeDisplacement[i],
					atomforce/0.6,
					pullforce/0.6);
				if((step%logOutputFreq - outputFreq == 0) && traj == 0)
				{
					printf("Pulling-plane CV output:\n");
					printf("%*s%*s%*s%*s%*s%*s%*s%*s\n",
							ENERGY_OUTPUT_WIDTH, "Step",
							ENERGY_OUTPUT_WIDTH, "Plane",
							ENERGY_OUTPUT_WIDTH, "X",
							ENERGY_OUTPUT_WIDTH, "Y",
							ENERGY_OUTPUT_WIDTH, "Z",
							ENERGY_OUTPUT_WIDTH, "Displacement",
							ENERGY_OUTPUT_WIDTH, "Atom Force",
							ENERGY_OUTPUT_WIDTH, "Pull Force");
				}
				printf("%*lld%*d%*f%*f%*f%*f%*f%*f\n",
						ENERGY_OUTPUT_WIDTH, step,
						ENERGY_OUTPUT_WIDTH, i,
						ENERGY_OUTPUT_WIDTH, potentialData.planePoint.x + potentialData.planeNorm.x*potentialData.h_planeDisplacement[i],
						ENERGY_OUTPUT_WIDTH, potentialData.planePoint.y + potentialData.planeNorm.y*potentialData.h_planeDisplacement[i],
						ENERGY_OUTPUT_WIDTH, potentialData.planePoint.z + potentialData.planeNorm.z*potentialData.h_planeDisplacement[i],
						ENERGY_OUTPUT_WIDTH, potentialData.h_planeDisplacement[i],
						ENERGY_OUTPUT_WIDTH, atomforce/0.6,
						ENERGY_OUTPUT_WIDTH, pullforce/0.6);
			}
		}
		fclose(outputFile);
		traj++;
	}
	cudaMemcpy(potentialData.d_force, potentialData.h_force, potentialData.Ntotal*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(potentialData.d_planeDisplacement, potentialData.h_planeDisplacement, potentialData.Nplane*sizeof(float),cudaMemcpyHostToDevice);
}

void destroyPlaneLocationUpdater()
{

}

void destroy()
{

}

#undef LOG

} //namespace pulling_plane_potential
