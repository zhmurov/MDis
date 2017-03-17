/*
 * PushingPlanePotential.cu
 *
 *  Created on: Apr 20, 2016
 *      Author: kir_min
 */
#include "../Core/global.h"
#include "../Util/Log.h"
#include "DrumPotential.cuh"

namespace drum_potential
{

class Log: public ILog
{
	virtual void Write(const char* message) const
	{
		std::cout << makeTimePrefix() << "<drum_potential> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void create()
{
	if(getYesNoParameter(PARAMETER_DRUM_POTENTIAL, DEFAULT_DRUM_POTENTIAL))
	{
		potential.compute = &computeDrumPotential;
		potential.destroy = &destroy;
		sprintf(potential.name, "Drum potential");
		potentials[potentialsCount] = &potential;
		potentialsCount++;
		
		init();
		LOG << "Done initializing Drum potential.";
	}
	else
	{
		LOG << "Drum didn't include.";
	}
}

void init()
{
	LOG << "Initializing Pulling-Plane potential...";
	
	blockSize = BLOCK_SIZE;
	blockCount = gsystem.Ntot/BLOCK_SIZE + 1;

	potentialData.sigma = getFloatParameter(PARAMETER_DRUM_POTENTIAL_SIGMA, 1.0, 1);
	potentialData.epsilon = getFloatParameter(PARAMETER_DRUM_POTENTIAL_EPSILON, 1.0, 1);
	potentialData.R = getFloatParameter(PARAMETER_DRUM_POTENTIAL_RADIUS);

	getVectorParameter(PARAMETER_DRUM_POTENTIAL_POINT, &potentialData.Point.x, &potentialData.Point.y, &potentialData.Point.z);
	getVectorParameter(PARAMETER_DRUM_POTENTIAL_NORM, &potentialData.Norm.x, &potentialData.Norm.y, &potentialData.Norm.z);	
	
	float4 h_n = potentialData.Norm;
	float h_n2 = h_n.x*h_n.x + h_n.y*h_n.y + h_n.z*h_n.z;
	float h_n_mod = sqrt(h_n2);
	h_n.x /= h_n_mod;
	h_n.y /= h_n_mod;
	h_n.z /= h_n_mod;
	potentialData.Norm = h_n;

	XYZ c_xyz;
	c_xyz.atomCount = 36*70;
	c_xyz.atoms = (XYZAtom*)calloc(c_xyz.atomCount, sizeof(XYZAtom));
	int alp;
	float h = 0;
	int count = 0;
	while (count < c_xyz.atomCount)
	{
		for(alp = 0; alp < 360; alp += 10)
		{
			c_xyz.atoms[count].name = 'C';
			c_xyz.atoms[count].x = h_n.x*h + 10.0*potentialData.Point.x;
			c_xyz.atoms[count].y = 10.0*cos(alp*M_PI/180.0)*potentialData.R + 10.0*potentialData.Point.y;
			c_xyz.atoms[count].z = 10.0*sin(alp*M_PI/180.0)*potentialData.R + 10.0*potentialData.Point.z;
			count ++;	
		}
		h ++;
	}
	getMaskedParameter(outputfilename, PARAMETER_DRUM_POTENTIAL_XYZ_OUTPUT_FILE);
	writeXYZ(outputfilename, &c_xyz);

	cudaMemcpyToSymbol(c_potentialData, &potentialData, sizeof(PotentialData), 0, cudaMemcpyHostToDevice);



}

__global__ void computeDrum_kernel()
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < c_gsystem.Ntot){
		float4 ri = c_gsystem.d_coord[i];
		float4 f = c_gsystem.d_forces[i];
		float4 r0 = c_potentialData.Point;
		float4 n = c_potentialData.Norm;
		float R = c_potentialData.R;
		float4 ri0;
		ri0.x = ri.x - r0.x;
		ri0.y = ri.y - r0.y;
		ri0.z = ri.z - r0.z;
		float scal_n_ri0 = n.x*ri0.x + n.y*ri0.y + n.z*ri0.z;
		float4 vec;
		vec.x = scal_n_ri0*n.x - ri0.x; 
		vec.y = scal_n_ri0*n.y - ri0.y;
		vec.z = scal_n_ri0*n.z - ri0.z;
		float vec2 = vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
		float mod_vec = sqrt(vec2);
		float4 r;
		r.x = (R - mod_vec)*vec.x/mod_vec;
		r.y = (R - mod_vec)*vec.y/mod_vec;
		r.z = (R - mod_vec)*vec.z/mod_vec;
		float r2 = r.x*r.x + r.y*r.y + r.z*r.z;
		float sor2 = c_potentialData.sigma*c_potentialData.sigma/r2;
		float sor6 = sor2*sor2*sor2;
		float mul = 6.0*c_potentialData.epsilon*sor6/r2;
		f.x += mul*r.x;
		f.y += mul*r.y;
		f.z += mul*r.z;
		c_gsystem.d_forces[i] = f;
	}
}

inline void computeDrumPotential()
{
	computeDrum_kernel<<<blockSize, blockCount>>>();
}

void destroy()
{

}

}
