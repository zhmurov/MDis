/*
 * ImproperPotential.cu
 *
 *  Created on: Aug 5, 2010
 *      Author: zhmurov
 */

#include "../Core/global.h"
#include "../Core/md.cuh"
#include "../Util/Log.h"
#include "ImproperPotential.cuh"

namespace improper_potential {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<harmonic_potential> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void create(){
	potential.compute = &compute;
	potential.destroy = &destroy;
	sprintf(potential.name, "Improper potentail");
	potentials[potentialsCount] = &potential;
	potentialsCount ++;
	energyOutput.computeValues = &computeEnergy;
	allocateCPU((void**)&energyOutput.values, parameters.Ntr*sizeof(float));
	strcpy(energyOutput.name, ENERGY_OUTPUT_NAME_IMPROPER);
	energyOutputs[energyOutputsCount] = &energyOutput;
	energyOutputsCount ++;
	init();
}

void init(){
	LOG << "Initializing improper potential...";
	improperData.I = topology.improperCount;
	improperData.Itot = topology.improperCount*parameters.Ntr;
	improperBlockSize = BLOCK_SIZE;
	improperBlockCount = improperData.Itot/BLOCK_SIZE + 1;
	improperSummBlockSize = BLOCK_SIZE;
	improperSummBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;

	if(improperData.Itot > 0){

		allocateCPU((void**)&improperData.h_improperCount, gsystem.Ntot*sizeof(int));
		allocateGPU((void**)&improperData.d_improperCount, gsystem.Ntot*sizeof(int));

		Improper improper;
		int i, d;

		for(i = 0; i < gsystem.N; i++){
			improperData.h_improperCount[i] = 0;
		}

		for(d = 0; d < topology.improperCount; d++){
			improper = topology.impropers[d];
			improperData.h_improperCount[improper.i] ++;
			improperData.h_improperCount[improper.j] ++;
			improperData.h_improperCount[improper.k] ++;
			improperData.h_improperCount[improper.l] ++;
		}

		improperData.maxImpropersPerAtom = 0;
		for(i = 0; i < gsystem.N; i++){
			if(improperData.h_improperCount[i] > improperData.maxImpropersPerAtom){
				improperData.maxImpropersPerAtom = improperData.h_improperCount[i];
			}
		}
		LOG << "Maximum impropers per atom is " << improperData.maxImpropersPerAtom;


		allocateCPU((void**)&improperData.h_impropers, improperData.Itot*sizeof(int4));
		allocateGPU((void**)&improperData.d_impropers, improperData.Itot*sizeof(int4));
		allocateCPU((void**)&improperData.h_improperRefs, improperData.Itot*sizeof(int4));
		allocateGPU((void**)&improperData.d_improperRefs, improperData.Itot*sizeof(int4));
		allocateCPU((void**)&improperData.h_improperParameters, improperData.Itot*sizeof(GImproperParameters));
		allocateGPU((void**)&improperData.d_improperParameters, improperData.Itot*sizeof(GImproperParameters));
		allocateCPU((void**)&improperData.h_improperForces,
				gsystem.widthTot*improperData.maxImpropersPerAtom*sizeof(float4));
		allocateGPU((void**)&improperData.d_improperForces,
				gsystem.widthTot*improperData.maxImpropersPerAtom*sizeof(float4));

		allocateCPU((void**)&improperData.h_improperEnergies, improperData.Itot*sizeof(float));
		allocateGPU((void**)&improperData.d_improperEnergies, improperData.Itot*sizeof(float));

		int multiplicityCount = 0;
		for(d = 0; d < topology.improperCount; d++){
			improper = topology.impropers[d];
			if(improper.multiplicity > 1){
				multiplicityCount += improper.multiplicity - 1;
			}
		}
		LOG << "Found " << multiplicityCount << " multiplicity values.";
		allocateCPU((void**)&improperData.h_multiplicityParameters,
				parameters.Ntr*multiplicityCount*sizeof(GImproperParameters));
		allocateGPU((void**)&improperData.d_multiplicityParameters,
				parameters.Ntr*multiplicityCount*sizeof(GImproperParameters));


		for(i = 0; i < gsystem.N; i++){
			improperData.h_improperCount[i] = 0;
		}

		multiplicityCount = 0;
		for(d = 0; d < topology.improperCount; d++){
			improper = topology.impropers[d];
			improperData.h_impropers[d].x = improper.i;
			improperData.h_impropers[d].y = improper.j;
			improperData.h_impropers[d].z = improper.k;
			improperData.h_impropers[d].w = improper.l;
			improperData.h_improperParameters[d].kpsi = improper.kpsi[0];
			improperData.h_improperParameters[d].psi0 = improper.psi0[0];
			improperData.h_improperParameters[d].n = (float)improper.n[0];
			if(improper.multiplicity == 1){
				improperData.h_improperParameters[d].multiplicityRef = -1;
			} else {
				DPRINTF("Adding multiplicity values (multiplicity = %d for improper #%d)\n",
						improper.multiplicity, d);
				improperData.h_improperParameters[d].multiplicityRef = multiplicityCount;
				for(i = 1; i < improper.multiplicity; i++){
					improperData.h_multiplicityParameters[multiplicityCount].kpsi = improper.kpsi[0];
					improperData.h_multiplicityParameters[multiplicityCount].psi0 = improper.psi0[0];
					improperData.h_multiplicityParameters[multiplicityCount].n = (float)improper.n[0];
					improperData.h_multiplicityParameters[multiplicityCount].multiplicityRef = multiplicityCount + 1;
					multiplicityCount ++;
				}
				improperData.h_multiplicityParameters[multiplicityCount - 1].multiplicityRef = -1;
			}
			improperData.h_improperRefs[d].x = improperData.h_improperCount[improper.i];
			improperData.h_improperRefs[d].y = improperData.h_improperCount[improper.j];
			improperData.h_improperRefs[d].z = improperData.h_improperCount[improper.k];
			improperData.h_improperRefs[d].w = improperData.h_improperCount[improper.l];
			improperData.h_improperCount[improper.i] ++;
			improperData.h_improperCount[improper.j] ++;
			improperData.h_improperCount[improper.k] ++;
			improperData.h_improperCount[improper.l] ++;
		}
		/*for(d = 0; d < topology.improperCount; d++){
			printf("%d: (%d-%d-%d-%d): K = %f, delta = %f, n = %f\n ", d,
					improperData.h_impropers[d].x,
					improperData.h_impropers[d].y,
					improperData.h_impropers[d].z,
					improperData.h_impropers[d].w,
					improperData.h_improperParameters[d].kpsi,
					improperData.h_improperParameters[d].psi0,
					improperData.h_improperParameters[d].n);
		}*/
		for(i = 0; i < gsystem.N; i++){
			if(improperData.h_improperCount[i] > improperData.maxImpropersPerAtom){
				DIE("Maximum impropers per atom exceeded the limit of %d on atom %d",
						improperData.maxImpropersPerAtom, i);
			}
		}

		int traj, itot, m, mtot;
		for(traj = 1; traj < parameters.Ntr; traj++){
			for(i = 0; i < improperData.I; i++){
				itot = traj*improperData.I + i;
				improperData.h_impropers[itot].x = improperData.h_impropers[i].x + traj*gsystem.N;
				improperData.h_impropers[itot].y = improperData.h_impropers[i].y + traj*gsystem.N;
				improperData.h_impropers[itot].z = improperData.h_impropers[i].z + traj*gsystem.N;
				improperData.h_impropers[itot].w = improperData.h_impropers[i].w + traj*gsystem.N;
				improperData.h_improperParameters[itot].kpsi = improperData.h_improperParameters[i].kpsi;
				improperData.h_improperParameters[itot].psi0 = improperData.h_improperParameters[i].psi0;
				improperData.h_improperParameters[itot].n = improperData.h_improperParameters[i].n;
				if(improperData.h_improperParameters[i].multiplicityRef == -1){
					improperData.h_improperParameters[itot].multiplicityRef = -1;
				} else {
					improperData.h_improperParameters[itot].multiplicityRef =
							improperData.h_improperParameters[i].multiplicityRef + multiplicityCount*traj;
				}

				improperData.h_improperRefs[itot].x = improperData.h_improperRefs[i].x;
				improperData.h_improperRefs[itot].y = improperData.h_improperRefs[i].y;
				improperData.h_improperRefs[itot].z = improperData.h_improperRefs[i].z;
				improperData.h_improperRefs[itot].w = improperData.h_improperRefs[i].w;
			}
			for(i = 0; i < gsystem.N; i++){
				itot = traj*gsystem.N + i;
				improperData.h_improperCount[itot] = improperData.h_improperCount[i];
			}
			for(m = 0; m < multiplicityCount; m++){
				mtot = traj*multiplicityCount + m;
				improperData.h_multiplicityParameters[mtot].kpsi = improperData.h_multiplicityParameters[m].kpsi;
				improperData.h_multiplicityParameters[mtot].psi0 = improperData.h_multiplicityParameters[m].psi0;
				improperData.h_multiplicityParameters[mtot].n = improperData.h_multiplicityParameters[m].n;
				if(improperData.h_multiplicityParameters[m].multiplicityRef == -1){
					improperData.h_multiplicityParameters[mtot].multiplicityRef = -1;
				} else {
					improperData.h_multiplicityParameters[mtot].multiplicityRef =
						improperData.h_multiplicityParameters[m].multiplicityRef + traj*multiplicityCount;
				}
			}
		}

		cudaMemcpy(improperData.d_improperCount, improperData.h_improperCount,
				gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(improperData.d_impropers, improperData.h_impropers,
				improperData.Itot*sizeof(int4), cudaMemcpyHostToDevice);
		cudaMemcpy(improperData.d_improperRefs, improperData.h_improperRefs,
				improperData.Itot*sizeof(int4), cudaMemcpyHostToDevice);
		cudaMemcpy(improperData.d_improperParameters, improperData.h_improperParameters,
				improperData.Itot*sizeof(GImproperParameters), cudaMemcpyHostToDevice);
		cudaMemcpy(improperData.d_multiplicityParameters, improperData.h_multiplicityParameters,
				parameters.Ntr*multiplicityCount*sizeof(GImproperParameters), cudaMemcpyHostToDevice);
	}

	cudaMemcpyToSymbol(c_improperData, &improperData,
				sizeof(GImproperData), 0, cudaMemcpyHostToDevice);

	LOG << "Done initializing improper potential";
}

__global__ void improperPotential_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_improperData.Itot){

		int4 improper = c_improperData.d_impropers[d_i];
		int4 ref = c_improperData.d_improperRefs[d_i];
		GImproperParameters par = c_improperData.d_improperParameters[d_i];

		float4 r1 = tex1Dfetch(t_coord, improper.x);
		float4 r2 = tex1Dfetch(t_coord, improper.y);
		float4 r3 = tex1Dfetch(t_coord, improper.z);
		float4 r4 = tex1Dfetch(t_coord, improper.w);

		float3 dr12, dr23, dr34;

		dr12.x = r1.x - r2.x;
		dr12.y = r1.y - r2.y;
		dr12.z = r1.z - r2.z;
		DO_PBC(dr12);

		dr23.x = r2.x - r3.x;
		dr23.y = r2.y - r3.y;
		dr23.z = r2.z - r3.z;
		DO_PBC(dr23);

		dr34.x = r3.x - r4.x;
		dr34.y = r3.y - r4.y;
		dr34.z = r3.z - r4.z;
		DO_PBC(dr34);

		float4 a, b, c;

		a.x = dr12.y*dr23.z - dr12.z*dr23.y;
		a.y = dr12.z*dr23.x - dr12.x*dr23.z;
		a.z = dr12.x*dr23.y - dr12.y*dr23.x;
		a.w = 1.0f/sqrtf(a.x*a.x + a.y*a.y + a.z*a.z);

		b.x = dr23.y*dr34.z - dr23.z*dr34.y;
		b.y = dr23.z*dr34.x - dr23.x*dr34.z;
		b.z = dr23.x*dr34.y - dr23.y*dr34.x;
		b.w = 1.0f/sqrtf(b.x*b.x + b.y*b.y + b.z*b.z);

		c.x = dr23.y*a.z - dr23.z*a.y;
		c.y = dr23.z*a.x - dr23.x*a.z;
		c.z = dr23.x*a.y - dr23.y*a.x;
		c.w = 1.0f/sqrtf(c.x*c.x + c.y*c.y + c.z*c.z);

		float cospsi = (a.x*b.x + a.y*b.y + a.z*b.z)*a.w*b.w;
		float sinpsi = (c.x*b.x + c.y*b.y + c.z*b.z)*c.w*b.w;

		float psi = -atan2(sinpsi, cospsi);

		float mult = 0.0f;
		do{
			if(par.n > 0.0f){
				mult += -par.n*par.kpsi*sinf(par.n*psi - par.psi0);
			} else {
				float diff = psi - par.psi0;
				if(diff < -M_PI){
					diff += 2.0f*M_PI;
				} else
				if(diff > M_PI){
					diff -= 2.0f*M_PI;
				}
				mult += 2.0f*par.kpsi*diff;
			}
			if(par.multiplicityRef != -1){
				par = c_improperData.d_multiplicityParameters[par.multiplicityRef];
			}
		} while(par.multiplicityRef != -1);

		float4 f1, f2, f3;

		b.x *= b.w;
		b.y *= b.w;
		b.z *= b.w;

		if(fabs(sinpsi) > 0.1f){

			a.x *= a.w;
			a.y *= a.w;
			a.z *= a.w;

			float3 dcosda, dcosdb;

			dcosda.x = a.w*(cospsi*a.x - b.x);
			dcosda.y = a.w*(cospsi*a.y - b.y);
			dcosda.z = a.w*(cospsi*a.z - b.z);

			dcosdb.x = b.w*(cospsi*b.x - a.x);
			dcosdb.y = b.w*(cospsi*b.y - a.y);
			dcosdb.z = b.w*(cospsi*b.z - a.z);

			mult = mult/sinpsi;

			f1.x = mult*(dr23.y*dcosda.z - dr23.z*dcosda.y);
			f1.y = mult*(dr23.z*dcosda.x - dr23.x*dcosda.z);
			f1.z = mult*(dr23.x*dcosda.y - dr23.y*dcosda.x);

			f3.x = mult*(dr23.z*dcosdb.y - dr23.y*dcosdb.z);
			f3.y = mult*(dr23.x*dcosdb.z - dr23.z*dcosdb.x);
			f3.z = mult*(dr23.y*dcosdb.x - dr23.x*dcosdb.y);

			f2.x = mult*(dr12.z*dcosda.y - dr12.y*dcosda.z + dr34.y*dcosdb.z - dr34.z*dcosdb.y);
			f2.y = mult*(dr12.x*dcosda.z - dr12.z*dcosda.x + dr34.z*dcosdb.x - dr34.x*dcosdb.z);
			f2.z = mult*(dr12.y*dcosda.x - dr12.x*dcosda.y + dr34.x*dcosdb.y - dr34.y*dcosdb.x);

		} else {

			c.x *= c.w;
			c.y *= c.w;
			c.z *= c.w;

			float3 dsindc, dsindb;

			dsindc.x = c.w*(sinpsi*c.x - b.x);
			dsindc.y = c.w*(sinpsi*c.y - b.y);
			dsindc.z = c.w*(sinpsi*c.z - b.z);

			dsindb.x = b.w*(sinpsi*b.x - c.x);
			dsindb.y = b.w*(sinpsi*b.y - c.y);
			dsindb.z = b.w*(sinpsi*b.z - c.z);

			mult = -mult/cospsi;

			f1.x = mult*((dr23.y*dr23.y + dr23.z*dr23.z)*dsindc.x - dr23.x*dr23.y*dsindc.y - dr23.x*dr23.z*dsindc.z);
			f1.y = mult*((dr23.z*dr23.z + dr23.x*dr23.x)*dsindc.y - dr23.y*dr23.z*dsindc.z - dr23.y*dr23.x*dsindc.x);
			f1.z = mult*((dr23.x*dr23.x + dr23.y*dr23.y)*dsindc.z - dr23.z*dr23.x*dsindc.x - dr23.z*dr23.y*dsindc.y);

			f3.x = mult*(dsindb.y*dr23.z - dsindb.z*dr23.y);
			f3.y = mult*(dsindb.z*dr23.x - dsindb.x*dr23.z);
			f3.z = mult*(dsindb.x*dr23.y - dsindb.y*dr23.x);

			f2.x = mult*(-(dr23.y*dr12.y + dr23.z*dr12.z)*dsindc.x + (2.0f*dr23.x*dr12.y - dr12.x*dr23.y)*dsindc.y
					+ (2.0f*dr23.x*dr12.z - dr12.x*dr23.z)*dsindc.z + dsindb.z*dr34.y - dsindb.y*dr34.z);
			f2.y = mult*(-(dr23.z*dr12.z + dr23.x*dr12.x)*dsindc.y + (2.0f*dr23.y*dr12.z - dr12.y*dr23.z)*dsindc.z
					+ (2.0f*dr23.y*dr12.x - dr12.y*dr23.x)*dsindc.x + dsindb.x*dr34.z - dsindb.z*dr34.x);
			f2.z = mult*(-(dr23.x*dr12.x + dr23.y*dr12.y)*dsindc.z + (2.0f*dr23.z*dr12.x - dr12.z*dr23.x)*dsindc.x
					+ (2.0f*dr23.z*dr12.y - dr12.z*dr23.y)*dsindc.y + dsindb.y*dr34.x - dsindb.x*dr34.y);
		}

		c_improperData.d_improperForces[c_gsystem.widthTot*ref.x + improper.x] = f1;
		f1.x = f2.x - f1.x;
		f1.y = f2.y - f1.y;
		f1.z = f2.z - f1.z;
		c_improperData.d_improperForces[c_gsystem.widthTot*ref.y + improper.y] = f1;
		f2.x = f3.x - f2.x;
		f2.y = f3.y - f2.y;
		f2.z = f3.z - f2.z;
		c_improperData.d_improperForces[c_gsystem.widthTot*ref.z + improper.z] = f2;
		f3.x = -f3.x;
		f3.y = -f3.y;
		f3.z = -f3.z;
		c_improperData.d_improperForces[c_gsystem.widthTot*ref.w + improper.w] = f3;
	}
}

__global__ void summImproperForces_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 f = c_gsystem.d_forces[d_i];
		float4 df;
		int i;
		for(i = 0; i < c_improperData.d_improperCount[d_i]; i++){
			df = c_improperData.d_improperForces[c_gsystem.widthTot*i + d_i];
			f.x += df.x;
			f.y += df.y;
			f.z += df.z;
		}
		c_gsystem.d_forces[d_i] = f;
	}
}

inline void compute(){
	improperPotential_kernel<<<improperBlockCount, improperBlockSize>>>();
	cudaThreadSynchronize();
	summImproperForces_kernel<<<improperSummBlockCount, improperSummBlockSize>>>();

	/*cudaMemcpy(gsystem.h_forces, gsystem.d_forces, atomCount*sizeof(float4), cudaMemcpyDeviceToHost);
	int i;
	float3 force = make_float3(0.0f, 0.0f, 0.0f);
	for(i = 0; i < atomCount; i++){
		force.x += gsystem.h_forces[i].x;
		force.y += gsystem.h_forces[i].y;
		force.z += gsystem.h_forces[i].z;
	}
	printf("Net force (impropers): (%f, %f, %f) %f\n", force.x, force.y, force.z,
			sqrtf(force.x*force.x + force.y*force.y + force.z*force.z));*/
}

__global__ void improperPotentialEnergy_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_improperData.Itot){

		int4 improper = c_improperData.d_impropers[d_i];
		GImproperParameters par = c_improperData.d_improperParameters[d_i];

		float4 r1 = tex1Dfetch(t_coord, improper.x);
		float4 r2 = tex1Dfetch(t_coord, improper.y);
		float4 r3 = tex1Dfetch(t_coord, improper.z);
		float4 r4 = tex1Dfetch(t_coord, improper.w);

		float3 dr12, dr23, dr34;

		dr12.x = r1.x - r2.x;
		dr12.y = r1.y - r2.y;
		dr12.z = r1.z - r2.z;
		DO_PBC(dr12);

		dr23.x = r2.x - r3.x;
		dr23.y = r2.y - r3.y;
		dr23.z = r2.z - r3.z;
		DO_PBC(dr23);

		dr34.x = r3.x - r4.x;
		dr34.y = r3.y - r4.y;
		dr34.z = r3.z - r4.z;
		DO_PBC(dr34);

		float4 a, b, c;

		a.x = dr12.y*dr23.z - dr12.z*dr23.y;
		a.y = dr12.z*dr23.x - dr12.x*dr23.z;
		a.z = dr12.x*dr23.y - dr12.y*dr23.x;
		a.w = 1.0f/sqrtf(a.x*a.x + a.y*a.y + a.z*a.z);

		b.x = dr23.y*dr34.z - dr23.z*dr34.y;
		b.y = dr23.z*dr34.x - dr23.x*dr34.z;
		b.z = dr23.x*dr34.y - dr23.y*dr34.x;
		b.w = 1.0f/sqrtf(b.x*b.x + b.y*b.y + b.z*b.z);

		c.x = dr23.y*a.z - dr23.z*a.y;
		c.y = dr23.z*a.x - dr23.x*a.z;
		c.z = dr23.x*a.y - dr23.y*a.x;
		c.w = 1.0f/sqrtf(c.x*c.x + c.y*c.y + c.z*c.z);

		float cospsi = (a.x*b.x + a.y*b.y + a.z*b.z)*a.w*b.w;
		float sinpsi = (c.x*b.x + c.y*b.y + c.z*b.z)*c.w*b.w;

		float psi = -atan2(sinpsi, cospsi);

		float pot = 0.0f;
		do{
			if(par.n > 0.0f){
				pot += par.kpsi*(1.0f + cosf(par.n*psi - par.psi0));
			} else {
				float diff = psi - par.psi0;
				if(diff < -M_PI){
					diff += 2.0f*M_PI;
				} else
				if(diff > M_PI){
					diff -= 2.0f*M_PI;
				}
				pot += par.kpsi*diff*diff;
			}
			if(par.multiplicityRef != -1){
				par = c_improperData.d_multiplicityParameters[par.multiplicityRef];
			}
		} while(par.multiplicityRef != -1);
		c_improperData.d_improperEnergies[d_i] = pot;
	}
}

inline void computeEnergy(){
	improperPotentialEnergy_kernel<<<improperBlockCount, improperBlockSize>>>();
	cudaMemcpy(improperData.h_improperEnergies, improperData.d_improperEnergies,
				improperData.Itot*sizeof(float), cudaMemcpyDeviceToHost);
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		float pot = 0.0f;
		for(i = 0; i < improperData.I; i++){
			pot += improperData.h_improperEnergies[i + traj*improperData.I];
		}
		energyOutput.values[traj] = pot;
	}
	checkCUDAError("improper energy");
}

void destroy(){

}

#undef LOG

} // namespace improper_potential
