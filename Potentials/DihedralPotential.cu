/*
 * DihedralPotential.cu
 *
 *  Created on: Aug 4, 2010
 *      Author: zhmurov
 */

#include "../Core/global.h"
#include "../Core/md.cuh"
#include "../Util/Log.h"
#include "DihedralPotential.cuh"

namespace dihedral_potential {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<dihedral_potential> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)
	
void create(){
	potential.compute = &compute;
	potential.destroy = &destroy;
	sprintf(potential.name, "Dihedral potential");
	potentials[potentialsCount] = &potential;
	potentialsCount ++;
	energyOutput.computeValues = &computeEnergy;
	allocateCPU((void**)&energyOutput.values, parameters.Ntr*sizeof(float));
	strcpy(energyOutput.name, ENERGY_OUTPUT_NAME_DIHEDRAL);
	energyOutputs[energyOutputsCount] = &energyOutput;
	energyOutputsCount ++;
	init();
}

void init(){
	LOG << "Initializing dihedral potential...";
	dihedralData.D = topology.dihedralCount;
	dihedralData.Dtot = topology.dihedralCount*parameters.Ntr;
	dihedralBlockSize = BLOCK_SIZE;
	dihedralBlockCount = dihedralData.Dtot/BLOCK_SIZE + 1;
	dihedralSummBlockSize = BLOCK_SIZE;
	dihedralSummBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;

	if(dihedralData.Dtot > 0){

		int i, d;
		Dihedral dihedral;

		allocateCPU((void**)&dihedralData.h_dihedralCount, gsystem.Ntot*sizeof(int));
		allocateGPU((void**)&dihedralData.d_dihedralCount, gsystem.Ntot*sizeof(int));

		for(i = 0; i < gsystem.N; i++){
			dihedralData.h_dihedralCount[i] = 0;
		}

		for(d = 0; d < topology.dihedralCount; d++){
			dihedral = topology.dihedrals[d];
			dihedralData.h_dihedralCount[dihedral.i] ++;
			dihedralData.h_dihedralCount[dihedral.j] ++;
			dihedralData.h_dihedralCount[dihedral.k] ++;
			dihedralData.h_dihedralCount[dihedral.l] ++;
		}

		dihedralData.maxDihedralsPerAtom = 0;
		for(i = 0; i < gsystem.N; i++){
			if(dihedralData.h_dihedralCount[i] > dihedralData.maxDihedralsPerAtom){
				dihedralData.maxDihedralsPerAtom = dihedralData.h_dihedralCount[i];
			}
		}
		LOG << "Maximum dihedrals per atom is " << dihedralData.maxDihedralsPerAtom;

		allocateCPU((void**)&dihedralData.h_dihedrals, dihedralData.Dtot*sizeof(int4));
		allocateGPU((void**)&dihedralData.d_dihedrals, dihedralData.Dtot*sizeof(int4));
		allocateCPU((void**)&dihedralData.h_dihedralRefs, dihedralData.Dtot*sizeof(int4));
		allocateGPU((void**)&dihedralData.d_dihedralRefs, dihedralData.Dtot*sizeof(int4));
		allocateCPU((void**)&dihedralData.h_dihedralParameters, dihedralData.Dtot*sizeof(GDihedralParameters));
		allocateGPU((void**)&dihedralData.d_dihedralParameters, dihedralData.Dtot*sizeof(GDihedralParameters));
		allocateCPU((void**)&dihedralData.h_dihedralForces,
				gsystem.widthTot*dihedralData.maxDihedralsPerAtom*sizeof(float4));
		allocateGPU((void**)&dihedralData.d_dihedralForces,
				gsystem.widthTot*dihedralData.maxDihedralsPerAtom*sizeof(float4));

		allocateCPU((void**)&dihedralData.h_dihedralEnergies, dihedralData.Dtot*sizeof(float));
		allocateGPU((void**)&dihedralData.d_dihedralEnergies, dihedralData.Dtot*sizeof(float));


		int multiplicityCount = 0;
		for(d = 0; d < dihedralData.D; d++){
			dihedral = topology.dihedrals[d];
			if(dihedral.multiplicity > 1){
				multiplicityCount += dihedral.multiplicity - 1;
			}
		}
		LOG << "Found " << multiplicityCount << " multiplicity values.";
		allocateCPU((void**)&dihedralData.h_multiplicityParameters,
				parameters.Ntr*multiplicityCount*sizeof(GDihedralParameters));
		allocateGPU((void**)&dihedralData.d_multiplicityParameters,
				parameters.Ntr*multiplicityCount*sizeof(GDihedralParameters));


		for(i = 0; i < gsystem.N; i++){
			dihedralData.h_dihedralCount[i] = 0;
		}

		multiplicityCount = 0;
		for(d = 0; d < topology.dihedralCount; d++){
			dihedral = topology.dihedrals[d];
			dihedralData.h_dihedrals[d].x = dihedral.i;
			dihedralData.h_dihedrals[d].y = dihedral.j;
			dihedralData.h_dihedrals[d].z = dihedral.k;
			dihedralData.h_dihedrals[d].w = dihedral.l;
			dihedralData.h_dihedralParameters[d].kchi = dihedral.kchi[0];
			dihedralData.h_dihedralParameters[d].delta = dihedral.delta[0];
			dihedralData.h_dihedralParameters[d].n = (float)dihedral.n[0];
			if(dihedral.multiplicity == 1){
				dihedralData.h_dihedralParameters[d].multiplicityRef = -1;
			} else {
				/*printf("Adding multiplicity values (multiplicity = %d for dihedral #%d)\n",
						dihedral.multiplicity, d);*/
				dihedralData.h_dihedralParameters[d].multiplicityRef = multiplicityCount;
				for(i = 1; i < dihedral.multiplicity; i++){
					dihedralData.h_multiplicityParameters[multiplicityCount].kchi = dihedral.kchi[i];
					dihedralData.h_multiplicityParameters[multiplicityCount].delta = dihedral.delta[i];
					dihedralData.h_multiplicityParameters[multiplicityCount].n = (float)dihedral.n[i];
					dihedralData.h_multiplicityParameters[multiplicityCount].multiplicityRef = multiplicityCount + 1;
					multiplicityCount ++;
				}
				dihedralData.h_multiplicityParameters[multiplicityCount - 1].multiplicityRef = -1;
			}
			dihedralData.h_dihedralRefs[d].x = dihedralData.h_dihedralCount[dihedral.i];
			dihedralData.h_dihedralRefs[d].y = dihedralData.h_dihedralCount[dihedral.j];
			dihedralData.h_dihedralRefs[d].z = dihedralData.h_dihedralCount[dihedral.k];
			dihedralData.h_dihedralRefs[d].w = dihedralData.h_dihedralCount[dihedral.l];
			dihedralData.h_dihedralCount[dihedral.i] ++;
			dihedralData.h_dihedralCount[dihedral.j] ++;
			dihedralData.h_dihedralCount[dihedral.k] ++;
			dihedralData.h_dihedralCount[dihedral.l] ++;
		}
		//int multiplicityRef;
		/*for(d = 0; d < topology.dihedralCount; d++){
			if(dihedralData.h_dihedralParameters[d].multiplicityRef != -1){
				printf("%d: (%d-%d-%d-%d) - (%s-%s-%s-%s): K = %f, delta = %f, n = %f, Multiplicity ref = %d\n", d,
						dihedralData.h_dihedrals[d].x,
						dihedralData.h_dihedrals[d].y,
						dihedralData.h_dihedrals[d].z,
						dihedralData.h_dihedrals[d].w,
						atomTypes[atoms[dihedralData.h_dihedrals[d].x].typeId].name,
						atomTypes[atoms[dihedralData.h_dihedrals[d].y].typeId].name,
						atomTypes[atoms[dihedralData.h_dihedrals[d].z].typeId].name,
						atomTypes[atoms[dihedralData.h_dihedrals[d].w].typeId].name,
						dihedralData.h_dihedralParameters[d].kchi,
						dihedralData.h_dihedralParameters[d].delta,
						dihedralData.h_dihedralParameters[d].n,
						dihedralData.h_dihedralParameters[d].multiplicityRef);
				multiplicityRef = dihedralData.h_dihedralParameters[d].multiplicityRef;
				while(multiplicityRef != -1){
					printf("\t\t K = %f, delta = %f, n = %f, Multiplicity ref = %d\n",
							dihedralData.h_multiplicityParameters[multiplicityRef].kchi,
							dihedralData.h_multiplicityParameters[multiplicityRef].delta,
							dihedralData.h_multiplicityParameters[multiplicityRef].n,
							dihedralData.h_multiplicityParameters[multiplicityRef].multiplicityRef);
					multiplicityRef = dihedralData.h_multiplicityParameters[multiplicityRef].multiplicityRef;
				}
			}
		}*/
		for(i = 0; i < gsystem.N; i++){
			if(dihedralData.h_dihedralCount[i] > dihedralData.maxDihedralsPerAtom){
				DIE("Maximum dihedrals per atom exceeded the limit of %d on atom %d\n",
						dihedralData.maxDihedralsPerAtom, i);
			}
		}

		int traj, dtot, itot, m, mtot;
		for(traj = 1; traj < parameters.Ntr; traj++){
			for(d = 0; d < dihedralData.D; d++){
				dtot = traj*dihedralData.D + d;
				dihedralData.h_dihedrals[dtot].x = dihedralData.h_dihedrals[d].x + traj*gsystem.N;
				dihedralData.h_dihedrals[dtot].y = dihedralData.h_dihedrals[d].y + traj*gsystem.N;
				dihedralData.h_dihedrals[dtot].z = dihedralData.h_dihedrals[d].z + traj*gsystem.N;
				dihedralData.h_dihedrals[dtot].w = dihedralData.h_dihedrals[d].w + traj*gsystem.N;
				dihedralData.h_dihedralParameters[dtot].kchi = dihedralData.h_dihedralParameters[d].kchi;
				dihedralData.h_dihedralParameters[dtot].delta = dihedralData.h_dihedralParameters[d].delta;
				dihedralData.h_dihedralParameters[dtot].n = dihedralData.h_dihedralParameters[d].n;
				if(dihedralData.h_dihedralParameters[d].multiplicityRef == -1){
					dihedralData.h_dihedralParameters[dtot].multiplicityRef = -1;
				} else {
					dihedralData.h_dihedralParameters[dtot].multiplicityRef =
							dihedralData.h_dihedralParameters[d].multiplicityRef + multiplicityCount*traj;
				}

				dihedralData.h_dihedralRefs[dtot].x = dihedralData.h_dihedralRefs[d].x;
				dihedralData.h_dihedralRefs[dtot].y = dihedralData.h_dihedralRefs[d].y;
				dihedralData.h_dihedralRefs[dtot].z = dihedralData.h_dihedralRefs[d].z;
				dihedralData.h_dihedralRefs[dtot].w = dihedralData.h_dihedralRefs[d].w;
			}
			for(i = 0; i < gsystem.N; i++){
				itot = traj*gsystem.N + i;
				dihedralData.h_dihedralCount[itot] = dihedralData.h_dihedralCount[i];
			}
			for(m = 0; m < multiplicityCount; m++){
				mtot = traj*multiplicityCount + m;
				dihedralData.h_multiplicityParameters[mtot].kchi = dihedralData.h_multiplicityParameters[m].kchi;
				dihedralData.h_multiplicityParameters[mtot].delta = dihedralData.h_multiplicityParameters[m].delta;
				dihedralData.h_multiplicityParameters[mtot].n = dihedralData.h_multiplicityParameters[m].n;
				if(dihedralData.h_multiplicityParameters[m].multiplicityRef == -1){
					dihedralData.h_multiplicityParameters[mtot].multiplicityRef = -1;
				} else {
					dihedralData.h_multiplicityParameters[mtot].multiplicityRef =
						dihedralData.h_multiplicityParameters[m].multiplicityRef + traj*multiplicityCount;
				}
			}
		}


		cudaMemcpy(dihedralData.d_dihedralCount, dihedralData.h_dihedralCount,
				gsystem.Ntot*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(dihedralData.d_dihedrals, dihedralData.h_dihedrals,
				dihedralData.Dtot*sizeof(int4), cudaMemcpyHostToDevice);
		cudaMemcpy(dihedralData.d_dihedralRefs, dihedralData.h_dihedralRefs,
				dihedralData.Dtot*sizeof(int4), cudaMemcpyHostToDevice);
		cudaMemcpy(dihedralData.d_dihedralParameters, dihedralData.h_dihedralParameters,
				dihedralData.Dtot*sizeof(GDihedralParameters), cudaMemcpyHostToDevice);
		cudaMemcpy(dihedralData.d_multiplicityParameters, dihedralData.h_multiplicityParameters,
				parameters.Ntr*multiplicityCount*sizeof(GDihedralParameters), cudaMemcpyHostToDevice);
	}

	cudaMemcpyToSymbol(c_dihedralData, &dihedralData,
				sizeof(GDihedralData), 0, cudaMemcpyHostToDevice);

	LOG << "Done initializing dihedral potential.";
}

__global__ void dihedralPotentialCHARMM_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_dihedralData.Dtot){

		int4 dihedral = c_dihedralData.d_dihedrals[d_i];
		int4 ref = c_dihedralData.d_dihedralRefs[d_i];
		GDihedralParameters par = c_dihedralData.d_dihedralParameters[d_i];

		float4 r1 = tex1Dfetch(t_coord, dihedral.x);
		float4 r2 = tex1Dfetch(t_coord, dihedral.y);
		float4 r3 = tex1Dfetch(t_coord, dihedral.z);
		float4 r4 = tex1Dfetch(t_coord, dihedral.w);

		float3 dr12, dr23, dr34;

		dr12.x = r1.x - r2.x;
		dr12.y = r1.y - r2.y;
		dr12.z = r1.z - r2.z;
		DO_PBC(dr12);

		dr23.x = r2.x - r3.x;
		dr23.y = r2.y - r3.y;
		dr23.z = r2.z - r3.z;
		DO_PBC(dr23);
		float r232 = dr23.x*dr23.x + dr23.y*dr23.y + dr23.z*dr23.z;

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

		float coschi = (a.x*b.x + a.y*b.y + a.z*b.z)*a.w*b.w;
		float sinchi = -sqrtf(r232)*a.w*b.w*(a.x*dr34.x + a.y*dr34.y + a.z*dr34.z);//(c.x*b.x + c.y*b.y + c.z*b.z)*c.w*b.w;

		float chi = -atan2(sinchi, coschi);

		float mult = 0.0f;

		/*do{
			if(par.n > 0.0f){
				mult += -par.n*par.kchi*sinf(par.n*chi - par.delta);
			} else {
				float diff = chi - par.delta;
				if(diff < -M_PI){
					diff += 2.0f*M_PI;
				} else
				if(diff > M_PI){
					diff -= 2.0f*M_PI;
				}
				mult += 2.0f*par.kchi*diff;
			}
			if(par.multiplicityRef != -1){
				par = c_dihedralData.d_multiplicityParameters[par.multiplicityRef];
			}
		} while(par.multiplicityRef != -1);*/
		int out = 0;
		do{
			if(par.n > 0.0f){
				int i;
				int n = (int)par.n;
				float E1 = 1.0f;
				float df = 0.0f;
				float ddf;

				for(i = 0; i < n; i++){
					ddf = E1*coschi - df*sinchi;
					df = E1*sinchi + df*coschi;
					E1 = ddf;
				}
				mult += par.n*par.kchi*(df*cosf(par.delta) - ddf*sinf(par.delta));//sinf(par.n*chi - par.delta);
			} else {
				float diff = chi - par.delta;
				if(diff < -M_PI){
					diff += 2.0f*M_PI;
				} else
				if(diff > M_PI){
					diff -= 2.0f*M_PI;
				}
				mult += 2.0f*par.kchi*diff;
			}
			if(par.multiplicityRef == -1){
				out = 1;
			} else {
				par = c_dihedralData.d_multiplicityParameters[par.multiplicityRef];
			}
		} while(out != 1);

		float4 f1, f2, f3;

		f1.w = -mult*a.w*a.w*sqrtf(r232);
		f1.x = f1.w*a.x;
		f1.y = f1.w*a.y;
		f1.z = f1.w*a.z;

		f2.w = mult/sqrtf(r232);
		float r1223 = dr12.x*dr23.x + dr12.y*dr23.y + dr12.z*dr23.z;
		float r2334 = dr23.x*dr34.x + dr23.y*dr34.y + dr23.z*dr34.z;
		float m1 = f2.w*r1223*a.w*a.w;
		float m2 = f2.w*r2334*b.w*b.w;
		f2.x = m1*a.x + m2*b.x;
		f2.y = m1*a.y + m2*b.y;
		f2.z = m1*a.z + m2*b.z;

		f3.w = -mult*b.w*b.w*sqrtf(r232);
		f3.x = f3.w*b.x;
		f3.y = f3.w*b.y;
		f3.z = f3.w*b.z;


		c_dihedralData.d_dihedralForces[c_gsystem.widthTot*ref.x + dihedral.x] = f1;
		f1.x = f2.x - f1.x;
		f1.y = f2.y - f1.y;
		f1.z = f2.z - f1.z;
		c_dihedralData.d_dihedralForces[c_gsystem.widthTot*ref.y + dihedral.y] = f1;
		f2.x = f3.x - f2.x;
		f2.y = f3.y - f2.y;
		f2.z = f3.z - f2.z;
		c_dihedralData.d_dihedralForces[c_gsystem.widthTot*ref.z + dihedral.z] = f2;
		f3.x = -f3.x;
		f3.y = -f3.y;
		f3.z = -f3.z;
		c_dihedralData.d_dihedralForces[c_gsystem.widthTot*ref.w + dihedral.w] = f3;
	}
}

__global__ void dihedralPotential_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_dihedralData.Dtot){

		int4 dihedral = c_dihedralData.d_dihedrals[d_i];
		int4 ref = c_dihedralData.d_dihedralRefs[d_i];
		GDihedralParameters par = c_dihedralData.d_dihedralParameters[d_i];

		float4 r1 = tex1Dfetch(t_coord, dihedral.x);
		float4 r2 = tex1Dfetch(t_coord, dihedral.y);
		float4 r3 = tex1Dfetch(t_coord, dihedral.z);
		float4 r4 = tex1Dfetch(t_coord, dihedral.w);

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

		float coschi = (a.x*b.x + a.y*b.y + a.z*b.z)*a.w*b.w;
		float sinchi = (c.x*b.x + c.y*b.y + c.z*b.z)*c.w*b.w;

		float chi = -atan2(sinchi, coschi);

		float mult = 0.0f;
		if(par.n > 0.0f){
			mult += -par.n*par.kchi*sinf(par.n*chi - par.delta);
		} else {
			float diff = chi - par.delta;
			if(diff < -M_PI){
				diff += 2.0f*M_PI;
			} else
			if(diff > M_PI){
				diff -= 2.0f*M_PI;
			}
			mult += 2.0f*par.kchi*diff;
		}
		while(par.multiplicityRef != -1){
			par = c_dihedralData.d_multiplicityParameters[par.multiplicityRef];
			if(par.n > 0.0f){
				mult += -par.n*par.kchi*sinf(par.n*chi - par.delta);
			} else {
				float diff = chi - par.delta;
				if(diff < -M_PI){
					diff += 2.0f*M_PI;
				} else
				if(diff > M_PI){
					diff -= 2.0f*M_PI;
				}
				mult += 2.0f*par.kchi*diff;
			}
		}

		float4 f1, f2, f3;

		b.x *= b.w;
		b.y *= b.w;
		b.z *= b.w;

		if(fabs(sinchi) > 0.1f){

			a.x *= a.w;
			a.y *= a.w;
			a.z *= a.w;

			float3 dcosda, dcosdb;

			dcosda.x = a.w*(coschi*a.x - b.x);
			dcosda.y = a.w*(coschi*a.y - b.y);
			dcosda.z = a.w*(coschi*a.z - b.z);

			dcosdb.x = b.w*(coschi*b.x - a.x);
			dcosdb.y = b.w*(coschi*b.y - a.y);
			dcosdb.z = b.w*(coschi*b.z - a.z);

			mult = mult/sinchi;

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

			dsindc.x = c.w*(sinchi*c.x - b.x);
			dsindc.y = c.w*(sinchi*c.y - b.y);
			dsindc.z = c.w*(sinchi*c.z - b.z);

			dsindb.x = b.w*(sinchi*b.x - c.x);
			dsindb.y = b.w*(sinchi*b.y - c.y);
			dsindb.z = b.w*(sinchi*b.z - c.z);

			mult = -mult/coschi;

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

		c_dihedralData.d_dihedralForces[c_gsystem.widthTot*ref.x + dihedral.x] = f1;
		f1.x = f2.x - f1.x;
		f1.y = f2.y - f1.y;
		f1.z = f2.z - f1.z;
		c_dihedralData.d_dihedralForces[c_gsystem.widthTot*ref.y + dihedral.y] = f1;
		f2.x = f3.x - f2.x;
		f2.y = f3.y - f2.y;
		f2.z = f3.z - f2.z;
		c_dihedralData.d_dihedralForces[c_gsystem.widthTot*ref.z + dihedral.z] = f2;
		f3.x = -f3.x;
		f3.y = -f3.y;
		f3.z = -f3.z;
		c_dihedralData.d_dihedralForces[c_gsystem.widthTot*ref.w + dihedral.w] = f3;
	}
}

__global__ void summDihedralForces_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 f = c_gsystem.d_forces[d_i];
		float4 df;
		int i;
		for(i = 0; i < c_dihedralData.d_dihedralCount[d_i]; i++){
			df = c_dihedralData.d_dihedralForces[c_gsystem.widthTot*i + d_i];
			f.x += df.x;
			f.y += df.y;
			f.z += df.z;
		}
		c_gsystem.d_forces[d_i] = f;
	}
}

inline void compute(){
	dihedralPotentialCHARMM_kernel<<<dihedralBlockCount, dihedralBlockSize>>>();
	cudaThreadSynchronize();
	summDihedralForces_kernel<<<dihedralSummBlockCount, dihedralSummBlockSize>>>();

	/*cudaMemcpy(gsystem.h_forces, gsystem.d_forces, atomCount*sizeof(float4), cudaMemcpyDeviceToHost);
	int i;
	float3 force = make_float3(0.0f, 0.0f, 0.0f);
	for(i = 0; i < atomCount; i++){
		force.x += gsystem.h_forces[i].x;
		force.y += gsystem.h_forces[i].y;
		force.z += gsystem.h_forces[i].z;
	}
	printf("Net force (dihedrals): (%f, %f, %f) %f\n", force.x, force.y, force.z,
			sqrtf(force.x*force.x + force.y*force.y + force.z*force.z));*/
}

__global__ void dihedralPotentialEnergy_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_dihedralData.Dtot){

		int4 dihedral = c_dihedralData.d_dihedrals[d_i];
		//int4 ref = c_dihedralData.d_dihedralRefs[d_i];
		GDihedralParameters par = c_dihedralData.d_dihedralParameters[d_i];

		float4 r1 = tex1Dfetch(t_coord, dihedral.x);
		float4 r2 = tex1Dfetch(t_coord, dihedral.y);
		float4 r3 = tex1Dfetch(t_coord, dihedral.z);
		float4 r4 = tex1Dfetch(t_coord, dihedral.w);

		float3 dr12, dr23, dr34;

		dr12.x = r1.x - r2.x;
		dr12.y = r1.y - r2.y;
		dr12.z = r1.z - r2.z;
		DO_PBC(dr12);

		dr23.x = r2.x - r3.x;
		dr23.y = r2.y - r3.y;
		dr23.z = r2.z - r3.z;
		DO_PBC(dr23);
		float r232 = dr23.x*dr23.x + dr23.y*dr23.y + dr23.z*dr23.z;

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

		float coschi = (a.x*b.x + a.y*b.y + a.z*b.z)*a.w*b.w;
		float sinchi = -sqrtf(r232)*a.w*b.w*(a.x*dr34.x + a.y*dr34.y + a.z*dr34.z);//(c.x*b.x + c.y*b.y + c.z*b.z)*c.w*b.w;

		float chi = -atan2(sinchi, coschi);

		float pot = 0.0f;
		pot += par.kchi*(1.0f + cosf(par.n*chi - par.delta));
		while(par.multiplicityRef != -1){
			par = c_dihedralData.d_multiplicityParameters[par.multiplicityRef];
			//if(par.n > 0.0f){
			pot += par.kchi*(1.0f + cosf(par.n*chi - par.delta));
			/*} else {
				float diff = chi - par.delta;
				if(diff < -M_PI){
					diff += 2.0f*M_PI;
				} else
				if(diff > M_PI){
					diff -= 2.0f*M_PI;
				}
				pot += par.kchi*diff*diff;
			}*/
			/*if(par.multiplicityRef != -1){
				par = c_dihedralData.d_multiplicityParameters[par.multiplicityRef];
			}*/
		}
		c_dihedralData.d_dihedralEnergies[d_i] = pot;
	}
}

inline void computeEnergy(){
	dihedralPotentialEnergy_kernel<<<dihedralBlockCount, dihedralBlockSize>>>();
	cudaMemcpy(dihedralData.h_dihedralEnergies, dihedralData.d_dihedralEnergies,
				dihedralData.Dtot*sizeof(float), cudaMemcpyDeviceToHost);
	int i, traj;
	for(traj = 0; traj < parameters.Ntr; traj++){
		float pot = 0.0f;
		for(i = 0; i < dihedralData.D; i++){
			pot += dihedralData.h_dihedralEnergies[i + traj*dihedralData.D];
		}
		energyOutput.values[traj] = pot;
	}
	checkCUDAError("dihedral energy");

}

void destroy(){

}

#undef LOG

} //namespace dihedral_potential
