/*
 * RepulsiveBoundaryPotential.cu
 *
 *  Created on: Feb 28, 2011
 *      Author: zhmurov
 */

#include "../Core/global.h"
#include "../Core/md.cuh"
#include "RepulsiveBoundaryPotential.cuh"

namespace repulsive_boundary {

class Log: public ILog {
	virtual void Write(const char* message) const {
		std::cout << makeTimePrefix() << "<repulsive_boundary> " << message << std::endl;
	}
} log;

#define LOG LogStream(log)

void create(){
	if (!getYesNoParameter(PARAMETER_REPULSIVE_BOUNDARY_ON, 0))
		return;
	if (getYesNoParameter(PARAMETER_REPULSIVE_BOUNDARY_FOR_MD, 1)) {
		LOG << "Repulsive boundary for MD enabled";
		potential.compute = &compute;
		potential.destroy = &destroy;
		sprintf(potential.name, "Repulsive boundary potential");
		potentials[potentialsCount] = &potential;
		potentialsCount ++;
	}
	init();
}

void init(){
	LOG << "Creating repulsive boundary...";
	repulsiveBoundaryBlockSize = BLOCK_SIZE;
	repulsiveBoundaryBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;
	float r = getFloatParameter(PARAMETER_REPULSIVE_BOUNDARY_RADIUS);
	repulsiveBoundaryData.geometry.w = r;
	getVectorParameter(PARAMETER_REPULSIVE_BOUNDARY_LOCATION,
			&repulsiveBoundaryData.geometry.x,
			&repulsiveBoundaryData.geometry.y,
			&repulsiveBoundaryData.geometry.z);
	repulsiveBoundaryData.epsilon = getFloatParameter(PARAMETER_REPULSIVE_BOUNDARY_EPSILON);
	repulsiveBoundaryData.sigma = getFloatParameter(PARAMETER_REPULSIVE_BOUNDARY_SIGMA);
	repulsiveBoundaryData.sixEpsilonSigma6 = repulsiveBoundaryData.sigma*repulsiveBoundaryData.sigma;
	repulsiveBoundaryData.sixEpsilonSigma6 =
			repulsiveBoundaryData.sixEpsilonSigma6*
			repulsiveBoundaryData.sixEpsilonSigma6*
			repulsiveBoundaryData.sixEpsilonSigma6;
	repulsiveBoundaryData.sixEpsilonSigma6 *= 6.0f*repulsiveBoundaryData.epsilon;
	cudaMemcpyToSymbol(c_repulsiveBoundaryData, &repulsiveBoundaryData,
					sizeof(RepulsiveBoundary), 0, cudaMemcpyHostToDevice);

	char filename[256];
	getMaskedParameter(filename, PARAMETER_REPULSIVE_BOUNDARY_PDB, PARAMETER_STRING_UNDEFINED);
	if(strcmp(filename, PARAMETER_STRING_UNDEFINED) != 0){
		LOG << "Creating PDB file with sphere representation";
		PDB pdbSphere;
		pdbSphere.ssCount = 0;
		int size = 36;
		pdbSphere.atomCount = size*size;
		pdbSphere.atoms = (PDBAtom*)calloc(pdbSphere.atomCount, sizeof(PDBAtom));
		pdbSphere.connections.connectCount = (int*)calloc(pdbSphere.atomCount, sizeof(int));
		pdbSphere.connections.connectMap = (int*)calloc(pdbSphere.atomCount*size, sizeof(int));
		float phi, theta;
		int i = 0;
		int j, k;
		for(j = 0; j < size; j++){
			for(k = 0; k < size; k++){
				i = j*size + k;
				/*if(j < size - 2){
					int jp0 = j*size + k;
					int jp1 = (j+1)*size + k;
					int jp2 = (j+2)*size + k;
					int count = pdbSphere.connections.connectCount[jp1];
					pdbSphere.connections.connectMap[size*size*count + jp1] = jp0;
					count ++;
					pdbSphere.connections.connectMap[size*size*count + jp1] = jp2;
					count ++;
					pdbSphere.connections.connectCount[jp1] = count;
				}*/
				//if(k < size - 2){
				int count = pdbSphere.connections.connectCount[i];
				if(k > 0){
					pdbSphere.connections.connectMap[size*size*count + i] = i-1;
					count ++;
				}
				if(k < size - 1){
					pdbSphere.connections.connectMap[size*size*count + i] = i+1;
					count ++;
				}
				if(j > 0){
					pdbSphere.connections.connectMap[size*size*count + i] = i - size;
					count ++;
				}
				if(j < size - 1){
					pdbSphere.connections.connectMap[size*size*count + i] = i + size;
					count ++;
				}
				pdbSphere.connections.connectCount[i] = count;
				//}
				phi = 2.0f*M_PI*((float)j)/((float)size)+0.000001f;
				theta = 2.0f*M_PI*((float)k)/((float)size)+0.000001f;
				pdbSphere.atoms[i].id = i+1;
				pdbSphere.atoms[i].resid = i+1;
				sprintf(pdbSphere.atoms[i].name, "CA");
				pdbSphere.atoms[i].altLoc = ' ';
				pdbSphere.atoms[i].chain = 'X';
				sprintf(pdbSphere.atoms[i].resName, "BON");
				pdbSphere.atoms[i].x = 10.0f*sinf(phi)*sinf(theta)*r;
				pdbSphere.atoms[i].y = 10.0f*sinf(phi)*cosf(theta)*r;
				pdbSphere.atoms[i].z = 10.0f*cosf(phi)*r;
			}
		}
		//printf("%d\n", i);
		writePDB(filename, &pdbSphere, 1);
	}
	LOG << "Done creating repulsive boundary.";
}

__global__ void repulsiveBoundarySpherePotential_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsystem.Ntot){
		float4 coord = tex1Dfetch(t_coord, d_i);
		float4 sphere = c_repulsiveBoundaryData.geometry;
		float4 f = c_gsystem.d_forces[d_i];
		float mult;
		coord.x -= sphere.x;
		coord.y -= sphere.y;
		coord.z -= sphere.z;
		coord.w = sqrtf(coord.x*coord.x + coord.y*coord.y + coord.z*coord.z);
		sphere.w -= coord.w;
		sphere.w -= c_repulsiveBoundaryData.sigma;
		if(sphere.w < 0){
			/*sphere.w = 1.0f/sphere.w;
			mult = sphere.w*sphere.w;
			mult = mult*mult*mult*sphere.w;
			mult *= c_repulsiveBoundaryData.sixEpsilonSigma6;
			mult /= coord.w;*/
			mult = 2.0f*c_repulsiveBoundaryData.epsilon*sphere.w/coord.w;
			f.x += mult*coord.x;
			f.y += mult*coord.y;
			f.z += mult*coord.z;
			c_gsystem.d_forces[d_i] = f;
		}

	}
}

inline void compute(){
	repulsiveBoundarySpherePotential_kernel<<<repulsiveBoundaryBlockCount, repulsiveBoundaryBlockSize>>>();
}

void destroy(){

}

#undef LOG

} // namespace repulsive_boundary

