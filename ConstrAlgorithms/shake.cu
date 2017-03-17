#include "../Core/global.h"
#include "../Core/md.cuh"
#include "../Util/Log.h"
#include "shake.cuh"

namespace shake_constrAlg {

class Log: public ILog {
    virtual void Write(const char* message) const {
        std::cout << makeTimePrefix() << "<shake_constrAlg> " << message << std::endl;
    }
} log;

#define LOG LogStream(log)

void create(){
	constrAlg.compute 	= &compute;
	constrAlg.destroy 	= &destroy;
	sprintf(constrAlg.name, "SHAKE constraint algorithm");
//	constrAlgs[constrAlgsCount] = &constrAlg;
//	constrAlgsCount ++;
	init();
	constrCount+=constrAlg.Nconstr;
	constrAlgs[constrAlgsCount] = &constrAlg;
	constrAlgsCount ++;
}


void assignHeavyAndHydr(int* shakeBonds, int k, int* heavyId, int* hydrId){
	if(topology.atoms[topology.bonds[shakeBonds[k]].i].type[0]=='H'){
		*hydrId=topology.bonds[shakeBonds[k]].i;
		*heavyId=topology.bonds[shakeBonds[k]].j;
	}
	else{
		*hydrId=topology.bonds[shakeBonds[k]].j;
		*heavyId=topology.bonds[shakeBonds[k]].i;
	}
}

void init(){
	shakeConstrData.tol = getFloatParameter(PARAMETER_RIGIDTOL, 0.00001f);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int shakeBondsCount=0;
	for(int k=0; k<topology.bondCount;k++){
		if ( (topology.atoms[topology.bonds[k].i].type[0]=='H')||(topology.atoms[topology.bonds[k].j].type[0]=='H') ){
			shakeBondsCount++;
//				topology.bonds[k].kb=0.0f;
		}
	}
	constrAlg.Nconstr=shakeBondsCount;

	shakeConstrData.h_shakeClusters 		= (  int*)calloc(shakeBondsCount, sizeof(  int)*4);//allocate more than we actually need
	shakeConstrData.h_shakeClustersParam 	= (float*)calloc(shakeBondsCount, sizeof(float)*4);

	int * shakeBonds = (int*)calloc(shakeBondsCount, sizeof(int));
	shakeBondsCount=0;
	for (int k=0; k<topology.bondCount; k++){
		if( (topology.atoms[topology.bonds[k].i].type[0]=='H')||(topology.atoms[topology.bonds[k].j].type[0]=='H') ){
			shakeBonds[shakeBondsCount]=k;
			shakeBondsCount++;
		}
	}

	shakeConstrData.shakeClustersCount=0;
	int heavyId, hydrId;
	int flag=0;
	for (int k=0; k< shakeBondsCount; k++){
		assignHeavyAndHydr(shakeBonds, k, &heavyId, &hydrId);
		flag=0;
		for(int l=0; l< k; l++){
			if(heavyId==shakeConstrData.h_shakeClusters[4*l+0]){
				if(shakeConstrData.h_shakeClusters[4*l+2]==-1){
					shakeConstrData.h_shakeClusters[4*l+2] = hydrId;
					flag=1;
					break;
				}
				else if(shakeConstrData.h_shakeClusters[4*l+3]==-1){
					shakeConstrData.h_shakeClusters[4*l+3] = hydrId;
					flag=1;
					break;
				}
				else{
					die("More than 3 hydrogens connected to one atom! Exiting...\n");
				}
			}
		}
		if (flag==0){
			shakeConstrData.h_shakeClusters[4*shakeConstrData.shakeClustersCount+0]=heavyId;
			shakeConstrData.h_shakeClusters[4*shakeConstrData.shakeClustersCount+1]=hydrId;
			shakeConstrData.h_shakeClusters[4*shakeConstrData.shakeClustersCount+2]=-1;
			shakeConstrData.h_shakeClusters[4*shakeConstrData.shakeClustersCount+3]=-1;
			shakeConstrData.h_shakeClustersParam[4*shakeConstrData.shakeClustersCount+0]=topology.atoms[heavyId].mass;
			shakeConstrData.h_shakeClustersParam[4*shakeConstrData.shakeClustersCount+1]=topology.atoms[hydrId].mass;
			shakeConstrData.h_shakeClustersParam[4*shakeConstrData.shakeClustersCount+2]=1.0f/topology.atoms[heavyId].mass + 1.0f/topology.atoms[hydrId].mass;
			shakeConstrData.h_shakeClustersParam[4*shakeConstrData.shakeClustersCount+3]=topology.bonds[shakeBonds[k]].b0;
				shakeConstrData.shakeClustersCount++;
		}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	shakeConstrData.shakeClustersTot=shakeConstrData.shakeClustersCount*parameters.Ntr;

	shakeConstrData.h_shakeClusters 		= (  int*)realloc(shakeConstrData.h_shakeClusters, parameters.Ntr*shakeConstrData.shakeClustersCount*sizeof(  int)*4);
	shakeConstrData.h_shakeClustersParam 	= (float*)realloc(shakeConstrData.h_shakeClustersParam, parameters.Ntr*shakeConstrData.shakeClustersCount*sizeof(float)*4);

	allocateGPU((void**)&shakeConstrData.d_shakeClusters, shakeConstrData.shakeClustersCount*parameters.Ntr*sizeof(int4));
	allocateGPU((void**)&shakeConstrData.d_shakeClustersParam, shakeConstrData.shakeClustersCount*parameters.Ntr*sizeof(float4));
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (int traj=1; traj < parameters.Ntr; traj++){
		for (int i = 0; i< shakeConstrData.shakeClustersCount; i++){
			for(int j = 0; j< 4; j++){
				if (shakeConstrData.h_shakeClusters[i*4+j] != -1){
					shakeConstrData.h_shakeClusters[(traj*shakeConstrData.shakeClustersCount + i)*4 + j]=shakeConstrData.h_shakeClusters[i*4+j] + traj*gsystem.N;
					shakeConstrData.h_shakeClustersParam[(traj*shakeConstrData.shakeClustersCount + i)*4 + j]=shakeConstrData.h_shakeClustersParam[i*4+j];
				}
				else{
					shakeConstrData.h_shakeClusters[(traj*shakeConstrData.shakeClustersCount + i)*4 + j] = -1;
					shakeConstrData.h_shakeClustersParam[(traj*shakeConstrData.shakeClustersCount + i)*4 + j]=shakeConstrData.h_shakeClustersParam[i*4+j];
				}
			}
		}
	}
	cudaMemcpy(shakeConstrData.d_shakeClusters,      shakeConstrData.h_shakeClusters,      shakeConstrData.shakeClustersCount*parameters.Ntr*sizeof(int)*4,   cudaMemcpyHostToDevice);
	cudaMemcpy(shakeConstrData.d_shakeClustersParam, shakeConstrData.h_shakeClustersParam, shakeConstrData.shakeClustersCount*parameters.Ntr*sizeof(float)*4, cudaMemcpyHostToDevice);


	shakeBlockCount = gsystem.Ntot/BLOCK_SIZE + 1;
	shakeBlockSize = BLOCK_SIZE;

	cudaMemcpyToSymbol(c_shakeConstrData, &shakeConstrData,
				sizeof(ShakeConstrData), 0, cudaMemcpyHostToDevice);
		LOG << "Done initializing shake algorithm.";
}




__global__ void shakeKernel(){
	int index=threadIdx.x+blockIdx.x*blockDim.x;
	if (index < c_shakeConstrData.shakeClustersTot){
		int4 cluster = c_shakeConstrData.d_shakeClusters[index];
		float4 param = c_shakeConstrData.d_shakeClustersParam[index];

		float mass_j    = param.x;
		float mass_i    = param.y;
		float invMasses = param.z;
		float d         = param.w;
		float d2		= d*d;

		float4 old_coord_j  = tex1Dfetch(t_coord, cluster.x);
		float4 old_coord_i0 = tex1Dfetch(t_coord, cluster.y);
		float4 old_coord_i1;
		float4 old_coord_i2;

		float4 delta_coord_j  = c_gsystem.d_midcoord[cluster.x];
		float4 delta_coord_i0 = c_gsystem.d_midcoord[cluster.y];
		float4 delta_coord_i1;
		float4 delta_coord_i2;


		if(cluster.z!=-1){
			old_coord_i1		= tex1Dfetch(t_coord, cluster.z);
			delta_coord_i1		= c_gsystem.d_midcoord[cluster.z];
		}
		if(cluster.w!=-1){
			old_coord_i2 		= tex1Dfetch(t_coord, cluster.w);
			delta_coord_i2		= c_gsystem.d_midcoord[cluster.w];
		}



        bool converged = false;
        int iteration = 0;
        while (iteration < 15 && !converged){
        	converged = true;
        	float4 oldDiff;
        	oldDiff.x = old_coord_i0.x - old_coord_j.x;
        	oldDiff.y = old_coord_i0.y - old_coord_j.y;
        	oldDiff.z = old_coord_i0.z - old_coord_j.z;
        	oldDiff.w = oldDiff.x*oldDiff.x + oldDiff.y*oldDiff.y + oldDiff.z*oldDiff.z;
        	float4 deltaDiff;
        	deltaDiff.x = delta_coord_i0.x - delta_coord_j.x;
        	deltaDiff.y = delta_coord_i0.y - delta_coord_j.y;
        	deltaDiff.z = delta_coord_i0.z - delta_coord_j.z;
        	deltaDiff.w = deltaDiff.x*deltaDiff.x + deltaDiff.y*deltaDiff.y + deltaDiff.z*deltaDiff.z;
        	float scalMult = oldDiff.x*deltaDiff.x + oldDiff.y*deltaDiff.y + oldDiff.z*deltaDiff.z;
        	float sigma = oldDiff.w + 2*scalMult + deltaDiff.w - d2;

        	float mult = 0.5f/invMasses;
        	if(fabs(sigma) > 2.0f*c_shakeConstrData.tol*d2){
        		converged = false;
        		sigma *= mult/(scalMult+oldDiff.w);
        		delta_coord_i0.x += -oldDiff.x*sigma/mass_i;
        		delta_coord_i0.y += -oldDiff.y*sigma/mass_i;
        		delta_coord_i0.z += -oldDiff.z*sigma/mass_i;

        		delta_coord_j.x  +=  oldDiff.x*sigma/mass_j;
        		delta_coord_j.y  +=  oldDiff.y*sigma/mass_j;
        		delta_coord_j.z  +=  oldDiff.z*sigma/mass_j;
        	}

        	if (cluster.z != -1){
        		oldDiff.x = old_coord_i1.x - old_coord_j.x;
        		oldDiff.y = old_coord_i1.y - old_coord_j.y;
        		oldDiff.z = old_coord_i1.z - old_coord_j.z;
        		oldDiff.w = oldDiff.x*oldDiff.x + oldDiff.y*oldDiff.y + oldDiff.z*oldDiff.z;

        		deltaDiff.x = delta_coord_i1.x - delta_coord_j.x;
        		deltaDiff.y = delta_coord_i1.y - delta_coord_j.y;
        		deltaDiff.z = delta_coord_i1.z - delta_coord_j.z;
        		deltaDiff.w = deltaDiff.x*deltaDiff.x + deltaDiff.y*deltaDiff.y + deltaDiff.z*deltaDiff.z;
        		scalMult = oldDiff.x*deltaDiff.x + oldDiff.y*deltaDiff.y + oldDiff.z*deltaDiff.z;
        		sigma = oldDiff.w + 2*scalMult + deltaDiff.w - d2;


        		if (fabs(sigma) > 2.0f*c_shakeConstrData.tol*d2){
        			converged = false;
        			sigma *= mult/(scalMult+oldDiff.w);
        			delta_coord_i1.x += -oldDiff.x*sigma/mass_i;
        			delta_coord_i1.y += -oldDiff.y*sigma/mass_i;
        			delta_coord_i1.z += -oldDiff.z*sigma/mass_i;

  	        		delta_coord_j.x  +=  oldDiff.x*sigma/mass_j;
  	        		delta_coord_j.y  +=  oldDiff.y*sigma/mass_j;
        			delta_coord_j.z  +=  oldDiff.z*sigma/mass_j;
        		}
        	}

        	if (cluster.w != -1){
        		oldDiff.x = old_coord_i2.x - old_coord_j.x;
        		oldDiff.y = old_coord_i2.y - old_coord_j.y;
        		oldDiff.z = old_coord_i2.z - old_coord_j.z;
        		oldDiff.w = oldDiff.x*oldDiff.x + oldDiff.y*oldDiff.y + oldDiff.z*oldDiff.z;

        		deltaDiff.x = delta_coord_i2.x - delta_coord_j.x;
        		deltaDiff.y = delta_coord_i2.y - delta_coord_j.y;
        		deltaDiff.z = delta_coord_i2.z - delta_coord_j.z;
        		deltaDiff.w = deltaDiff.x*deltaDiff.x + deltaDiff.y*deltaDiff.y + deltaDiff.z*deltaDiff.z;
        		scalMult = oldDiff.x*deltaDiff.x + oldDiff.y*deltaDiff.y + oldDiff.z*deltaDiff.z;
        		sigma = oldDiff.w + 2*scalMult + deltaDiff.w - d2;


        		if (fabs(sigma) > 2.0f*c_shakeConstrData.tol*d2){
        			converged = false;
        			sigma *= mult/(scalMult+oldDiff.w);
        			delta_coord_i2.x += -oldDiff.x*sigma/mass_i;
        			delta_coord_i2.y += -oldDiff.y*sigma/mass_i;
        			delta_coord_i2.z += -oldDiff.z*sigma/mass_i;

  	        		delta_coord_j.x  +=  oldDiff.x*sigma/mass_j;
  	        		delta_coord_j.y  +=  oldDiff.y*sigma/mass_j;
        			delta_coord_j.z  +=  oldDiff.z*sigma/mass_j;
        		}
        	}
            iteration++;
        }
        if(iteration > 10){
       //	printf("unconverged block\n");
        }
        c_gsystem.d_midcoord[cluster.x]=delta_coord_j;
        c_gsystem.d_midcoord[cluster.y]=delta_coord_i0;

        if(param.z != -1){
        	c_gsystem.d_midcoord[cluster.z]=delta_coord_i1;
        }
        if(param.w != -1){
        	c_gsystem.d_midcoord[cluster.w]=delta_coord_i2;

        }
	}
}

/*
__global__ void empty_kernel(){
	int index=threadIdx.x+blockIdx.x*blockDim.x;
}
*/

inline void compute(){
//	printf("NumSHAKEclusters = %d\n", shakeConstrData.shakeClustersCount);
	if ( shakeConstrData.shakeClustersCount > 0){
		shakeKernel<<<shakeBlockCount, BLOCK_SIZE>>>();
		/*
		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord(start, 0);
		empty_kernel<<<shakeBlockCount, BLOCK_SIZE>>>();
	    cudaEventRecord(stop, 0);
	    cudaEventSynchronize(stop);
	    float time;
	    cudaEventElapsedTime(&time, start, stop);
	    printf ("Time for the kernel: %f ms\n", time);
	    */
	}


}

void destroy(){

}

#undef LOG

} // namespace shake_constrAlg
