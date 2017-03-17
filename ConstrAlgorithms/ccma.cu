#include "../Core/global.h"
#include "../Core/md.cuh"
#include "../Util/Log.h"
#include "./quern/quern.h"
#include "ccma.cuh"


namespace ccma_constrAlg {

class Log: public ILog {
    virtual void Write(const char* message) const {
        std::cout << makeTimePrefix() << "<ccma_constrAlg> " << message << std::endl;
    }
} log;

#define LOG LogStream(log)

void create(){
	constrAlg.compute 	= &compute;
	constrAlg.destroy 	= &destroy;
	sprintf(constrAlg.name, "CCMA constraint algorithm");
	init();
	constrCount+=constrAlg.Nconstr;
	constrAlgs[constrAlgsCount] = &constrAlg;
	constrAlgsCount ++;
}

float findAngle(int i, int j){
	for(int k=0; k< topology.angleCount; k++){
		if ( ((i==topology.angles[k].i)&&(j==topology.angles[k].k)) || ((j==topology.angles[k].i)&&(i==topology.angles[k].k))  ){
			return topology.angles[k].theta0;
		}
	}
	die("couldn't find an appropriate angle for two ccma bonds");
}

void init(){
	ccmaConstrData.tol = getFloatParameter(PARAMETER_RIGIDTOL, 0.00001f);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	ccmaConstrData.ccmaAtomsCount = 0;
	ccmaConstrData.ccmaConstrCount= 0;//!to be used this way only with implicit solvent
//get number of bonds handled by ccma
	for (int i=0; i< topology.bondCount; i++){
		if (true){//check if settle bond
			ccmaConstrData.ccmaConstrCount++;
		}
	}
	constrAlg.Nconstr=ccmaConstrData.ccmaConstrCount;
	ccmaConstrData.ccmaConstrTot =ccmaConstrData.ccmaConstrCount*parameters.Ntr;
//

	allocateCPU((void**)&ccmaConstrData.h_atomPair, ccmaConstrData.ccmaConstrTot*sizeof(int2));
	allocateGPU((void**)&ccmaConstrData.d_atomPair, ccmaConstrData.ccmaConstrTot*sizeof(int2));

	allocateCPU((void**)&ccmaConstrData.h_M, ccmaConstrData.ccmaConstrTot*sizeof(float));
	allocateGPU((void**)&ccmaConstrData.d_M, ccmaConstrData.ccmaConstrTot*sizeof(float));

	allocateCPU((void**)&ccmaConstrData.h_constrDir, ccmaConstrData.ccmaConstrTot*sizeof(float4));
	allocateGPU((void**)&ccmaConstrData.d_constrDir, ccmaConstrData.ccmaConstrTot*sizeof(float4));

	//allocate some more gpu memory
	allocateGPU((void**)&ccmaConstrData.d_constrSigma, ccmaConstrData.ccmaConstrTot*sizeof(float));
	allocateGPU((void**)&ccmaConstrData.d_constrLambda, ccmaConstrData.ccmaConstrTot*sizeof(float));

	int index=0;
	for (int i=0; i< topology.bondCount; i++){
		if (true){//check if settle bond
			ccmaConstrData.h_atomPair[index].x=topology.bonds[i].i;
			ccmaConstrData.h_atomPair[index].y=topology.bonds[i].j;
			ccmaConstrData.h_M[index]=1.0f/topology.atoms[ccmaConstrData.h_atomPair[index].x].mass+1.0f/topology.atoms[ccmaConstrData.h_atomPair[index].y].mass;
			ccmaConstrData.h_constrDir[index].w=topology.bonds[i].b0;
			index++;
		}
	}

//currently all the atoms are considered to be ccma atoms
//some of them just don't have any ccma bonds attached
//this is ok while we have no water, but may become a problem afterwards
	ccmaConstrData.ccmaAtomsCount	= topology.atomCount;
	ccmaConstrData.ccmaAtomsTot	= ccmaConstrData.ccmaAtomsCount*parameters.Ntr;

	//ccmaConstrData.h_attachedBonds 	= (int*)calloc( ccmaConstrData.ccmaAtomsCount ,sizeof(int)*4);//calloc is used to zero the memory
	allocateCPU((void**)&ccmaConstrData.h_attachedBonds, ccmaConstrData.ccmaAtomsTot*sizeof(int)*4);
	allocateGPU((void**)&ccmaConstrData.d_attachedBonds, ccmaConstrData.ccmaAtomsTot*sizeof(int)*4);

	for(int i=0; i<ccmaConstrData.ccmaAtomsTot; i++){
		for(int j=0; j<4; j++)
		ccmaConstrData.h_attachedBonds[i*4+j]=0;
	}


	int * nAttachedConstr			= (int*)calloc( ccmaConstrData.ccmaAtomsCount ,sizeof(int));
	for(int i=0; i<ccmaConstrData.ccmaConstrCount; i++){
		ccmaConstrData.h_attachedBonds[4*ccmaConstrData.h_atomPair[i].x+nAttachedConstr[ccmaConstrData.h_atomPair[i].x]]=(i+1);//the sign of constraint index means the position of atom in it (i or j),
		ccmaConstrData.h_attachedBonds[4*ccmaConstrData.h_atomPair[i].y+nAttachedConstr[ccmaConstrData.h_atomPair[i].y]]=-(i+1);//so we add 1 to each to deal with the zero constraint.
		nAttachedConstr[ccmaConstrData.h_atomPair[i].x]++;
		nAttachedConstr[ccmaConstrData.h_atomPair[i].y]++;
	}


	//copy everything that's left on GPU
	for(int i=0; i< ccmaConstrData.ccmaConstrCount; i++){
		for(int j=1; j< parameters.Ntr; j++){
			ccmaConstrData.h_atomPair[i+j*ccmaConstrData.ccmaConstrCount].x=ccmaConstrData.h_atomPair[i].x+j*gsystem.N;
			ccmaConstrData.h_atomPair[i+j*ccmaConstrData.ccmaConstrCount].y=ccmaConstrData.h_atomPair[i].y+j*gsystem.N;
			ccmaConstrData.h_M[i+j*ccmaConstrData.ccmaConstrCount]=ccmaConstrData.h_M[i];
			ccmaConstrData.h_constrDir[i+j*ccmaConstrData.ccmaConstrCount].w=ccmaConstrData.h_constrDir[i].w;
		}
	}

	for(int i=0; i< 4*ccmaConstrData.ccmaAtomsCount; i++){
		for(int j=1; j< parameters.Ntr; j++){
			if(ccmaConstrData.h_attachedBonds[i]!=0){
				float sign = (ccmaConstrData.h_attachedBonds[i]>0) ? 1.0 : -1.0;
				ccmaConstrData.h_attachedBonds[i+j*4*ccmaConstrData.ccmaAtomsCount]=ccmaConstrData.h_attachedBonds[i]+sign*j*ccmaConstrData.ccmaConstrCount;
			}
		}
	}

	cudaMemcpy(ccmaConstrData.d_M, ccmaConstrData.h_M,ccmaConstrData.ccmaConstrTot*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy(ccmaConstrData.d_atomPair, ccmaConstrData.h_atomPair,ccmaConstrData.ccmaConstrTot*sizeof(int2), cudaMemcpyHostToDevice );
	cudaMemcpy(ccmaConstrData.d_constrDir, ccmaConstrData.h_constrDir,ccmaConstrData.ccmaConstrTot*sizeof(float4), cudaMemcpyHostToDevice );
	cudaMemcpy(ccmaConstrData.d_attachedBonds, ccmaConstrData.h_attachedBonds,ccmaConstrData.ccmaAtomsTot*4*sizeof(int), cudaMemcpyHostToDevice );


//fill the coupling matrix
	int* matrixRowStart = (int*)calloc(ccmaConstrData.ccmaConstrCount+1, sizeof(int));
	int* matrixColIndex = (int*)calloc(ccmaConstrData.ccmaConstrCount*6, sizeof(int));
	double* matrixValue	= (double*)calloc(ccmaConstrData.ccmaConstrCount*6, sizeof(double));

	int counter=0;
	float invmasses;
	float mi, mj, angle;
//	printf("constrCount=%d\n", ccmaConstrData.ccmaConstrCount);
	for(int k=0; k< ccmaConstrData.ccmaConstrCount; k++){
		matrixRowStart[k]=counter;
		//printf("matrixRowStart[%d] = %d\n", k, matrixRowStart[k]);
		mi = topology.atoms[ccmaConstrData.h_atomPair[k].x].mass;
		mj = topology.atoms[ccmaConstrData.h_atomPair[k].y].mass;
		invmasses=1.0f/mi+1.0f/mj;
		//ccmaConstrData.h_M[k]=invmasses;
		//printf("mi = %f, mj=%f, invm=%f", mi, mj, invmasses);
		for(int l=0; l< ccmaConstrData.ccmaConstrCount; l++){
			if(k==l){
				matrixValue[counter]=1.0f;
				matrixColIndex[counter]=l;
				counter++;
			}
			else{
				if(ccmaConstrData.h_atomPair[l].x==ccmaConstrData.h_atomPair[k].x){
					angle = findAngle(ccmaConstrData.h_atomPair[l].y, ccmaConstrData.h_atomPair[k].y);
					//printf("angle=%f\n", angle);
					matrixValue[counter]=(double)cos(angle)/(mi*invmasses);
					matrixColIndex[counter]=l;
					counter++;
				}
				else if(ccmaConstrData.h_atomPair[l].y==ccmaConstrData.h_atomPair[k].x){
					angle = findAngle(ccmaConstrData.h_atomPair[l].x, ccmaConstrData.h_atomPair[k].y);
					matrixValue[counter]=(double)cos(angle)/(mi*invmasses);
					matrixColIndex[counter]=l;
					counter++;
				}
				else if(ccmaConstrData.h_atomPair[l].x==ccmaConstrData.h_atomPair[k].y){
					angle = findAngle(ccmaConstrData.h_atomPair[l].y, ccmaConstrData.h_atomPair[k].x);
					matrixValue[counter]=(double)cos(angle)/(mj*invmasses);
					matrixColIndex[counter]=l;
					counter++;
				}
				else if(ccmaConstrData.h_atomPair[l].y==ccmaConstrData.h_atomPair[k].y){
					angle = findAngle(ccmaConstrData.h_atomPair[l].x, ccmaConstrData.h_atomPair[k].x);
					matrixValue[counter]=(double)cos(angle)/(mj*invmasses);
					matrixColIndex[counter]=l;
					counter++;
				}
			}
		}

	}
	matrixRowStart[ccmaConstrData.ccmaConstrCount]=counter;


//invert matrix using QR-decomposition
	int *qRowStart, *qColIndex, *rRowStart, *rColIndex;
	double *qValue, *rValue;
	int result = QUERN_compute_qr(ccmaConstrData.ccmaConstrCount, ccmaConstrData.ccmaConstrCount, &matrixRowStart[0], &matrixColIndex[0], &matrixValue[0], NULL,
	                &qRowStart, &qColIndex, &qValue, &rRowStart, &rColIndex, &rValue);
	double* rhs = (double*)calloc(ccmaConstrData.ccmaConstrCount, sizeof(double));
	double* res = (double*)calloc(ccmaConstrData.ccmaConstrCount, sizeof(double));

	free(matrixValue);
	free(matrixColIndex);

	float** invMval = (float**)calloc(ccmaConstrData.ccmaConstrCount, sizeof(float*));
	int** 	invMcol = (int**)calloc(ccmaConstrData.ccmaConstrCount, sizeof(int*));
	int* 	Nused= (int*)calloc(ccmaConstrData.ccmaConstrCount, sizeof(int));
	for(int i=0; i< ccmaConstrData.ccmaConstrCount; i++){
		invMval[i] = (float*)calloc(32, sizeof(float));
		invMcol[i] = (int*)calloc(32, sizeof(int));
	}

	for (int i = 0; i < ccmaConstrData.ccmaConstrCount; i++) {
	// Extract column i of the inverse matrix.
	//	printf("did it, %d\n",i);
		for (int j = 0; j < ccmaConstrData.ccmaConstrCount; j++)
			rhs[j] = (i == j ? 1.0 : 0.0);
	    result = QUERN_multiply_with_q_transpose(ccmaConstrData.ccmaConstrCount, qRowStart, qColIndex, qValue, rhs);
	    result = QUERN_solve_with_r(ccmaConstrData.ccmaConstrCount, rRowStart, rColIndex, rValue, rhs, res);

	   	for (int j = 0; j < ccmaConstrData.ccmaConstrCount; j++) {
	   		double value = res[j];//*ccmaConstrData.h_constrDir[i].w/ccmaConstrData.h_constrDir[j].w;
	   		if (abs(value) > 0.01f){
	   			invMval[j][Nused[j]]=(float)value;
	   			invMcol[j][Nused[j]]=i;
	   			Nused[j]++;
	   			if(Nused[j]%32==0){
	   				printf("ups\n");
	   				invMval[j]=(float*)realloc(invMval[j],(Nused[j]+32)*sizeof(float));
	   				invMcol[j]=(int*)realloc(invMcol[j],(Nused[j]+32)*sizeof(int));
	   			}
	   		}
	   	}
	}
	QUERN_free_result(qRowStart, qColIndex, qValue);
	QUERN_free_result(rRowStart, rColIndex, rValue);




	//int maxRowElements = 0;
	for (unsigned i = 0; i < ccmaConstrData.ccmaConstrCount; i++)
		maxRowElements = (Nused[i]>maxRowElements) ? Nused[i] : maxRowElements;
//	maxRowElements++;//why?
	if(maxRowElements%4!=0)
		maxRowElements=((int)(maxRowElements/4)+1)*4;//align16
	printf("maxRowElements = %d\n",maxRowElements);
	/*

	 // Sort the constraints.

	*/

//prepare matrix for copying to gpu

	allocateCPU((void**)&ccmaConstrData.h_matrEls, ccmaConstrData.ccmaConstrCount*sizeof(float)*maxRowElements);
	allocateCPU((void**)&ccmaConstrData.h_matrElIds, ccmaConstrData.ccmaConstrCount*sizeof(int)*maxRowElements);
	allocateCPU((void**)&ccmaConstrData.h_numValidEls, ccmaConstrData.ccmaConstrCount*sizeof(int));

	allocateGPU((void**)&ccmaConstrData.d_matrEls, ccmaConstrData.ccmaConstrCount*sizeof(float)*maxRowElements);
	allocateGPU((void**)&ccmaConstrData.d_matrElIds, ccmaConstrData.ccmaConstrCount*sizeof(int)*maxRowElements);
	allocateGPU((void**)&ccmaConstrData.d_numValidEls, ccmaConstrData.ccmaConstrCount*sizeof(int));

	for(int i=0; i< ccmaConstrData.ccmaConstrCount; i++){
		ccmaConstrData.h_numValidEls[i]=Nused[i];
		for(int j=0; j< Nused[i]; j++){
			ccmaConstrData.h_matrEls[maxRowElements*i+j]=invMval[i][j];
			ccmaConstrData.h_matrElIds[maxRowElements*i+j]=invMcol[i][j];
		}
	}
	cudaMemcpy(ccmaConstrData.d_matrEls, ccmaConstrData.h_matrEls,ccmaConstrData.ccmaConstrCount*maxRowElements*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy(ccmaConstrData.d_matrElIds, ccmaConstrData.h_matrElIds,ccmaConstrData.ccmaConstrCount*maxRowElements*sizeof(int), cudaMemcpyHostToDevice );
	cudaMemcpy(ccmaConstrData.d_numValidEls, ccmaConstrData.h_numValidEls,ccmaConstrData.ccmaConstrCount*sizeof(int), cudaMemcpyHostToDevice );

	cudaHostAlloc((void**)&ccmaConstrData.h_ccmaConvergedFlag, sizeof(int), cudaHostAllocMapped);
	cudaHostGetDevicePointer((void**) &ccmaConstrData.d_ccmaConvergedFlag, (void*) ccmaConstrData.h_ccmaConvergedFlag, 0);

	cudaMemcpyToSymbol(c_ccmaConstrData, &ccmaConstrData,
				sizeof(CcmaConstrData), 0, cudaMemcpyHostToDevice);

//	cudaBindTexture(0, t_matrEl, ccmaConstrData.d_matrEls, ccmaConstrData.ccmaConstrCount*maxRowElements*sizeof(float));
//	cudaBindTexture(0, t_matrElid, ccmaConstrData.d_matrElIds, ccmaConstrData.ccmaConstrCount*maxRowElements*sizeof(int));
//	cudaBindTexture(0, t_sigma, ccmaConstrData.d_constrSigma, ccmaConstrData.ccmaConstrTot*sizeof(float));

	ccmaBlockSize = BLOCK_SIZE;
	ccmaBlockCount1 = ccmaConstrData.ccmaConstrTot/ccmaBlockSize + 1;
	ccmaBlockCount2 = ccmaConstrData.ccmaAtomsTot/ccmaBlockSize + 1;


		LOG << "Done initializing ccma algorithm.";


}




__global__ void calculateDirections_0(){
	int index = threadIdx.x+blockIdx.x*blockDim.x;

	__shared__ int converged;
	if(threadIdx.x==0){
		converged=1;
	}
	if(index<c_ccmaConstrData.ccmaConstrTot){

		float4 oldDiff = c_ccmaConstrData.d_constrDir[index];
		float4 deltaDiff;
		int2 pair = c_ccmaConstrData.d_atomPair[index];
		float4 oldPosi = c_gsystem.d_coord[pair.x];
		float4 oldPosj = c_gsystem.d_coord[pair.y];
		float4 deltaCoordi = c_gsystem.d_midcoord[pair.x];
		float4 deltaCoordj = c_gsystem.d_midcoord[pair.y];

		oldDiff.x=oldPosi.x-oldPosj.x;
		oldDiff.y=oldPosi.y-oldPosj.y;
		oldDiff.z=oldPosi.z-oldPosj.z;
		c_ccmaConstrData.d_constrDir[index]=oldDiff;
		float d2 = oldDiff.w*oldDiff.w;
		oldDiff.w=oldDiff.x*oldDiff.x+oldDiff.y*oldDiff.y+oldDiff.z*oldDiff.z;
		deltaDiff.x=deltaCoordi.x-deltaCoordj.x;
		deltaDiff.y=deltaCoordi.y-deltaCoordj.y;
		deltaDiff.z=deltaCoordi.z-deltaCoordj.z;
		deltaDiff.w=deltaDiff.x*deltaDiff.x+deltaDiff.y*deltaDiff.y+deltaDiff.z*deltaDiff.z;
		float scalMult = oldDiff.x*deltaDiff.x+oldDiff.y*deltaDiff.y+oldDiff.z*deltaDiff.z;
		float sigma	=	oldDiff.w + 2*scalMult + deltaDiff.w - d2;
		c_ccmaConstrData.d_constrSigma[index]=sigma/((oldDiff.w+scalMult)*c_ccmaConstrData.d_M[index]);
		if ((converged==1) && (fabs(sigma) > 2*c_ccmaConstrData.tol*d2)){
			converged=0;
			*c_ccmaConstrData.d_ccmaConvergedFlag = 0;
		}
	}
}


__global__ void calculateDirections_1(){
	int index = threadIdx.x+blockIdx.x*blockDim.x;

	__shared__ int converged;
	if(threadIdx.x==0){
		converged=1;
	}
	if(index<c_ccmaConstrData.ccmaConstrTot){

		float4 oldDiff = c_ccmaConstrData.d_constrDir[index];
		int2 pair = c_ccmaConstrData.d_atomPair[index];
		float4 deltaDiff;
		float d2 = oldDiff.w*oldDiff.w;
		oldDiff.w=oldDiff.x*oldDiff.x+oldDiff.y*oldDiff.y+oldDiff.z*oldDiff.z;
		float4 deltaCoordi = c_gsystem.d_midcoord[pair.x];
		float4 deltaCoordj = c_gsystem.d_midcoord[pair.y];
		deltaDiff.x=deltaCoordi.x-deltaCoordj.x;
		deltaDiff.y=deltaCoordi.y-deltaCoordj.y;
		deltaDiff.z=deltaCoordi.z-deltaCoordj.z;
		deltaDiff.w=deltaDiff.x*deltaDiff.x+deltaDiff.y*deltaDiff.y+deltaDiff.z*deltaDiff.z;
		float scalMult = oldDiff.x*deltaDiff.x+oldDiff.y*deltaDiff.y+oldDiff.z*deltaDiff.z;
		float sigma	=	oldDiff.w + 2*scalMult + deltaDiff.w - d2;
		c_ccmaConstrData.d_constrSigma[index]=sigma/((oldDiff.w+scalMult)*c_ccmaConstrData.d_M[index]);
		if ((converged==1) && (fabs(sigma) > 2*c_ccmaConstrData.tol*d2)){
			converged=0;
			*c_ccmaConstrData.d_ccmaConvergedFlag = 0;
		}
	}
}



__global__ void multByCcmaMatrix(int maxRowElements, int ntr ){
	int index = threadIdx.x+blockIdx.x*blockDim.x;
	if (index < c_ccmaConstrData.ccmaConstrCount){
		for(int n=0; n< ntr; n++){
			float sum=0.0f;
			for (int i=0; i< c_ccmaConstrData.d_numValidEls[index]; i++){
				int j=c_ccmaConstrData.d_matrElIds[index*maxRowElements+i];//c_ccmaConstrData.d_matrElIds[index*maxRowElements+i];//tex1Dfetch(t_matrElid, index*maxRowElements+i);
				float mel = c_ccmaConstrData.d_matrEls[index*maxRowElements+i];//c_ccmaConstrData.d_matrEls[index*maxRowElements+i];//tex1Dfetch(t_matrEl, index*maxRowElements+i);
				float sigma = c_ccmaConstrData.d_constrSigma[j+n*c_ccmaConstrData.ccmaConstrCount];//c_ccmaConstrData.d_constrSigma[j+n*c_ccmaConstrData.ccmaConstrCount];//tex1Dfetch(t_sigma, j+n*c_ccmaConstrData.ccmaConstrCount);
				sum+=mel*sigma;
			}
			c_ccmaConstrData.d_constrLambda[index+n*c_ccmaConstrData.ccmaConstrCount]=sum;
		}
	}
}


__global__ void moveCcmaAtoms(){
	int index=threadIdx.x+blockIdx.x*blockDim.x;
	if (index<c_ccmaConstrData.ccmaAtomsTot){
		//int atomId=c_gsystem.d_CCMAatoms[index];
		float4 position = c_gsystem.d_midcoord[index];
		float mass = tex1Dfetch(t_m, (int)position.w);
		union constrIds{
			int4 i4;
			int i[4];
		};
		union constrIds constrId;
		constrId.i4=c_ccmaConstrData.d_attachedBonds[index];//here  is the reason of the ugliness below

		float sign;
		float lambda;
		float4 dir;
		float mult;
		for (int j=0; j< 4; j++){
			if(constrId.i[j]==0)
				break;
			sign = (constrId.i[j] > 0)? 1.0f : -1.0f;
			constrId.i[j]*=(int)sign;
			constrId.i[j]-=1;
			dir = c_ccmaConstrData.d_constrDir[constrId.i[j]];
			lambda = c_ccmaConstrData.d_constrLambda[constrId.i[j]];
			mult = -0.5 * sign * lambda/mass;

			position.x+= mult * dir.x;
			position.y+= mult * dir.y;
			position.z+= mult * dir.z;
		}

		c_gsystem.d_midcoord[index]=position;
	}
}


inline void compute(){
	int blocks = deviceProps.multiProcessorCount;
	int threads= ccmaConstrData.ccmaConstrTot/blocks + 1;


	if ( ccmaConstrData.ccmaConstrCount > 0 ){
		calculateDirections_0<<< ccmaBlockCount1, ccmaBlockSize>>>();
		const int checkInt = 1;
		for(int i=0; i<100; i++){
			//printf("%d\n", i);
			multByCcmaMatrix<<<ccmaConstrData.ccmaConstrCount/ccmaBlockSize + 1, ccmaBlockSize>>>(maxRowElements, parameters.Ntr);
			moveCcmaAtoms<<<ccmaBlockCount2, ccmaBlockSize>>>();
			if((i+1)%checkInt==0){
				cudaThreadSynchronize();
				*ccmaConstrData.h_ccmaConvergedFlag = 1;
			}
			calculateDirections_1<<< ccmaBlockCount1, ccmaBlockSize>>>();
			if ((i+1)%checkInt==0){
				cudaThreadSynchronize();
				if(*ccmaConstrData.h_ccmaConvergedFlag==1){
					break;
				}
			}
		}
	}
}


void destroy(){

}

#undef LOG

} // namespace leapfrog_integrator
