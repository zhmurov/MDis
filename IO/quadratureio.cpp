/*
 * quadratureio.cpp
 *
 *  Created on: Jan 4, 2011
 *      Author: zhmurov
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "quadratureio.h"

#define buf_size 80

void readAngularQuadrature(const char* filename, AngularQuadrature* quadr){
	printf("Reading angular quadrature from %s.\n", filename);
	FILE* file = fopen(filename, "r");
	char buffer[buf_size];
	int pointCount = 0;
	if(file != NULL){
		while(fgets(buffer, buf_size, file) != NULL){
			pointCount++;
		}
		printf("Found %d points in quadrature.\n", pointCount);
		rewind(file);
		quadr->pointCount = pointCount;
		quadr->points = (AngularQuadraturePoint*)calloc(pointCount, sizeof(AngularQuadraturePoint));
		pointCount = 0;
		while(fgets(buffer, buf_size, file) != NULL){
			char* pch = strtok(buffer, " \t");
			quadr->points[pointCount].theta = atof(pch)*M_PI/180.0f;
			pch = strtok(NULL, " \t");
			quadr->points[pointCount].phi = atof(pch)*M_PI/180.0f;
			pch = strtok(NULL, " \t");
			quadr->points[pointCount].w = atof(pch);
			quadr->points[pointCount].x = cosf(quadr->points[pointCount].theta)*sinf(quadr->points[pointCount].phi);
			quadr->points[pointCount].y = sinf(quadr->points[pointCount].theta)*sinf(quadr->points[pointCount].phi);
			quadr->points[pointCount].z = cosf(quadr->points[pointCount].phi);
			pointCount++;
		}
		printf("Quadrature read:\n");
		printf("     Theta:      Phi:        Weight:     x:          y:          z:\n");
		int i;
		for(i = 0; i < quadr->pointCount; i++){
			printf("%-*d%-*f%-*f%-*f%-*f%-*f%-*f\n",
					5, i + 1,
					12, quadr->points[i].theta,
					12, quadr->points[i].phi,
					12, quadr->points[i].w,
					12, quadr->points[i].x,
					12, quadr->points[i].y,
					12, quadr->points[i].z);
		}
	} else {
		printf("File '%s' not found.\n", filename);
	}
	fclose(file);
	printf("Done reading angular quadrature.\n");
}

/*
 * Taken from Numerical Recipes: The Art of Scientific Computing
 * by W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery
 * http://www.nr.com/
 *
 * Given the lower and upper limits of integration x1 and x2, this routine returns arrays x[0..n-1]
 * and w[0..n-1] of length n, containing the abscissas and weights of the Gauss-Legendre n-point
 * quadrature formula.*/

void createGaussianLegendreRadialQuadrature(float x1, float x2, int pointCount, RadialQuadrature* quadr){
	printf("Computing weights/coordinates for Gaussian-Legandre quadrature for [%5.3f,%5.3f] with %d integration points.\n", x1, x2, pointCount);
	const double EPS = 3.0e-14; // EPS is the relative precision.
	quadr->pointCount = pointCount;
	quadr->points = (RadialQuadraturePoint*)calloc(pointCount, sizeof(RadialQuadraturePoint));
	quadr->x1 = x1;
	quadr->x2 = x2;
	double z1, z, xm, xl, pp, p3, p2, p1;
	int n = pointCount;
	int m = (n + 1)/2; //The roots are symmetric in the interval, so we only have to find half of them.
	xm = 0.5*(x2 + x1);
	xl = 0.5*(x2 - x1);
	int i, j;
	for(i = 0; i < m; i++){ //Loop over the desired roots.
		z = cos(M_PI*((double)i + 0.75)/((double)n + 0.5)); //Starting with this approximation to the ith root,
											//we enter the main loop of refinement by Newton’s method.
		do {
			p1 = 1.0;
			p2 = 0.0;
			for(j = 0; j < n; j++){ //Loop up the recurrence relation to get the
				p3 = p2; //Legendre polynomial evaluated at z.
				p2 = p1;
				p1 = ((2.0*((double)j) + 1.0)*z*p2 - ((double)j)*p3)/((double)j + 1.0);
			}   //p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
			    //by a standard relation involving also p2, the polynomial of one lower order.
			pp = ((double)n)*(z*p1 - p2)/(z*z - 1.0);
			z1 = z;
			z = z1 - p1/pp; //Newton’s method.
		} while(fabs(z - z1) > EPS);
		quadr->points[i].r = xm - xl*z; //Scale the root to the desired interval,
		quadr->points[n-1-i].r = xm + xl*z; //and put in its symmetric counterpart.
		quadr->points[i].w = 2.0*xl/((1.0 - z*z)*pp*pp); //Compute the weight
		quadr->points[n-1-i].w = quadr->points[i].w; //and its symmetric counterpart.
	}
	printf("Will use following quadrature data:\n");
	printf("r, nm       w\n");
	for(i = 0; i < quadr->pointCount; i++){
		printf("%-*f%-*f\n", 12, quadr->points[i].r, 12, quadr->points[i].w);
	}

	printf("Done computing weights/coordinates for Gaussian-Legandre quadrature.\n");
}


float integrateAngularQuadrature(float (*func)(float, float, float), AngularQuadrature quadr, float r){
	int i;
	float val = 0.0f;
	for(i = 0; i < quadr.pointCount; i++){
		val += quadr.points[i].w*(*func)(r*quadr.points[i].x, r*quadr.points[i].y, r*quadr.points[i].z);
	}
	val *= 4.0f*M_PI;
	return val;
}

float integrateRadialQuadrature(float (*func)(float), RadialQuadrature quadr){
	int i;
	float val = 0.0f;
	for(i = 0; i < quadr.pointCount; i++){
		val += quadr.points[i].w*(*func)(quadr.points[i].r);
	}
	return val;
}
