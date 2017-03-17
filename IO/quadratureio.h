/*
 * quadratureio.h
 *
 *  Created on: Jan 4, 2011
 *      Author: zhmurov
 */

#pragma once

typedef struct {
	float theta;
	float phi;
	float x;
	float y;
	float z;
	float w;
} AngularQuadraturePoint;

typedef struct {
	int pointCount;
	AngularQuadraturePoint* points;
} AngularQuadrature;

typedef struct {
	float r;
	float w;
} RadialQuadraturePoint;

typedef struct {
	int pointCount;
	RadialQuadraturePoint* points;
	float x1;
	float x2;
} RadialQuadrature;

void readAngularQuadrature(const char* filename, AngularQuadrature* quadr);
void createGaussianLegendreRadialQuadrature(float x1, float x2, int pointCount, RadialQuadrature* quadr);

float integrateAngularQuadrature(float (*func)(float, float, float), AngularQuadrature quadr, float r);
float integrateRadialQuadrature(float (*func)(float), RadialQuadrature quadr);

