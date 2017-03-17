/*
 * mdunitsconverter.c
 *
 *  Created on: Aug 9, 2009
 *      Author: zhmurov
 */
#include <math.h>

float convertTimeCHARMMtoMD(float fs){
	return fs*1.0e-3;
}

float convertKbCHARMMtoMD(float kcalpermolA2){
	return kcalpermolA2*100.0f*4.184f;
}

float convertDistanceCHARMMtoMD(float A){
	return A/10.0f;
}

float convertKthetaCHARMMtoMD(float kcalpermolrad2){
	return kcalpermolrad2*4.184f;
}

float convertKpsiCHARMMtoMD(float kcalpermolrad2){
	return kcalpermolrad2*4.184f;
}

float convertAngleCHARMMtoMD(float deg){
	return deg*M_PI/180.0f;
}

float convertKubCHARMMtoMD(float kcalpermolA2){
	return kcalpermolA2*100.0f*4.184f;
}

float convertKchiCHARMMtoMD(float kcalpermol){
	return kcalpermol*4.184f;
}

float convertEpsilonCHARMMtoMD(float kcalpermol){
	return kcalpermol*4.184f;
}
