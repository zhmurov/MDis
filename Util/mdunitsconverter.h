/*
 * mdunitsconverter.h
 *
 *  Created on: Aug 9, 2009
 *      Author: zhmurov
 */

#pragma once

#define FEL_MD 1.38935456e2
#define Na 6.02214179e23

float convertTimeCHARMMtoMD(float fs);
float convertKbCHARMMtoMD(float kcalpermolA2);
float convertDistanceCHARMMtoMD(float A);
float convertKthetaCHARMMtoMD(float kcalpermolrad2);
float convertKpsiCHARMMtoMD(float kcalpermolrad2);
float convertAngleCHARMMtoMD(float deg);
float convertKubCHARMMtoMD(float kcalpermolA2);
float convertKchiCHARMMtoMD(float kcalpermol);
float convertEpsilonCHARMMtoMD(float kcalpermol);
