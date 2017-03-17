/*
 * ffutils.h
 *
 *  Created on: Jul 5, 2009
 *      Author: zhmurov
 */

#pragma once

FFBondType findBondType(char atomType1[5], char atomType2[5], ForceField* forceField);
FFAngleType findAngleType(char atomType1[5], char atomType2[5], char atomType3[5], ForceField* forceField);
int findDihedralType(char atomType1[5], char atomType2[5], char atomType3[5], char atomType4[5],
		ForceField* forceField, FFDihedralType* dt);
int findImproperType(char atomType1[5], char atomType2[5], char atomType3[5], char atomType4[5],
		ForceField* forceField, FFImproperType* it);
FFCMAPType findCMAPType(
		char atomType1[5], char atomType2[5], char atomType3[5], char atomType4[5],
		char atomType5[5], char atomType6[5], char atomType7[5], char atomType8[5],
		ForceField* forceField);
FFNonbondedType findNonbondedType(char atomType[5], ForceField* forceField);
