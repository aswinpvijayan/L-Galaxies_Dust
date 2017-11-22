/*
 * Created on: Oct2016
 * Author: scottclay
 *
 * Various functions to help when coding with dust. 
 * 
 * Note (Nov2017) - Many of these functions are no longer in use! 
 * Many of them contain references to dust structures that no longer exist. 
 * So be careful!
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"


#ifdef DETAILED_DUST

struct DustMass DustMass_init() {
	struct DustMass dust;
	dust.AGB.SiC    = 0.0;
	dust.AGB.Sil    = 0.0;
	dust.AGB.Cb     = 0.0;
	dust.AGB.Fe     = 0.0;
	dust.SNII.SiC   = 0.0;
	dust.SNII.Sil   = 0.0;
	dust.SNII.Cb    = 0.0;
	dust.SNII.Fe    = 0.0;
	dust.SNIa.SiC   = 0.0;
	dust.SNIa.Sil   = 0.0;
	dust.SNIa.Cb    = 0.0;
	dust.SNIa.Fe    = 0.0;
	dust.Growth.SiC = 0.0;
	dust.Growth.Sil = 0.0;
	dust.Growth.Cb  = 0.0;
	dust.Growth.Fe  = 0.0;
	dust.Destruction_SNe.SiC = 0.0;
	dust.Destruction_SNe.Sil = 0.0;
	dust.Destruction_SNe.Cb = 0.0;
	dust.Destruction_SNe.Fe = 0.0;
	dust.Destruction_SF.SiC = 0.0;
	dust.Destruction_SF.Sil = 0.0;
	dust.Destruction_SF.Cb = 0.0;
	dust.Destruction_SF.Fe = 0.0;
	return (dust);
}

float DustMass_Total(struct DustMass dust) {
	float Sum = 0.0;
	Sum += dust.AGB.SiC;
	Sum += dust.AGB.Sil;
	Sum += dust.AGB.Fe;
	Sum += dust.AGB.Cb;
	
	Sum += dust.SNII.SiC;
	Sum += dust.SNII.Sil;
	Sum += dust.SNII.Fe;
	Sum += dust.SNII.Cb;
	
	Sum += dust.SNIa.Fe;
	
	Sum += dust.Growth.Fe;
	
	return Sum;
}

void Print_DustMass(struct DustMass dust) {
	printf("Start of galaxy dust output\n");
	printf("Dust.AGB.SiC = %g\n",dust.AGB.SiC);
	printf("Dust.AGB.Sil = %g\n",dust.AGB.Sil);
	printf("Dust.AGB.Fe = %g\n",dust.AGB.Fe);
	printf("Dust.AGB.Cb = %g\n",dust.AGB.Cb);
	printf("Dust.SNII.SiC = %g\n",dust.SNII.SiC);
	printf("Dust.SNII.Sil = %g\n",dust.SNII.Sil);
	printf("Dust.SNII.Fe = %g\n",dust.SNII.Fe);
	printf("Dust.SNII.Cb = %g\n",dust.SNII.Cb);
	printf("Dust.SNIa.Fe = %g\n",dust.SNIa.Fe);
	printf("Dust.Growth.Fe = %g\n",dust.Growth.Fe);
	}
float DustMass_Total_SiC(struct DustMass dust) {
	return(dust.AGB.SiC+dust.SNII.SiC+dust.SNIa.SiC+dust.Growth.SiC-dust.Destruction_SNe.SiC-dust.Destruction_SF.SiC);
}

float DustMass_Total_Sil(struct DustMass dust) {
	return(dust.AGB.Sil+dust.SNII.Sil+dust.SNIa.Sil+dust.Growth.Sil-dust.Destruction_SNe.Sil-dust.Destruction_SF.Sil);
}

float DustMass_Total_Cb(struct DustMass dust) {
	return(dust.AGB.Cb+dust.SNII.Cb+dust.SNIa.Cb+dust.Growth.Cb-dust.Destruction_SNe.Cb-dust.Destruction_SF.Cb);
}

float DustMass_Total_Fe(struct DustMass dust) {
	return(dust.AGB.Fe+dust.SNII.Fe+dust.SNIa.Fe+dust.Growth.Fe-dust.Destruction_SNe.Fe-dust.Destruction_SF.Fe);
}


float DustMass_AGB_Total(struct DustMass dust) {
	return(dust.AGB.SiC+dust.AGB.Sil+dust.AGB.Cb+dust.AGB.Fe);
}

float DustMass_SNII_Total(struct DustMass dust) {
	return(dust.SNII.SiC+dust.SNII.Sil+dust.SNII.Cb+dust.SNII.Fe);
}

float DustMass_SNIa_Total(struct DustMass dust) {
	return(dust.SNIa.SiC+dust.SNIa.Sil+dust.SNIa.Cb+dust.SNIa.Fe);
}

float DustMass_Growth_Total(struct DustMass dust) {
	return(dust.Growth.SiC+dust.Growth.Sil+dust.Growth.Cb+dust.Growth.Fe);
}

float DustMass_Destruction_Total(struct DustMass dust) {
	return(dust.Destruction_SNe.SiC + dust.Destruction_SNe.Sil + dust.Destruction_SNe.Cb + dust.Destruction_SNe.Fe + dust.Destruction_SF.SiC + dust.Destruction_SF.Sil + dust.Destruction_SF.Cb + dust.Destruction_SF.Fe);
}

float DustMass_Total_Created(struct DustMass dust) {
	float AGB,SNII,SNIa,Growth;
	AGB = DustMass_AGB_Total(dust);
	SNII= DustMass_SNII_Total(dust);
	SNIa= DustMass_SNIa_Total(dust);
	Growth = DustMass_Growth_Total(dust);
	return (AGB + SNII + SNIa + Growth);
}
	
float DustMass_Total_Created_Minus_Destroy(struct DustMass dust) {
	float AGB,SNII,SNIa,Growth,Destruction;
	AGB = DustMass_AGB_Total(dust);
	SNII= DustMass_SNII_Total(dust);
	SNIa= DustMass_SNIa_Total(dust);
	Growth = DustMass_Growth_Total(dust);
	Destruction = DustMass_Destruction_Total(dust);
	return (AGB + SNII + SNIa + Growth - Destruction);
}

struct DustMass DustMass_add(struct DustMass dust1, struct DustMass dust2, float fraction)
{
	struct DustMass dust;
	dust.AGB.SiC = dust1.AGB.SiC + fraction*dust2.AGB.SiC;
	dust.AGB.Sil = dust1.AGB.Sil + fraction*dust2.AGB.Sil;
	dust.AGB.Cb  = dust1.AGB.Cb  + fraction*dust2.AGB.Cb;
	dust.AGB.Fe  = dust1.AGB.Fe  + fraction*dust2.AGB.Fe;
	
	dust.SNII.SiC = dust1.SNII.SiC + fraction*dust2.SNII.SiC;
	dust.SNII.Sil = dust1.SNII.Sil + fraction*dust2.SNII.Sil;
	dust.SNII.Cb  = dust1.SNII.Cb  + fraction*dust2.SNII.Cb;
	dust.SNII.Fe  = dust1.SNII.Fe  + fraction*dust2.SNII.Fe;
	
	dust.SNIa.SiC = dust1.SNIa.SiC + fraction*dust2.SNIa.SiC;
	dust.SNIa.Sil = dust1.SNIa.Sil + fraction*dust2.SNIa.Sil;
	dust.SNIa.Cb  = dust1.SNIa.Cb  + fraction*dust2.SNIa.Cb;
	dust.SNIa.Fe  = dust1.SNIa.Fe  + fraction*dust2.SNIa.Fe;

	dust.Growth.SiC = dust1.Growth.SiC + fraction*dust2.Growth.SiC;
	dust.Growth.Sil = dust1.Growth.Sil + fraction*dust2.Growth.Sil;
	dust.Growth.Cb  = dust1.Growth.Cb  + fraction*dust2.Growth.Cb;
	dust.Growth.Fe  = dust1.Growth.Fe  + fraction*dust2.Growth.Fe;
	
	//Not sure how to deal with destruction...
	dust.Destruction_SNe.SiC = dust1.Destruction_SNe.SiC + fraction*dust2.Destruction_SNe.SiC;
	dust.Destruction_SNe.Sil = dust1.Destruction_SNe.Sil + fraction*dust2.Destruction_SNe.Sil;
	dust.Destruction_SNe.Cb  = dust1.Destruction_SNe.Cb  + fraction*dust2.Destruction_SNe.Cb;
	dust.Destruction_SNe.Fe  = dust1.Destruction_SNe.Fe  + fraction*dust2.Destruction_SNe.Fe;
	
	dust.Destruction_SF.SiC = dust1.Destruction_SF.SiC + fraction*dust2.Destruction_SF.SiC;
	dust.Destruction_SF.Sil = dust1.Destruction_SF.Sil + fraction*dust2.Destruction_SF.Sil;
	dust.Destruction_SF.Cb  = dust1.Destruction_SF.Cb  + fraction*dust2.Destruction_SF.Cb;
	dust.Destruction_SF.Fe  = dust1.Destruction_SF.Fe  + fraction*dust2.Destruction_SF.Fe;

	
	return(dust);
}
#endif //DETAILED_DUST
