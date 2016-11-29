/*
 * model_dustyields.c
 *
 *  Created on: Oct2016
 *      Author: scottclay
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

#ifdef DETAILED_DUST
void update_dust_mass(int p, int centralgal, double dt, int nstep)
{
	int Zi;
	double timestep_width; //Width of current timestep in CODE UNITS
	int TimeBin; //Bin in Yield arrays corresponding to current timestep
	double Zi_disp, NormAGBDustYieldRate_actual[AGB_DUST_TYPE_NUM];
	double MassDiff;
	double timet, sfh_time;
	//double time_to_ts; //Time from high-z (upper) edge of SFH bin to middle of current timestep (used for massive SNII to hot) [in Myrs]
	//double tcut; //Maximum lifetime of stars that have their ejected put straight into the HotGas [in Myrs]
	double ColdGasSurfaceDensity, fwind, SNIIEjectaToHot; //Required for metal-rich wind implementation
	double DiskSFR, step_width_times_DiskSFR, DiskSFR_physical_units, step_width_times_DiskSFR_physical_units, inverse_DiskMass_physical_units;
	double BulgeSFR, step_width_times_BulgeSFR, BulgeSFR_physical_units, step_width_times_BulgeSFR_physical_units, inverse_BulgeMass_physical_units;
	double ICMSFR, step_width_times_ICMSFR, ICMSFR_physical_units, step_width_times_ICMSFR_physical_units, inverse_ICM_physical_units;
	double Disk_total_metallicity, Bulge_total_metallicity, ICM_total_metallicity;
	double NormMassEjecRateSumAllTypes;
	double TotalMassReturnedToColdDiskGas, TotalMassReturnedToHotGas;
	int n; //Iterator used for loop over NOUT when updating MassWeightedAge
	double AgeCorrectionDisk[NOUT];
	double AgeCorrectionBulge[NOUT];
	double CarOxyRatio;
	double SumAGBDust;
	int dust_scenario = 0;
	int dust_check = 0;

	TotalMassReturnedToColdDiskGas=0.0;
	TotalMassReturnedToHotGas=0.0;
	
	
	timestep_width = dt; //Width of current timestep in CODE UNITS (units cancel out when dividing by SFH bin width, sfh_dt) (12-04-12)
	TimeBin = (STEPS*(Halo[Gal[p].HaloNr].SnapNum-1.0))+nstep; //TimeBin = (STEPS*Gal[p].SnapNum)+nstep; //Bin in Yield tables corresponding to current timestep //TEST!: BRUNO: Snapnum would be +1 too low for a 'jumping' galaxy (14-11-13)
	timet = NumToTime((Halo[Gal[p].HaloNr].SnapNum-1.0)) - (nstep + 0.5) * dt; //Time from middle of the current timestep to z=0 (used here for MassWeightAge corrections)
	//NB: NumToTime(Gal[p].SnapNum) is the time to z=0 from start of current snapshot
	//    nstep is the number of the current timestep (0-19)
	//    dt is the width of one timestep within current snapshot
#ifdef METALRICHWIND			//if dust follows the metals we should leave this in
#ifdef GASDENSITYFWIND
	ColdGasSurfaceDensity = ((Gal[p].ColdGas*(1.0e10/Hubble_h))/(4.0*3.14159265*Gal[p].GasDiskRadius*Gal[p].GasDiskRadius/Hubble_h));
	fwind = min(1.0, (1.0/(ColdGasSurfaceDensity/5.0e12))); //1.0e13 //Fraction of SN-II ejecta put directly into HotGas
	if (Gal[p].ColdGas != (float)Gal[p].ColdGas) {fwind = 1.0;}
#endif
#ifndef GASDENSITYFWIND
	fwind = fwind_value;
#endif
#endif
#ifndef METALRICHWIND
	fwind = 0.0; //For all stellar ejecta (from disk) to ColdGas
#endif

    int i,j,k;
    for (i=0;i<=Gal[p].sfh_ibin;i++) //LOOP OVER SFH BINS
    {
    	sfh_time=Gal[p].sfh_t[i]+(0.5*Gal[p].sfh_dt[i]);
    	//time_to_ts = ((sfh_time+(0.5*Gal[p].sfh_dt[i])) - timet)*(UnitTime_in_years/Hubble_h)/1.0e6; //Time from high-z (upper) edge of SFH bin to middle of current timestep [in Myrs]
    	//tcut = 2.0*((Gal[p].Rvir/Gal[p].Vvir)/0.0001); //Maximum lifetime of stars that have their ejected put straight into the HotGas [in Myrs]



//*****************************************
//AGB ENRICHMENT FROM DISK STARS INTO COLD GAS:
//*****************************************
#ifdef DUST_AGB		
    if (Gal[p].sfh_DiskMass[i] > 0.0) {
     	//pre-calculations to speed up the code
    	DiskSFR = Gal[p].sfh_DiskMass[i]/Gal[p].sfh_dt[i];
    	step_width_times_DiskSFR = timestep_width * DiskSFR;
    	DiskSFR_physical_units = DiskSFR * (1.0e10/Hubble_h); //Note: This is NOT in physical units (i.e. NOT in Msun/yr, but in Msun/[code_time_units]). But this is ok, as code-time-units cancel out when multiplying by timestep_width to get 'step_width_times_DiskSFR_physical_units' on the line below ('DiskSFR_physical_units' is never used itself).
    	step_width_times_DiskSFR_physical_units = timestep_width * DiskSFR_physical_units;
    	inverse_DiskMass_physical_units=Hubble_h/(Gal[p].sfh_DiskMass[i]*1.0e10);
    	Disk_total_metallicity=metals_total(Gal[p].sfh_MetalsDiskMass[i])/Gal[p].sfh_DiskMass[i];

    	Zi = find_initial_metallicity_dust(p, i, 1, 1);
    	Zi_disp = (Disk_total_metallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.

			

 		//interpolates yields from lookup tables we produced in dust_yield_integrals.c
	    for (k=0;k<AGB_DUST_TYPE_NUM;k++)	//add numdustypes to allvars =4		//these are what we need for each type of dust
	    {
	    	NormAGBDustYieldRate_actual[k] = NormAGBDustYieldRate[TimeBin][i][Zi][k] + ((NormAGBDustYieldRate[TimeBin][i][Zi+1][k] - NormAGBDustYieldRate[TimeBin][i][Zi][k])*Zi_disp);	    	
	    }
	    
		CarOxyRatio = Gal[p].ColdGas_elements.Cb/Gal[p].ColdGas_elements.O;
		if (CarOxyRatio < 0.85) {
			dust_check = 1;
			//printf("1\n");
			}
		else if((CarOxyRatio >= 0.85) && (CarOxyRatio < 1.00)) {
			dust_check = 2;
			//printf("2\n");
			}
		else if(CarOxyRatio >= 1.00) {
			dust_check = 3;
			//printf("3\n");
			}
		else {
			dust_check = 3;
			//printf("3out\n");
		}
		/* 
		C/O ratios taken from Ferrarotti2006
		The transition from the silicate dominated mineral composition
		of M stars to the peculiar mixture of solids of S stars
		thus can be expected to occur around C/O = 0.85 while the
		carbon dominated dust mixture starts above C/O ≈ 1.00 
		*/

		switch (dust_check)
		{			
			case 1: //M stars
							
				//Create dust
				Gal[p].DustISM.AGB.Sil += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[0])); //M_forsterite
				Gal[p].DustISM.AGB.Sil += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[1])); //M_fayalite
				Gal[p].DustISM.AGB.Sil += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[2])); //M_enstatite
				Gal[p].DustISM.AGB.Sil += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[3])); //M_ferrosilite
				Gal[p].DustISM.AGB.Sil += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[4])); //M_quartz
		
				Gal[p].DustISM.AGB.Fe  += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[5])); //M_iron
		
				Gal[p].DustISM.AGB.SiC += 0.0; 
				Gal[p].DustISM.AGB.Cb   += 0.0; 
				
				SumAGBDust = 0.0;
				for (j=0; j<6; j++){
				SumAGBDust += (step_width_times_DiskSFR * NormAGBDustYieldRate_actual[j]);
				}			

				Gal[p].MetalsColdGas.agb -= max(0.0,SumAGBDust);
									
				//Remove new dust from specific elements in Coldgas
				Gal[p].ColdGas_elements.Si -= max(0.0,SIL_SI_DUST_FRACTION*(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[0]));
				Gal[p].ColdGas_elements.Si -= max(0.0,SIL_SI_DUST_FRACTION*(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[1]));
				Gal[p].ColdGas_elements.Si -= max(0.0,SIL_SI_DUST_FRACTION*(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[2]));
				Gal[p].ColdGas_elements.Si -= max(0.0,SIL_SI_DUST_FRACTION*(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[3]));
				Gal[p].ColdGas_elements.Si -= max(0.0,SIL_SI_DUST_FRACTION*(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[4]));

				Gal[p].ColdGas_elements.Mg -= max(0.0,SIL_MG_DUST_FRACTION*(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[0]));
				Gal[p].ColdGas_elements.Mg -= max(0.0,SIL_MG_DUST_FRACTION*(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[1]));
				Gal[p].ColdGas_elements.Mg -= max(0.0,SIL_MG_DUST_FRACTION*(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[2]));
				Gal[p].ColdGas_elements.Mg -= max(0.0,SIL_MG_DUST_FRACTION*(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[3]));
				Gal[p].ColdGas_elements.Mg -= max(0.0,SIL_MG_DUST_FRACTION*(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[4]));
		
				Gal[p].ColdGas_elements.Cb -= 0.0;
				Gal[p].ColdGas_elements.Fe -= 0.0;

				break;
		
			case 2: //S stars

				Gal[p].DustISM.AGB.Sil += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[6])); //S_quartz
		
				Gal[p].DustISM.AGB.Fe  += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[7])); //S_iron
		
				Gal[p].DustISM.AGB.SiC += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[8])); //S_SiC
		
				Gal[p].DustISM.AGB.Cb   += 0.0; 
				
				SumAGBDust = 0.0;
				for (j=6; j<9; j++){
				SumAGBDust += (step_width_times_DiskSFR * NormAGBDustYieldRate_actual[j]);
				}
				Gal[p].MetalsColdGas.agb -= max(0.0,SumAGBDust);
				
				Gal[p].ColdGas_elements.Si -= max(0.0,SIL_SI_DUST_FRACTION*(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[6]));
				Gal[p].ColdGas_elements.Mg -= max(0.0,SIL_MG_DUST_FRACTION*(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[6]));
		
				Gal[p].ColdGas_elements.Si -= max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[8]));

				Gal[p].ColdGas_elements.Cb -= 0.0;
		
				Gal[p].ColdGas_elements.Fe -= max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[7]));

				break;
				
			case 3: //C stars

				Gal[p].DustISM.AGB.Sil += 0.0; //C_quartz = none
		
				Gal[p].DustISM.AGB.Fe  += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[9])); //C_iron
		
				Gal[p].DustISM.AGB.SiC += 0.0; //C_SiC = none
		
				Gal[p].DustISM.AGB.Cb   += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[10])); //C_carbon

				SumAGBDust = 0.0;
				for (j=9; j<11; j++){
				SumAGBDust += (step_width_times_DiskSFR * NormAGBDustYieldRate_actual[j]);
				}
				Gal[p].MetalsColdGas.agb -= max(0.0,SumAGBDust);
		
				Gal[p].ColdGas_elements.Cb -= max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[10]));
		
				Gal[p].ColdGas_elements.Fe -= max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[9]));
				
				break;
		}
		} //if sfh_DM >0
    

#endif //DUST_AGB

//*****************************************
//AGB ENRICHMENT FROM ICL STARS INTO HOT GAS:			
//*****************************************
//Not yet implemented


//*****************************************
//SNII dust enrichment 			
//*****************************************


#ifdef DUST_SNII
if (Gal[p].sfh_DiskMass[i] > 0.0) {
		
	float eta_SNII_Sil = 0.00035;
	float eta_SNII_Fe  = 0.001;
	float eta_SNII_SiC = 0.0003;
	float eta_SNII_Cb  = 0.15;

	float A_Sil_dust = 121.41;
	float A_Fe_dust  = 55.85;
	float A_SiC_dust = 40.10; //SiC dust does not form in the ISM
	float A_Cb_dust  = 12.01;

	float A_Si = 28.0855;
	float A_Cb = 12.01;
	float A_Fe = 55.85;

	
	//These look like rates...but the dt is taken care of in recipe_yields and 
	//incorperated into SNII_prevstep........etc.
	
	Gal[p].DustISM.SNII.Sil += SNII_prevstep_Cold_Si[i] * eta_SNII_Sil * A_Sil_dust/A_Si;
	Gal[p].DustISM.SNII.Fe  += SNII_prevstep_Cold_Fe[i] * eta_SNII_Fe  * A_Fe_dust/A_Fe;
	Gal[p].DustISM.SNII.SiC += SNII_prevstep_Cold_Si[i] * eta_SNII_SiC * A_SiC_dust/A_Si;
	Gal[p].DustISM.SNII.Cb  += SNII_prevstep_Cold_Cb[i] * eta_SNII_Cb  * A_Cb_dust/A_Cb;
		
	Gal[p].ColdGas_elements.Si -= SNII_prevstep_Cold_Si[i] * eta_SNII_Sil * A_Sil_dust/A_Si;
	Gal[p].ColdGas_elements.Fe -= SNII_prevstep_Cold_Fe[i] * eta_SNII_Fe  * A_Fe_dust/A_Fe;
	Gal[p].ColdGas_elements.Si -= SNII_prevstep_Cold_Si[i] * eta_SNII_SiC * A_SiC_dust/A_Si;
	Gal[p].ColdGas_elements.Cb -= SNII_prevstep_Cold_Cb[i] * eta_SNII_Cb  * A_Cb_dust/A_Cb;
		
	Gal[p].MetalsColdGas.type2 -= (SNII_prevstep_Cold_Si[i] * eta_SNII_Sil * A_Sil_dust/A_Si)/(1.0e10/Hubble_h);
	Gal[p].MetalsColdGas.type2 -= (SNII_prevstep_Cold_Fe[i] * eta_SNII_Fe  * A_Fe_dust/A_Fe)/(1.0e10/Hubble_h);
	Gal[p].MetalsColdGas.type2 -= (SNII_prevstep_Cold_Si[i] * eta_SNII_SiC * A_SiC_dust/A_Si)/(1.0e10/Hubble_h);
	Gal[p].MetalsColdGas.type2 -= (SNII_prevstep_Cold_Cb[i] * eta_SNII_Cb  * A_Cb_dust/A_Cb)/(1.0e10/Hubble_h);
	
}
#endif //DUST_SNII
	
	
//*****************************************
//SNIa dust enrichment 			
//*****************************************

	// SNIa only produces Iron dust
	// ^^Need to check this! 
	// THIS NEEDS MULTIPLYING BY THE SNIa RATE
	// SEE ZHUKOVSKA PAPER
	// TURN OFF FOR NOW
	
#ifdef DUST_SNIA		
if (Gal[p].sfh_DiskMass[i] > 0.0) {
	float eta_SNIa_Sil = 0.0;
	float eta_SNIa_Fe  = 0.005;
	float eta_SNIa_SiC = 0.0;
	float eta_SNIa_Cb  = 0.0;

	float A_Fe_dust  = 55.85;
	float A_Fe = 55.85;

	Gal[p].DustISM.SNIa.Fe  += SNIa_prevstep_Cold_Fe[i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe;
	Gal[p].ColdGas_elements.Fe -= SNIa_prevstep_Cold_Fe[i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe;
	Gal[p].MetalsColdGas.type1a -= (SNIa_prevstep_Cold_Fe[i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe)/(1.0e10/Hubble_h);

}
#endif //DUST_SNIA


//*****************************************
//Growth of dust inside MC //maybe go last??			
//*****************************************

#ifdef DUST_GROWTH
    if ( (Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas.type2>0.0) && (Gal[p].MetalsColdGas.agb>0.0) ) {
			float t_acc_0, Xc, Z_sun, Z_coldgas, Z_fraction;
		
			Z_sun = 0.02;
			t_acc_0 = 2.0/UnitTime_in_Megayears;

			if (Gal[p].ColdGas > 0.00) {
				Z_coldgas = (Gal[p].MetalsColdGas.agb+Gal[p].MetalsColdGas.type2+Gal[p].MetalsColdGas.type1a)/Gal[p].ColdGas;
				Z_fraction = Z_coldgas/Z_sun;
				}
			else {
				Z_fraction = 0.0;
				}

			Xc = 0.5;	// Coldgas fraction
			//This should be the actual fraction of gas that is cold, but for now, 50%
			
			//**************************
			// Assumes 50/50 split between agb and type2. This needs fixing in the future
			Gal[p].DustISM.Growth.Sil += max(0.0,(Xc*Gal[p].DustISM.Growth.Sil*(dt/t_acc_0)*Z_fraction));			
			Gal[p].ColdGas_elements.Si -= max(0.0,(Xc*Gal[p].DustISM.Growth.Sil*(dt/t_acc_0)*Z_fraction));
			Gal[p].MetalsColdGas.agb -= max(0.0,0.5*(Xc*Gal[p].DustISM.Growth.Sil*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));
			Gal[p].MetalsColdGas.type2 -= max(0.0,0.5*(Xc*Gal[p].DustISM.Growth.Sil*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));
			
			Gal[p].DustISM.Growth.SiC += max(0.0,(Xc*Gal[p].DustISM.Growth.SiC*(dt/t_acc_0)*Z_fraction));			
			Gal[p].ColdGas_elements.Si -= max(0.0,(Xc*Gal[p].DustISM.Growth.SiC*(dt/t_acc_0)*Z_fraction));
			Gal[p].MetalsColdGas.agb -= max(0.0,0.5*(Xc*Gal[p].DustISM.Growth.SiC*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));
			Gal[p].MetalsColdGas.type2 -= max(0.0,0.5*(Xc*Gal[p].DustISM.Growth.SiC*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));

			Gal[p].DustISM.Growth.Fe += max(0.0,(Xc*Gal[p].DustISM.Growth.Fe*(dt/t_acc_0)*Z_fraction));			
			Gal[p].ColdGas_elements.Fe -= max(0.0,(Xc*Gal[p].DustISM.Growth.Fe*(dt/t_acc_0)*Z_fraction));
			Gal[p].MetalsColdGas.agb -= max(0.0,0.5*(Xc*Gal[p].DustISM.Growth.Fe*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));
			Gal[p].MetalsColdGas.type2 -= max(0.0,0.5*(Xc*Gal[p].DustISM.Growth.Fe*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));

			Gal[p].DustISM.Growth.Cb += max(0.0,(Xc*Gal[p].DustISM.Growth.Cb*(dt/t_acc_0)*Z_fraction));			
			Gal[p].ColdGas_elements.Cb -= max(0.0,(Xc*Gal[p].DustISM.Growth.Cb*(dt/t_acc_0)*Z_fraction));
			Gal[p].MetalsColdGas.agb -= max(0.0,0.5*(Xc*Gal[p].DustISM.Growth.Cb*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));
			Gal[p].MetalsColdGas.type2 -= max(0.0,0.5*(Xc*Gal[p].DustISM.Growth.Cb*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));
			
			//**************************
			
			
			Gal[p].DustISM.Growth.Sil += max(0.0,(Xc*Gal[p].DustISM.AGB.Sil*(dt/t_acc_0)*Z_fraction));			
			Gal[p].ColdGas_elements.Si -= max(0.0,(Xc*Gal[p].DustISM.AGB.Sil*(dt/t_acc_0)*Z_fraction));
			Gal[p].MetalsColdGas.agb -= max(0.0,(Xc*Gal[p].DustISM.AGB.Sil*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));

			Gal[p].DustISM.Growth.Sil += max(0.0,(Xc*Gal[p].DustISM.SNII.Sil*(dt/t_acc_0)*Z_fraction));			
			Gal[p].ColdGas_elements.Si -= max(0.0,(Xc*Gal[p].DustISM.SNII.Sil*(dt/t_acc_0)*Z_fraction));
			Gal[p].MetalsColdGas.type2 -= max(0.0,(Xc*Gal[p].DustISM.SNII.Sil*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));

			Gal[p].DustISM.Growth.SiC += max(0.0,(Xc*Gal[p].DustISM.AGB.SiC*(dt/t_acc_0)*Z_fraction));			
			Gal[p].ColdGas_elements.Si -= max(0.0,(Xc*Gal[p].DustISM.AGB.SiC*(dt/t_acc_0)*Z_fraction));
			Gal[p].MetalsColdGas.agb -= max(0.0,(Xc*Gal[p].DustISM.AGB.SiC*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));

			Gal[p].DustISM.Growth.SiC += max(0.0,(Xc*Gal[p].DustISM.SNII.SiC*(dt/t_acc_0)*Z_fraction));			
			Gal[p].ColdGas_elements.Si -= max(0.0,(Xc*Gal[p].DustISM.SNII.SiC*(dt/t_acc_0)*Z_fraction));
			Gal[p].MetalsColdGas.type2 -= max(0.0,(Xc*Gal[p].DustISM.SNII.SiC*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));
			
			Gal[p].DustISM.Growth.Fe += max(0.0,(Xc*Gal[p].DustISM.AGB.Fe*(dt/t_acc_0)*Z_fraction));			
			Gal[p].ColdGas_elements.Fe -= max(0.0,(Xc*Gal[p].DustISM.AGB.Fe*(dt/t_acc_0)*Z_fraction));
			Gal[p].MetalsColdGas.agb -= max(0.0,(Xc*Gal[p].DustISM.AGB.Fe*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));
			
			Gal[p].DustISM.Growth.Fe += max(0.0,(Xc*Gal[p].DustISM.SNII.Fe*(dt/t_acc_0)*Z_fraction));			
			Gal[p].ColdGas_elements.Fe -= max(0.0,(Xc*Gal[p].DustISM.SNII.Fe*(dt/t_acc_0)*Z_fraction));
			Gal[p].MetalsColdGas.type2 -= max(0.0,(Xc*Gal[p].DustISM.SNII.Fe*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));

			Gal[p].DustISM.Growth.Cb += max(0.0,(Xc*Gal[p].DustISM.AGB.Cb*(dt/t_acc_0)*Z_fraction));			
			Gal[p].ColdGas_elements.Cb -= max(0.0,(Xc*Gal[p].DustISM.AGB.Cb*(dt/t_acc_0)*Z_fraction));
			Gal[p].MetalsColdGas.agb -= max(0.0,(Xc*Gal[p].DustISM.AGB.Cb*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));

			Gal[p].DustISM.Growth.Cb += max(0.0,(Xc*Gal[p].DustISM.SNII.Cb*(dt/t_acc_0)*Z_fraction));			
			Gal[p].ColdGas_elements.Cb -= max(0.0,(Xc*Gal[p].DustISM.SNII.Cb*(dt/t_acc_0)*Z_fraction));
			Gal[p].MetalsColdGas.type2 -= max(0.0,(Xc*Gal[p].DustISM.SNII.Cb*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));
			
#ifdef DUST_SNIA		
			Gal[p].DustISM.Growth.Fe += max(0.0,(Xc*Gal[p].DustISM.SNIa.Fe*(dt/t_acc_0)*Z_fraction));			
			Gal[p].ColdGas_elements.Fe -= max(0.0,(Xc*Gal[p].DustISM.SNIa.Fe*(dt/t_acc_0)*Z_fraction));
			Gal[p].MetalsColdGas.type1a -= max(0.0,(Xc*Gal[p].DustISM.SNIa.Fe*(dt/t_acc_0)*Z_fraction)/(1.0e10/Hubble_h));
#endif //DUST_SNIA
} 
#endif //DUST_GROWTH

//*****************************************
//Dust destruction			
//*****************************************

#ifdef DUST_DESTRUCTION
    if ( (Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas.type2>0.0) && (Gal[p].MetalsColdGas.agb>0.0) ) {
			float t_des, M_cleared, f_SN, agb_sne_ratio;
			float des_frac, R_SN; 
			M_cleared = 1000; //Msol
			f_SN = 0.36; //Dimensionless
    	    DiskSFR = Gal[p].sfh_DiskMass[i]/Gal[p].sfh_dt[i];
			
			if( (DiskSFR>0.0) && (Gal[p].ColdGas>0.0) ) {
				t_des = (Gal[p].ColdGas*(1.0e10/Hubble_h))/M_cleared * 15.14/(0.1233*f_SN) * (Hubble_h * UnitTime_in_years)/(DiskSFR*1.0e10);
				des_frac = dt*UnitTime_in_years/t_des;
			}
			else {
				t_des = 0.0;
				des_frac = 0.0;
			}				
							
			if (des_frac < 1.0) {	
    		
    			//This AGB_SNE_RATIO needs fixing in the future. For now its assuming a 
    			// 50/50 split return between AGB and SNE. 
    			
    			//Sil dust destruction 
    			
				if( (Gal[p].DustISM.AGB.Sil>0.0) && (Gal[p].DustISM.SNII.Sil>0.0)) {
					//agb_sne_ratio = Gal[p].DustISM.AGB.Sil / Gal[p].DustISM.SNII.Sil;
					agb_sne_ratio = 0.5;
					}
				else {
					agb_sne_ratio = 0.5;
					}
				
				Gal[p].ColdGas_elements.Si +=  max(0.0,Gal[p].DustISM.AGB.Sil *des_frac);
				Gal[p].MetalsColdGas.agb += max(0.0,(Gal[p].DustISM.AGB.Sil *des_frac)/(1.0e10/Hubble_h));
				Gal[p].DustISM.AGB.Sil -= max(0.0,Gal[p].DustISM.AGB.Sil *des_frac);
			
				Gal[p].ColdGas_elements.Si +=  max(0.0,Gal[p].DustISM.SNII.Sil *des_frac);
				Gal[p].MetalsColdGas.type2 += max(0.0,(Gal[p].DustISM.SNII.Sil *des_frac)/(1.0e10/Hubble_h));
				Gal[p].DustISM.SNII.Sil -= max(0.0,Gal[p].DustISM.SNII.Sil *des_frac);

				Gal[p].ColdGas_elements.Si +=  max(0.0,Gal[p].DustISM.Growth.Sil *des_frac);
				Gal[p].MetalsColdGas.agb  += max(0.0,agb_sne_ratio*(Gal[p].DustISM.Growth.Sil *des_frac)/(1.0e10/Hubble_h));
				Gal[p].MetalsColdGas.type2+= max(0.0,(1.0-agb_sne_ratio)*(Gal[p].DustISM.Growth.Sil *des_frac)/(1.0e10/Hubble_h));			
				Gal[p].DustISM.Growth.Sil -= max(0.0,Gal[p].DustISM.Growth.Sil *des_frac);
 
 				//Fe dust
				if( (Gal[p].DustISM.AGB.Fe>0.0) && (Gal[p].DustISM.SNII.Fe>0.0)) {
					//agb_sne_ratio = Gal[p].DustISM.AGB.Fe / Gal[p].DustISM.SNII.Fe;
					agb_sne_ratio = 0.5;
					}
				else {
					agb_sne_ratio = 0.5;
					}

				Gal[p].ColdGas_elements.Fe +=  max(0.0,Gal[p].DustISM.AGB.Fe *des_frac);
				Gal[p].MetalsColdGas.agb += max(0.0,(Gal[p].DustISM.AGB.Fe *des_frac)/(1.0e10/Hubble_h));
				Gal[p].DustISM.AGB.Fe -= max(0.0,Gal[p].DustISM.AGB.Fe *des_frac);

				Gal[p].ColdGas_elements.Fe += max(0.0, Gal[p].DustISM.SNII.Fe *des_frac);
				Gal[p].MetalsColdGas.type2 += max(0.0,(Gal[p].DustISM.SNII.Fe *des_frac)/(1.0e10/Hubble_h));
				Gal[p].DustISM.SNII.Fe -= max(0.0,Gal[p].DustISM.SNII.Fe *des_frac);

				Gal[p].ColdGas_elements.Fe +=  max(0.0,Gal[p].DustISM.Growth.Fe *des_frac);
				Gal[p].MetalsColdGas.agb  += max(0.0,agb_sne_ratio*(Gal[p].DustISM.Growth.Fe *des_frac)/(1.0e10/Hubble_h));
				Gal[p].MetalsColdGas.type2+= max(0.0,(1.0-agb_sne_ratio)*(Gal[p].DustISM.Growth.Fe *des_frac)/(1.0e10/Hubble_h));			
				Gal[p].DustISM.Growth.Fe -= max(0.0,Gal[p].DustISM.Growth.Fe *des_frac);

#ifdef DUST_SNIA		
				Gal[p].ColdGas_elements.Fe += max(0.0, Gal[p].DustISM.SNIa.Fe *des_frac);
				Gal[p].MetalsColdGas.type1a += max(0.0,(Gal[p].DustISM.SNIa.Fe *des_frac)/(1.0e10/Hubble_h));
				Gal[p].DustISM.SNIa.Fe -= max(0.0,Gal[p].DustISM.SNIa.Fe *des_frac);
#endif //DUST_SNIA		
				
			
				//SiC dust
			
				if( (Gal[p].DustISM.AGB.SiC>0.0) && (Gal[p].DustISM.SNII.SiC>0.0)) {
					//agb_sne_ratio = Gal[p].DustISM.AGB.SiC / Gal[p].DustISM.SNII.SiC;
					agb_sne_ratio = 0.5;
					}
				else {
					agb_sne_ratio = 0.5;
					}

				Gal[p].ColdGas_elements.Si +=  max(0.0,Gal[p].DustISM.AGB.SiC *des_frac);
				Gal[p].MetalsColdGas.agb += max(0.0,(Gal[p].DustISM.AGB.SiC *des_frac)/(1.0e10/Hubble_h));
				Gal[p].DustISM.AGB.SiC -= max(0.0,Gal[p].DustISM.AGB.SiC *des_frac);

				Gal[p].ColdGas_elements.Si +=  max(0.0,Gal[p].DustISM.SNII.SiC *des_frac);
				Gal[p].MetalsColdGas.type2 += max(0.0,(Gal[p].DustISM.SNII.SiC *des_frac)/(1.0e10/Hubble_h));
				Gal[p].DustISM.SNII.SiC -= max(0.0,Gal[p].DustISM.SNII.SiC *des_frac);

				Gal[p].ColdGas_elements.Si +=  max(0.0,Gal[p].DustISM.Growth.SiC *des_frac);
				Gal[p].MetalsColdGas.agb  += max(0.0,agb_sne_ratio*(Gal[p].DustISM.Growth.SiC *des_frac)/(1.0e10/Hubble_h));
				Gal[p].MetalsColdGas.type2+= max(0.0,(1.0-agb_sne_ratio)*(Gal[p].DustISM.Growth.SiC *des_frac)/(1.0e10/Hubble_h));			
				Gal[p].DustISM.Growth.SiC -= max(0.0,Gal[p].DustISM.Growth.SiC *des_frac);


				//Cb dust
			
				if( (Gal[p].DustISM.AGB.Cb>0.0) && (Gal[p].DustISM.SNII.Cb>0.0)) {
					//agb_sne_ratio = Gal[p].DustISM.AGB.Cb / Gal[p].DustISM.SNII.Cb;
					agb_sne_ratio = 0.5;
					}
				else {
					agb_sne_ratio = 0.5;
					}

				Gal[p].ColdGas_elements.Cb +=  max(0.0,Gal[p].DustISM.AGB.Cb *des_frac);
				Gal[p].MetalsColdGas.agb += max(0.0,(Gal[p].DustISM.AGB.Cb *des_frac)/(1.0e10/Hubble_h));
				Gal[p].DustISM.AGB.Cb -= max(0.0,Gal[p].DustISM.AGB.Cb *des_frac);

				Gal[p].ColdGas_elements.Cb +=  max(0.0,Gal[p].DustISM.SNII.Cb *des_frac);
				Gal[p].MetalsColdGas.type2 += max(0.0,(Gal[p].DustISM.SNII.Cb *des_frac)/(1.0e10/Hubble_h));
				Gal[p].DustISM.SNII.Cb -= max(0.0,Gal[p].DustISM.SNII.Cb *des_frac);

				Gal[p].ColdGas_elements.Cb +=  max(0.0,Gal[p].DustISM.Growth.Cb *des_frac);
				Gal[p].MetalsColdGas.agb  += max(0.0,agb_sne_ratio*(Gal[p].DustISM.Growth.Cb *des_frac)/(1.0e10/Hubble_h));
				Gal[p].MetalsColdGas.type2+= max(0.0,(1.0-agb_sne_ratio)*(Gal[p].DustISM.Growth.Cb *des_frac)/(1.0e10/Hubble_h));			
				Gal[p].DustISM.Growth.Cb -= max(0.0,Gal[p].DustISM.Growth.Cb *des_frac);

			
			} //des_Frac > 1.0
		}		
#endif //DUST_DESTRUCTION

	
//} //if coldgas > 1.0e7
	
    } //for (i=0;i<=Gal[p].sfh_ibin;i++) //MAIN LOOP OVER SFH BINS
    			//AGB NEED to be inside SFH bin loop as it depends on current SFR
    			//SNII and Ia NEED to be inside SFH bin loop as it uses something from recipe_yields which
    			//depends on the SFH bin. 
    			//Growth and Destruction_SNe should NOT BE INSIDE THIS LOOP (or should it??????????????)
    			//Destruction_SF NEEDS to be inside this loop. 
		//printf("Finishing dust yield code\n");
		//printf("METALS.AGB = %g\n",Gal[p].MetalsColdGas.agb);

} //update dust mass 


int find_initial_metallicity_SNe_rates(int Zi, int sfh_bin, int table_type)
{
	int i, Zi_bin;
	double Z_in;

	Zi_bin = -1;
	i = 0;
	Z_in = lifetimeMetallicities[Zi];

	switch (table_type)
	{
		case 2: //SN-II metallicity table
			while (Zi_bin == -1)
			{
				if (SNIIMetallicities[i] < Z_in)
				{
					i++;
					if (i == SNII_Z_NUM) Zi_bin = i-1;
				}
				else Zi_bin = i;
			}
			break;

	}
	return Zi_bin;
}



int find_initial_metallicity_dust(int p, int sfh_bin, int table_type, int component_type)
{
	if (component_type == 1) //Disk stars
	{
	int i, Zi_bin;
	double initMetals, Z_disk;

	initMetals = metals_total(Gal[p].sfh_MetalsDiskMass[sfh_bin]); //IN [10^10/h Msun]
	Zi_bin = -1;
	i = 0;
	if (initMetals == 0.0 || Gal[p].sfh_DiskMass[sfh_bin] == 0.0)
	{
		Z_disk = 0.0;
	}
	else Z_disk = initMetals/Gal[p].sfh_DiskMass[sfh_bin]; //Dimensionless

	switch (table_type)
	{
		case 1: //Lifetime metallicity table
			while (Zi_bin == -1)
			{
				if (lifetimeMetallicities[i] < Z_disk)
				{
					i++;
					if (i == LIFETIME_Z_NUM) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
				}
				else Zi_bin = i;
			}
			break;
		//case 3 //SNIa yields are NOT metallicity dependent
		case 4: //Dust metallicity table
			while (Zi_bin == -1)
			{
				if (AGBDustMetallicities[i] < Z_disk)
				{
					i++;
					if (i == 3) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
				}
				else Zi_bin = i;
			}
			break;
	}

	if (Zi_bin == 0 ) return Zi_bin;
	else return Zi_bin-1;
	}
	else if (component_type == 2) //Bulge stars
	{
		int i, Zi_bin;
		double initMetals, Z_bulge;

		initMetals = metals_total(Gal[p].sfh_MetalsBulgeMass[sfh_bin]); //IN [10^10/h Msun]
		Zi_bin = -1;
		i = 0;
		if (initMetals == 0.0 || Gal[p].sfh_BulgeMass[sfh_bin] == 0.0)
		{
			Z_bulge = 0.0;
		}
		else Z_bulge = initMetals/Gal[p].sfh_BulgeMass[sfh_bin];

		switch (table_type)
		{
			case 1: //Lifetime metallicity table
				while (Zi_bin == -1)
				{
					if (lifetimeMetallicities[i] < Z_bulge) //Gal[p].sfh_MetalsDiskMass[sfh_bin].type2/Gal[p].sfh_DiskMass[sfh_bin])
					{
						i++;
						if (i == LIFETIME_Z_NUM) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
					}
					else Zi_bin = i;
				}
				break;
			//case 3 //SNIa yields are NOT metallicity dependent
			case 4: //Dust metallicity table
				while (Zi_bin == -1)
				{
					if (AGBDustMetallicities[i] < Z_bulge)
					{
						i++;
						if (i == 3) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
					}
					else Zi_bin = i;
				}
				break;
		}
		if (Zi_bin == 0 ) return Zi_bin;
		else return Zi_bin-1;
	}
	else if (component_type == 3) //ICL stars
		{
			int i, Zi_bin;
			double initMetals, Z_ICM;

			initMetals = metals_total(Gal[p].sfh_MetalsICM[sfh_bin]); //IN [10^10/h Msun]
			Zi_bin = -1;
			i = 0;
			if (initMetals == 0.0 || Gal[p].sfh_ICM[sfh_bin] == 0.0)
			{
				Z_ICM = 0.0;
			}
			else Z_ICM = initMetals/Gal[p].sfh_ICM[sfh_bin];

			switch (table_type)
			{
				case 1: //Lifetime metallicity table
					while (Zi_bin == -1)
					{
						if (lifetimeMetallicities[i] < Z_ICM)
						{
							i++;
							if (i == LIFETIME_Z_NUM) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
						}
						else Zi_bin = i;
					}
					break;
				//case 3 //SNIa yields are NOT metallicity dependent
				case 4: //AGB metallicity table
					while (Zi_bin == -1)
					{
						if (AGBDustMetallicities[i] < Z_ICM)
						{
							i++;
							if (i == 3) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
						}
						else Zi_bin = i;
					}
					break;
			}
			if (Zi_bin == 0 ) return Zi_bin;
			else return Zi_bin-1;
		}
	else { printf("Wrong stellar component type for Z_init calculation: Use either 1 (disk), 2 (bulge) or 3 (ICL)"); exit(1);}
}

#endif ///DDust

