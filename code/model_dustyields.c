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
	double agb_ratio, type2_ratio, type1a_ratio;

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

//if( (Gal[p].MetalsColdGas.agb > 0.0) && (Gal[p].MetalsColdGas.type2 > 0.0) && (Gal[p].MetalsColdGas.type1a > 0.0) ) { 

//*****************************************
//AGB ENRICHMENT FROM DISK STARS INTO COLD GAS:
//*****************************************

printf("Total Dust = %g\n",elements_total(Gal[p].Dust_elements));

#ifdef DUST_AGB		
    if ( (Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas.agb >0.0) ) {
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
	    	//NormAGBDustYieldRate_actual[k] = NormAGBDustYieldRate[TimeBin][i][Zi][k] + ((NormAGBDustYieldRate[TimeBin][i][Zi+1][k] - NormAGBDustYieldRate[TimeBin][i][Zi][k])*Zi_disp);	    	
	    	NormAGBDustYieldRate_actual[k] = NormAGBDustYieldRate[TimeBin][i][Zi_saved][k] + ((NormAGBDustYieldRate[TimeBin][i][Zi_saved+1][k] - NormAGBDustYieldRate[TimeBin][i][Zi_saved][k])*Zi_disp_saved);	    	
	    }
	    
							
#ifdef FULL_DUST							
		Gal[p].DustISM.AGB.Sil += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[0])); //M_forsterite
		Gal[p].DustISM.AGB.Sil += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[1])); //M_fayalite
		Gal[p].DustISM.AGB.Sil += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[2])); //M_enstatite
		Gal[p].DustISM.AGB.Sil += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[3])); //M_ferrosilite
		Gal[p].DustISM.AGB.Sil += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[4])); //M_quartz
		Gal[p].DustISM.AGB.Fe  += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[5])); //M_iron
		Gal[p].DustISM.AGB.Sil += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[6])); //S_quartz
        Gal[p].DustISM.AGB.Fe  += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[7])); //S_iron
		Gal[p].DustISM.AGB.SiC += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[8])); //C_SiC
        Gal[p].DustISM.AGB.Fe  += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[9])); //C_iron
        Gal[p].DustISM.AGB.Cb  += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[10])); //C_carbon
#endif				
		//Calculate the amount of dust CREATED ----------------------------------------------------------------------
		
		double Dust_Forsterite  = max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[0])); //M_forsterite
		double Dust_Fayalite    = max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[1])); //M_fayalite
		double Dust_Enstatite   = max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[2])); //M_enstatite
		double Dust_Ferrosilite = max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[3])); //M_ferrosilite
		double Dust_Quartz      = max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[4])); //M_quartz
		double Dust_Iron        = max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[5])); //M_iron
		double Dust_SiC		    = max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[8])); //C_SiC
		double Dust_Carbon		= max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[10])); //C_carbon

		Dust_Quartz     += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[6])); //S_quartz
		Dust_Iron		+= max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[7])); //S_iron
		Dust_Iron		+= max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[9])); //C_iron
		
		//Remove total dust created from metallicity-----------------------------------------------------------------
		SumAGBDust = 0.0;
		for (j=0; j<11; j++){
		SumAGBDust += (step_width_times_DiskSFR * NormAGBDustYieldRate_actual[j]);
		}			
		Gal[p].MetalsColdGas.agb -= max(0.0,SumAGBDust);			
		
		
		//Ferrosilite Mg2SiO4 ----------------------------------------
		Gal[p].Dust_elements.Mg += Dust_Forsterite * 0.345504;
		Gal[p].Dust_elements.Si += Dust_Forsterite * 0.199622;
		Gal[p].Dust_elements.O  += Dust_Forsterite * 0.454874;
		
		Gal[p].ColdGas_elements.Mg -= Dust_Forsterite * 0.345504;
		Gal[p].ColdGas_elements.Si -= Dust_Forsterite * 0.199622;
		Gal[p].ColdGas_elements.O  -= Dust_Forsterite * 0.454874;
		
		//Fayalite Fe2SiO4 --------------------------------------------
		Gal[p].Dust_elements.Fe += Dust_Fayalite * 0.548110;
		Gal[p].Dust_elements.Si += Dust_Fayalite * 0.137827;
		Gal[p].Dust_elements.O  += Dust_Fayalite * 0.314063;

		Gal[p].ColdGas_elements.Fe -= Dust_Fayalite * 0.548110;
		Gal[p].ColdGas_elements.Si -= Dust_Fayalite * 0.137827;
		Gal[p].ColdGas_elements.O  -= Dust_Fayalite * 0.314063;
		
		//Enstatite MgSi03 --------------------------------------------
		Gal[p].Dust_elements.Mg += Dust_Enstatite * 0.243050;
		Gal[p].Dust_elements.Si += Dust_Enstatite * 0.279768;
		Gal[p].Dust_elements.O  += Dust_Enstatite * 0.478124;
		
		Gal[p].ColdGas_elements.Mg -= Dust_Enstatite * 0.243050;
		Gal[p].ColdGas_elements.Si -= Dust_Enstatite * 0.279768;
		Gal[p].ColdGas_elements.O  -= Dust_Enstatite * 0.478124;
		
		//Ferrosilite Fe2Si206 ----------------------------------------
		Gal[p].Dust_elements.Fe += Dust_Ferrosilite * 0.423297;
		Gal[p].Dust_elements.Si += Dust_Ferrosilite * 0.212884;
		Gal[p].Dust_elements.O  += Dust_Ferrosilite * 0.363819;

		Gal[p].ColdGas_elements.Fe -= Dust_Ferrosilite * 0.423297;
		Gal[p].ColdGas_elements.Si -= Dust_Ferrosilite * 0.212884;
		Gal[p].ColdGas_elements.O  -= Dust_Ferrosilite * 0.363819;
		
		//Quartz SiO4 -------------------------------------------------
		Gal[p].Dust_elements.Si += Dust_Quartz * 0.305002;
		Gal[p].Dust_elements.O  += Dust_Quartz * 0.694998;
		
		Gal[p].ColdGas_elements.Si -= Dust_Quartz * 0.305002;
		Gal[p].ColdGas_elements.O  -= Dust_Quartz * 0.694998;
		
		//SiC SiC -----------------------------------------------------
		Gal[p].Dust_elements.Si += Dust_SiC * 0.305002;
		Gal[p].Dust_elements.O  += Dust_SiC * 0.694998;

		Gal[p].ColdGas_elements.Si -= Dust_SiC * 0.305002;
		Gal[p].ColdGas_elements.O  -= Dust_SiC * 0.694998;
		
		//Iron Fe -----------------------------------------------------
		Gal[p].Dust_elements.Fe += Dust_Iron * 1.0;
		Gal[p].ColdGas_elements.Fe -= Dust_Iron * 1.0;
		
		//Carbon C ----------------------------------------------------
		Gal[p].Dust_elements.Cb += Dust_Carbon * 1.0;
		Gal[p].ColdGas_elements.Cb -= Dust_Carbon * 1.0;


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
if ((Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas.type2 >0.0)) {
		
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
	
#ifdef FULL_DUST
		Gal[p].DustISM.SNII.Sil += SNII_prevstep_Cold_Si[i] * eta_SNII_Sil * A_Sil_dust/A_Si;
		Gal[p].DustISM.SNII.Fe  += SNII_prevstep_Cold_Fe[i] * eta_SNII_Fe  * A_Fe_dust/A_Fe;
		Gal[p].DustISM.SNII.SiC += SNII_prevstep_Cold_Si[i] * eta_SNII_SiC * A_SiC_dust/A_Si;
		Gal[p].DustISM.SNII.Cb  += SNII_prevstep_Cold_Cb[i] * eta_SNII_Cb  * A_Cb_dust/A_Cb;	
#endif

		//Create dust--------------------------------------------------------------------------------

		double Dust_Silicates = SNII_prevstep_Cold_Si[i] * eta_SNII_Sil * A_Sil_dust/A_Si;
		double Dust_Iron      = SNII_prevstep_Cold_Fe[i] * eta_SNII_Fe  * A_Fe_dust/A_Fe;
		double Dust_SiC	      = SNII_prevstep_Cold_Si[i] * eta_SNII_SiC * A_SiC_dust/A_Si;
		double Dust_Carbon    = SNII_prevstep_Cold_Cb[i] * eta_SNII_Cb  * A_Cb_dust/A_Cb;	


		//Remove total dust created from metallicity-----------------------------------------------------------------

		Gal[p].MetalsColdGas.type2 -= Dust_Silicates/(1.0e10/Hubble_h);
		Gal[p].MetalsColdGas.type2 -= Dust_Iron/(1.0e10/Hubble_h);
		Gal[p].MetalsColdGas.type2 -= Dust_SiC/(1.0e10/Hubble_h);
		Gal[p].MetalsColdGas.type2 -= Dust_Carbon/(1.0e10/Hubble_h);


		//SNII Silicates ---------------------------------------------------------------------
		
		Gal[p].Dust_elements.Si += Dust_Silicates * 0.210432;
		Gal[p].Dust_elements.Mg += Dust_Silicates * 0.091053;
		Gal[p].Dust_elements.Fe += Dust_Silicates * 0.278948;
		Gal[p].Dust_elements.O  += Dust_Silicates * 0.419567;

		Gal[p].ColdGas_elements.Si -= Dust_Silicates * 0.210432;
		Gal[p].ColdGas_elements.Mg -= Dust_Silicates * 0.091053;
		Gal[p].ColdGas_elements.Fe -= Dust_Silicates * 0.278948;
		Gal[p].ColdGas_elements.O  -= Dust_Silicates * 0.419567;
		
		//SNII SiC ---------------------------------------------------------------------------

		Gal[p].Dust_elements.Si += Dust_SiC * 0.305002;
		Gal[p].Dust_elements.O  += Dust_SiC * 0.694998;

		Gal[p].ColdGas_elements.Si -= Dust_SiC * 0.305002;
		Gal[p].ColdGas_elements.O  -= Dust_SiC * 0.694998;

		//SNII Fe ---------------------------------------------------------------------------

		Gal[p].Dust_elements.Fe += Dust_Iron * 1.0;
		Gal[p].ColdGas_elements.Fe -= Dust_Iron * 1.0;

		//SNII Cb ---------------------------------------------------------------------------

		Gal[p].Dust_elements.Cb += Dust_Carbon * 1.0;
		Gal[p].ColdGas_elements.Cb -= Dust_Carbon * 1.0;

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
if ((Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas.type1a >0.0)) {
	float eta_SNIa_Fe  = 0.005;
	float A_Fe_dust  = 55.85;
	float A_Fe = 55.85;

	double Dust_Iron = SNIa_prevstep_Cold_Fe[i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe;
	
#ifdef FULL_DUST
	Gal[p].DustISM.SNIa.Fe  += SNIa_prevstep_Cold_Fe[i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe;
#endif	
	
	Gal[p].Dust_elements.Fe += Dust_Iron * 1.0;
	Gal[p].ColdGas_elements.Fe -= Dust_Iron * 1.0;
	
	Gal[p].MetalsColdGas.type1a -= Dust_Iron/(1.0e10/Hubble_h);

}
#endif //DUST_SNIA

//*****************************************
//Growth of dust inside MC //maybe go last??			
//*****************************************

#ifdef DUST_GROWTH
    if ( (Gal[p].sfh_DiskMass[i] > 0.0) && (metals_total(Gal[p].MetalsColdGas)>0.0) ) {//){// && (Gal[p].MetalsColdGas.type2>0.0) && (Gal[p].MetalsColdGas.agb>0.0) ) {
		float t_acc_0, Xc, Z_sun, Z_coldgas, Z_fraction;
		float Current_Dust, Growth_Fraction;
	
		Z_sun = 0.02;
		t_acc_0 = 15.0/UnitTime_in_Megayears;

		if (Gal[p].ColdGas > 0.00) {
			Z_coldgas = metals_total(Gal[p].MetalsColdGas)/Gal[p].ColdGas;
			Z_fraction = Z_coldgas/Z_sun;
			}
		else {
			Z_fraction = 0.0;
			}

		Xc = 0.5;	// Coldgas fraction

		//Grow Dust
		//Current_Dust = DustMass_Total(Gal[p].DustISM)*0.08;
		//Growth_Fraction = Xc*(dt/t_acc_0)*Z_fraction;
		
		//////////////////////////////////
		
		float K=4.926E-5;   // (units pc^4) / (M_solar ^2)    
		float rmol=pow((K*pow((((Gal[p].GasDiskRadius*1E6)/Hubble_h)/3.0),-4)*(1E10 * (Gal[p].ColdGas/Hubble_h))*((1E10*Gal[p].ColdGas/Hubble_h)+0.4*(1E10*(Gal[p].DiskMass/Hubble_h)))),0.8);
      	float rmolgal=pow((3.44*pow(rmol,-0.506)+4.82*pow(rmol,-1.054)),-1);
      	float H2MassNoh = ((1E10*Gal[p].ColdGas*0.74/Hubble_h)*rmolgal)/(1+rmolgal);
      	float H2Gas = (H2MassNoh/1E10) * Hubble_h;  // Units M_solar/h
      	float H2MassNoh2 = ((Gal[p].ColdGas_elements.H)*rmolgal)/(1+rmolgal);
      	float H2Gas2 = (H2MassNoh2/1E10) * Hubble_h;  // Units M_solar/h
     
     
#ifndef DUST_GROWTH_H2_FRACTION
		Current_Dust = DustMass_Total(Gal[p].DustISM);
		Growth_Fraction = Xc*(dt/t_acc_0)*Z_fraction;   
#else
		Current_Dust = DustMass_Total(Gal[p].DustISM);
		Growth_Fraction = Xc*(dt/t_acc_0)*Z_fraction*(H2Gas2/Gal[p].ColdGas);   

#endif

		
		
		
		//printf("Code = %g\t0.75CG = %g\tRobs = %g\n",Gal[p].ColdGas*1.0E10/Hubble_h, H2MassNoh, H2MassNoh2);
		
		
		
		
		////////////////////////////////
		
		
			//	printf("start ***********\n");
			//	printf("Before Growth\n");
			//	Print_DustMass(Gal[p].DustISM);
			//	printf("***********\n");
		
		
		
		
		
	//	printf("i = %d\tstart *********\n",i);
	//	printf("Current = %g\tGrowth = %g\n", Current_Dust, Growth_Fraction);
	//	printf("%g\t%g\t%g\t%g\n",Xc,dt,t_acc_0,Z_fraction);
	//	Print_DustMass(Gal[p].DustISM);		
	
		Gal[p].DustISM.Growth.Fe += (Current_Dust * Growth_Fraction);
	
	//	Print_DustMass(Gal[p].DustISM);		

		//Remove from metals
		float current_metals_lol = metals_total(Gal[p].MetalsColdGas) * 1.0E10/Hubble_h;
	//	printf("grow 1 %g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",Gal[p].ColdGas_elements.Cb,Gal[p].ColdGas_elements.N ,Gal[p].ColdGas_elements.O,Gal[p].ColdGas_elements.Ne,Gal[p].ColdGas_elements.Mg,Gal[p].ColdGas_elements.Si,Gal[p].ColdGas_elements.S ,Gal[p].ColdGas_elements.Ca,Gal[p].ColdGas_elements.Fe);
		Gal[p].ColdGas_elements.Cb -= ((Gal[p].ColdGas_elements.Cb/current_metals_lol) * Current_Dust * Growth_Fraction);
		Gal[p].ColdGas_elements.N  -= ((Gal[p].ColdGas_elements.N/ current_metals_lol) * Current_Dust * Growth_Fraction);
		Gal[p].ColdGas_elements.O  -= ((Gal[p].ColdGas_elements.O/ current_metals_lol) * Current_Dust * Growth_Fraction);
		Gal[p].ColdGas_elements.Ne -= ((Gal[p].ColdGas_elements.Ne/current_metals_lol) * Current_Dust * Growth_Fraction);
		Gal[p].ColdGas_elements.Mg -= ((Gal[p].ColdGas_elements.Mg/current_metals_lol) * Current_Dust * Growth_Fraction);
		Gal[p].ColdGas_elements.Si -= ((Gal[p].ColdGas_elements.Si/current_metals_lol) * Current_Dust * Growth_Fraction);
		Gal[p].ColdGas_elements.S  -= ((Gal[p].ColdGas_elements.S/ current_metals_lol) * Current_Dust * Growth_Fraction);
		Gal[p].ColdGas_elements.Ca -= ((Gal[p].ColdGas_elements.Ca/current_metals_lol) * Current_Dust * Growth_Fraction);
		Gal[p].ColdGas_elements.Fe -= ((Gal[p].ColdGas_elements.Fe/current_metals_lol) * Current_Dust * Growth_Fraction);
		//printf("grow 2 %g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",Gal[p].ColdGas_elements.Cb,Gal[p].ColdGas_elements.N ,Gal[p].ColdGas_elements.O,Gal[p].ColdGas_elements.Ne,Gal[p].ColdGas_elements.Mg,Gal[p].ColdGas_elements.Si,Gal[p].ColdGas_elements.S ,Gal[p].ColdGas_elements.Ca,Gal[p].ColdGas_elements.Fe);


		//Remove from metallicities
		//printf("1 %g\t%g\t%g\n",Gal[p].MetalsColdGas.agb,Gal[p].MetalsColdGas.type2,Gal[p].MetalsColdGas.type1a);
		
		agb_ratio    = Gal[p].MetalsColdGas.agb/metals_total(Gal[p].MetalsColdGas);
		type2_ratio  = Gal[p].MetalsColdGas.type2/metals_total(Gal[p].MetalsColdGas);
		type1a_ratio = Gal[p].MetalsColdGas.type1a/metals_total(Gal[p].MetalsColdGas);
		
		Gal[p].MetalsColdGas.agb    -= (agb_ratio   * Current_Dust * Growth_Fraction)/(1.0E10/Hubble_h);
		Gal[p].MetalsColdGas.type2  -= (type2_ratio * Current_Dust * Growth_Fraction)/(1.0E10/Hubble_h);
		Gal[p].MetalsColdGas.type1a -= (type1a_ratio* Current_Dust * Growth_Fraction)/(1.0E10/Hubble_h);
		//printf("2 %g\t%g\t%g\n",Gal[p].MetalsColdGas.agb,Gal[p].MetalsColdGas.type2,Gal[p].MetalsColdGas.type1a);

		//printf("end *********\n");
			//	printf("After Growth\n");
			//	Print_DustMass(Gal[p].DustISM);
			//	printf("***********\n");




} 
#endif //DUST_GROWTH

//*****************************************
//Dust destruction			
//*****************************************

#ifdef DUST_DESTRUCTION
    if ( (Gal[p].sfh_DiskMass[i] > 0.0) ){ // && (Gal[p].MetalsColdGas.type2>0.0) && (Gal[p].MetalsColdGas.agb>0.0) ) {
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
			
				//Destroy Dust
				float Current_Dust = DustMass_Total(Gal[p].DustISM);
			
				agb_ratio    = Gal[p].MetalsColdGas.agb/metals_total(Gal[p].MetalsColdGas);
				type2_ratio  = Gal[p].MetalsColdGas.type2/metals_total(Gal[p].MetalsColdGas);
				type1a_ratio = Gal[p].MetalsColdGas.type1a/metals_total(Gal[p].MetalsColdGas);


				Gal[p].DustISM.AGB.Sil -= (Gal[p].DustISM.AGB.Sil)*des_frac;
				Gal[p].DustISM.AGB.SiC -= (Gal[p].DustISM.AGB.SiC)*des_frac;
				Gal[p].DustISM.AGB.Fe  -= (Gal[p].DustISM.AGB.Fe )*des_frac;
				Gal[p].DustISM.AGB.Cb  -= (Gal[p].DustISM.AGB.Cb )*des_frac;			

				Gal[p].DustISM.SNII.Sil -= (Gal[p].DustISM.SNII.Sil)*des_frac;
				Gal[p].DustISM.SNII.SiC -= (Gal[p].DustISM.SNII.SiC)*des_frac;
				Gal[p].DustISM.SNII.Fe  -= (Gal[p].DustISM.SNII.Fe )*des_frac;
				Gal[p].DustISM.SNII.Cb  -= (Gal[p].DustISM.SNII.Cb )*des_frac;			

				Gal[p].DustISM.SNIa.Fe  -= (Gal[p].DustISM.SNIa.Fe )*des_frac;

				Gal[p].DustISM.Growth.Sil -= (Gal[p].DustISM.Growth.Sil)*des_frac;
				Gal[p].DustISM.Growth.SiC -= (Gal[p].DustISM.Growth.SiC)*des_frac;
				Gal[p].DustISM.Growth.Fe  -= (Gal[p].DustISM.Growth.Fe )*des_frac;
				Gal[p].DustISM.Growth.Cb  -= (Gal[p].DustISM.Growth.Cb )*des_frac;			

				//Update metals
				float current_metals_lol = metals_total(Gal[p].MetalsColdGas) * 1.0E10/Hubble_h;
			//	printf("grow 1 %g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",Gal[p].ColdGas_elements.Cb,Gal[p].ColdGas_elements.N ,Gal[p].ColdGas_elements.O,Gal[p].ColdGas_elements.Ne,Gal[p].ColdGas_elements.Mg,Gal[p].ColdGas_elements.Si,Gal[p].ColdGas_elements.S ,Gal[p].ColdGas_elements.Ca,Gal[p].ColdGas_elements.Fe);
				Gal[p].ColdGas_elements.Cb += ((Gal[p].ColdGas_elements.Cb/current_metals_lol) * Current_Dust * des_frac);
				Gal[p].ColdGas_elements.N  += ((Gal[p].ColdGas_elements.N/ current_metals_lol) * Current_Dust * des_frac);
				Gal[p].ColdGas_elements.O  += ((Gal[p].ColdGas_elements.O/ current_metals_lol) * Current_Dust * des_frac);
				Gal[p].ColdGas_elements.Ne += ((Gal[p].ColdGas_elements.Ne/current_metals_lol) * Current_Dust * des_frac);
				Gal[p].ColdGas_elements.Mg += ((Gal[p].ColdGas_elements.Mg/current_metals_lol) * Current_Dust * des_frac);
				Gal[p].ColdGas_elements.Si += ((Gal[p].ColdGas_elements.Si/current_metals_lol) * Current_Dust * des_frac);
				Gal[p].ColdGas_elements.S  += ((Gal[p].ColdGas_elements.S/ current_metals_lol) * Current_Dust * des_frac);
				Gal[p].ColdGas_elements.Ca += ((Gal[p].ColdGas_elements.Ca/current_metals_lol) * Current_Dust * des_frac);
				Gal[p].ColdGas_elements.Fe += ((Gal[p].ColdGas_elements.Fe/current_metals_lol) * Current_Dust * des_frac);

				//Update metallicities
				Gal[p].MetalsColdGas.agb    += (agb_ratio   * Current_Dust * des_frac)/(1.0E10/Hubble_h);
				Gal[p].MetalsColdGas.type2  += (type2_ratio * Current_Dust * des_frac)/(1.0E10/Hubble_h);
				Gal[p].MetalsColdGas.type1a += (type1a_ratio* Current_Dust * des_frac)/(1.0E10/Hubble_h);

			//	printf("After Dest\n");
			//	Print_DustMass(Gal[p].DustISM);
			//	printf("end***********\n");


			} //des_Frac > 1.0
		}		
#endif //DUST_DESTRUCTION

	
//} //if coldgas > 1.0e7
//} //metals > 0.0
	
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

