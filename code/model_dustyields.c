/*
 * model_dustyields.c
 *
 *  Created on: Oct2016
 *  Last modified: Nov 2017
 *      Author: scottclay
 * 
 *  Adds a model of dust production (via AGB stars, SNe remnants and grain growth
 *  molecular clouds) and dust destruction (via SNe shock waves). 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

#ifdef DETAILED_DUST

double calc_dust_growth(double dt, float element, float dust_element, float growthfac) {
    double grown_dust;
    double f_cond0, f_cond;

/* This partially implements PSG17.
 * There are some dubious approximations made in that paper.
 * Our implmentation here is even worse as it assumes a fixed density for the molecular gas
 * (or maybe not actually as Popping's model is completely arbitrary). */
    grown_dust=0.;
    if (element > 0.0) {
        f_cond0 = dust_element / element;
        f_cond = pow((pow((f_cond0 * growthfac),-2.0) + 1.0),-0.5);
        grown_dust = max(f_cond * element - dust_element, 0.);
    }
        
    return grown_dust;

}

#endif

void update_dust_mass(int p, int centralgal, double dt, int nstep, int halonr)
{
	int Zi;
	double timestep_width; //Width of current timestep in CODE UNITS
	int TimeBin; //Bin in Yield arrays corresponding to current timestep
	double Zi_disp, NormAGBDustYieldRate_actual[AGB_DUST_TYPE_NUM];
	double timet, sfh_time;
	double fwind; //Required for metal-rich wind implementation
	double DiskSFR, step_width_times_DiskSFR, DiskSFR_physical_units, step_width_times_DiskSFR_physical_units, inverse_DiskMass_physical_units;
	double Disk_total_metallicity;//, Bulge_total_metallicity, ICM_total_metallicity;
	double TotalMassReturnedToColdDiskGas, TotalMassReturnedToHotGas;

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


//*****************************************
//DUST ENRICHMENT FROM AGB DISK STARS INTO COLD GAS:
//*****************************************

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
	    for (k=0;k<AGB_DUST_TYPE_NUM;k++)	
	    {
	    	NormAGBDustYieldRate_actual[k] = NormAGBDustYieldRate[TimeBin][i][Zi_saved][k] + ((NormAGBDustYieldRate[TimeBin][i][Zi_saved+1][k] - NormAGBDustYieldRate[TimeBin][i][Zi_saved][k])*Zi_disp_saved);	    	
	    }
	    
#ifdef FULL_DUST_RATES        
        Gal[p].DustColdGasRates.AGB += max(0.0,(DiskSFR_physical_units * NormAGBDustYieldRate_actual[0])) /(UnitTime_in_years);
        Gal[p].DustColdGasRates.AGB += max(0.0,(DiskSFR_physical_units * NormAGBDustYieldRate_actual[1])) /(UnitTime_in_years);
        Gal[p].DustColdGasRates.AGB += max(0.0,(DiskSFR_physical_units * NormAGBDustYieldRate_actual[2])) /(UnitTime_in_years);
        Gal[p].DustColdGasRates.AGB += max(0.0,(DiskSFR_physical_units * NormAGBDustYieldRate_actual[3])) /(UnitTime_in_years);
        Gal[p].DustColdGasRates.AGB += max(0.0,(DiskSFR_physical_units * NormAGBDustYieldRate_actual[4])) /(UnitTime_in_years);
        Gal[p].DustColdGasRates.AGB += max(0.0,(DiskSFR_physical_units * NormAGBDustYieldRate_actual[5])) /(UnitTime_in_years);
        Gal[p].DustColdGasRates.AGB += max(0.0,(DiskSFR_physical_units * NormAGBDustYieldRate_actual[6])) /(UnitTime_in_years);
        Gal[p].DustColdGasRates.AGB += max(0.0,(DiskSFR_physical_units * NormAGBDustYieldRate_actual[7])) /(UnitTime_in_years);
        Gal[p].DustColdGasRates.AGB += max(0.0,(DiskSFR_physical_units * NormAGBDustYieldRate_actual[8])) /(UnitTime_in_years);
        Gal[p].DustColdGasRates.AGB += max(0.0,(DiskSFR_physical_units * NormAGBDustYieldRate_actual[9])) /(UnitTime_in_years);
        Gal[p].DustColdGasRates.AGB += max(0.0,(DiskSFR_physical_units * NormAGBDustYieldRate_actual[10]))/(UnitTime_in_years);
#endif				
		//Calculate the amount of dust CREATED ----------------------------------------------------------------------
		//These are calculated based on pre-code calculations in dustyield_integrals.c and then multiplied
		//by the SFR here to get the amount of dust created for each specific type(quartz, iron, carbon etc.)
		//and for 3 types of star (M,C,S). 
		
		double Dust_Forsterite  = max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[0])); //M_forsterite
		double Dust_Fayalite    = max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[1])); //M_fayalite
		double Dust_Enstatite   = max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[2])); //M_enstatite
		double Dust_Ferrosilite = max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[3])); //M_ferrosilite
		double Dust_Quartz      = max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[4])); //M_quartz
		double Dust_Iron        = max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[5])); //M_iron
		double Dust_SiC	        = max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[8])); //C_SiC
		double Dust_Carbon	= max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[10])); //C_carbon

		Dust_Quartz     += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[6])); //S_quartz
		Dust_Iron		+= max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[7])); //S_iron
		Dust_Iron		+= max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[9])); //C_iron
		
		
		//Element Conversion -----------------------------------------------------------------------------------
		//Conversion of dust species (i.e. Ferrosilite) into Actual elements to store
		//in correct arrays (i.e. Ferrosilite -> Mg/Si/O)
		//All the following conversions are done by mass fraction
		
		//Ferrosilite Mg2SiO4 ----------------------------------------
		Gal[p].DustColdGas_elements.Mg += Dust_Forsterite * 0.345504;
		Gal[p].DustColdGas_elements.Si += Dust_Forsterite * 0.199622;
		Gal[p].DustColdGas_elements.O  += Dust_Forsterite * 0.454874;
				
		//Fayalite Fe2SiO4 --------------------------------------------
		Gal[p].DustColdGas_elements.Fe += Dust_Fayalite * 0.548110;
		Gal[p].DustColdGas_elements.Si += Dust_Fayalite * 0.137827;
		Gal[p].DustColdGas_elements.O  += Dust_Fayalite * 0.314063;
		
		//Enstatite MgSi03 --------------------------------------------
		Gal[p].DustColdGas_elements.Mg += Dust_Enstatite * 0.243050;
		Gal[p].DustColdGas_elements.Si += Dust_Enstatite * 0.279768;
		Gal[p].DustColdGas_elements.O  += Dust_Enstatite * 0.478124;
		
		//Ferrosilite Fe2Si206 ----------------------------------------
		Gal[p].DustColdGas_elements.Fe += Dust_Ferrosilite * 0.423297;
		Gal[p].DustColdGas_elements.Si += Dust_Ferrosilite * 0.212884;
		Gal[p].DustColdGas_elements.O  += Dust_Ferrosilite * 0.363819;
		
		//Quartz SiO4 -------------------------------------------------
		Gal[p].DustColdGas_elements.Si += Dust_Quartz * 0.305002;
		Gal[p].DustColdGas_elements.O  += Dust_Quartz * 0.694998;
		
		//SiC SiC -----------------------------------------------------
		Gal[p].DustColdGas_elements.Si += Dust_SiC * 0.700453;
		Gal[p].DustColdGas_elements.Cb += Dust_SiC * 0.299547;

		//Iron Fe -----------------------------------------------------
		Gal[p].DustColdGas_elements.Fe += Dust_Iron * 1.0;
		
		//Carbon C ----------------------------------------------------
		Gal[p].DustColdGas_elements.Cb += Dust_Carbon * 1.0;
		
	} //if sfh_DM >0
    

#endif //DUST_AGB

//*****************************************
//DUST ENRICHMENT FROM SNII FROM DISK STARS INTO COLD GAS:
//*****************************************

#ifdef DUST_SNII
	if ((Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas.type2 >0.0)) {
	
		//eta (dust condensation eff.) and the atomic weights A_x are taken from Zhukovska2008	
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

	
	#ifdef FULL_DUST_RATES
			Gal[p].DustColdGasRates.SNII += (SNII_prevstep_Cold_Si[i] * eta_SNII_Sil * A_Sil_dust/A_Si)/(dt * UnitTime_in_years);
			Gal[p].DustColdGasRates.SNII += (SNII_prevstep_Cold_Fe[i] * eta_SNII_Fe  * A_Fe_dust/A_Fe )/(dt * UnitTime_in_years);
			Gal[p].DustColdGasRates.SNII += (SNII_prevstep_Cold_Si[i] * eta_SNII_SiC * A_SiC_dust/A_Si)/(dt * UnitTime_in_years);
			Gal[p].DustColdGasRates.SNII += (SNII_prevstep_Cold_Cb[i] * eta_SNII_Cb  * A_Cb_dust/A_Cb) /(dt * UnitTime_in_years);	
	#endif

			//Create dust (based on the prescription of Zhukovska2008)---------------------------
			//SNII_prevstep_x is calculated in model_yields.c 
			//It is the amount of a specific metal (i.e. Si) produced in the last timestep

			double Dust_Silicates = SNII_prevstep_Cold_Si[i] * eta_SNII_Sil * A_Sil_dust/A_Si;
			double Dust_Iron      = SNII_prevstep_Cold_Fe[i] * eta_SNII_Fe  * A_Fe_dust/A_Fe;
			double Dust_SiC	      = SNII_prevstep_Cold_Si[i] * eta_SNII_SiC * A_SiC_dust/A_Si;
			double Dust_Carbon    = SNII_prevstep_Cold_Cb[i] * eta_SNII_Cb  * A_Cb_dust/A_Cb;	

			//Element conversion -----------------------------------------------------------------
			//Conversion of dust species (i.e. Silicates) into Actual elements to store
			//in correct arrays (i.e. Silicates -> Mg/Si/Fe/O)
			//All the following conversions are done by mass fraction

			//SNII Silicates -------------------
		
			Gal[p].DustColdGas_elements.Si += Dust_Silicates * 0.210432;
			Gal[p].DustColdGas_elements.Mg += Dust_Silicates * 0.091053;
			Gal[p].DustColdGas_elements.Fe += Dust_Silicates * 0.278948;
			Gal[p].DustColdGas_elements.O  += Dust_Silicates * 0.419567;

			//SNII SiC --------------------------

			Gal[p].DustColdGas_elements.Si += Dust_SiC * 0.299547;
			Gal[p].DustColdGas_elements.Cb  += Dust_SiC * 0.700453;

			//SNII Fe -------------------------

			Gal[p].DustColdGas_elements.Fe += Dust_Iron * 1.0;

			//SNII Cb ------------------------

			Gal[p].DustColdGas_elements.Cb += Dust_Carbon * 1.0;

	}//if sfh_DM >0
#endif //DUST_SNII
	
	
//*****************************************
//DUST ENRICHMENT FROM SNIA FROM DISK STARS INTO COLD GAS:
//*****************************************

#ifdef DUST_SNIA		
	if ((Gal[p].sfh_DiskMass[i] > 0.0) && (Gal[p].MetalsColdGas.type1a >0.0)) {
		
		//eta and A_x taken from Zhukovska2008
		float eta_SNIa_Fe  = 0.005; //dust condensation eff.
		float A_Fe_dust  = 55.85; //atomic weight
		float A_Fe = 55.85; //atomic weight

		double Dust_Iron = SNIa_prevstep_Cold_Fe[i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe;
	
#ifdef FULL_DUST_RATES		
		Gal[p].DustColdGasRates.SNIA  += (SNIa_prevstep_Cold_Fe[i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe)/(dt * UnitTime_in_years);
#endif	
	
		Gal[p].DustColdGas_elements.Fe += Dust_Iron * 1.0;
	}//if sfh_DM >0
#endif //DUST_SNIA

} //loop over SFH bins


//*****************************************
//Dust grain growth inside molecular clouds 
//*****************************************

#ifdef DUST_GROWTH
    if ((metals_total(Gal[p].MetalsColdGas)>0.0) && (Gal[p].ColdGas > 0.00)) {
		float Z_sun, Z_coldgas, Z_fraction;
	
		Z_sun = 0.02; //Solar metallicity
		
		//Calculate cold gas metallicity
		if (Gal[p].ColdGas > 0.00) {
			Z_coldgas = metals_total(Gal[p].MetalsColdGas)/Gal[p].ColdGas;
			Z_fraction = Z_coldgas/Z_sun;
			}
		else {
			Z_fraction = 0.0;
			}
			
		//Dust growth inside MCs requires an approximation of H2 gas content --------------
		//This is taken from the Martindale2017 paper 
		//(with code provided by Hazel)

		float K=4.926E-5;   // (units pc^4) / (M_solar ^2)    
		float rmol=pow((K*pow((((Gal[p].GasDiskRadius*1E6)/Hubble_h)/3.0),-4)*(1E10 * (Gal[p].ColdGas/Hubble_h))*((1E10*Gal[p].ColdGas/Hubble_h)+0.4*(1E10*(Gal[p].DiskMass/Hubble_h)))),0.8);
		float rmolgal=pow((3.44*pow(rmol,-0.506)+4.82*pow(rmol,-1.054)),-1);
		float H2MassNoh2 = ((Gal[p].ColdGas_elements.H)*rmolgal)/(1+rmolgal);
		float H2Gas2 = (H2MassNoh2/1E10) * Hubble_h;  // Units M_solar/h

		//------------------------------------------------------
		
		//calculate exchange timescales and molecular gas fractions
		float t_exchange, t_exchange_eff;
		float f_mol;
		
		t_exchange = DUST_TEXCH; //Exchange timescale
		f_mol = H2Gas2/Gal[p].ColdGas; //Molecular gas fraction 
		t_exchange_eff = t_exchange * ((1 - f_mol)/f_mol); //Effective exchange timescale


		// Calculate amount of dust grown of each element --------------------------------------------------------------------
		double Dust_Total_old = Gal[p].ColdGas_elements.Cb + Gal[p].ColdGas_elements.O  + Gal[p].ColdGas_elements.N  +
					Gal[p].ColdGas_elements.Ne + Gal[p].ColdGas_elements.Mg + Gal[p].ColdGas_elements.Si +
		                        Gal[p].ColdGas_elements.S  + Gal[p].ColdGas_elements.Ca + Gal[p].ColdGas_elements.Fe;

		if (Gal[p].ColdGas > 0.0) {
		    // The following units conversion should be given by parameters, but for now hard-wire in
		    double coldgas = Gal[p].ColdGas*1e10/Hubble_h;
		    double t_acc_0 = DUST_TACC0;
		    double t_exchange = DUST_TEXCH;
		    double t_acc = t_acc_0 * coldgas / (
			Gal[p].ColdGas_elements.Cb + 
			Gal[p].ColdGas_elements.O + 
			Gal[p].ColdGas_elements.N + 
			// Gal[p].ColdGas_elements.Ne + 
			Gal[p].ColdGas_elements.Mg + 
			Gal[p].ColdGas_elements.Si + 
			Gal[p].ColdGas_elements.S + 
			Gal[p].ColdGas_elements.Ca +
			Gal[p].ColdGas_elements.Fe
			) ;
		    double growthfac = (1.0 + t_exchange/t_acc);
		    double ratefac = dt*UnitTime_in_years/t_exchange_eff;
		    Gal[p].DustColdGas_elements.Cb += calc_dust_growth(dt, Gal[p].ColdGas_elements.Cb, Gal[p].DustColdGas_elements.Cb, growthfac) * ratefac;
		    Gal[p].DustColdGas_elements.O  += calc_dust_growth(dt, Gal[p].ColdGas_elements.O,  Gal[p].DustColdGas_elements.O,  growthfac) * ratefac;
		    Gal[p].DustColdGas_elements.N  += calc_dust_growth(dt, Gal[p].ColdGas_elements.N,  Gal[p].DustColdGas_elements.N,  growthfac) * ratefac;
		    //Gal[p].DustColdGas_elements.Ne += calc_dust_growth(dt, Gal[p].ColdGas_elements.Ne, Gal[p].DustColdGas_elements.Ne, growthfac) * ratefac;
		    Gal[p].DustColdGas_elements.Mg += calc_dust_growth(dt, Gal[p].ColdGas_elements.Mg, Gal[p].DustColdGas_elements.Mg, growthfac) * ratefac;
		    Gal[p].DustColdGas_elements.Si += calc_dust_growth(dt, Gal[p].ColdGas_elements.Si, Gal[p].DustColdGas_elements.Si, growthfac) * ratefac;
		    Gal[p].DustColdGas_elements.S  += calc_dust_growth(dt, Gal[p].ColdGas_elements.S,  Gal[p].DustColdGas_elements.S,  growthfac) * ratefac;
		    Gal[p].DustColdGas_elements.Ca += calc_dust_growth(dt, Gal[p].ColdGas_elements.Ca, Gal[p].DustColdGas_elements.Ca, growthfac) * ratefac;
		    Gal[p].DustColdGas_elements.Fe += calc_dust_growth(dt, Gal[p].ColdGas_elements.Fe, Gal[p].DustColdGas_elements.Fe, growthfac) * ratefac;
		}

		// Check for over-production of dust
		int dustFlag = 0;
		if (Gal[p].DustColdGas_elements.Cb > Gal[p].ColdGas_elements.Cb) {
		    Gal[p].DustColdGas_elements.Cb = Gal[p].ColdGas_elements.Cb;
		    dustFlag=1;
		}
		if (Gal[p].DustColdGas_elements.N > Gal[p].ColdGas_elements.N) {
		    Gal[p].DustColdGas_elements.N = Gal[p].ColdGas_elements.N;
		    dustFlag=1;
		}
		if (Gal[p].DustColdGas_elements.O > Gal[p].ColdGas_elements.O) {
		    Gal[p].DustColdGas_elements.O = Gal[p].ColdGas_elements.O;
		    dustFlag=1;
		}
		if (Gal[p].DustColdGas_elements.Ne > Gal[p].ColdGas_elements.Ne) {
		    Gal[p].DustColdGas_elements.Ne = Gal[p].ColdGas_elements.Ne;
		    dustFlag=1;
		}
		if (Gal[p].DustColdGas_elements.Mg > Gal[p].ColdGas_elements.Mg) {
		    Gal[p].DustColdGas_elements.Mg = Gal[p].ColdGas_elements.Mg;
		    dustFlag=1;
		}
		if (Gal[p].DustColdGas_elements.Si > Gal[p].ColdGas_elements.Si) {
		    Gal[p].DustColdGas_elements.Si = Gal[p].ColdGas_elements.Si;
		    dustFlag=1;
		}
		if (Gal[p].DustColdGas_elements.S > Gal[p].ColdGas_elements.S) {
		    Gal[p].DustColdGas_elements.S = Gal[p].ColdGas_elements.S;
		    dustFlag=1;
		}
		if (Gal[p].DustColdGas_elements.Ca > Gal[p].ColdGas_elements.Ca) {
		    Gal[p].DustColdGas_elements.Ca = Gal[p].ColdGas_elements.Ca;
		    dustFlag=1;
		}
		if (Gal[p].DustColdGas_elements.Fe > Gal[p].ColdGas_elements.Fe) {
		    Gal[p].DustColdGas_elements.Fe = Gal[p].ColdGas_elements.Fe;
		    dustFlag=1;
		}
		if (dustFlag) printf("Excess dust produced\n");
		
		
#ifdef FULL_DUST_RATES
		// Looks like we are accumulating to average over a snap?  WOuld instantaneous be better?
		Gal[p].DustColdGasRates.GROW += (Gal[p].ColdGas_elements.Cb + Gal[p].ColdGas_elements.O  + Gal[p].ColdGas_elements.N  +
					     Gal[p].ColdGas_elements.Ne + Gal[p].ColdGas_elements.Mg + Gal[p].ColdGas_elements.Si +
					     Gal[p].ColdGas_elements.S  + Gal[p].ColdGas_elements.Ca + Gal[p].ColdGas_elements.Fe -
					     Dust_Total_old) / (dt * UnitTime_in_years);
#endif

}
#endif //DUST_GROWTH

//**************************************************
//Dust grain destruction from supernovae shock waves
//**************************************************

#ifdef DUST_DESTRUCTION
    if ((metals_total(Gal[p].MetalsColdGas)>0.0) ) {
		//For dust destruction we follow the prescription of McKee1989.
		float t_des, M_cleared, f_SN;
		float des_frac; 
		M_cleared = 1000.0; //Msol //Mass of gas that is cleared of dust by an average SNe
		f_SN = 0.36; //Dimensionless
		DiskSFR = Gal[p].Sfr;
		float R_SN_IMF = 0.2545/19.87; //Rate of supernovae as determined from the IMF
		float R_SN = R_SN_IMF * DiskSFR * (1.0E10/Hubble_h) * (1/UnitTime_in_years); //Actual rate of SNe

		//Calculate destruction timescale and destruction fraction
		if(R_SN>0.0) {
			t_des = (Gal[p].ColdGas*(1.0e10/Hubble_h))/(M_cleared * f_SN * R_SN);
			des_frac = (dt*UnitTime_in_years/t_des);
		}
		else {
			t_des = 0.0;
			des_frac = 0.0;
		}				
		
		//Calculate destroyed dust ---------------------------------------------------------------------------
		double Dust_Cb =(Gal[p].DustColdGas_elements.Cb * des_frac);
		double Dust_N  =(Gal[p].DustColdGas_elements.N  * des_frac);
		double Dust_O  =(Gal[p].DustColdGas_elements.O  * des_frac);
		double Dust_Ne =(Gal[p].DustColdGas_elements.Ne * des_frac);
		double Dust_Mg =(Gal[p].DustColdGas_elements.Mg * des_frac);
		double Dust_Si =(Gal[p].DustColdGas_elements.Si * des_frac);
		double Dust_S  =(Gal[p].DustColdGas_elements.S  * des_frac);
		double Dust_Ca =(Gal[p].DustColdGas_elements.Ca * des_frac);
		double Dust_Fe =(Gal[p].DustColdGas_elements.Fe * des_frac);
		double Dust_Total = Dust_Cb+Dust_N+Dust_O+Dust_Ne+Dust_Mg+Dust_Si+Dust_S+Dust_Ca+Dust_Fe;
		//double Dust_Total2 = min(Dust_Cb,Gal[p].DustColdGas_elements.Cb)+min(Dust_N,Gal[p].DustColdGas_elements.N)+min(Dust_O,Gal[p].DustColdGas_elements.O)+min(Dust_Ne,Gal[p].DustColdGas_elements.Ne)+min(Dust_Mg,Gal[p].DustColdGas_elements.Mg)+min(Dust_Si,Gal[p].DustColdGas_elements.Si)+min(Dust_S,Gal[p].DustColdGas_elements.S)+min(Dust_Ca,Gal[p].DustColdGas_elements.Ca)+min(Dust_Fe,Gal[p].DustColdGas_elements.Fe);
		
		//Remove destroyed dust to array---------------------------------------------------------------
		
		Gal[p].DustColdGas_elements.Cb -= min(Dust_Cb,Gal[p].DustColdGas_elements.Cb);
		Gal[p].DustColdGas_elements.N  -= min(Dust_N,Gal[p].DustColdGas_elements.N);
		Gal[p].DustColdGas_elements.O  -= min(Dust_O,Gal[p].DustColdGas_elements.O);
		Gal[p].DustColdGas_elements.Ne -= min(Dust_Ne,Gal[p].DustColdGas_elements.Ne);
		Gal[p].DustColdGas_elements.Mg -= min(Dust_Mg,Gal[p].DustColdGas_elements.Mg);
		Gal[p].DustColdGas_elements.Si -= min(Dust_Si,Gal[p].DustColdGas_elements.Si);
		Gal[p].DustColdGas_elements.S  -= min(Dust_S,Gal[p].DustColdGas_elements.S); 
		Gal[p].DustColdGas_elements.Ca -= min(Dust_Ca,Gal[p].DustColdGas_elements.Ca);
		Gal[p].DustColdGas_elements.Fe -= min(Dust_Fe,Gal[p].DustColdGas_elements.Fe);
		
#ifdef FULL_DUST_RATES		
		Gal[p].DustColdGasRates.DEST += Dust_Total/(dt * UnitTime_in_years);
#endif

		}	
#endif //DUST_DESTRUCTION

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
