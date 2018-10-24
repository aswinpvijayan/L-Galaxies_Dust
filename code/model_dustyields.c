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

/*

Note: 
* Dust in contained in metals is contained in gas.
* ColdGasDiff + ColdGasClouds = ColdGas
* mu = ColdGasClouds/ColdGas
* fI = fraction of metals locked up in dust in the diffuse medium
     = DustColdGasDiff/MetalsColdGasDiff
* fC = fraction of metals locked up in dust in cold clouds
     = DustColdGasClouds/MetalsColdGasClouds
* Any dust production mechanism except for growth in molecular clouds enriches the diffuse medium. 
* When SNe destroys dust, a fraction of dust in both the diffuse and the molecular medium is destroyed.
* There is also a slow sputtering of dust (by cosmic rays and hot gas) that occurs in the diffuse medium.
* Also we donot include Nitrogen, Neon and Sulphur in our dust growth model.
*/

// Differential equations for dust growth inside and outside molecular clouds.
float dfC_dt(float fC, float fI, float tacc, float texch, float mu) {
     /* fC, fI - dust fraction inside and outside.
     * tacc - growth timescale (depends upon dust fraction, i.e. changes as fC, fI grow)
     * texch - exchnage timescale of clouds with surroundings
     * mu - mass fraction in molecular clouds. */
    
    return (1.-fC)/tacc + (fI-fC)/texch;
}
float dfI_dt(float fC, float fI, float tacc, float texch, float mu) {
    /* Differential equations for dust growth inside and outside molecular clouds.
     * fC, fI - dust fraction inside and outside.
     * tacc - growth timescale (depends upon dust fraction, i.e. changes as fC, fI grow)
     * texch - exchnage timescale of clouds with surroundings
     * mu - mass fraction in molecular clouds. */
    
    return (fC-fI)*mu/(texch*(1.-mu));
}


//Approximations for the differential equations presented above for the growth of 
//dust. 
double fitfc(float fC0, float fI0, float t, float mu, float t_e) {
    
    float t_e_eff=(1.0-mu)*t_e/mu;
    float t_a_eff=(mu+t_e)/(1+t_e);
    float accFac=1.0-exp(-t/t_a_eff);
    float fCequil=(fI0+t_e)/(1.0+t_e);
    float fCexchFac=1.0-exp(-t/(1.0/mu+t_e_eff));
    float fC=fC0+(fCequil-fC0)*accFac+(1.0-fCequil)*fCexchFac;
    return fC;
}

double fitfi(float fC0, float fI0, float t, float mu, float t_e) {
    
    float t_e_eff=(1.0-mu)*t_e/mu;
    float fIexchFac=1.0-exp(-t/(1.0-mu+t_e_eff));
    float fI=fI0+(fC0-fI0)*fIexchFac;
    return fI;

}

void update_fractions(float dt_bin, float t_acc, float t_exch, float mu, int p, int n) {
    
    float fc, acc_fac, exch_fac, fmean, fdiff;
    
    acc_fac=1.-exp(-dt_bin/t_acc);
    exch_fac=exp(-dt_bin/((1.-mu)*t_exch));
    fc = Gal[p].f_c[n] + (Gal[p].f_cmax[n] - Gal[p].f_c[n])*acc_fac;
    Gal[p].f_c[n] = fc;
    fmean = mu*Gal[p].f_c[n]+(1. - mu)*Gal[p].f_i[n];
    fdiff = (Gal[p].f_c[n] - Gal[p].f_i[n])*exch_fac;
    Gal[p].f_c[n] = fmean + (1. - mu)*fdiff;
    Gal[p].f_i[n] = fmean - mu*fdiff;
    
    return;
}

void drop_fnan(int n, int p) {

    if (isnan(Gal[p].f_i[n])) {Gal[p].f_i[n] = 0.;}
    if (isnan(Gal[p].f_c[n])) {Gal[p].f_c[n] = 0.;}

    return;
}

void drop_dustnan(int p) {

    if (isnan(Gal[p].DustColdGasClouds_elements.Cb)) {Gal[p].DustColdGasClouds_elements.Cb=0.;}
    if (isnan(Gal[p].DustColdGasClouds_elements.O)) {Gal[p].DustColdGasClouds_elements.O=0.;}
    if (isnan(Gal[p].DustColdGasClouds_elements.Mg)) {Gal[p].DustColdGasClouds_elements.Mg=0.;}
    if (isnan(Gal[p].DustColdGasClouds_elements.Si)) {Gal[p].DustColdGasClouds_elements.Si=0.;}
    if (isnan(Gal[p].DustColdGasClouds_elements.S)) {Gal[p].DustColdGasClouds_elements.S=0.;}
    if (isnan(Gal[p].DustColdGasClouds_elements.Ca)) {Gal[p].DustColdGasClouds_elements.Ca=0.;}
    if (isnan(Gal[p].DustColdGasClouds_elements.Fe)) {Gal[p].DustColdGasClouds_elements.Fe=0.;}
    
    if (isnan(Gal[p].DustColdGasDiff_elements.Cb)) {Gal[p].DustColdGasDiff_elements.Cb=0.;}
    if (isnan(Gal[p].DustColdGasDiff_elements.O)) {Gal[p].DustColdGasDiff_elements.O=0.;}
    if (isnan(Gal[p].DustColdGasDiff_elements.Mg)) {Gal[p].DustColdGasDiff_elements.Mg=0.;}
    if (isnan(Gal[p].DustColdGasDiff_elements.Si)) {Gal[p].DustColdGasDiff_elements.Si=0.;}
    if (isnan(Gal[p].DustColdGasDiff_elements.S)) {Gal[p].DustColdGasDiff_elements.S=0.;}
    if (isnan(Gal[p].DustColdGasDiff_elements.Ca)) {Gal[p].DustColdGasDiff_elements.Ca=0.;}
    if (isnan(Gal[p].DustColdGasDiff_elements.Fe)) {Gal[p].DustColdGasDiff_elements.Fe=0.;}
    
    return;
}

/*
double fitfc_new(float fC0, float fI0, float t, float mu, float t_e) {

    float t_a_eff=min(1.0+1.0/t_e,1.0/mu);
    float accFac=fC0/(fC0+(1.0-fC0)*exp(-t/t_a_eff));
    float t_e_eff=(1.0-mu)*t_e/mu;
    float fCexchFac=1.0-exp(-t/((1.0+mu)/mu+t_e_eff));
    float fCequil=max(1.0-(1.0-mu)*t_e,(t_e-1.0+sqrt((t_e-1.0)*(t_e-1.0)+4.0*fI0*t_e))/(2.0*t_e));
    float X=1.0-fC0;
    float Y=fCequil-fC0;
    float xX=(accFac-fC0)/X;
    float fC=fC0+Y*xX +(1-fCequil)*fCexchFac;
    return fC;
}


double fitfi_new(float fC0, float fI0, float t, float mu, float t_e) {

    float t_e_eff=(1.0-mu)*t_e/mu;
    float fIexchFac=1.0-exp(-t/(1.0/mu+t_e_eff));
    float fI=fI0+(fC0-fI0)*fIexchFac;
    return fI;
}
*/

void update_dust(int p) {

    Gal[p].DustColdGasClouds_elements.Cb = Gal[p].f_c[0]*Gal[p].ColdGasClouds_elements.Cb;
    //Gal[p].DustColdGasClouds_elements.N = Gal[p].f_c[1]*Gal[p].ColdGasClouds_elements.N;  
    Gal[p].DustColdGasClouds_elements.O = Gal[p].f_c[2]*Gal[p].ColdGasClouds_elements.O;
    //Gal[p].DustColdGasClouds_elements.Ne = Gal[p].f_c[3]*Gal[p].ColdGasClouds_elements.Ne;  
    Gal[p].DustColdGasClouds_elements.Mg = Gal[p].f_c[4]*Gal[p].ColdGasClouds_elements.Mg; 
    Gal[p].DustColdGasClouds_elements.Si = Gal[p].f_c[5]*Gal[p].ColdGasClouds_elements.Si;
    Gal[p].DustColdGasClouds_elements.S = Gal[p].f_c[6]*Gal[p].ColdGasClouds_elements.S;
    Gal[p].DustColdGasClouds_elements.Ca = Gal[p].f_c[7]*Gal[p].ColdGasClouds_elements.Ca;
    Gal[p].DustColdGasClouds_elements.Fe = Gal[p].f_c[8]*Gal[p].ColdGasClouds_elements.Fe;
    
    return;
}


#ifdef DETAILED_DUST

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
    
    float mu_gas;
	
	double previoustime, newtime, deltaT;
	
	previoustime = NumToTime(Gal[p].SnapNum - 1);
	newtime = NumToTime(Gal[p].SnapNum);
	deltaT = previoustime - newtime;
	
	
	mu_gas = Gal[p].mu_gas;
    
	TotalMassReturnedToColdDiskGas=0.0;
	TotalMassReturnedToHotGas=0.0;
	
	mass_checks("Dust at beginning of dust growth",p);
	
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
		t_des = (Gal[p].ColdGas*(1.0e10/Hubble_h))/(M_cleared * f_SN * R_SN);
		float survive_frac = exp(-dt*UnitTime_in_years/t_des);
        
        #ifdef FULL_DUST_RATES	
	        //Total dust before dust destruction
	        double dust_total = elements_total(Gal[p].DustColdGasDiff_elements) + elements_total(Gal[p].DustColdGasClouds_elements);
        #endif
		//We assume that the SNR will destroy equal amounts of dust in cold clouds and 
		//the diffuse medium, but all those will end up as diffuse gas, I guess.  
		//Then some will be reaccreted onto cold clouds.
	    //Simplest approximation is just to destroy the same fraction in each.
		
		Gal[p].DustColdGasDiff_elements=elements_add(elements_init(),Gal[p].DustColdGasDiff_elements,survive_frac);
	    Gal[p].DustColdGasClouds_elements=elements_add(elements_init(),Gal[p].DustColdGasClouds_elements,survive_frac);
		
		#ifdef DCR_Dest
		    //****************************
		    //Destruction of dust from miscellaneous processes, only done in the diffuse medium
		    //****************************
		    //Timescale of destruction chosen to be 1Gyr (arbitrary)
		    float CR_timescale = 1e9;
	        survive_frac = exp(-(dt*UnitTime_in_years/CR_timescale));
		    
		    Gal[p].DustColdGasDiff_elements=elements_add(elements_init(),Gal[p].DustColdGasDiff_elements,survive_frac);
	    #endif    
		
		#ifdef DDestHIIregion
		    //****************************
		    //Dust destruction in the diffused medium by evaporation in HII regions
		    //****************************
		    //This is a rough implementation based on the mass of HII regions created by OB stars
		    //as presented in Tielens' book on ISM. Need to get better values of the photon flux for
		    //the star of average mass produced in the particular time step. Current implementation
		    //has only negligible effect on the dust mass.
		    float MetClouds_tot = elements_total(Gal[p].ColdGasClouds_elements);
		    
	        float MetDiff_tot = elements_total(Gal[p].ColdGasDiff_elements);
		    
		    double sfr = Gal[p].Sfr * UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
		    double gas_cleared = (0.0496173*sfr*dt*UnitTime_in_years/13.7054)*2.24;
		    double frac_diff = gas_cleared/MetDiff_tot;
		    double frac_clouds = 0.;
		    
		    if (frac_diff > 1.) {
		    
		        frac_diff = 1.;
		        frac_clouds = (gas_cleared-MetDiff_tot)/MetClouds_tot;
	        }
	        
            Gal[p].DustColdGasDiff_elements=elements_add(Gal[p].DustColdGasDiff_elements,Gal[p].DustColdGasDiff_elements,-frac_diff);
		    Gal[p].DustColdGasClouds_elements=elements_add(Gal[p].DustColdGasDiff_elements,Gal[p].DustColdGasDiff_elements,-frac_clouds);
		#endif    
		
		
        #ifdef FULL_DUST_RATES	
	    Gal[p].DustColdGasRates.DEST += (dust_total - elements_total(Gal[p].DustColdGasDiff_elements)					                    - elements_total(Gal[p].DustColdGasClouds_elements))/(deltaT * UnitTime_in_years);
        #endif
        mass_checks("Dust from destruction",p);
		}	
#endif //DUST_DESTRUCTION



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
	    	NormAGBDustYieldRate_actual[k] = NormAGBDustYieldRate[TimeBin][i][Zi_saved[i]][k] + ((NormAGBDustYieldRate[TimeBin][i][Zi_saved[i]+1][k] - NormAGBDustYieldRate[TimeBin][i][Zi_saved[i]][k])*Zi_disp_saved[i]);	    	
	    }

	    
    #ifdef FULL_DUST_RATES        
        for (int l=0; l<11; l++) {
		    //Gal[p].DustColdGasRates.AGB +=  NormAGBDustYieldRate_actual[l] * DiskSFR_physical_units / (UnitTime_in_years);
		    Gal[p].DustColdGasRates.AGB +=  NormAGBDustYieldRate_actual[l] * DiskSFR_physical_units*dt / (deltaT*UnitTime_in_years);
	    }
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
		double Dust_SiC		    = max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[8])); //C_SiC
		double Dust_Carbon		= max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[10])); //C_carbon

		Dust_Quartz     += max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[6])); //S_quartz
		Dust_Iron		+= max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[7])); //S_iron
		Dust_Iron		+= max(0.0,(step_width_times_DiskSFR_physical_units * NormAGBDustYieldRate_actual[9])); //C_iron
		
		//Element Conversion -----------------------------------------------------------------------------------
		//Conversion of dust species (i.e. Ferrosilite) into Actual elements to store
		//in correct arrays (i.e. Ferrosilite -> Mg/Si/O)
		//All the following conversions are done by mass fraction
		//Corrections added at places to avoid dust masses to go beyond their mass in metals
		
		//Ferrosilite Mg2SiO4 ----------------------------------------
		Gal[p].DustColdGasDiff_elements.Mg += min(Dust_Forsterite * 0.345504 * (1.0 - mu_gas), Gal[p].ColdGasDiff_elements.Mg - Gal[p].DustColdGasDiff_elements.Mg);
		Gal[p].DustColdGasDiff_elements.Si += Dust_Forsterite * 0.199622 * (1.0 - mu_gas);
		Gal[p].DustColdGasDiff_elements.O  += Dust_Forsterite * 0.454874 * (1.0 - mu_gas);
		
		Gal[p].DustColdGasClouds_elements.Mg += min(Dust_Forsterite * 0.345504 * mu_gas, Gal[p].ColdGasClouds_elements.Mg - Gal[p].DustColdGasClouds_elements.Mg);
		Gal[p].DustColdGasClouds_elements.Si += Dust_Forsterite * 0.199622 * mu_gas;
		Gal[p].DustColdGasClouds_elements.O  += Dust_Forsterite * 0.454874 * mu_gas;
				
		//Fayalite Fe2SiO4 --------------------------------------------
		Gal[p].DustColdGasDiff_elements.Fe += min(Dust_Fayalite * 0.548110 * (1.0 - mu_gas), Gal[p].ColdGasDiff_elements.Fe - Gal[p].DustColdGasDiff_elements.Fe);
		Gal[p].DustColdGasDiff_elements.Si += Dust_Fayalite * 0.137827 * (1.0 - mu_gas);
		Gal[p].DustColdGasDiff_elements.O  += Dust_Fayalite * 0.314063 * (1.0 - mu_gas);
		
		Gal[p].DustColdGasClouds_elements.Fe += min(Dust_Fayalite * 0.548110 * mu_gas, Gal[p].ColdGasClouds_elements.Fe - Gal[p].DustColdGasClouds_elements.Fe);
		Gal[p].DustColdGasClouds_elements.Si += Dust_Fayalite * 0.137827 * mu_gas;
		Gal[p].DustColdGasClouds_elements.O  += Dust_Fayalite * 0.314063 * mu_gas;
		
		//Enstatite MgSi03 --------------------------------------------
		Gal[p].DustColdGasDiff_elements.Mg += min(Dust_Enstatite * 0.243050 * (1.0 - mu_gas), Gal[p].ColdGasDiff_elements.Mg - Gal[p].DustColdGasDiff_elements.Mg);
		Gal[p].DustColdGasDiff_elements.Si += Dust_Enstatite * 0.279768 * (1.0 - mu_gas);
		Gal[p].DustColdGasDiff_elements.O  += Dust_Enstatite * 0.478124 * (1.0 - mu_gas);
		
		Gal[p].DustColdGasClouds_elements.Mg += min(Dust_Enstatite * 0.243050 * mu_gas, Gal[p].ColdGasClouds_elements.Mg - Gal[p].DustColdGasClouds_elements.Mg);
		Gal[p].DustColdGasClouds_elements.Si += Dust_Enstatite * 0.279768 * mu_gas;
		Gal[p].DustColdGasClouds_elements.O  += Dust_Enstatite * 0.478124 * mu_gas;
		
		//Ferrosilite Fe2Si206 ----------------------------------------
		Gal[p].DustColdGasDiff_elements.Fe += min(Dust_Ferrosilite * 0.423297 * (1.0 - mu_gas), Gal[p].ColdGasDiff_elements.Fe - Gal[p].DustColdGasDiff_elements.Fe);
		Gal[p].DustColdGasDiff_elements.Si += Dust_Ferrosilite * 0.212884 * (1.0 - mu_gas);
		Gal[p].DustColdGasDiff_elements.O  += Dust_Ferrosilite * 0.363819 * (1.0 - mu_gas);
		
		Gal[p].DustColdGasClouds_elements.Fe += min(Dust_Ferrosilite * 0.423297 * mu_gas, Gal[p].ColdGasClouds_elements.Fe - Gal[p].DustColdGasClouds_elements.Fe);
		Gal[p].DustColdGasClouds_elements.Si += Dust_Ferrosilite * 0.212884 * mu_gas;
		Gal[p].DustColdGasClouds_elements.O  += Dust_Ferrosilite * 0.363819 * mu_gas;
		
		//Quartz SiO4 -------------------------------------------------
		Gal[p].DustColdGasDiff_elements.Si += Dust_Quartz * 0.305002 * (1.0 - mu_gas);
		Gal[p].DustColdGasDiff_elements.O  += Dust_Quartz * 0.694998 * (1.0 - mu_gas);
		
		Gal[p].DustColdGasClouds_elements.Si += Dust_Quartz * 0.305002 * mu_gas;
		Gal[p].DustColdGasClouds_elements.O  += Dust_Quartz * 0.694998 * mu_gas;
		
		//SiC SiC -----------------------------------------------------
		Gal[p].DustColdGasDiff_elements.Si += Dust_SiC * 0.700453 * (1.0 - mu_gas);
	    Gal[p].DustColdGasDiff_elements.Cb  += Dust_SiC * 0.299547 * (1.0 - mu_gas);
		
		Gal[p].DustColdGasClouds_elements.Si += Dust_SiC * 0.700453 * mu_gas;
	    Gal[p].DustColdGasClouds_elements.Cb  += Dust_SiC * 0.299547 * mu_gas;

		//Iron Fe -----------------------------------------------------
		Gal[p].DustColdGasDiff_elements.Fe += min(Dust_Iron * 1.0 * (1.0 - mu_gas), Gal[p].ColdGasDiff_elements.Fe - Gal[p].DustColdGasDiff_elements.Fe);
		
		Gal[p].DustColdGasClouds_elements.Fe += min(Dust_Iron * 1.0 * mu_gas, Gal[p].ColdGasClouds_elements.Fe - Gal[p].DustColdGasClouds_elements.Fe);
		
		//Carbon C ----------------------------------------------------
		Gal[p].DustColdGasDiff_elements.Cb += Dust_Carbon * 1.0 * (1.0 - mu_gas);
		
		Gal[p].DustColdGasClouds_elements.Cb += Dust_Carbon * 1.0 * mu_gas;
		
	} //if sfh_DM >0
    
    mass_checks("Dust from AGB",p);
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
        //Mg = 24.305
        //O = 15.999
	
	    #ifdef FULL_DUST_RATES
	        //This is estimating the dust in various compounds (say silicates) using the amount of the particular element (say silicon) produced in a process
			Gal[p].DustColdGasRates.SNII += (SNII_prevstep_Cold_Si[i] * eta_SNII_Sil * A_Sil_dust/A_Si)/(deltaT * UnitTime_in_years);
			Gal[p].DustColdGasRates.SNII += (SNII_prevstep_Cold_Fe[i] * eta_SNII_Fe  * A_Fe_dust/A_Fe )/(deltaT * UnitTime_in_years);
			Gal[p].DustColdGasRates.SNII += (SNII_prevstep_Cold_Si[i] * eta_SNII_SiC * A_SiC_dust/A_Si)/(deltaT * UnitTime_in_years);
			Gal[p].DustColdGasRates.SNII += (SNII_prevstep_Cold_Cb[i] * eta_SNII_Cb  * A_Cb_dust/A_Cb) /(deltaT * UnitTime_in_years);	
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
		
			Gal[p].DustColdGasDiff_elements.Si += Dust_Silicates * 0.210432 * (1.0 - mu_gas);
			Gal[p].DustColdGasDiff_elements.Mg += min(Dust_Silicates * 0.091053 * (1.0 - mu_gas), Gal[p].ColdGasDiff_elements.Mg - Gal[p].DustColdGasDiff_elements.Mg);
			Gal[p].DustColdGasDiff_elements.Fe += min(Dust_Silicates * 0.278948 * (1.0 - mu_gas), Gal[p].ColdGasDiff_elements.Fe - Gal[p].DustColdGasDiff_elements.Fe);
			Gal[p].DustColdGasDiff_elements.O  += Dust_Silicates * 0.419567 * (1.0 - mu_gas);
			
			Gal[p].DustColdGasClouds_elements.Si += Dust_Silicates * 0.210432 * mu_gas;
			Gal[p].DustColdGasClouds_elements.Mg += min(Dust_Silicates * 0.091053 * mu_gas, Gal[p].ColdGasClouds_elements.Mg - Gal[p].DustColdGasClouds_elements.Mg);
			Gal[p].DustColdGasClouds_elements.Fe += min(Dust_Silicates * 0.278948 * mu_gas, Gal[p].ColdGasClouds_elements.Fe - Gal[p].DustColdGasClouds_elements.Fe);
			Gal[p].DustColdGasClouds_elements.O  += Dust_Silicates * 0.419567 * mu_gas;

			//SNII SiC --------------------------

			Gal[p].DustColdGasDiff_elements.Si += Dust_SiC * 0.700453 * (1.0 - mu_gas);
			Gal[p].DustColdGasDiff_elements.Cb  += Dust_SiC * 0.299547 * (1.0 - mu_gas);
			
			Gal[p].DustColdGasClouds_elements.Si += Dust_SiC * 0.700453 * mu_gas;
			Gal[p].DustColdGasClouds_elements.Cb  += Dust_SiC * 0.299547 * mu_gas;

			//SNII Fe -------------------------

			Gal[p].DustColdGasDiff_elements.Fe += min(Dust_Iron * 1.0 * (1.0 - mu_gas), Gal[p].ColdGasDiff_elements.Fe - Gal[p].DustColdGasDiff_elements.Fe);
			
			Gal[p].DustColdGasClouds_elements.Fe += Dust_Iron * 1.0  * mu_gas;

			//SNII Cb ------------------------

			Gal[p].DustColdGasDiff_elements.Cb += Dust_Carbon * 1.0 * (1.0 - mu_gas);
			
			Gal[p].DustColdGasClouds_elements.Cb += Dust_Carbon * 1.0 * mu_gas;
            
            mass_checks("Dust from SNII",p);
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
		    Gal[p].DustColdGasRates.SNIA  += (SNIa_prevstep_Cold_Fe[i] * eta_SNIa_Fe  * A_Fe_dust/A_Fe)/(deltaT * UnitTime_in_years);
        #endif	
	
		Gal[p].DustColdGasDiff_elements.Fe += min(Dust_Iron * 1.0 * (1.0 - mu_gas), Gal[p].ColdGasDiff_elements.Fe - Gal[p].DustColdGasDiff_elements.Fe);
		
		Gal[p].DustColdGasClouds_elements.Fe += min(Dust_Iron * 1.0 * mu_gas, Gal[p].ColdGasClouds_elements.Fe - Gal[p].DustColdGasClouds_elements.Fe);
		
		mass_checks("Dust from SNIA",p);
		
	}//if sfh_DM >0
#endif //DUST_SNIA

} //loop over SFH bins


//*****************************************
//Dust grain growth inside molecular clouds 
//*****************************************

#ifdef DUST_GROWTH
    if ((metals_total(Gal[p].MetalsColdGas)>0.0) && (Gal[p].ColdGas > 0.00)) {
			
        //Growth implementation requires the fraction of elements in both the media
        //to be distributed according to the method of molecular hydrogen computation 
        //method used, since mu_gas is passed into the dust growth equations.
        shuffle_ISM(p);
    
        
        float A_H = 1.008;
        float A_He = 4.003;
        float A_Cb = 12.01; //0
        float A_N = 14.007; //1
        float A_O = 15.999; //2
        float A_Ne = 20.1797; //3
        float A_Mg = 24.305; //4
        float A_Si = 28.085; //5
        float A_S = 32.06; //6
        float A_Ca = 40.078; //7
        float A_Fe = 55.845; //8
        
        
        //The number of CO molecules that can be produced from available Carbon and Oxygen
        //in the clouds. xxAssuming only 30% C is locked up as CO in cloudsxx
        float num_CO = min((Gal[p].ColdGasClouds_elements.Cb-Gal[p].DustColdGasClouds_elements.Cb)/A_Cb, 
                     min(Gal[p].ColdGasClouds_elements.Cb/A_Cb*0.3,
	                 (Gal[p].ColdGasClouds_elements.O-Gal[p].DustColdGasClouds_elements.O)/A_O));
        
        //Carbon and Oxygen available in the clouds for grain growth
        double Cb_clouds = Gal[p].ColdGasClouds_elements.Cb - num_CO*A_Cb;
        double O_clouds = Gal[p].ColdGasClouds_elements.O - num_CO*A_O;
        
        /*
        We assume that oxygen is volatile such that it can exist in the dust phase only in the 
        form of compounds. We adopt from Zhukovska et al.(2008) that the dominant silicate species
        are olivine and pyroxene in the ratio 32:68. 
        Olivine: [Mg_x Fe_{1-x}]_2 Si O_4
        Pyroxene: Mg_x Fe_{1-x} Si O_3 
        x is assumed to be 0.8 in their paper, this value doesn't seem to affect the dust masses 
        that much. It is because this doesn't significantly affect the amount of oxygen left in 
        the cold gas.
        */
        float f_ol = 0.32;
        float x = 0.8;
        //float M_olivine = 2.*x*A_Mg + 2.*(1-x)*A_Fe + A_Si + 4.*A_O;
        //float M_pyroxene = x*A_Mg + (1.-x)*A_Fe + A_Si + 3.*A_O;
        float N_O = f_ol*4. + (1.-f_ol)*3.;  
        float N_Mg = f_ol*2.*x + (1.-f_ol)*x;
        float N_Si = f_ol*1. + (1.-f_ol)*1.;
        float N_Fe = f_ol*2.*(1-x) + (1.-f_ol)*(1.-x);
        
        float num_Silicates = min((O_clouds-Gal[p].DustColdGasClouds_elements.O)/(N_O*A_O), min((Gal[p].ColdGasClouds_elements.Mg-Gal[p].DustColdGasClouds_elements.Mg)/(N_Mg*A_Mg), min((Gal[p].ColdGasClouds_elements.Si-Gal[p].DustColdGasClouds_elements.Si)/(N_Si*A_Si), (Gal[p].ColdGasClouds_elements.Fe-Gal[p].DustColdGasClouds_elements.Fe)/(N_Fe*A_Fe))));
        
        //Major iron oxide species: Fe3O4 and Fe2O3. Oxygen combining with iron.  
        // 5 Fe for 7 oxygen. N_Oxide = 1.4
        float N_Oxide = 1.4;
        float num_iron_oxide = min((O_clouds-num_Silicates*N_O*A_O-Gal[p].DustColdGasClouds_elements.O)/(N_Oxide*A_O), (Gal[p].ColdGasClouds_elements.Fe-num_Silicates*N_Fe*A_Fe-Gal[p].DustColdGasClouds_elements.Fe)/A_Fe);
        
        
        float f_Cb_max =  1.0 - num_CO*A_Cb/Gal[p].ColdGasClouds_elements.Cb;
        float f_O_max =  (num_Silicates*N_O*A_O + num_iron_oxide*N_Oxide*A_O + Gal[p].DustColdGasClouds_elements.O)/Gal[p].ColdGasClouds_elements.O;
        /*
        float f_Mg_max = (num_Silicates*N_Mg*A_Mg + Gal[p].DustColdGasClouds_elements.Mg)/Gal[p].ColdGasClouds_elements.Mg;
        float f_Fe_max = (num_Silicates*N_Fe*A_Fe + Gal[p].DustColdGasClouds_elements.Fe)/Gal[p].ColdGasClouds_elements.Fe;
        */
        
        float Clouds_tot = elements_total(Gal[p].ColdGasClouds_elements);
        
        float DustDiff_init = elements_total(Gal[p].DustColdGasDiff_elements); 
        float DustClouds_init = elements_total(Gal[p].DustColdGasClouds_elements);    
        
        
        //--------------------------------------------------------
        //These are some numbers that should actually go into the calculation of
        //t_acc_0 of a particular galaxy, but is not implemented at the moment
        //Ideally the equation should be t_acc = t_acc_0*(Mass in clouds/Mass of dust in clouds)
        //                                       *(1/(mean molecular weight * number density))
        float met_mean = (Gal[p].ColdGasClouds_elements.Cb/A_Cb + 
                          Gal[p].ColdGasClouds_elements.N/A_N   + 
                          Gal[p].ColdGasClouds_elements.O/A_O   + 
                          Gal[p].ColdGasClouds_elements.Ne/A_Ne + 
                          Gal[p].ColdGasClouds_elements.Mg/A_Mg + 
                          Gal[p].ColdGasClouds_elements.Si/A_Si + 
                          Gal[p].ColdGasClouds_elements.S/A_S   + 
                          Gal[p].ColdGasClouds_elements.Ca/A_Ca + 
                          Gal[p].ColdGasClouds_elements.Fe/A_Fe)/(Clouds_tot-Gal[p].ColdGasClouds_elements.H - Gal[p].ColdGasClouds_elements.He);
        float oneovermeanmu = Gal[p].ColdGasClouds_elements.H/Clouds_tot +
                              Gal[p].ColdGasClouds_elements.He/(4*Clouds_tot) +
                              met_mean*((Clouds_tot-Gal[p].ColdGasClouds_elements.H - Gal[p].ColdGasClouds_elements.He)/Clouds_tot) ;
        float meanmu = 1./oneovermeanmu;
        
        //40.76 is the conversion factor for solar mass/pc^3 to g/cm^3 multiplied by avagadro's number
        //to get the number of atoms to get an approximate number density.
        //double naprox = 40.76*3.0*Clouds_tot/(4.0*M_PI*pow((((Gal[p].GasDiskRadius*1E6)/Hubble_h)/3.0),3));  
		//------------------------------------------------------
		
		
		//calculate exchange timescales and molecular gas fractions
		float rho, rho_crit, t_exchange, tdyn;
		
		if(Gal[p].Type == 0)
        {
            tdyn = (((Gal[p].GasDiskRadius*1E6)/Hubble_h)/3.0) / Gal[p].Vmax;
        }
        else
        {
            tdyn = (((Gal[p].GasDiskRadius*1E6)/Hubble_h)/3.0) / Gal[p].InfallVmax;
        }
		
		tdyn *= 3.086e13;  //converting from units of pc s/km to s
		float Grav = 6.674e-11; //in SI units
		rho =  1./(Grav*tdyn*tdyn);  //in units of kg/m^3
		double num_den = (6.022e26)*(1e-6)*rho/(meanmu);  //in units of cm^-3
		
        float Z_coldgas, Z_fraction;
        float Z_sun = 0.0134;
        if (Gal[p].ColdGas > 0.00) 
        {
		    Z_coldgas = metals_total(Gal[p].MetalsColdGas)/Gal[p].ColdGas;
		    Z_fraction = Z_coldgas/Z_sun;
	    }
	    else 
	    {
		    Z_fraction = 0.0;
	    }
	    
	    float Sigma_gas = (1.0e10 * (Gal[p].ColdGas/Hubble_h))/(M_PI*((((Gal[p].GasDiskRadius*1.0e6)/Hubble_h)/3.0)*(((Gal[p].GasDiskRadius*1.0e6)/Hubble_h)/3.0)));
        
	    float c_f;
	    if (Z_fraction < 1)
	    {
	        c_f = pow(Z_fraction, -0.7);
	    }
	    else
	    {
	        c_f = 1.0;
        }
	    
	    float Sigma_comp = c_f*Sigma_gas;
        
        		
		//rho = (1.0e10 * (Gal[p].ColdGas/Hubble_h))/((4./3.)*M_PI*(((Gal[p].GasDiskRadius*1.0e6)/Hubble_h)/3.0)*(((Gal[p].GasDiskRadius*1.0e6)/Hubble_h)/3.0)*(((Gal[p].GasDiskRadius*1.0e6)/Hubble_h)/3.0)); // in Msol/pc^3
		//rho *= 6.768e-20; //in kg/m^3
		//rho_crit = 9.9e-27; // in kg/m^3
		
		double t_acc_0 = 15e3; //At the moment value chosen for no reason other than the 
		                      //fact that it seems to give good results
		t_exchange = 2e7 ; //Exchange timescale, approximate timescale of cloud lifetime
		                  //from literature
        
        
        if (DustClouds_init == 0) {
            
            Gal[p].t_acc = 1e15;
        
        }
        else {
            
            Gal[p].t_acc = t_acc_0*(Clouds_tot/DustClouds_init);
        
        }
        
               
        //Here we use an approximation for calculating the fraction of dust produced in
        //clouds and the diffused ISM which is valid for a constant accretion timescale
        //across the timestep. Since that is not the case as it depends on the amount of 
        //dust in the clouds, we divide the calculation in constant bins of the time step
        //in between snaps. The dust mass in clouds is calculated again to update t_acc.
        
        //Note: Dust molecules injected into the ISM has elements bonded with itself and 
        //only certain other elements. 
        // * Cb: Cb, Si
        // * O: Mg, Si, Fe
        // * Mg: Si, O
        // * Si: with all
        // * Fe: Si, Fe 
        // * Maybe bring this into play like implemented by Zhukovska
        
        Gal[p].f_cmax[0] = f_Cb_max;
        if (isnan(Gal[p].f_cmax[0])) {Gal[p].f_cmax[0] = 1.;}
        Gal[p].f_cmax[2] = f_O_max;
        if (isnan(Gal[p].f_cmax[2])) {Gal[p].f_cmax[2] = 1.;}
        
        /*
        Gal[p].f_cmax[4] = f_Mg_max;
        if (isnan(Gal[p].f_cmax[4])) {Gal[p].f_cmax[4] = 1.;}
        Gal[p].f_cmax[8] = f_Fe_max;
        if (isnan(Gal[p].f_cmax[8])) {Gal[p].f_cmax[4] = 1.;}
        */
        
        float t_dt = (dt*UnitTime_in_years);
        int l;
        int N = 3;
        float DustClouds_tot;
        float dt_bins = t_dt/N;
        /*if (t_acc*5 < t_dt)
        {
            N = t_dt/t_acc+1;
            dt_bins = t_dt;
        }*/ 
        for(l=0; l<N; l++) {
            
            Gal[p].f_i[0] = Gal[p].DustColdGasDiff_elements.Cb/Gal[p].ColdGasDiff_elements.Cb;    
            Gal[p].f_c[0] = Gal[p].DustColdGasClouds_elements.Cb/Gal[p].ColdGasClouds_elements.Cb;
            drop_fnan(0, p);
            update_fractions(dt_bins, Gal[p].t_acc, t_exchange, mu_gas, p, 0);
            drop_fnan(0, p);
            
            Gal[p].f_i[2] = Gal[p].DustColdGasDiff_elements.O/Gal[p].ColdGasDiff_elements.O;
            Gal[p].f_c[2] = Gal[p].DustColdGasClouds_elements.O/Gal[p].ColdGasClouds_elements.O;
            drop_fnan(2, p);
            update_fractions(dt_bins, Gal[p].t_acc, t_exchange, mu_gas, p, 2);
            drop_fnan(2, p);
            
            Gal[p].f_i[4] = Gal[p].DustColdGasDiff_elements.Mg/Gal[p].ColdGasDiff_elements.Mg;
            Gal[p].f_c[4] = Gal[p].DustColdGasClouds_elements.Mg/Gal[p].ColdGasClouds_elements.Mg;
            drop_fnan(4, p);
            update_fractions(dt_bins, Gal[p].t_acc, t_exchange, mu_gas, p, 4);
            drop_fnan(4, p);
            
            Gal[p].f_i[5] = Gal[p].DustColdGasDiff_elements.Si/Gal[p].ColdGasDiff_elements.Si;
            Gal[p].f_c[5] = Gal[p].DustColdGasClouds_elements.Si/Gal[p].ColdGasClouds_elements.Si;
            drop_fnan(5, p);
            update_fractions(dt_bins, Gal[p].t_acc, t_exchange, mu_gas, p, 5);
            drop_fnan(5, p);
            
            Gal[p].f_i[6] = Gal[p].DustColdGasDiff_elements.S/Gal[p].ColdGasDiff_elements.S;
            Gal[p].f_c[6] = Gal[p].DustColdGasClouds_elements.S/Gal[p].ColdGasClouds_elements.S;
            drop_fnan(6, p);
            update_fractions(dt_bins, Gal[p].t_acc, t_exchange, mu_gas, p, 6);
            drop_fnan(6, p);
            
            Gal[p].f_i[7] = Gal[p].DustColdGasDiff_elements.Ca/Gal[p].ColdGasDiff_elements.Ca;
            Gal[p].f_c[7] = Gal[p].DustColdGasClouds_elements.Ca/Gal[p].ColdGasClouds_elements.Ca;
            drop_fnan(7, p);
            update_fractions(dt_bins, Gal[p].t_acc, t_exchange, mu_gas, p, 7);
            drop_fnan(7, p);
            
            Gal[p].f_i[8] = Gal[p].DustColdGasDiff_elements.Fe/Gal[p].ColdGasDiff_elements.Fe;
            Gal[p].f_c[8] = Gal[p].DustColdGasClouds_elements.Fe/Gal[p].ColdGasClouds_elements.Fe;
            drop_fnan(8, p);
            update_fractions(dt_bins, Gal[p].t_acc, t_exchange, mu_gas, p, 8);
            drop_fnan(8, p);
            
            //drop_dustnan(p);
            //Updating the amount of dust after this mini-step.
            update_dust(p); 
            
            num_CO = min((Gal[p].ColdGasClouds_elements.Cb-Gal[p].DustColdGasClouds_elements.Cb)/A_Cb, 
                     min(Gal[p].ColdGasClouds_elements.Cb/A_Cb*0.3,
	                 (Gal[p].ColdGasClouds_elements.O-Gal[p].DustColdGasClouds_elements.O)/A_O));
            
            O_clouds = Gal[p].ColdGasClouds_elements.O - num_CO*A_O;
            
            f_Cb_max =  1.0 - num_CO*A_Cb/Gal[p].ColdGasClouds_elements.Cb;
            
            num_Silicates = min((O_clouds-Gal[p].DustColdGasClouds_elements.O)/(N_O*A_O), min((Gal[p].ColdGasClouds_elements.Mg-Gal[p].DustColdGasClouds_elements.Mg)/(N_Mg*A_Mg), min((Gal[p].ColdGasClouds_elements.Si-Gal[p].DustColdGasClouds_elements.Si)/(N_Si*A_Si), (Gal[p].ColdGasClouds_elements.Fe-Gal[p].DustColdGasClouds_elements.Fe)/(N_Fe*A_Fe))));
            
            num_iron_oxide = min((O_clouds-num_Silicates*N_O*A_O-Gal[p].DustColdGasClouds_elements.O)/(N_Oxide*A_O), (Gal[p].ColdGasClouds_elements.Fe-num_Silicates*N_Fe*A_Fe-Gal[p].DustColdGasClouds_elements.Fe)/A_Fe);
            
            //f_O_max =  (num_Silicates*N_O*A_O + Gal[p].DustColdGasClouds_elements.O)/Gal[p].ColdGasClouds_elements.O;
            f_O_max =  (num_Silicates*N_O*A_O + num_iron_oxide*N_Oxide*A_O + Gal[p].DustColdGasClouds_elements.O)/Gal[p].ColdGasClouds_elements.O;
            //f_Mg_max = (num_Silicates*N_Mg*A_Mg + Gal[p].DustColdGasClouds_elements.Mg)/Gal[p].ColdGasClouds_elements.Mg;
            //f_Fe_max = (num_Silicates*N_Fe*A_Fe + Gal[p].DustColdGasClouds_elements.Fe)/Gal[p].ColdGasClouds_elements.Fe;
            
            Gal[p].f_cmax[0] = f_Cb_max;
            if (isnan(Gal[p].f_cmax[0])) {Gal[p].f_cmax[0] = 1.;}
            Gal[p].f_cmax[2] = f_O_max;
            if (isnan(Gal[p].f_cmax[2])) {Gal[p].f_cmax[2] = 1.;}
            
            /*
            Gal[p].f_cmax[4] = f_Mg_max;
            if (isnan(Gal[p].f_cmax[4])) {Gal[p].f_cmax[4] = 0.;}
            Gal[p].f_cmax[8] = f_Fe_max;
            if (isnan(Gal[p].f_cmax[8])) {Gal[p].f_cmax[4] = 0.;}
            */
                        
            DustClouds_tot = elements_total(Gal[p].DustColdGasClouds_elements);  
            
            Gal[p].t_acc = t_acc_0*(Clouds_tot/DustClouds_tot);
        
        }

        //Gal[p].t_acc = surf_rho;        
        		
		Gal[p].DustColdGasDiff_elements.Cb = Gal[p].f_i[0] * Gal[p].ColdGasDiff_elements.Cb ;
		//Gal[p].DustColdGasDiff_elements.N  = Gal[p].f_i[1] * Gal[p].ColdGasDiff_elements.N ;
		Gal[p].DustColdGasDiff_elements.O  = Gal[p].f_i[2] * Gal[p].ColdGasDiff_elements.O ;
		//Gal[p].DustColdGasDiff_elements.Ne = Gal[p].f_i[3] * Gal[p].ColdGasDiff_elements.Ne;
		Gal[p].DustColdGasDiff_elements.Mg = Gal[p].f_i[4] * Gal[p].ColdGasDiff_elements.Mg ;
		Gal[p].DustColdGasDiff_elements.Si = Gal[p].f_i[5] * Gal[p].ColdGasDiff_elements.Si;
		Gal[p].DustColdGasDiff_elements.S  = Gal[p].f_i[6] * Gal[p].ColdGasDiff_elements.S ;
		Gal[p].DustColdGasDiff_elements.Ca = Gal[p].f_i[7] * Gal[p].ColdGasDiff_elements.Ca ;
		Gal[p].DustColdGasDiff_elements.Fe = Gal[p].f_i[8] * Gal[p].ColdGasDiff_elements.Fe ;
		
		
		Gal[p].DustColdGasClouds_elements.Cb = Gal[p].f_c[0] * Gal[p].ColdGasClouds_elements.Cb;
		//Gal[p].DustColdGasClouds_elements.N  = Gal[p].f_c[1] * Gal[p].ColdGasClouds_elements.N;
		Gal[p].DustColdGasClouds_elements.O  = Gal[p].f_c[2] * Gal[p].ColdGasClouds_elements.O;
		//Gal[p].DustColdGasClouds_elements.Ne = Gal[p].f_c[3] * Gal[p].ColdGasClouds_elements.Ne;
		Gal[p].DustColdGasClouds_elements.Mg = Gal[p].f_c[4] * Gal[p].ColdGasClouds_elements.Mg;
		Gal[p].DustColdGasClouds_elements.Si = Gal[p].f_c[5] * Gal[p].ColdGasClouds_elements.Si;
		Gal[p].DustColdGasClouds_elements.S  = Gal[p].f_c[6] * Gal[p].ColdGasClouds_elements.S;
		Gal[p].DustColdGasClouds_elements.Ca = Gal[p].f_c[7] * Gal[p].ColdGasClouds_elements.Ca;
		Gal[p].DustColdGasClouds_elements.Fe = Gal[p].f_c[8] * Gal[p].ColdGasClouds_elements.Fe;
		
		float Dust_DiffGrowth = elements_total(Gal[p].DustColdGasDiff_elements) - DustDiff_init;
		float Dust_CloudsGrowth = elements_total(Gal[p].DustColdGasClouds_elements) - DustClouds_init;
		
        #ifdef FULL_DUST_RATES
		Gal[p].DustColdGasRates.GROW += (Dust_DiffGrowth + Dust_CloudsGrowth)/(deltaT * UnitTime_in_years);
        #endif
        
        mass_checks("Dust from growth",p);
        
}
#endif //DUST_GROWTH
    
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
