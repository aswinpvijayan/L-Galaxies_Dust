/*  Copyright (C) <2016>  <L-Galaxies>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/> */

/*
 * recipe_yields.c
 *
 *  Created on: 18.11.2011
 *      Author: robyates
 *
 *
 * NOTE: 


 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"


double mu_Obr(int p) {
    
    float K = 4.949e-5;   // 4.949E-5(units pc^4) / (M_solar ^2)    Wolfram alpha gave me a value of 4.949E-5 
    //float factor = (((Gal[p].GasDiskRadius*1.0e6)/Hubble_h)/3.0)*(((Gal[p].GasDiskRadius*1.0e6)/Hubble_h)/3.0)*(((Gal[p].GasDiskRadius*1.0e6)/Hubble_h)/3.0)*(((Gal[p].GasDiskRadius*1.0e6)/Hubble_h)/3.0);
	//float rmol=pow((K*(1.0/factor)*(1.0e10 * (Gal[p].ColdGas/Hubble_h))*((1.0e10*Gal[p].ColdGas/Hubble_h)+0.4*(1E10*(Gal[p].DiskMass/Hubble_h)))),0.8);  //Obreshkow et al. (2009) uses the scale radius for his calculation that is why there is a division by 3 for the GasDisk radius since an exponenetial profile is assumed
	float rmol=pow((K*pow((((Gal[p].GasDiskRadius*1E6)/Hubble_h)/3.0),-4)*(1E10 * (Gal[p].ColdGas/Hubble_h))*((1E10*Gal[p].ColdGas/Hubble_h)+0.4*(1E10*(Gal[p].DiskMass/Hubble_h)))),0.8);  //Obreshkow et al. (2009) uses the scale radius for his calculation that is why there is a division by 3 for the GasDisk radius since an exponenetial profile is assumed
  	float rmolgal=pow((3.44*pow(rmol,-0.506)+4.82*pow(rmol,-1.054)),-1);
    float mu = rmolgal/(1+rmolgal); //Molecular gas fraction
    
    return mu;
}


double mu_GK11(int p) {

    float pi = 3.14159265358979323846;
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

    
	float SurfHIH2 = (Gal[p].ColdGas_elements.H)/(pi*(((Gal[p].GasDiskRadius*1.0e6)/Hubble_h)/3.0)*(((Gal[p].GasDiskRadius*1.0e6)/Hubble_h)/3.0));	
	float Umw = (Gal[p].Sfr * UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS)/1.0; //From Gnedin and Kratsow (2011) UV photon proportional to the SFR/SFRmw where SFRmw = 1 Msol/yr
	float Dmw = Z_fraction;
	
	float Dstar = (1.5e-3)*log(1+pow(3*Umw, 1.7));	
	float alpha = 5*(Umw/2)/(1+Umw*Umw/4.0);
	float s = 0.04/(Dstar + Dmw);
	float g = (1+alpha*s+s*s)/(1+s);
	float Lambda = log(1 + g*pow(Dmw, 3/7)*(pow(Umw/15, 4/7)));
	float SurfDbar = 20 * (pow(Lambda, 4/7)*Dmw)*pow(1+Umw*pow(Dmw, 2), -0.5);
	
	float f_h2 = 1./((1 + (SurfDbar/SurfHIH2))*(1 + (SurfDbar/SurfHIH2)));
	
    return f_h2;

}


double mu_Krumholz(int p) {

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
	
	if (Z_fraction < 0.01) {Z_fraction = 0.01;}
	
	float chi = 3.1*(1.0 + 3.1*pow(Z_fraction, 0.365))/4.1;
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
	
	float tau_c = 0.066*Sigma_comp*Z_fraction;
	float s = log(1.0 + 0.6*chi + 0.01*chi*chi)/(0.6*tau_c);
	float f_h2;
	if (s < 2.)
	{
	    f_h2 = 1. - (0.75)*(s)/(1. + 0.25*s);
    }
	else 
	{
	    f_h2 = 0.;
    }
	
	return f_h2; 

}


void update_yields_and_return_mass(int p, int centralgal, double dt, int nstep)
{
	int Zi, igal;
	double timestep_width; //Width of current timestep in CODE UNITS
	int TimeBin; //Bin in Yield arrays corresponding to current timestep
	double Zi_disp, NormSNIIMassEjecRate_actual, NormSNIaMassEjecRate_actual, NormAGBMassEjecRate_actual, NormSNIIMetalEjecRate_actual, NormSNIaMetalEjecRate_actual, NormAGBMetalEjecRate_actual;
	
	struct elements coldgas_init = elements_init();
    struct elements coldgas_add = elements_init();	
	
	#ifdef Obreshkow
	    Gal[p].mu_gas = mu_Obr(p);
    #endif
	
	#ifdef GK11
	    Gal[p].mu_gas = mu_GK11(p);
	#endif

	#ifdef Krumholz
        Gal[p].mu_gas = mu_Krumholz(p);
    #endif
	
	if (Gal[p].mu_gas > 1.)
	{
	    printf("\n*** mu_gas greater than 1 ***\n");
	    terminate("");
	}
	else if (Gal[p].mu_gas < 0.)
	{
	    printf("\n*** mu_gas less than 0 ***\n");
	    terminate("");
	}
	
#ifdef INDIVIDUAL_ELEMENTS
	double NormSNIIYieldRate_actual[NUM_ELEMENTS], NormSNIaYieldRate_actual[NUM_ELEMENTS], NormAGBYieldRate_actual[NUM_ELEMENTS];
#endif
	double timet, sfh_time;
	double fwind, SNIIEjectaToHot; //Required for metal-rich wind implementation
	double DiskSFR, step_width_times_DiskSFR, DiskSFR_physical_units, step_width_times_DiskSFR_physical_units, inverse_DiskMass_physical_units;
	double BulgeSFR, step_width_times_BulgeSFR, BulgeSFR_physical_units, step_width_times_BulgeSFR_physical_units, inverse_BulgeMass_physical_units;
	double ICMSFR, step_width_times_ICMSFR, ICMSFR_physical_units, step_width_times_ICMSFR_physical_units, inverse_ICM_physical_units;
	double Disk_total_metallicity, Bulge_total_metallicity, ICM_total_metallicity;
	double NormMassEjecRateSumAllTypes;

	int n; //Iterator used for loop over NOUT when updating MassWeightedAge
	double AgeCorrectionDisk[NOUT];
	double AgeCorrectionBulge[NOUT];

	TotalMassReturnedToColdDiskGas=0.0;
	TotalMassReturnedToHotGas=0.0;
    
    SN2_rate = 0.;
    SN1_rate = 0.;
    
	for(n=0;n<NOUT;n++)
	{
		AgeCorrectionDisk[n] = 0.0;
		AgeCorrectionBulge[n] = 0.0;
	}

	timestep_width = dt; //Width of current timestep in CODE UNITS (units cancel out when dividing by SFH bin width, sfh_dt) (12-04-12)
	TimeBin = (STEPS*Gal[p].SnapNum)+nstep;//Bin in Yield tables corresponding to current timestep
	timet = NumToTime(Gal[p].SnapNum) - (nstep + 0.5) * dt; //Time from middle of the current timestep to z=0 (used here for MassWeightAge corrections)
	//NB: NumToTime(Gal[p].SnapNum) is the time to z=0 from start of current snapshot
	//    nstep is the number of the current timestep (0-19)
	//    dt is the width of one timestep within current snapshot
	
	
#ifdef METALRICHWIND
	ColdGasSurfaceDensity = max(0.0, (Gal[p].ColdGas*(1.0e10/Hubble_h))/(4.0*3.14159265*Gal[p].GasDiskRadius*Gal[p].GasDiskRadius/Hubble_h));
	fwind = min(1.0, max(0.0, 1.0/(ColdGasSurfaceDensity/5.0e12))); //Fraction of SN-II ejecta put directly into HotGas
	if (Gal[p].ColdGas != (float)Gal[p].ColdGas) {fwind = 1.0;}
#endif
#ifndef METALRICHWIND
	fwind = 0.0; //For all stellar ejecta (from disk) to ColdGas
#endif

    //for stars dying that enrich the Hot gas directly
    if(Gal[p].Type==2)
      igal=Gal[p].CentralGal;
    else
      igal=p;

#ifdef DETAILED_DUST
    coldgas_init = elements_add(elements_init(), Gal[p].ColdGas_elements, 1.);
#endif  

    int i;
    for (i=0;i<=Gal[p].sfh_ibin;i++) //LOOP OVER SFH BINS
    {
    	sfh_time=Gal[p].sfh_t[i]+(0.5*Gal[p].sfh_dt[i]);

#ifdef DETAILED_DUST
// These variables store the amount of metals produced by SNe in a timestep for use in
// model_dustyields.c
		SNII_prevstep_Cold_Si[i] = 0.0;
		SNII_prevstep_Cold_Fe[i] = 0.0;
		SNII_prevstep_Cold_Cb[i] = 0.0;
		SNIa_prevstep_Cold_Fe[i] = 0.0;
#endif //DETAILED_DUST


    //*****************************************
    //ENRICHMENT FROM DISK STARS INTO COLD GAS:
    //*****************************************
    if (Gal[p].DiskMass > 0.0) //Only calculate enrichment from the disc if there is a disc at the current timestep.
    if (Gal[p].sfh_DiskMass[i] > 0.0)
    {
     	//pre-calculations to speed up the code
    	DiskSFR = Gal[p].sfh_DiskMass[i]/Gal[p].sfh_dt[i];
    	step_width_times_DiskSFR = timestep_width * DiskSFR;
    	DiskSFR_physical_units = DiskSFR * (1.0e10/Hubble_h); //Note: This is NOT in physical units (i.e. NOT in Msun/yr, but in Msun/[code_time_units]). But this is ok, as code-time-units cancel out when multiplying by timestep_width to get 'step_width_times_DiskSFR_physical_units' on the line below ('DiskSFR_physical_units' is never used itself).
    	step_width_times_DiskSFR_physical_units = timestep_width * DiskSFR_physical_units;
    	inverse_DiskMass_physical_units=Hubble_h/(Gal[p].sfh_DiskMass[i]*1.0e10);
    	Disk_total_metallicity=metals_total(Gal[p].sfh_MetalsDiskMass[i])/Gal[p].sfh_DiskMass[i];


    	Zi = find_initial_metallicity(p, i, 1, 1);
    	//Interpolate the disk metallicity on the lifetimeMetallicities tables:
    	Zi_disp = (Disk_total_metallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.
#ifdef DETAILED_DUST
		Zi_saved[i] = Zi;
		Zi_disp_saved[i] = Zi_disp;		
#endif    	
    	NormSNIIMassEjecRate_actual = NormSNIIMassEjecRate[TimeBin][i][Zi] + ((NormSNIIMassEjecRate[TimeBin][i][Zi+1] - NormSNIIMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIaMassEjecRate_actual = NormSNIaMassEjecRate[TimeBin][i][Zi] + ((NormSNIaMassEjecRate[TimeBin][i][Zi+1] - NormSNIaMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormAGBMassEjecRate_actual = NormAGBMassEjecRate[TimeBin][i][Zi] + ((NormAGBMassEjecRate[TimeBin][i][Zi+1] - NormAGBMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIIMetalEjecRate_actual = NormSNIIMetalEjecRate[TimeBin][i][Zi] + ((NormSNIIMetalEjecRate[TimeBin][i][Zi+1] - NormSNIIMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIaMetalEjecRate_actual = NormSNIaMetalEjecRate[TimeBin][i][Zi] + ((NormSNIaMetalEjecRate[TimeBin][i][Zi+1] - NormSNIaMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormAGBMetalEjecRate_actual = NormAGBMetalEjecRate[TimeBin][i][Zi] + ((NormAGBMetalEjecRate[TimeBin][i][Zi+1] - NormAGBMetalEjecRate[TimeBin][i][Zi])*Zi_disp);

        //SNe rate calculations for getting dust destruction
    	SN2_rate += max(0., DiskSFR_physical_units*(SNIIRate[TimeBin][i][Zi] + (SNIIRate[TimeBin][i][Zi+1] - SNIIRate[TimeBin][i][Zi])*Zi_disp))/(UnitTime_in_years);
        SN1_rate += max(0., DiskSFR_physical_units*(SNIaRate[TimeBin][i][Zi] + (SNIaRate[TimeBin][i][Zi+1] - SNIaRate[TimeBin][i][Zi])*Zi_disp))/(UnitTime_in_years);

#ifdef INSTANTANEOUS_RECYCLE //to recover results from instantaneous recycling approximation
    	reset_ejection_rates(i, Gal[p].sfh_ibin,
    			&NormSNIIMassEjecRate_actual, &NormSNIIMetalEjecRate_actual,
    			&NormSNIaMassEjecRate_actual, &NormAGBMassEjecRate_actual,
    			&NormSNIaMetalEjecRate_actual, &NormAGBMetalEjecRate_actual);
#endif //INSTANTANEOUS_RECYCLE

    	//pre-calculations to speed up the code
     	NormMassEjecRateSumAllTypes = NormSNIIMassEjecRate_actual + NormSNIaMassEjecRate_actual + NormAGBMassEjecRate_actual;
#ifdef INDIVIDUAL_ELEMENTS
    	int k;
	    for (k=0;k<NUM_ELEMENTS;k++)
	    {
	    	NormSNIIYieldRate_actual[k] = NormSNIIYieldRate[TimeBin][i][Zi][k] + ((NormSNIIYieldRate[TimeBin][i][Zi+1][k] - NormSNIIYieldRate[TimeBin][i][Zi][k])*Zi_disp);
	    	NormSNIaYieldRate_actual[k] = NormSNIaYieldRate[TimeBin][i][Zi][k] + ((NormSNIaYieldRate[TimeBin][i][Zi+1][k] - NormSNIaYieldRate[TimeBin][i][Zi][k])*Zi_disp);
	    	NormAGBYieldRate_actual[k] = NormAGBYieldRate[TimeBin][i][Zi][k] + ((NormAGBYieldRate[TimeBin][i][Zi+1][k] - NormAGBYieldRate[TimeBin][i][Zi][k])*Zi_disp);
	    }
#endif

#ifdef PORTINARI
	    SNIIEjectaToHot = max(0.0, fwind * step_width_times_DiskSFR * (NormSNIIMetalEjecRate_actual + (Disk_total_metallicity * NormSNIIMassEjecRate_actual)));
	    Gal[p].MetalsHotGas.type2 += SNIIEjectaToHot;
	    Gal[p].MetalsColdGas.type2 += max(0.0, (1.0-fwind) * step_width_times_DiskSFR * (NormSNIIMetalEjecRate_actual + (Disk_total_metallicity * NormSNIIMassEjecRate_actual)));
#endif
#ifdef CHIEFFI
	    SNIIEjectaToHot = max(0.0, fwind * step_width_times_DiskSFR * NormSNIIMetalEjecRate_actual);
	    Gal[p].MetalsHotGas.type2 += SNIIEjectaToHot;
	    Gal[p].MetalsColdGas.type2 += max(0.0, (1.0-fwind) * step_width_times_DiskSFR * NormSNIIMetalEjecRate_actual);
#endif

#ifndef SNIATOHOT
	    Gal[p].HotGas += SNIIEjectaToHot;
    	Gal[p].ColdGas += max(0.0, (step_width_times_DiskSFR * NormMassEjecRateSumAllTypes)-SNIIEjectaToHot);
    	TotalMassReturnedToColdDiskGas += max(0.0, (step_width_times_DiskSFR * NormMassEjecRateSumAllTypes)-SNIIEjectaToHot); //Only use energy from SNe that eject into ColdGas to reheat
	    TotalMassReturnedToHotGas += SNIIEjectaToHot;
	    //TotalMassReturnedToColdDiskGas += max(0.0, step_width_times_DiskSFR * NormMassEjecRateSumAllTypes); //Use energy from ALL SNe (that eject into ColdGas and HotGas) to reheat
#else
	    Gal[p].HotGas += max(0.0, step_width_times_DiskSFR * NormSNIaMassEjecRate_actual) + SNIIEjectaToHot;
	    Gal[p].ColdGas += max(0.0, step_width_times_DiskSFR * (NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)-SNIIEjectaToHot);
	    TotalMassReturnedToColdDiskGas += max(0.0, step_width_times_DiskSFR * (NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)-SNIIEjectaToHot); //Only use energy from SNe that eject into ColdGas to reheat
	    //TotalMassReturnedToColdDiskGas += max(0.0, step_width_times_DiskSFR * (NormMassEjecRateSumAllTypes)); //Use energy from ALL SNe (that eject into ColdGas and HotGas) to reheat
	    TotalMassReturnedToHotGas += max(0.0, step_width_times_DiskSFR * NormSNIaMassEjecRate_actual) + SNIIEjectaToHot;
#endif

#ifndef SNIATOHOT
	    Gal[p].MetalsColdGas.type1a += max(0.0, step_width_times_DiskSFR * NormSNIaMetalEjecRate_actual);
#else
    	Gal[p].MetalsHotGas.type1a += max(0.0, step_width_times_DiskSFR * NormSNIaMetalEjecRate_actual);
#endif
	    Gal[p].MetalsColdGas.agb += max(0.0, step_width_times_DiskSFR * (NormAGBMetalEjecRate_actual + (Disk_total_metallicity * NormAGBMassEjecRate_actual)));

#ifdef INDIVIDUAL_ELEMENTS
#ifdef PORTINARI
#ifndef SNIATOHOT
    		Gal[p].HotGas_elements.H += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[0] + (Gal[p].sfh_ElementsDiskMass[i].H * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual)); //SN-II ejecta to HotGas in metal-rich wind (fwind)
    		Gal[p].ColdGas_elements.H += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[0] + (Gal[p].sfh_ElementsDiskMass[i].H * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual)); //SN-II ejecta to ColdGas (1.0-fwind)
    		Gal[p].ColdGas_elements.H += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[0] + NormAGBYieldRate_actual[0]) + (Gal[p].sfh_ElementsDiskMass[i].H*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual)); //SN-Ia and AGB ejecta to ColdGas
    		Gal[p].HotGas_elements.He += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[1] + (Gal[p].sfh_ElementsDiskMass[i].He * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.He += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[1] + (Gal[p].sfh_ElementsDiskMass[i].He * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.He += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[1] + NormAGBYieldRate_actual[1]) + (Gal[p].sfh_ElementsDiskMass[i].He*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
#ifndef MAINELEMENTS
    		Gal[p].HotGas_elements.Cb += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[2] + (Gal[p].sfh_ElementsDiskMass[i].Cb * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Cb += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[2] + (Gal[p].sfh_ElementsDiskMass[i].Cb * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Cb += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsDiskMass[i].Cb*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.N += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[3] + (Gal[p].sfh_ElementsDiskMass[i].N * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.N += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[3] + (Gal[p].sfh_ElementsDiskMass[i].N * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.N += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsDiskMass[i].N*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.O += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[4] + (Gal[p].sfh_ElementsDiskMass[i].O * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.O += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[4] + (Gal[p].sfh_ElementsDiskMass[i].O * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsDiskMass[i].O*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Ne += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[5] + (Gal[p].sfh_ElementsDiskMass[i].Ne * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Ne += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[5] + (Gal[p].sfh_ElementsDiskMass[i].Ne * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Ne += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[5] + NormAGBYieldRate_actual[5]) + (Gal[p].sfh_ElementsDiskMass[i].Ne*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Mg += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[6] + (Gal[p].sfh_ElementsDiskMass[i].Mg * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Mg += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[6] + (Gal[p].sfh_ElementsDiskMass[i].Mg * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[6] + NormAGBYieldRate_actual[6]) + (Gal[p].sfh_ElementsDiskMass[i].Mg*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Si += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[7] + (Gal[p].sfh_ElementsDiskMass[i].Si * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Si += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[7] + (Gal[p].sfh_ElementsDiskMass[i].Si * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Si += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[7] + NormAGBYieldRate_actual[7]) + (Gal[p].sfh_ElementsDiskMass[i].Si*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.S += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[8] + (Gal[p].sfh_ElementsDiskMass[i].S * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.S += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[8] + (Gal[p].sfh_ElementsDiskMass[i].S * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.S += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[8] + NormAGBYieldRate_actual[8]) + (Gal[p].sfh_ElementsDiskMass[i].S*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Ca += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[9] + (Gal[p].sfh_ElementsDiskMass[i].Ca * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Ca += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[9] + (Gal[p].sfh_ElementsDiskMass[i].Ca * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Ca += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[9] + NormAGBYieldRate_actual[9]) + (Gal[p].sfh_ElementsDiskMass[i].Ca*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Fe += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[10] + (Gal[p].sfh_ElementsDiskMass[i].Fe * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Fe += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[10] + (Gal[p].sfh_ElementsDiskMass[i].Fe * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[10] + NormAGBYieldRate_actual[10]) + (Gal[p].sfh_ElementsDiskMass[i].Fe*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
#else
    		Gal[p].HotGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsDiskMass[i].O*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    		Gal[p].HotGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsDiskMass[i].Mg*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    		Gal[p].HotGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsDiskMass[i].Fe*inverse_DiskMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
#endif //MAINELEMENTS
#ifdef DETAILED_DUST //SNIA TO COLD! //PORTINARI
			SNII_prevstep_Cold_Si[i] += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[7] + (Gal[p].sfh_ElementsDiskMass[i].Si * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
			SNII_prevstep_Cold_Fe[i] += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[10] + (Gal[p].sfh_ElementsDiskMass[i].Fe * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
			SNII_prevstep_Cold_Cb[i] += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[2] + (Gal[p].sfh_ElementsDiskMass[i].Cb * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
			SNIa_prevstep_Cold_Fe[i] += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[10])));
#endif //DETAILED_DUST
#endif //NOT SNIATOHOT
#ifdef SNIATOHOT
    		Gal[p].HotGas_elements.H += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[0] + (Gal[p].sfh_ElementsDiskMass[i].H * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual)); //SN-II ejecta to HotGas in metal-rich wind (fwind)
    		Gal[p].HotGas_elements.H += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[0]); //SN-Ia ejecta to HotGas
    		Gal[p].ColdGas_elements.H += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[0] + (Gal[p].sfh_ElementsDiskMass[i].H * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual)); //SN-II ejecta to ColdGas (1.0-fwind)
    		Gal[p].ColdGas_elements.H += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[0] + (Gal[p].sfh_ElementsDiskMass[i].H*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual)); //AGB ejecta to ColdGas
    		Gal[p].HotGas_elements.He += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[1] + (Gal[p].sfh_ElementsDiskMass[i].He * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].HotGas_elements.He += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[1]);
    		Gal[p].ColdGas_elements.He += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[1] + (Gal[p].sfh_ElementsDiskMass[i].He * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.He += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[1] + (Gal[p].sfh_ElementsDiskMass[i].He*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
#ifndef MAINELEMENTS
    		Gal[p].HotGas_elements.Cb += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[2] + (Gal[p].sfh_ElementsDiskMass[i].Cb * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].HotGas_elements.Cb += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[2]);
    		Gal[p].ColdGas_elements.Cb += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[2] + (Gal[p].sfh_ElementsDiskMass[i].Cb * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Cb += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[2] + (Gal[p].sfh_ElementsDiskMass[i].Cb*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.N += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[3] + (Gal[p].sfh_ElementsDiskMass[i].N * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].HotGas_elements.N += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[3]);
    		Gal[p].ColdGas_elements.N += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[3] + (Gal[p].sfh_ElementsDiskMass[i].N * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.N += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[3] + (Gal[p].sfh_ElementsDiskMass[i].N*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.O += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[4] + (Gal[p].sfh_ElementsDiskMass[i].O * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].HotGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[4]);
    		Gal[p].ColdGas_elements.O += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[4] + (Gal[p].sfh_ElementsDiskMass[i].O * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[4] + (Gal[p].sfh_ElementsDiskMass[i].O*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Ne += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[5] + (Gal[p].sfh_ElementsDiskMass[i].Ne * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].HotGas_elements.Ne += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[5]);
    		Gal[p].ColdGas_elements.Ne += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[5] + (Gal[p].sfh_ElementsDiskMass[i].Ne * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Ne += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[5] + (Gal[p].sfh_ElementsDiskMass[i].Ne*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Mg += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[6] + (Gal[p].sfh_ElementsDiskMass[i].Mg * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].HotGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[6]);
    		Gal[p].ColdGas_elements.Mg += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[6] + (Gal[p].sfh_ElementsDiskMass[i].Mg * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[6] + (Gal[p].sfh_ElementsDiskMass[i].Mg*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Si += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[7] + (Gal[p].sfh_ElementsDiskMass[i].Si * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].HotGas_elements.Si += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[7]);
    		Gal[p].ColdGas_elements.Si += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[7] + (Gal[p].sfh_ElementsDiskMass[i].Si * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Si += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[7] + (Gal[p].sfh_ElementsDiskMass[i].Si*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.S += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[8] + (Gal[p].sfh_ElementsDiskMass[i].S * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].HotGas_elements.S += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[8]);
    		Gal[p].ColdGas_elements.S += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[8] + (Gal[p].sfh_ElementsDiskMass[i].S * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.S += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[8] + (Gal[p].sfh_ElementsDiskMass[i].S*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Ca += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[9] + (Gal[p].sfh_ElementsDiskMass[i].Ca * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].HotGas_elements.Ca += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[9]);
    		Gal[p].ColdGas_elements.Ca += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[9] + (Gal[p].sfh_ElementsDiskMass[i].Ca * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Ca += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[9] + (Gal[p].sfh_ElementsDiskMass[i].Ca*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Fe += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[10] + (Gal[p].sfh_ElementsDiskMass[i].Fe * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].HotGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[10]);
    		Gal[p].ColdGas_elements.Fe += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[10] + (Gal[p].sfh_ElementsDiskMass[i].Fe * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
    		Gal[p].ColdGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[10] + (Gal[p].sfh_ElementsDiskMass[i].Fe*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
#else //MAINELEMENTS
       		Gal[p].HotGas_elements.O += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[2] + (Gal[p].sfh_ElementsDiskMass[i].O * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
        	Gal[p].HotGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[2]);
        	Gal[p].ColdGas_elements.O += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[2] + (Gal[p].sfh_ElementsDiskMass[i].O * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
        	Gal[p].ColdGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[2] + (Gal[p].sfh_ElementsDiskMass[i].O*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
       		Gal[p].HotGas_elements.Mg += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[3] + (Gal[p].sfh_ElementsDiskMass[i].Mg * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
        	Gal[p].HotGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[3]);
        	Gal[p].ColdGas_elements.Mg += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[3] + (Gal[p].sfh_ElementsDiskMass[i].Mg * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
        	Gal[p].ColdGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[3] + (Gal[p].sfh_ElementsDiskMass[i].Mg*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
       		Gal[p].HotGas_elements.Fe += max(0.0, fwind * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[4] + (Gal[p].sfh_ElementsDiskMass[i].Fe * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
        	Gal[p].HotGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[4]);
        	Gal[p].ColdGas_elements.Fe += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[4] + (Gal[p].sfh_ElementsDiskMass[i].Fe * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
        	Gal[p].ColdGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[4] + (Gal[p].sfh_ElementsDiskMass[i].Fe*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
#endif //MAINELEMENTS
#ifdef DETAILED_DUST //SNIA TO HOT! //PORTINARI
			SNII_prevstep_Cold_Si[i] += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[7] + (Gal[p].sfh_ElementsDiskMass[i].Si * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
			SNII_prevstep_Cold_Fe[i] += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[10] + (Gal[p].sfh_ElementsDiskMass[i].Fe * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
			SNII_prevstep_Cold_Cb[i] += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * (NormSNIIYieldRate_actual[2] + (Gal[p].sfh_ElementsDiskMass[i].Cb * inverse_DiskMass_physical_units) * NormSNIIMassEjecRate_actual));
			SNIa_prevstep_Cold_Fe[i] += 0.0; //SNIA has gone to hot gas
#endif //DETAILED_DUST
#endif //SNIATOHOT
#endif //PORTINARI

#ifdef CHIEFFI
#ifndef SNIATOHOT
    		Gal[p].HotGas_elements.H += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[0]); //SN-II ejecta to HotGas in metal-rich wind (fwind) //NB: No unsynth component required for SN-II ejecta when using the CL04 SN-II yields
    		Gal[p].ColdGas_elements.H += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[0]); //SN-II ejecta to ColdGas (1.0-fwind)
    		Gal[p].ColdGas_elements.H += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[0] + NormAGBYieldRate_actual[0]) + (Gal[p].sfh_ElementsDiskMass[i].H*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual)); //SN-Ia and AGB ejecta to ColdGas
    		Gal[p].HotGas_elements.He += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[1]);
    		Gal[p].ColdGas_elements.He += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[1]);
    		Gal[p].ColdGas_elements.He += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[1] + NormAGBYieldRate_actual[1]) + (Gal[p].sfh_ElementsDiskMass[i].He*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
#ifndef MAINELEMENTS
    		Gal[p].HotGas_elements.Cb += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[2]);
    		Gal[p].ColdGas_elements.Cb += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[2]);
    		Gal[p].ColdGas_elements.Cb += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsDiskMass[i].Cb*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.N += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[3]);
    		Gal[p].ColdGas_elements.N += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[3]);
    		Gal[p].ColdGas_elements.N += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsDiskMass[i].N*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.O += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[4]);
    		Gal[p].ColdGas_elements.O += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[4]);
    		Gal[p].ColdGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsDiskMass[i].O*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Ne += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[5]);
    		Gal[p].ColdGas_elements.Ne += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[5]);
    		Gal[p].ColdGas_elements.Ne += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[5] + NormAGBYieldRate_actual[5]) + (Gal[p].sfh_ElementsDiskMass[i].Ne*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Mg += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[6]);
    		Gal[p].ColdGas_elements.Mg += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[6]);
    		Gal[p].ColdGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[6] + NormAGBYieldRate_actual[6]) + (Gal[p].sfh_ElementsDiskMass[i].Mg*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Si += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[7]);
    		Gal[p].ColdGas_elements.Si += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[7]);
    		Gal[p].ColdGas_elements.Si += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[7] + NormAGBYieldRate_actual[7]) + (Gal[p].sfh_ElementsDiskMass[i].Si*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.S += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[8]);
    		Gal[p].ColdGas_elements.S += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[8]);
    		Gal[p].ColdGas_elements.S += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[8] + NormAGBYieldRate_actual[8]) + (Gal[p].sfh_ElementsDiskMass[i].S*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Ca += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[9]);
    		Gal[p].ColdGas_elements.Ca += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[9]);
    		Gal[p].ColdGas_elements.Ca += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[9] + NormAGBYieldRate_actual[9]) + (Gal[p].sfh_ElementsDiskMass[i].Ca*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Fe += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[10]);
    		Gal[p].ColdGas_elements.Fe += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[10]);
    		Gal[p].ColdGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[10] + NormAGBYieldRate_actual[10]) + (Gal[p].sfh_ElementsDiskMass[i].Fe*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
#else //MAINELEMENTS
    		Gal[p].HotGas_elements.O += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[2]);
    		Gal[p].ColdGas_elements.O += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[2]);
    		Gal[p].ColdGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsDiskMass[i].O*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Mg += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[3]);
    		Gal[p].ColdGas_elements.Mg += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[3]);
    		Gal[p].ColdGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsDiskMass[i].Mg*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Fe += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[4]);
    		Gal[p].ColdGas_elements.Fe += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[4]);
    		Gal[p].ColdGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units * ((NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsDiskMass[i].Fe*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
#endif //MAINELEMENTS
#ifdef DETAILED_DUST //SNIA TO COLD //CHIEFFI 
			SNII_prevstep_Cold_Si[i] += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[7]);
			SNII_prevstep_Cold_Fe[i] += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[10]);
			SNII_prevstep_Cold_Cb[i] += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[2]);
			SNIa_prevstep_Cold_Fe[i] += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[10]);
#endif //DDust
#endif //SNIATOHOT
#ifdef SNIATOHOT
    		Gal[p].HotGas_elements.H += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[0]); //SN-II ejecta to HotGas in metal-rich wind (fwind) //NB: No unsynth component required for SN-II ejecta when using the CL04 SN-II yields
    		Gal[p].HotGas_elements.H += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[0]); //SN-Ia ejecta to HotGas
    		Gal[p].ColdGas_elements.H += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[0]); //SN-II ejecta to ColdGas (1.0-fwind)
    		Gal[p].ColdGas_elements.H += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[0] + (Gal[p].sfh_ElementsDiskMass[i].H*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual)); //AGB ejecta to ColdGas
    		Gal[p].HotGas_elements.He += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[1]);
    		Gal[p].HotGas_elements.He += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[1]);
    		Gal[p].ColdGas_elements.He += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[1]);
    		Gal[p].ColdGas_elements.He += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[1] + (Gal[p].sfh_ElementsDiskMass[i].He*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
#ifndef MAINELEMENTS
    		Gal[p].HotGas_elements.Cb += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[2]);
        	Gal[p].HotGas_elements.Cb += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[2]);
        	Gal[p].ColdGas_elements.Cb += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[2]);
        	Gal[p].ColdGas_elements.Cb += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[2] + (Gal[p].sfh_ElementsDiskMass[i].Cb*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.N += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[3]);
    		Gal[p].HotGas_elements.N += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[3]);
    		Gal[p].ColdGas_elements.N += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[3]);
    		Gal[p].ColdGas_elements.N += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[3] + (Gal[p].sfh_ElementsDiskMass[i].N*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.O += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[4]);
    		Gal[p].HotGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[4]);
    		Gal[p].ColdGas_elements.O += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[4]);
    		Gal[p].ColdGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[4] + (Gal[p].sfh_ElementsDiskMass[i].O*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Ne += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[5]);
    		Gal[p].HotGas_elements.Ne += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[5]);
    		Gal[p].ColdGas_elements.Ne += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[5]);
    		Gal[p].ColdGas_elements.Ne += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[5] + (Gal[p].sfh_ElementsDiskMass[i].Ne*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Mg += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[6]);
    		Gal[p].HotGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[6]);
    		Gal[p].ColdGas_elements.Mg += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[6]);
    		Gal[p].ColdGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[6] + (Gal[p].sfh_ElementsDiskMass[i].Mg*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Si += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[7]);
    		Gal[p].HotGas_elements.Si += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[7]);
    		Gal[p].ColdGas_elements.Si += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[7]);
    		Gal[p].ColdGas_elements.Si += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[7] + (Gal[p].sfh_ElementsDiskMass[i].Si*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.S += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[8]);
    		Gal[p].HotGas_elements.S += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[8]);
    		Gal[p].ColdGas_elements.S += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[8]);
    		Gal[p].ColdGas_elements.S += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[8] + (Gal[p].sfh_ElementsDiskMass[i].S*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Ca += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[9]);
    		Gal[p].HotGas_elements.Ca += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[9]);
    		Gal[p].ColdGas_elements.Ca += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[9]);
    		Gal[p].ColdGas_elements.Ca += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[9] + (Gal[p].sfh_ElementsDiskMass[i].Ca*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Fe += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[10]);
    		Gal[p].HotGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[10]);
    		Gal[p].ColdGas_elements.Fe += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[10]);
    		Gal[p].ColdGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[10] + (Gal[p].sfh_ElementsDiskMass[i].Fe*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
#else //MAINELEMENTS
    		Gal[p].HotGas_elements.O += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[2]);
    		Gal[p].HotGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[2]);
    		Gal[p].ColdGas_elements.O += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[2]);
    		Gal[p].ColdGas_elements.O += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[2] + (Gal[p].sfh_ElementsDiskMass[i].O*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Mg += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[3]);
    		Gal[p].HotGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[3]);
    		Gal[p].ColdGas_elements.Mg += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[3]);
    		Gal[p].ColdGas_elements.Mg += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[3] + (Gal[p].sfh_ElementsDiskMass[i].Mg*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
    		Gal[p].HotGas_elements.Fe += max(0.0, fwind * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[4]);
    		Gal[p].HotGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units * NormSNIaYieldRate_actual[4]);
    		Gal[p].ColdGas_elements.Fe += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[4]);
    		Gal[p].ColdGas_elements.Fe += max(0.0, step_width_times_DiskSFR_physical_units * (NormAGBYieldRate_actual[4] + (Gal[p].sfh_ElementsDiskMass[i].Fe*inverse_DiskMass_physical_units)*NormAGBMassEjecRate_actual));
#endif //MAINELEMENTS
#ifdef DETAILED_DUST //Inside SNIA TO HOT //CHIEFFI 
			SNII_prevstep_Cold_Si[i] += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[7]);
			SNII_prevstep_Cold_Fe[i] += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[10]);
			SNII_prevstep_Cold_Cb[i] += max(0.0, (1.0-fwind) * step_width_times_DiskSFR_physical_units * NormSNIIYieldRate_actual[2]);
			SNIa_prevstep_Cold_Fe[i] += 0.0;
#endif //DDust
#endif //SNIATOHOT
#endif //CHIEFFI
#endif //INDIVIDUAL_ELEMENTS
    	//UPDATE DISK MASS COMPONENTS:
    	/*ROB (13-02-13): All the mass/metals/elements in the stars that die in this timestep are lost from the stellar component.
    	//i.e. All the mass/metals/elements in the stars at birth are removed...
    	//...Some goes to the gas (+ newly synthesised component), the rest goes into the 'stellar remnants' which are not tracked and do not contribute to the stellar component's mass/metals/elements budget.*/
    	Gal[p].DiskMass -= max(0.0, step_width_times_DiskSFR * NormMassEjecRateSumAllTypes);
    	Gal[p].MetalsDiskMass.type2 -= max(0.0, step_width_times_DiskSFR * (Disk_total_metallicity * NormSNIIMassEjecRate_actual));
    	Gal[p].MetalsDiskMass.type1a -= max(0.0, step_width_times_DiskSFR * (Disk_total_metallicity * NormSNIaMassEjecRate_actual));
	    Gal[p].MetalsDiskMass.agb -= max(0.0, step_width_times_DiskSFR * (Disk_total_metallicity * NormAGBMassEjecRate_actual));

#ifdef INDIVIDUAL_ELEMENTS
	    Gal[p].DiskMass_elements.H -= max(0.0, step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].H*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes));
	    Gal[p].DiskMass_elements.He -= max(0.0, step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].He*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes));
#ifndef MAINELEMENTS
	    Gal[p].DiskMass_elements.Cb -= max(0.0, step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].Cb*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes));
	    Gal[p].DiskMass_elements.N -= max(0.0, step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].N*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes));
#endif
	    Gal[p].DiskMass_elements.O -= max(0.0, step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].O*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes));
#ifndef MAINELEMENTS
	    Gal[p].DiskMass_elements.Ne -= max(0.0, step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].Ne*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes));
#endif
	    Gal[p].DiskMass_elements.Mg -= max(0.0, step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].Mg*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes));
#ifndef MAINELEMENTS
	    Gal[p].DiskMass_elements.Si -= max(0.0, step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].Si*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes));
	    Gal[p].DiskMass_elements.S -= max(0.0, step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].S*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes));
	    Gal[p].DiskMass_elements.Ca -= max(0.0, step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].Ca*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes));
#endif
	    Gal[p].DiskMass_elements.Fe -= max(0.0, step_width_times_DiskSFR_physical_units * ((Gal[p].sfh_ElementsDiskMass[i].Fe*inverse_DiskMass_physical_units)*NormMassEjecRateSumAllTypes));
#endif //INDIVIDUAL_ELEMENTS
	    //Update ages:
    	for(n=0;n<NOUT;n++)
    	{
    		AgeCorrectionDisk[n] += max(0.0, (sfh_time-NumToTime(ListOutputSnaps[n]))*(step_width_times_DiskSFR * NormMassEjecRateSumAllTypes));
    		if (AgeCorrectionDisk[n] < 0.0) AgeCorrectionDisk[n] = 0.0;
    	}
    } //if (Gal[p].sfh_DiskMass[i] > 0.0)

    //*****************************************
    //ENRICHMENT FROM BULGE STARS INTO HOT GAS:
    //*****************************************
    if (Gal[p].BulgeMass > 0.0)
    if (Gal[p].sfh_BulgeMass[i] > 0.0)
    {
    	//pre-calculations to speed up the code
    	BulgeSFR = Gal[p].sfh_BulgeMass[i]/Gal[p].sfh_dt[i];
    	step_width_times_BulgeSFR = timestep_width * BulgeSFR;
    	BulgeSFR_physical_units = BulgeSFR * (1.0e10/Hubble_h);
    	step_width_times_BulgeSFR_physical_units = timestep_width * BulgeSFR_physical_units;
    	inverse_BulgeMass_physical_units=Hubble_h/(Gal[p].sfh_BulgeMass[i]*1.0e10);
    	Bulge_total_metallicity=metals_total(Gal[p].sfh_MetalsBulgeMass[i])/Gal[p].sfh_BulgeMass[i];
    	Zi = find_initial_metallicity(p, i, 1, 2);
    	//Interpolate the bulge luminosity on the lifetimeMetallicities tables:
    	Zi_disp = (Bulge_total_metallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.

    	NormSNIIMassEjecRate_actual = NormSNIIMassEjecRate[TimeBin][i][Zi] + ((NormSNIIMassEjecRate[TimeBin][i][Zi+1] - NormSNIIMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIaMassEjecRate_actual = NormSNIaMassEjecRate[TimeBin][i][Zi] + ((NormSNIaMassEjecRate[TimeBin][i][Zi+1] - NormSNIaMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormAGBMassEjecRate_actual = NormAGBMassEjecRate[TimeBin][i][Zi] + ((NormAGBMassEjecRate[TimeBin][i][Zi+1] - NormAGBMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIIMetalEjecRate_actual = NormSNIIMetalEjecRate[TimeBin][i][Zi] + ((NormSNIIMetalEjecRate[TimeBin][i][Zi+1] - NormSNIIMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIaMetalEjecRate_actual = NormSNIaMetalEjecRate[TimeBin][i][Zi] + ((NormSNIaMetalEjecRate[TimeBin][i][Zi+1] - NormSNIaMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormAGBMetalEjecRate_actual = NormAGBMetalEjecRate[TimeBin][i][Zi] + ((NormAGBMetalEjecRate[TimeBin][i][Zi+1] - NormAGBMetalEjecRate[TimeBin][i][Zi])*Zi_disp);

    	//pre-calculations to speed up the code
    	NormMassEjecRateSumAllTypes = NormSNIIMassEjecRate_actual + NormSNIaMassEjecRate_actual + NormAGBMassEjecRate_actual;
#ifdef INDIVIDUAL_ELEMENTS
    	int k;
	    for (k=0;k<NUM_ELEMENTS;k++)
	    {
	    	NormSNIIYieldRate_actual[k] = NormSNIIYieldRate[TimeBin][i][Zi][k] + ((NormSNIIYieldRate[TimeBin][i][Zi+1][k] - NormSNIIYieldRate[TimeBin][i][Zi][k])*Zi_disp);
	    	NormSNIaYieldRate_actual[k] = NormSNIaYieldRate[TimeBin][i][Zi][k] + ((NormSNIaYieldRate[TimeBin][i][Zi+1][k] - NormSNIaYieldRate[TimeBin][i][Zi][k])*Zi_disp);
	    	NormAGBYieldRate_actual[k] = NormAGBYieldRate[TimeBin][i][Zi][k] + ((NormAGBYieldRate[TimeBin][i][Zi+1][k] - NormAGBYieldRate[TimeBin][i][Zi][k])*Zi_disp);
	    }
#endif //INDIVIDUAL_ELEMENTS

#ifdef INSTANTANEOUS_RECYCLE //to recover results from instantaneous recycling approximation
    	reset_ejection_rates(i, Gal[p].sfh_ibin,
    			&NormSNIIMassEjecRate_actual, &NormSNIIMetalEjecRate_actual,
    			&NormSNIaMassEjecRate_actual, &NormAGBMassEjecRate_actual,
    			&NormSNIaMetalEjecRate_actual, &NormAGBMetalEjecRate_actual);
#endif //INSTANTANEOUS_RECYCLE

    	//UPDATE HOT GAS COMPONENTS:
#ifndef BULGE_TO_COLD
    	Gal[p].HotGas += max(0.0, step_width_times_BulgeSFR * NormMassEjecRateSumAllTypes);
    	TotalMassReturnedToHotGas += max(0.0, step_width_times_BulgeSFR * NormMassEjecRateSumAllTypes);
#ifdef PORTINARI
    	Gal[p].MetalsHotGas.type2 += max(0.0, step_width_times_BulgeSFR * (NormSNIIMetalEjecRate_actual + (Bulge_total_metallicity * NormSNIIMassEjecRate_actual)));
#endif
#ifdef CHIEFFI
    	Gal[p].MetalsHotGas.type2 += max(0.0, step_width_times_BulgeSFR * NormSNIIMetalEjecRate_actual);
#endif
    	//Gal[p].MetalsHotGas.type1a += step_width_times_BulgeSFR * (NormSNIaMetalEjecRate_actual + (Bulge_total_metallicity * NormSNIaMassEjecRate_actual));
    	Gal[p].MetalsHotGas.type1a += max(0.0, step_width_times_BulgeSFR * NormSNIaMetalEjecRate_actual);
    	Gal[p].MetalsHotGas.agb += max(0.0, step_width_times_BulgeSFR * (NormAGBMetalEjecRate_actual + (Bulge_total_metallicity * NormAGBMassEjecRate_actual)));

#ifdef INDIVIDUAL_ELEMENTS
#ifdef PORTINARI
    	Gal[p].HotGas_elements.H += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[0] + NormSNIaYieldRate_actual[0] + NormAGBYieldRate_actual[0]) + (Gal[p].sfh_ElementsBulgeMass[i].H*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.He += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[1] + NormSNIaYieldRate_actual[1] + NormAGBYieldRate_actual[1]) + (Gal[p].sfh_ElementsBulgeMass[i].He*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
#ifndef MAINELEMENTS
    	Gal[p].HotGas_elements.Cb += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsBulgeMass[i].Cb*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.N += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsBulgeMass[i].N*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsBulgeMass[i].O*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Ne += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[5] + NormSNIaYieldRate_actual[5] + NormAGBYieldRate_actual[5]) + (Gal[p].sfh_ElementsBulgeMass[i].Ne*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[6] + NormSNIaYieldRate_actual[6] + NormAGBYieldRate_actual[6]) + (Gal[p].sfh_ElementsBulgeMass[i].Mg*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Si += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[7] + NormSNIaYieldRate_actual[7] + NormAGBYieldRate_actual[7]) + (Gal[p].sfh_ElementsBulgeMass[i].Si*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.S += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[8] + NormSNIaYieldRate_actual[8] + NormAGBYieldRate_actual[8]) + (Gal[p].sfh_ElementsBulgeMass[i].S*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Ca += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[9] + NormSNIaYieldRate_actual[9] + NormAGBYieldRate_actual[9]) + (Gal[p].sfh_ElementsBulgeMass[i].Ca*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[10] + NormSNIaYieldRate_actual[10] + NormAGBYieldRate_actual[10]) + (Gal[p].sfh_ElementsBulgeMass[i].Fe*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
#else
    	Gal[p].HotGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsBulgeMass[i].O*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsBulgeMass[i].Mg*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsBulgeMass[i].Fe*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
#endif //MAINELEMENTS
#endif //PORTINARI
#ifdef CHIEFFI
    	Gal[p].HotGas_elements.H += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[0] + NormSNIaYieldRate_actual[0] + NormAGBYieldRate_actual[0]) + (Gal[p].sfh_ElementsBulgeMass[i].H*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));  //ROB: No unsynth component required for SN-II ejecta, when using the Chieffi & Limongi 92007) yield tables/
    	Gal[p].HotGas_elements.He += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[1] + NormSNIaYieldRate_actual[1] + NormAGBYieldRate_actual[1]) + (Gal[p].sfh_ElementsBulgeMass[i].He*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
#ifndef MAINELEMENTS
    	Gal[p].HotGas_elements.Cb += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsBulgeMass[i].Cb*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.N += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsBulgeMass[i].N*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsBulgeMass[i].O*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Ne += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[5] + NormSNIaYieldRate_actual[5] + NormAGBYieldRate_actual[5]) + (Gal[p].sfh_ElementsBulgeMass[i].Ne*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[6] + NormSNIaYieldRate_actual[6] + NormAGBYieldRate_actual[6]) + (Gal[p].sfh_ElementsBulgeMass[i].Mg*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Si += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[7] + NormSNIaYieldRate_actual[7] + NormAGBYieldRate_actual[7]) + (Gal[p].sfh_ElementsBulgeMass[i].Si*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.S += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[8] + NormSNIaYieldRate_actual[8] + NormAGBYieldRate_actual[8]) + (Gal[p].sfh_ElementsBulgeMass[i].S*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Ca += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[9] + NormSNIaYieldRate_actual[9] + NormAGBYieldRate_actual[9]) + (Gal[p].sfh_ElementsBulgeMass[i].Ca*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[10] + NormSNIaYieldRate_actual[10] + NormAGBYieldRate_actual[10]) + (Gal[p].sfh_ElementsBulgeMass[i].Fe*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
#else
    	Gal[p].HotGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsBulgeMass[i].O*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsBulgeMass[i].Mg*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsBulgeMass[i].Fe*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
#endif //MAINELEMENTS
#endif //CHIEFFI
#endif //INDIVIDUAL_ELEMENTS
#else //BULGE_TO_COLD
    	Gal[p].ColdGas += max(0.0, step_width_times_BulgeSFR * NormMassEjecRateSumAllTypes);
    	TotalMassReturnedToColdDiskGas += max(0.0, step_width_times_BulgeSFR * NormMassEjecRateSumAllTypes);
#ifdef PORTINARI
    	Gal[p].MetalsColdGas.type2 += max(0.0, step_width_times_BulgeSFR * (NormSNIIMetalEjecRate_actual + (Bulge_total_metallicity * NormSNIIMassEjecRate_actual)));
#endif
#ifdef CHIEFFI
    	Gal[p].MetalsColdGas.type2 += max(0.0, step_width_times_BulgeSFR * NormSNIIMetalEjecRate_actual);
#endif
    	//Gal[p].MetalsColdGas.type1a += step_width_times_BulgeSFR * (NormSNIaMetalEjecRate_actual + (Bulge_total_metallicity * NormSNIaMassEjecRate_actual));
    	Gal[p].MetalsColdGas.type1a += max(0.0, step_width_times_BulgeSFR * NormSNIaMetalEjecRate_actual);
    	Gal[p].MetalsColdGas.agb += max(0.0, step_width_times_BulgeSFR * (NormAGBMetalEjecRate_actual + (Bulge_total_metallicity * NormAGBMassEjecRate_actual)));

#ifdef INDIVIDUAL_ELEMENTS
#ifdef PORTINARI
    	Gal[p].ColdGas_elements.H += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[0] + NormSNIaYieldRate_actual[0] + NormAGBYieldRate_actual[0]) + (Gal[p].sfh_ElementsBulgeMass[i].H*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.He += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[1] + NormSNIaYieldRate_actual[1] + NormAGBYieldRate_actual[1]) + (Gal[p].sfh_ElementsBulgeMass[i].He*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
#ifndef MAINELEMENTS
    	Gal[p].ColdGas_elements.Cb += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsBulgeMass[i].Cb*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.N += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsBulgeMass[i].N*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsBulgeMass[i].O*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.Ne += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[5] + NormSNIaYieldRate_actual[5] + NormAGBYieldRate_actual[5]) + (Gal[p].sfh_ElementsBulgeMass[i].Ne*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[6] + NormSNIaYieldRate_actual[6] + NormAGBYieldRate_actual[6]) + (Gal[p].sfh_ElementsBulgeMass[i].Mg*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.Si += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[7] + NormSNIaYieldRate_actual[7] + NormAGBYieldRate_actual[7]) + (Gal[p].sfh_ElementsBulgeMass[i].Si*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.S += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[8] + NormSNIaYieldRate_actual[8] + NormAGBYieldRate_actual[8]) + (Gal[p].sfh_ElementsBulgeMass[i].S*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.Ca += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[9] + NormSNIaYieldRate_actual[9] + NormAGBYieldRate_actual[9]) + (Gal[p].sfh_ElementsBulgeMass[i].Ca*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[10] + NormSNIaYieldRate_actual[10] + NormAGBYieldRate_actual[10]) + (Gal[p].sfh_ElementsBulgeMass[i].Fe*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
#else
    	Gal[p].ColdGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsBulgeMass[i].O*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsBulgeMass[i].Mg*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsBulgeMass[i].Fe*inverse_BulgeMass_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
#endif //MAINELEMENTS
#endif //PORTINARI
#ifdef CHIEFFI
    	Gal[p].ColdGas_elements.H += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[0] + NormSNIaYieldRate_actual[0] + NormAGBYieldRate_actual[0]) + (Gal[p].sfh_ElementsBulgeMass[i].H*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));  //ROB: No unsynth component required for SN-II ejecta, when using the Chieffi & Limongi 92007) yield tables/
    	Gal[p].ColdGas_elements.He += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[1] + NormSNIaYieldRate_actual[1] + NormAGBYieldRate_actual[1]) + (Gal[p].sfh_ElementsBulgeMass[i].He*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
#ifndef MAINELEMENTS
    	Gal[p].ColdGas_elements.Cb += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsBulgeMass[i].Cb*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.N += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsBulgeMass[i].N*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsBulgeMass[i].O*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.Ne += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[5] + NormSNIaYieldRate_actual[5] + NormAGBYieldRate_actual[5]) + (Gal[p].sfh_ElementsBulgeMass[i].Ne*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[6] + NormSNIaYieldRate_actual[6] + NormAGBYieldRate_actual[6]) + (Gal[p].sfh_ElementsBulgeMass[i].Mg*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.Si += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[7] + NormSNIaYieldRate_actual[7] + NormAGBYieldRate_actual[7]) + (Gal[p].sfh_ElementsBulgeMass[i].Si*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.S += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[8] + NormSNIaYieldRate_actual[8] + NormAGBYieldRate_actual[8]) + (Gal[p].sfh_ElementsBulgeMass[i].S*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.Ca += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[9] + NormSNIaYieldRate_actual[9] + NormAGBYieldRate_actual[9]) + (Gal[p].sfh_ElementsBulgeMass[i].Ca*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[10] + NormSNIaYieldRate_actual[10] + NormAGBYieldRate_actual[10]) + (Gal[p].sfh_ElementsBulgeMass[i].Fe*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
#else
    	Gal[p].ColdGas_elements.O += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsBulgeMass[i].O*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.Mg += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsBulgeMass[i].Mg*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].ColdGas_elements.Fe += max(0.0, step_width_times_BulgeSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsBulgeMass[i].Fe*inverse_BulgeMass_physical_units)*(NormAGBMassEjecRate_actual)));
#endif //MAINELEMENTS
#endif //CHIEFFI
#endif //INDIVIDUAL_ELEMENTS
#endif //BULGE_TO_COLD
    	//UPDATE BULGE MASS COMPONENTS:
    	Gal[p].BulgeMass -= max(0.0, step_width_times_BulgeSFR * NormMassEjecRateSumAllTypes);
    	Gal[p].MetalsBulgeMass.type2 -= max(0.0, step_width_times_BulgeSFR * (Bulge_total_metallicity * NormSNIIMassEjecRate_actual));
    	Gal[p].MetalsBulgeMass.type1a -= max(0.0, step_width_times_BulgeSFR * (Bulge_total_metallicity * NormSNIaMassEjecRate_actual));
    	Gal[p].MetalsBulgeMass.agb -= max(0.0, step_width_times_BulgeSFR * (Bulge_total_metallicity * NormAGBMassEjecRate_actual));

#ifdef INDIVIDUAL_ELEMENTS
    	Gal[p].BulgeMass_elements.H -= max(0.0, step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].H*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes));
    	Gal[p].BulgeMass_elements.He -= max(0.0, step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].He*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes));
#ifndef MAINELEMENTS
    	Gal[p].BulgeMass_elements.Cb -= max(0.0, step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].Cb*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes));
    	Gal[p].BulgeMass_elements.N -= max(0.0, step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].N*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes));
#endif
    	Gal[p].BulgeMass_elements.O -= max(0.0, step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].O*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes));
#ifndef MAINELEMENTS
    	Gal[p].BulgeMass_elements.Ne -= max(0.0, step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].Ne*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes));
#endif
    	Gal[p].BulgeMass_elements.Mg -= max(0.0, step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].Mg*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes));
#ifndef MAINELEMENTS
    	Gal[p].BulgeMass_elements.Si -= max(0.0, step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].Si*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes));
    	Gal[p].BulgeMass_elements.S -= max(0.0, step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].S*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes));
    	Gal[p].BulgeMass_elements.Ca -= max(0.0, step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].Ca*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes));
#endif
    	Gal[p].BulgeMass_elements.Fe -= max(0.0, step_width_times_BulgeSFR_physical_units * ((Gal[p].sfh_ElementsBulgeMass[i].Fe*inverse_BulgeMass_physical_units)*NormMassEjecRateSumAllTypes));
#endif //INDIVIDUAL_ELEMENTS
    	//Update ages:
        for(n=0;n<NOUT;n++)
        {
        	AgeCorrectionBulge[n] += max(0.0, (sfh_time-NumToTime(ListOutputSnaps[n]))*(step_width_times_BulgeSFR * NormMassEjecRateSumAllTypes));
        	if (AgeCorrectionBulge[n] < 0.0) AgeCorrectionBulge[n] = 0.0;
        }
    } //if (Gal[p].sfh_BulgeMass[i] > 0.0) //BULGE


    //*****************************************
    //ENRICHMENT FROM ICL STARS INTO HOT GAS:
    //*****************************************

    if (Gal[p].sfh_ICM[i] > 0.0)
    {
    	//pre-calculations to speed up the code
    	ICMSFR = Gal[p].sfh_ICM[i]/Gal[p].sfh_dt[i];
    	step_width_times_ICMSFR = timestep_width * ICMSFR;
    	ICMSFR_physical_units = ICMSFR * (1.0e10/Hubble_h);
    	step_width_times_ICMSFR_physical_units = timestep_width * ICMSFR_physical_units;
    	inverse_ICM_physical_units=Hubble_h/(Gal[p].sfh_ICM[i]*1.0e10);
    	ICM_total_metallicity=metals_total(Gal[p].sfh_MetalsICM[i])/Gal[p].sfh_ICM[i];

    	Zi = find_initial_metallicity(p, i, 1, 3);
    	//Interpolate the ICM metallicity on the lifetimeMetallicities tables:
    	Zi_disp = (ICM_total_metallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.

    	NormSNIIMassEjecRate_actual = NormSNIIMassEjecRate[TimeBin][i][Zi] + ((NormSNIIMassEjecRate[TimeBin][i][Zi+1] - NormSNIIMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIaMassEjecRate_actual = NormSNIaMassEjecRate[TimeBin][i][Zi] + ((NormSNIaMassEjecRate[TimeBin][i][Zi+1] - NormSNIaMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormAGBMassEjecRate_actual = NormAGBMassEjecRate[TimeBin][i][Zi] + ((NormAGBMassEjecRate[TimeBin][i][Zi+1] - NormAGBMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIIMetalEjecRate_actual = NormSNIIMetalEjecRate[TimeBin][i][Zi] + ((NormSNIIMetalEjecRate[TimeBin][i][Zi+1] - NormSNIIMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIaMetalEjecRate_actual = NormSNIaMetalEjecRate[TimeBin][i][Zi] + ((NormSNIaMetalEjecRate[TimeBin][i][Zi+1] - NormSNIaMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormAGBMetalEjecRate_actual = NormAGBMetalEjecRate[TimeBin][i][Zi] + ((NormAGBMetalEjecRate[TimeBin][i][Zi+1] - NormAGBMetalEjecRate[TimeBin][i][Zi])*Zi_disp);

    	//pre-calculations to speed up the code
    	NormMassEjecRateSumAllTypes = NormSNIIMassEjecRate_actual + NormSNIaMassEjecRate_actual + NormAGBMassEjecRate_actual;

#ifdef INDIVIDUAL_ELEMENTS
    	int k;
	    for (k=0;k<NUM_ELEMENTS;k++)
	    {
	    	NormSNIIYieldRate_actual[k] = NormSNIIYieldRate[TimeBin][i][Zi][k] + ((NormSNIIYieldRate[TimeBin][i][Zi+1][k] - NormSNIIYieldRate[TimeBin][i][Zi][k])*Zi_disp);
	    	NormSNIaYieldRate_actual[k] = NormSNIaYieldRate[TimeBin][i][Zi][k] + ((NormSNIaYieldRate[TimeBin][i][Zi+1][k] - NormSNIaYieldRate[TimeBin][i][Zi][k])*Zi_disp);
	    	NormAGBYieldRate_actual[k] = NormAGBYieldRate[TimeBin][i][Zi][k] + ((NormAGBYieldRate[TimeBin][i][Zi+1][k] - NormAGBYieldRate[TimeBin][i][Zi][k])*Zi_disp);
	    }
#endif //INDIVIDUAL_ELEMENTS

#ifdef INSTANTANEOUS_RECYCLE //to recover results from instantaneous recycling approximation
    	reset_ejection_rates(i, Gal[p].sfh_ibin,
    			&NormSNIIMassEjecRate_actual, &NormSNIIMetalEjecRate_actual,
    			&NormSNIaMassEjecRate_actual, &NormAGBMassEjecRate_actual,
    			&NormSNIaMetalEjecRate_actual, &NormAGBMetalEjecRate_actual);
#endif //INSTANTANEOUS_RECYCLE

    	//UPDATE HOT GAS COMPONENTS:
    	Gal[p].HotGas += max(0.0, step_width_times_ICMSFR * NormMassEjecRateSumAllTypes);
    	TotalMassReturnedToHotGas += max(0.0, step_width_times_ICMSFR * NormMassEjecRateSumAllTypes);
#ifdef PORTINARI
    	Gal[p].MetalsHotGas.type2 += max(0.0, step_width_times_ICMSFR * (NormSNIIMetalEjecRate_actual + (ICM_total_metallicity * NormSNIIMassEjecRate_actual)));
#endif
#ifdef CHIEFFI
    	Gal[p].MetalsHotGas.type2 += max(0.0, step_width_times_ICMSFR * NormSNIIMetalEjecRate_actual);
#endif
    	//Gal[p].MetalsHotGas.type1a += step_width_times_ICMSFR * (NormSNIaMetalEjecRate_actual + (ICM_total_metallicity * NormSNIaMassEjecRate_actual));
    	Gal[p].MetalsHotGas.type1a += max(0.0, step_width_times_ICMSFR * NormSNIaMetalEjecRate_actual);
    	Gal[p].MetalsHotGas.agb += max(0.0, step_width_times_ICMSFR * (NormAGBMetalEjecRate_actual + (ICM_total_metallicity * NormAGBMassEjecRate_actual)));

#ifdef INDIVIDUAL_ELEMENTS
#ifdef PORTINARI
    	Gal[p].HotGas_elements.H += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[0] + NormSNIaYieldRate_actual[0] + NormAGBYieldRate_actual[0]) + (Gal[p].sfh_ElementsICM[i].H*inverse_ICM_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.He += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[1] + NormSNIaYieldRate_actual[1] + NormAGBYieldRate_actual[1]) + (Gal[p].sfh_ElementsICM[i].He*inverse_ICM_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
#ifndef MAINELEMENTS
    	Gal[p].HotGas_elements.Cb += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsICM[i].Cb*inverse_ICM_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.N += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsICM[i].N*inverse_ICM_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.O += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsICM[i].O*inverse_ICM_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Ne += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[5] + NormSNIaYieldRate_actual[5] + NormAGBYieldRate_actual[5]) + (Gal[p].sfh_ElementsICM[i].Ne*inverse_ICM_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Mg += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[6] + NormSNIaYieldRate_actual[6] + NormAGBYieldRate_actual[6]) + (Gal[p].sfh_ElementsICM[i].Mg*inverse_ICM_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Si += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[7] + NormSNIaYieldRate_actual[7] + NormAGBYieldRate_actual[7]) + (Gal[p].sfh_ElementsICM[i].Si*inverse_ICM_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.S += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[8] + NormSNIaYieldRate_actual[8] + NormAGBYieldRate_actual[8]) + (Gal[p].sfh_ElementsICM[i].S*inverse_ICM_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Ca += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[9] + NormSNIaYieldRate_actual[9] + NormAGBYieldRate_actual[9]) + (Gal[p].sfh_ElementsICM[i].Ca*inverse_ICM_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Fe += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[10] + NormSNIaYieldRate_actual[10] + NormAGBYieldRate_actual[10]) + (Gal[p].sfh_ElementsICM[i].Fe*inverse_ICM_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
#else
    	Gal[p].HotGas_elements.O += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsICM[i].O*inverse_ICM_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Mg += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsICM[i].Mg*inverse_ICM_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Fe += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsICM[i].Fe*inverse_ICM_physical_units)*(NormSNIIMassEjecRate_actual + NormAGBMassEjecRate_actual)));
#endif //MAINELEMENTS
#endif //PORTINARI
#ifdef CHIEFFI
    	Gal[p].HotGas_elements.H += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[0] + NormSNIaYieldRate_actual[0] + NormAGBYieldRate_actual[0]) + (Gal[p].sfh_ElementsICM[i].H*inverse_ICM_physical_units)*(NormAGBMassEjecRate_actual)));  //ROB: No unsynth component required for SN-II ejecta, when using the Chieffi & Limongi 92007) yield tables/
    	Gal[p].HotGas_elements.He += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[1] + NormSNIaYieldRate_actual[1] + NormAGBYieldRate_actual[1]) + (Gal[p].sfh_ElementsICM[i].He*inverse_ICM_physical_units)*(NormAGBMassEjecRate_actual)));
#ifndef MAINELEMENTS
    	Gal[p].HotGas_elements.Cb += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsICM[i].Cb*inverse_ICM_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.N += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsICM[i].N*inverse_ICM_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.O += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsICM[i].O*inverse_ICM_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Ne += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[5] + NormSNIaYieldRate_actual[5] + NormAGBYieldRate_actual[5]) + (Gal[p].sfh_ElementsICM[i].Ne*inverse_ICM_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Mg += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[6] + NormSNIaYieldRate_actual[6] + NormAGBYieldRate_actual[6]) + (Gal[p].sfh_ElementsICM[i].Mg*inverse_ICM_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Si += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[7] + NormSNIaYieldRate_actual[7] + NormAGBYieldRate_actual[7]) + (Gal[p].sfh_ElementsICM[i].Si*inverse_ICM_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.S += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[8] + NormSNIaYieldRate_actual[8] + NormAGBYieldRate_actual[8]) + (Gal[p].sfh_ElementsICM[i].S*inverse_ICM_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Ca += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[9] + NormSNIaYieldRate_actual[9] + NormAGBYieldRate_actual[9]) + (Gal[p].sfh_ElementsICM[i].Ca*inverse_ICM_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Fe += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[10] + NormSNIaYieldRate_actual[10] + NormAGBYieldRate_actual[10]) + (Gal[p].sfh_ElementsICM[i].Fe*inverse_ICM_physical_units)*(NormAGBMassEjecRate_actual)));
#else
    	Gal[p].HotGas_elements.O += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[2] + NormSNIaYieldRate_actual[2] + NormAGBYieldRate_actual[2]) + (Gal[p].sfh_ElementsICM[i].O*inverse_ICM_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Mg += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[3] + NormSNIaYieldRate_actual[3] + NormAGBYieldRate_actual[3]) + (Gal[p].sfh_ElementsICM[i].Mg*inverse_ICM_physical_units)*(NormAGBMassEjecRate_actual)));
    	Gal[p].HotGas_elements.Fe += max(0.0, step_width_times_ICMSFR_physical_units * ((NormSNIIYieldRate_actual[4] + NormSNIaYieldRate_actual[4] + NormAGBYieldRate_actual[4]) + (Gal[p].sfh_ElementsICM[i].Fe*inverse_ICM_physical_units)*(NormAGBMassEjecRate_actual)));
#endif //MAINELEMENTS
#endif //CHIEFFI
#endif //INDIVIDUAL_ELEMENTS

    	//UPDATE ICL COMPONENTS:
    	Gal[p].ICM -= max(0.0, step_width_times_ICMSFR * NormMassEjecRateSumAllTypes);
    	Gal[p].MetalsICM.type2 -= max(0.0, step_width_times_ICMSFR * (Bulge_total_metallicity * NormSNIIMassEjecRate_actual));
    	Gal[p].MetalsICM.type1a -= max(0.0, step_width_times_ICMSFR * (Bulge_total_metallicity * NormSNIaMassEjecRate_actual));
    	Gal[p].MetalsICM.agb -= max(0.0, step_width_times_ICMSFR * (Bulge_total_metallicity * NormAGBMassEjecRate_actual));

#ifdef INDIVIDUAL_ELEMENTS
    	Gal[p].ICM_elements.H -= max(0.0, step_width_times_ICMSFR_physical_units * ((Gal[p].sfh_ElementsICM[i].H*inverse_ICM_physical_units)*NormMassEjecRateSumAllTypes));
    	Gal[p].ICM_elements.He -= max(0.0, step_width_times_ICMSFR_physical_units * ((Gal[p].sfh_ElementsICM[i].He*inverse_ICM_physical_units)*NormMassEjecRateSumAllTypes));
#ifndef MAINELEMENTS
    	Gal[p].ICM_elements.Cb -= max(0.0, step_width_times_ICMSFR_physical_units * ((Gal[p].sfh_ElementsICM[i].Cb*inverse_ICM_physical_units)*NormMassEjecRateSumAllTypes));
    	Gal[p].ICM_elements.N -= max(0.0, step_width_times_ICMSFR_physical_units * ((Gal[p].sfh_ElementsICM[i].N*inverse_ICM_physical_units)*NormMassEjecRateSumAllTypes));
#endif
    	Gal[p].ICM_elements.O -= max(0.0, step_width_times_ICMSFR_physical_units * ((Gal[p].sfh_ElementsICM[i].O*inverse_ICM_physical_units)*NormMassEjecRateSumAllTypes));
#ifndef MAINELEMENTS
    	Gal[p].ICM_elements.Ne -= max(0.0, step_width_times_ICMSFR_physical_units * ((Gal[p].sfh_ElementsICM[i].Ne*inverse_ICM_physical_units)*NormMassEjecRateSumAllTypes));
#endif
    	Gal[p].ICM_elements.Mg -= max(0.0, step_width_times_ICMSFR_physical_units * ((Gal[p].sfh_ElementsICM[i].Mg*inverse_ICM_physical_units)*NormMassEjecRateSumAllTypes));
#ifndef MAINELEMENTS
    	Gal[p].ICM_elements.Si -= max(0.0, step_width_times_ICMSFR_physical_units * ((Gal[p].sfh_ElementsICM[i].Si*inverse_ICM_physical_units)*NormMassEjecRateSumAllTypes));
    	Gal[p].ICM_elements.S -= max(0.0, step_width_times_ICMSFR_physical_units * ((Gal[p].sfh_ElementsICM[i].S*inverse_ICM_physical_units)*NormMassEjecRateSumAllTypes));
    	Gal[p].ICM_elements.Ca -= max(0.0, step_width_times_ICMSFR_physical_units * ((Gal[p].sfh_ElementsICM[i].Ca*inverse_ICM_physical_units)*NormMassEjecRateSumAllTypes));
#endif
    	Gal[p].ICM_elements.Fe -= max(0.0, step_width_times_ICMSFR_physical_units * ((Gal[p].sfh_ElementsICM[i].Fe*inverse_ICM_physical_units)*NormMassEjecRateSumAllTypes));
#endif //INDIVIDUAL_ELEMENTS

    } //if (Gal[p].sfh_ICM[i] > 0.0) //ICM

    } //for (i=0;i<=Gal[p].sfh_ibin;i++) //MAIN LOOP OVER SFH BINS

#ifdef DETAILED_DUST
    coldgas_add = elements_add(Gal[p].ColdGas_elements, coldgas_init, -1.);
    Gal[p].ColdGasClouds_elements = elements_add(Gal[p].ColdGasClouds_elements, coldgas_add, Gal[p].mu_gas);
    Gal[p].ColdGasDiff_elements = elements_add(Gal[p].ColdGasDiff_elements, coldgas_add, 1.0 - Gal[p].mu_gas);
#endif

#ifndef DETAILED_DUST
	// IF DETAILED_DUST is switched ON 
	// SN-Feedback must be called AFTER BOTH model_yields AND model_dustyields
	// and now appears inside main.c if this is the case! 

    //CALL SN-FEEDBACK RECIPE: Sending total mass returned to ColdGas to calculate FB energy:
    SN_feedback(p, centralgal, TotalMassReturnedToColdDiskGas);
#endif     

    //Update Mass-weighted ages:
    for(n=0;n<NOUT;n++)
    {
    	Gal[p].MassWeightAge[n] -= (AgeCorrectionDisk[n]+AgeCorrectionBulge[n]);
    }
    
#ifdef H2_AND_RINGS
    double TotalMassReturnedToColdDiskGasr[RNUM], TotalMassReturnedToHotGasr[RNUM];
    double Coldmetallicityr[RNUM], Hotmetallicity[RNUM];
    int ii;
    for(ii=0;ii<RNUM;ii++)
    {
    	TotalMassReturnedToColdDiskGasr[ii]= TotalMassReturnedToColdDiskGas/((float)RNUM);
    	TotalMassReturnedToHotGasr[ii]=TotalMassReturnedToHotGasr/((float)RNUM);
    	Coldmetallicityr[ii]=metals_total(Gal[p].MetalsColdGas)/Gal[p].ColdGas/((float)RNUM);
    	Hotmetallicity[ii]=metals_total(Gal[p].MetalsHotGas)/Gal[p].HotGas/((float)RNUM);
    }
#endif
}

int find_initial_metallicity(int p, int sfh_bin, int table_type, int component_type)
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
	else Z_disk = initMetals/Gal[p].sfh_DiskMass[sfh_bin];

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
		case 2: //SN-II metallicity table
			while (Zi_bin == -1)
			{
				if (SNIIMetallicities[i] < Z_disk)
				{
					i++;
					if (i == 5) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
				}
				else Zi_bin = i;
			}
			break;
		//case 3 //SNIa yields are NOT metallicity dependent
		case 4: //AGB metallicity table
			while (Zi_bin == -1)
			{
				if (AGBMetallicities[i] < Z_disk)
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
			case 2: //SN-II metallicity table
				while (Zi_bin == -1)
				{
					if (SNIIMetallicities[i] < Z_bulge)
					{
						i++;
						if (i == 5) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
					}
					else Zi_bin = i;
				}
				break;
			//case 3 //SNIa yields are NOT metallicity dependent
			case 4: //AGB metallicity table
				while (Zi_bin == -1)
				{
					if (AGBMetallicities[i] < Z_bulge)
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
				case 2: //SN-II metallicity table
					while (Zi_bin == -1)
					{
						if (SNIIMetallicities[i] < Z_ICM)
						{
							i++;
							if (i == 5) Zi_bin = i; //If galaxy's Z is higher than max Z from table, then just take max Z from table
						}
						else Zi_bin = i;
					}
					break;
				//case 3 //SNIa yields are NOT metallicity dependent
				case 4: //AGB metallicity table
					while (Zi_bin == -1)
					{
						if (AGBMetallicities[i] < Z_ICM)
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


#ifdef INSTANTANEOUS_RECYCLE //to recover results from instantaneous recycling approximation
 void reset_ejection_rates(int i, int sfh_ibin,
		 double *NormSNIIMassEjecRate_actual, double *NormSNIIMetalEjecRate_actual,
		 double *NormSNIaMassEjecRate_actual, double *NormAGBMassEjecRate_actual,
		 double *NormSNIaMetalEjecRate_actual, double *NormAGBMetalEjecRate_actual)
 {
    	if(i==sfh_ibin)
    	{
    		*NormSNIIMassEjecRate_actual = 0.43;
    		*NormSNIIMetalEjecRate_actual = 0.03;
    	}
    	else
    	{
    		*NormSNIIMassEjecRate_actual = 0.0;
    		*NormSNIIMetalEjecRate_actual = 0.0;
    	}
    	*NormSNIaMassEjecRate_actual = 0.0;
    	*NormAGBMassEjecRate_actual =  0.0;
    	*NormSNIaMetalEjecRate_actual = 0.0;
    	*NormAGBMetalEjecRate_actual =  0.0;
 }
#endif //INSTANTANEOUS_RECYCLE

