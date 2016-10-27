/*
 * dustyields_integrals.c
 *
 * Similar to yields_integrals which does....
 * Pre-calculates the normalised ejecta rates at every timestep, assuming 1 Msun populations.
 * Multiply by SFR from SFH bins (and interpolate between default metallicities) to obtain
 * true ejecta rates (done in recipe_yields.c).
 *
 * Created: Oct2016
 * Author: ScottClay
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

#ifdef DETAILED_DUST

void init_integrated_dust_yields()
{
	int ii, jj, kk, ll;

	for(ii=0;ii<STEPS*(LastDarkMatterSnapShot+1);ii++) {
		for(jj=0;jj<SFH_NBIN;jj++) {
			for(kk=0;kk<LIFETIME_Z_NUM;kk++) {
				for(ll=0;ll<AGB_DUST_TYPE_NUM;ll++) {
					NormAGBDustYieldRate[ii][jj][kk][ll]=0.0;									
}}}}}

void print_integrated_dust_yields()
{
	int ii, jj, kk, ll;

	for(ii=0;ii<STEPS*(LastDarkMatterSnapShot+1);ii++) {
		for(jj=0;jj<SFH_NBIN;jj++) {
			for(kk=0;kk<LIFETIME_Z_NUM;kk++) {
				for(ll=0;ll<AGB_DUST_TYPE_NUM;ll++) {
					//NormAGBDustYieldRate[ii][jj][kk][ll]=0.0;
					//if(NormAGBDustYieldRate[ii][jj][kk][ll] < 0.0) {
						printf("step = %d\t SFH bin = %d\t metallicity = %d\t dusttype = %d\n",ii,jj,kk,ll);
						printf("NormAGBDustYieldRate[ii][jj][kk][ll] = %f \n", NormAGBDustYieldRate[ii][jj][kk][ll]);
					//}
}}}/*printf("Normalise dust yields\t%d\n",ii);*/}}



void integrate_dust_yields()
{
	double previoustime, newtime, deltaT;
	int snap, step,i,mb;
	double timet;

	int Mi_lower, Mi_upper,Zi_AGB;
	int Mi_lower_AGB, Mi_upper_AGB, t_lower_lifetime, t_upper_lifetime;
	int width_in_timesteps, mbmax; //Number of current timesteps that fit in any given SFH bin, and the number of mini bins considered for any given SFH bin (max. = 30, for memory considerations)
	double dt, t_lower, t_upper, Mi_lower_actual, Mi_upper_actual;
	double AGBYields_lower_actual[AGB_DUST_TYPE_NUM], AGBYields_upper_actual[AGB_DUST_TYPE_NUM];
	
	
	for(snap=0;snap<(LastDarkMatterSnapShot+1)-1;snap++) //LOOP OVER SNAPSHOTS
	{
	    previoustime = NumToTime(snap); //Time to z=0 from start of current snapshot [in code units]
	    newtime = NumToTime(snap+1); //Time to z=0 from end of current snapshot [in code units]
	    deltaT = previoustime - newtime; //timestep

	    for(step=0;step<STEPS;step++) //LOOP OVER TIMESTEPS
	    {
	    	dt = deltaT/STEPS;  //Time-width of a timestep in current snapshot [in code units]
	    	timet = previoustime - (step + 0.5) * dt; //Time from middle of the current timestep to z=0 [in code units]
	        for (i=0;i<=SFH_ibin[snap][step];i++) //LOOP OVER SFH BINS
	        {
	   	  	//New method: sub-dividing SFH bins into SFH_NMINIBIN number of 'mini bins': Can later choose inside code which mini bin the characteristic SF time is in:
	        	width_in_timesteps = SFH_dt[snap][step][i]/dt; //Width of SFH bin in number of current timesteps [in code units] //NB: Typecasting a float to an integer here (width_in_timesteps is and integer).
	        	if(width_in_timesteps < 1) width_in_timesteps = 1;
	        	//if (width_in_timesteps > SFH_NMINIBIN) {mbmax = SFH_NMINIBIN;} else {mbmax = width_in_timesteps;}
	        	mbmax = width_in_timesteps;
	        	for (mb=1;mb<=mbmax;mb++) //LOOP OVER MINI BINS (New method)
	        	{
	        		//From lower/upper edge of mini-bin to middle of current timestep:
	        		t_lower = (SFH_t[snap][step][i] + (SFH_dt[snap][step][i]) - (mb*(SFH_dt[snap][step][i]/mbmax)) - timet) * UnitTime_in_years/Hubble_h; //IN YEARS //Time from low-z (lower) edge of SFH mini-bin j to middle of current timestep
	        		t_upper = (SFH_t[snap][step][i] + (SFH_dt[snap][step][i]) - (mb*(SFH_dt[snap][step][i]/mbmax)) + (SFH_dt[snap][step][i]/mbmax) - timet) * UnitTime_in_years/Hubble_h; //IN YEARS //Time from high-z (upper) edge of SFH mini-bin j to middle of current timestep

	        	  int Zi;
	        	  for (Zi=0;Zi<LIFETIME_Z_NUM;Zi++) //LOOP OVER POSSIBLE INITIAL METALLICITIES
	        	  {
	        	  	//Mi_lower_upper goes to get the bin number, so an integer, from the lifetime array
	        	  	//Mi_lower_upper_actual is the actual mass in Msol
	        		
	        		Mi_lower = find_initial_mass_dust(t_upper, Zi); //Mass bin (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep, from SFH bin i.
	        		Mi_upper = find_initial_mass_dust(t_lower, Zi); //Mass bin (lifetime arrays) corresp. to highest mass of star to 'die' in current timestep, from SFH bin i.
	        
	        		Mi_lower_actual = lifetimeMasses[Mi_lower] + ((lifetimeMasses[Mi_lower+1]-lifetimeMasses[Mi_lower]) * ((t_upper-lifetimes[Zi][Mi_lower])/(lifetimes[Zi][Mi_lower+1]-lifetimes[Zi][Mi_lower]))); //IN MSUN  //Lowest mass of star to 'die' in current timestep from SFH bin i.
	        		Mi_upper_actual = lifetimeMasses[Mi_upper] + ((lifetimeMasses[Mi_upper+1]-lifetimeMasses[Mi_upper]) * ((t_lower-lifetimes[Zi][Mi_upper])/(lifetimes[Zi][Mi_upper+1]-lifetimes[Zi][Mi_upper]))); //IN MSUN  //Highest mass of star to 'die' in current timestep from SFH bin i.


	        		if (Mi_upper_actual <= 0.0 || Mi_upper_actual > 120.0 || Mi_upper == LIFETIME_MASS_NUM-1)
	        			Mi_upper_actual = SNII_MAX_MASS; //Mi_upper_actual can be	< 0.0 for the highest Zi, as lifetimes[5][149+1] = 0.0 (i.e. doesn't exist!)

	        		if (Mi_lower_actual <= 0.85) Mi_lower_actual = AGB_MIN_MASS; //No stars below 0.85 Msun contribute to chemical enrichment.
	        		//if (Mi_upper == LIFETIME_MASS_NUM-1) Mi_upper_actual = SNII_MAX_MASS; //No stars above 120 Msun assumed to exist.
	        		if (lifetimeMasses[Mi_upper] >= SNII_MAX_MASS) Mi_upper_actual = SNII_MAX_MASS; //No stars of mass above max. SN-II progenitor assumed to exist.

	        		//if (Mi_upper_actual <= 0.0) {Mi_upper_actual = Mi_lower_actual;} //ROB: Just a condition added when artificially changing lifetimes. (08-01-13) //printf("USED!!\n");
					
					if (Mi_upper_actual <= 0.85) Mi_upper_actual = AGB_MIN_MASS;
					if (Mi_lower_actual >= 120.0) Mi_lower_actual = SNII_MAX_MASS;
					
					/*
					if ( Mi_lower_actual > Mi_upper_actual) {
						//printf("%f \t %f \n",Mi_lower_actual,Mi_upper_actual);
						printf("%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n",Mi_lower_actual,Mi_upper_actual,Mi_upper,Mi_upper+1,Zi,lifetimeMasses[Mi_upper],lifetimeMasses[Mi_upper+1],t_lower,lifetimes[Zi][Mi_upper]);
						}
					*/




	        		if (t_upper >= lifetimes[Zi][LIFETIME_MASS_NUM-1]) //If the longest time from SFH bin i to current timestep is shorter than the shortest possible lifetime, there's no enrichment, so skip calculations.
	        		{
	        			//*****************************************
	        			//AGB Winds (Disc and Bulge):
	        			//*****************************************
	        			//Zi = lifetime metallicity bun number, Zi_AGB goes to find the metallicity from the dust table that matches this lifetime metallicity. 
	        	    	Zi_AGB = find_initial_metallicity_comp_dust(Zi, i, 4); //finds the metallicity in dust table that corresponds to metallicity bin in the lifetimes table 
						
						//Mi_lower_upper = integer mass bin from lifetime arrays
	        			Mi_lower_AGB = max_Mi_lower(Mi_lower,4); //Mass bin (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep, from SFH bin i.
	        			Mi_upper_AGB = min_Mi_upper(Mi_upper,4); //Mass bin (lifetime arrays) corresp. to highest mass of star to 'die' in current timestep, from SFH bin i.

	        	    	if (Mi_lower_AGB <= Mi_upper_AGB)
	        	    	{
	        	    		//Mi_lower_upper_AGB = INTEGER = dust agb mass bin number 
	        	    		Mi_lower_AGB = find_agb_mass_bin_dust(lifetimeMasses[Mi_lower_AGB]);
	        	    		Mi_upper_AGB = find_agb_mass_bin_dust(lifetimeMasses[Mi_upper_AGB]);

	        	    		find_actual_ejecta_limits_dust(4, Mi_lower_actual, Mi_upper_actual, Mi_lower_AGB, Mi_upper_AGB, Zi_AGB,
	        	    				AGBYields_lower_actual, AGBYields_upper_actual);

							//Actual integration
	        	    		int j;
	        	    		for (j=Mi_lower_AGB;j<Mi_upper_AGB;j++)
	        	    		//for (j=Mi_lower_AGB;j<=Mi_upper_AGB;j++)

	        	    		{       
	        	    			int k;
								for (k=0;k<AGB_DUST_TYPE_NUM;k++) { 
									NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k] += (AGBDustMasses[j+1]-AGBDustMasses[j]) * ((AGBDustCreated[Zi_AGB][j][k] + AGBDustCreated[Zi_AGB][j+1][k])/2.0);
        	    					if (NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k] < 0.0) {
        	    						printf("0\t%d\t\t%f \t %f \t %f \t %f \t %f \n",j,NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k], AGBDustMasses[j+1], AGBDustMasses[j], AGBDustCreated[Zi_AGB][j][k], AGBDustCreated[Zi_AGB][j+1][k]);
        	    						}
        	    					}
	        	    		
	        	    		 	    				
        	    			/*	
        	    				if (j != Mi_lower_AGB && j != Mi_upper_AGB)
        	    				{
									int k;
									for (k=0;k<AGB_DUST_TYPE_NUM;k++) { 
										NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k] += (AGBDustMasses[j+1]-AGBDustMasses[j]) * ((AGBDustCreated[Zi_AGB][j][k] + AGBDustCreated[Zi_AGB][j+1][k])/2.0);
        	    						if (NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k] < 0.0) {
        	    							printf("1\t%f \t %f \t %f \t %f \t %f \n",NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k], AGBDustMasses[j+1], AGBDustMasses[j], AGBDustCreated[Zi_AGB][j][k], AGBDustCreated[Zi_AGB][j+1][k]);
        	    							}
        	    					}
        	    				} 
        	    				else {
									int k;
									for (k=0;k<AGB_DUST_TYPE_NUM;k++) { 
        	    					NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k] = 0.0;
        	
        	    				}*/
        	    				/*
        	    				else if (j == Mi_lower_AGB && j == Mi_upper_AGB && Mi_lower_actual <= AGB_MAX_MASS)
        	    				{
									int k;
									for (k=0;k<AGB_DUST_TYPE_NUM;k++) { 
										NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k] += (Mi_upper_actual-Mi_lower_actual) * ((AGBYields_lower_actual[k] + AGBYields_upper_actual[k])/2.0);
        	    						if (NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k] < 0.0) {
        	    							printf("2\t%f \t %f \t %f \t %f \t %f \n",NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k], Mi_upper_actual, Mi_lower_actual, AGBYields_lower_actual[k], AGBYields_upper_actual[k]);
        	    							}
        	    					
        	    					}
        	    				}
        	    				else if (j == Mi_lower_AGB && j != Mi_upper_AGB)
        	    				{
									int k;
									for (k=0;k<AGB_DUST_TYPE_NUM;k++) { 
										NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k] += (AGBDustMasses[j+1]-Mi_lower_actual) * ((AGBYields_lower_actual[k] + AGBDustCreated[Zi_AGB][j+1][k])/2.0);
        	    						if (NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k] < 0.0) {
        	    							printf("3\t%f \t %f \t %f \t %f \t %f \n",NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k], AGBDustMasses[j+1],Mi_lower_actual, AGBYields_lower_actual[k], AGBDustCreated[Zi_AGB][j+1][k]);
        	    							}
        	    					
        	    					}
        	    				}
        	    				else if (j == Mi_upper_AGB && j != Mi_lower_AGB)
        	    				{
        	    					int k;
									for (k=0;k<AGB_DUST_TYPE_NUM;k++) { 
										NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k] += (Mi_upper_actual-AGBDustMasses[j]) * ((AGBDustCreated[Zi_AGB][j][k] + AGBYields_upper_actual[k])/2.0);
        	    						if (NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k] < 0.0) {
        	    						
        	    							printf("4\t%f \t %f \t %f \t %f \t %f \n",NormAGBDustYieldRate[(STEPS*snap)+step][i][Zi][k], Mi_upper_actual, AGBDustMasses[j], AGBDustCreated[Zi_AGB][j][k], AGBYields_upper_actual[k]);
        	    							}
									
									
        	    					}
         	    				}*/
        	    			}	//for (j=Mi_lower_AGB;j<=Mi_upper_AGB;j++)
	        	    	} //if (Mi_lower_AGB <= Mi_upper_AGB)
	        	    	
	        	   } //if (t_upper >= lifetimes[Zi][Mi_lower+1])

	        	} //for (Zi=0;Zi<LIFETIME_Z_NUM;Zi++)
	          } //for (j=0;j<=width_in_timesteps;j++) //MINI_BINS	
	        } //for (i=0;i<=SFH_ibin_structure[(SFH_NBIN*snap)+step];i++)
	        	
	    } //for(step=0;step<STEPS;step++)
	} //for(snap=0;snap<(LastDarkMatterSnapShot+1);snap++)
	printf("Dust Yield integrals calculated.\n");
}

int find_initial_metallicity_comp_dust(int Zi, int sfh_bin, int table_type)
{
	int i, Zi_bin;
	double Z_in;

	Zi_bin = -1; //This is the bin we want. The closest in dust tables to lifetime metallicity
	i = 0;
	Z_in = lifetimeMetallicities[Zi]; //Actual lifetime metallicity 

	switch (table_type)
	{
		case 4: //AGB dust metallicity table
			while (Zi_bin == -1)
			{
				if (AGBDustMetallicities[i] < Z_in)
				{
					i++;
					if (i == AGB_DUST_METAL_NUM) Zi_bin = i-1;
				}
				else Zi_bin = i;
			}
			break;
	}
	return Zi_bin;
}

int find_initial_mass_dust(double lifetime, int Zi_bin)
{
    if (lifetime == 0.0) return LIFETIME_MASS_NUM-1; //If the bin 'touches now', then return max mass (ie: star of shortest lifetime) ie: bin for 120Msun
    else if (lifetime > lifetimes[Zi_bin][0]) return 0; //If true lifetime is longer than max lifetime in table (shouldn't be), then return element number 0
    else
    {
	int Mi_bin;

	Mi_bin = -1;
	int i;
	i = 0;
	while (Mi_bin == -1)
	{
		if (lifetimes[Zi_bin][i] > lifetime)
		{
			i++;
			if (i == LIFETIME_MASS_NUM) Mi_bin = i; //If lifetime is shorter than min lifetime from table, then just return max mass (120 Msun)
		}
		else Mi_bin = i;
	}
	return Mi_bin-1; //This returns element number i for lifetimeMasses[Zi][i] array BELOW the true initial mass corresponding to t_lower or t_upper.
    }
}

int max_Mi_lower_dust(int Mi_lower, int channel_type)
{
	switch (channel_type)
		{
			case 4: //AGB mass limits
					if (lifetimeMasses[Mi_lower] > AGB_MIN_MASS) return Mi_lower;
					else
					{
						int i;
						i = 0;
						do { i++; }
						while (lifetimeMasses[i] < AGB_MIN_MASS);
						return i;
					}
					break;
			default: printf("Wrong ejection mode chosen in max_Mi_lower: Use either 2 (SNe-II), 3 (SNe-Ia) or 4 (AGB winds)"); exit(1);
		}
}

int min_Mi_upper_dust(int Mi_upper, int channel_type)
{
	switch (channel_type)
		{
			case 4: //AGB mass limits
					if (lifetimeMasses[Mi_upper] < AGB_MAX_MASS) return Mi_upper;
					else
					{
						int i;
						i = LIFETIME_MASS_NUM-1;
						do { i--; }
						while (lifetimeMasses[i] > AGB_MAX_MASS);
						return i;
					}
					break;
			default: printf("Wrong ejection mode chosen in min_Mi_upper: Use either 2 (SNe-II), 3 (SNe-Ia) or 4 (AGB winds)"); exit(1);
		}
}


int find_agb_mass_bin_dust(double masslimit)
{
	if (masslimit == AGB_MAX_MASS) return AGB_DUST_MASS_NUM-1; //If AGB_MAX_MASS == mass limit = 7Msol, this is max for agb dust tables too. So this works. If you change AGB_Max_Mass - this will need changing
	else
	{
	int Mi_bin;

	Mi_bin = -1;
	int i;
	i = 0;
	while (Mi_bin == -1)
	{
		if (AGBDustMasses[i] < masslimit)
		{
			i++;
			if (i == AGB_DUST_MASS_NUM) Mi_bin = i-1; //If mass is greater than max mass for AGB winds (shouldn't be), then just return max mass (5.0 Msun)
		}
		else Mi_bin = i;
	}
	return Mi_bin;
	}
}

void find_actual_ejecta_limits_dust(int channel_type, double Mi_lower_actual, double Mi_upper_actual, int Mi_lower, int Mi_upper, int Zi,
double* Yields_lower_actual, double* Yields_upper_actual)
{
	switch (channel_type)
	{
	case 4: //AGB
		if (Mi_lower == 0)
		{
		    int k;
		    for (k=0;k<AGB_DUST_TYPE_NUM;k++)
		    {
		    	Yields_lower_actual[k] = AGBDustCreated[Zi][0][k];
		    }
		}
		else
		{
		    int k;
		    for (k=0;k<AGB_DUST_TYPE_NUM;k++)
		    {
				Yields_lower_actual[k] = AGBDustCreated[Zi][Mi_lower][k] + ((AGBDustCreated[Zi][Mi_lower+1][k]-AGBDustCreated[Zi][Mi_lower][k]) * ((Mi_lower_actual-AGBDustMasses[Mi_lower])/(AGBDustMasses[Mi_lower+1]-AGBDustMasses[Mi_lower])));
		    }
		}

		if (Mi_upper == AGB_MASS_NUM-1)
		{
		    int k;
		    for (k=0;k<AGB_DUST_TYPE_NUM;k++)
		    {
		    	Yields_upper_actual[k] = AGBDustCreated[Zi][AGB_DUST_MASS_NUM-1][k];
		    }
		}
		else
		{
		    int k;
		    for (k=0;k<AGB_DUST_TYPE_NUM;k++)
		    {
		    	Yields_upper_actual[k] = AGBDustCreated[Zi][Mi_upper][k] + ((AGBDustCreated[Zi][Mi_upper+1][k]-AGBDustCreated[Zi][Mi_upper][k]) * ((Mi_upper_actual-AGBDustMasses[Mi_upper])/(AGBDustMasses[Mi_upper+1]-AGBDustMasses[Mi_upper])));
		     }
		}
		break;
	}
}

#endif //DETAILED_EMRICHMENT

