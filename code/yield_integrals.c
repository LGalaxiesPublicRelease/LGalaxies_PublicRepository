/*
 * yield_integrals.c
 *
 * Pre-calculates the normalised ejecta rates at every timestep, assuming 1 Msun populations.
 * Multiply by SFR from SFH bins (and interpolate between default metallicities) to obtain
 * true ejecta rates (done in recipe_yields.c).
 *
 *  Created on: 10.05.2012
 *      Author: robyates
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

#ifdef DETAILED_METALS_AND_MASS_RETURN

void init_integrated_yields()
{
	int ii, jj, kk, ll;

	for(ii=0;ii<STEPS*(LastDarkMatterSnapShot+1);ii++)
		for(jj=0;jj<SFH_NBIN;jj++)
			for(kk=0;kk<LIFETIME_Z_NUM;kk++)
			{
				NormSNIIMassEjecRate[ii][jj][kk]=0.;
				NormSNIIMetalEjecRate[ii][jj][kk]=0.;
#ifdef INDIVIDUAL_ELEMENTS
				for(ll=0;ll<NUM_ELEMENTS;ll++)
					NormSNIIYieldRate[ii][jj][kk][ll]=0.;
#endif
				NormAGBMassEjecRate[ii][jj][kk]=0.;
				NormAGBMetalEjecRate[ii][jj][kk]=0.;
#ifdef INDIVIDUAL_ELEMENTS
				for(ll=0;ll<NUM_ELEMENTS;ll++)
					NormAGBYieldRate[ii][jj][kk][ll]=0.;
#endif
				NormSNIaMassEjecRate[ii][jj][kk]=0.;
				NormSNIaMetalEjecRate[ii][jj][kk]=0.;
#ifdef INDIVIDUAL_ELEMENTS
				for(ll=0;ll<NUM_ELEMENTS;ll++)
					NormSNIaYieldRate[ii][jj][kk][ll]=0.;
#endif
			}
}


void integrate_yields()
{
	double previoustime, newtime, deltaT;
	int snap, step,i,mb;
	double timet;

	int Mi_lower, Mi_upper, Mi_lower_SNII, Mi_upper_SNII, Zi_SNII, Zi_AGB;
	int Mi_lower_AGB, Mi_upper_AGB, t_lower_lifetime, t_upper_lifetime;
	int width_in_timesteps, mbmax; //Number of current timesteps that fit in any given SFH bin, and the number of mini bins considered for any given SFH bin (max. = 30, for memory considerations)
	double dt, t_lower, t_upper, DTD_lower, DTD_upper, Mi_lower_actual, Mi_upper_actual;
	double SNIIEjectedMasses_lower_actual, SNIIEjectedMasses_upper_actual, SNIITotalMetals_lower_actual, SNIITotalMetals_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
	double SNIIYields_lower_actual[NUM_ELEMENTS], SNIIYields_upper_actual[NUM_ELEMENTS];
#endif
	double AGBEjectedMasses_lower_actual, AGBEjectedMasses_upper_actual, AGBTotalMetals_lower_actual, AGBTotalMetals_upper_actual;
#ifdef INDIVIDUAL_ELEMENTS
	double AGBYields_lower_actual[NUM_ELEMENTS], AGBYields_upper_actual[NUM_ELEMENTS];
#endif

	FRACCOUNTA = 0;

	int L1a,L2a,L3a,L4a,L5a,L1b,L2b,L3b,L4b;
	L1a=0;L2a=0;L3a=0;L4a=0;L5a=0;L1b=0;L2b=0;L3b=0;L4b=0;

	int counta;
	TheSFH[0] = 1.0/(tau_dt[0]*UnitTime_in_years/Hubble_h);
	for(counta=1;counta<SFH_NBIN;counta++)
	{
		TheSFH[counta] = 1.0/(tau_dt[counta]*UnitTime_in_years/Hubble_h);
	}
	//TheSFH[0] = 0.0; TheSFH[1] = 3.8; TheSFH[2] = 5.0; TheSFH[3] = 5.5; TheSFH[4] = 5.2; TheSFH[5] = 5.0; TheSFH[6] = 4.8; TheSFH[7] = 4.2; TheSFH[8] = 4.0; TheSFH[9] = 3.8; TheSFH[10] = 3.2; TheSFH[11] = 3.0; TheSFH[12] = 2.8; TheSFH[13] = 2.5; TheSFH[14] = 2.2; TheSFH[15] = 2.2; TheSFH[16] = 2.0; TheSFH[17] = 1.8; TheSFH[18] = 1.8; TheSFH[19] = 1.5;

	//Set KALPHA //Integral of the IMF (by number) from 0.1 - MAX Msun:
	double KALPHA;
	double F316;
	if (SNII_MAX_MASS == 120.0)
	{
		KALPHA = 1.4765;
		F316 = 0.0385147;
	}
	else if (SNII_MAX_MASS == 100.0)
	{
		KALPHA = 1.49272;
		F316 = 0.0389188;
	}
	else if (SNII_MAX_MASS == 70.0)
	{
		KALPHA = 1.52679;
		F316 = 0.0398185;
	}
	else if (SNII_MAX_MASS == 60.0)
	{
		KALPHA = 1.54319;
		F316 = 0.0402507;
	}
	else if (SNII_MAX_MASS == 50.0)
	{
		KALPHA = 1.56408;
		F316 = 0.0408018;
	}
	else if (SNII_MAX_MASS == 40.0)
	{
		KALPHA = 1.59203;
		F316 = 0.0415416;
	}
	else if (SNII_MAX_MASS == 30.0)
	{
		KALPHA = 1.63252;
		F316 = 0.0426178;
	}
	else
	{
		KALPHA = 1.49272;
		F316 = 0.0389188;
		printf("****************\n");
		printf("In yield_integrals.c:\n");
		printf("\nSNII_MAX_MASS neither 30, 40, 50, 60, 70, 100, nor 120 Msun\n");
		printf("KALPHA set to 1.49272\n");
		printf("F316 set to 0.0389188\n");
		printf("(These are the values for SNII_MAX_MASS = 100.0)\n");
		printf("****************\n\n");
	}

	for(snap=0;snap<(LastDarkMatterSnapShot+1)-1;snap++) //LOOP OVER SNAPSHOTS
	{
	    previoustime = NumToTime(snap); //Time to z=0 from start of current snapshot [in code units]
	    newtime = NumToTime(snap+1); //Time to z=0 from end of current snapshot [in code units]
	    deltaT = previoustime - newtime;

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
	        		Mi_lower = find_initial_mass(t_upper, Zi); //Mass bin (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep, from SFH bin i.
	        		Mi_upper = find_initial_mass(t_lower, Zi); //Mass bin (lifetime arrays) corresp. to highest mass of star to 'die' in current timestep, from SFH bin i.
	        
	        		Mi_lower_actual = lifetimeMasses[Mi_lower] + ((lifetimeMasses[Mi_lower+1]-lifetimeMasses[Mi_lower]) * ((t_upper-lifetimes[Zi][Mi_lower])/(lifetimes[Zi][Mi_lower+1]-lifetimes[Zi][Mi_lower]))); //IN MSUN  //Lowest mass of star to 'die' in current timestep from SFH bin i.
	        		Mi_upper_actual = lifetimeMasses[Mi_upper] + ((lifetimeMasses[Mi_upper+1]-lifetimeMasses[Mi_upper]) * ((t_lower-lifetimes[Zi][Mi_upper])/(lifetimes[Zi][Mi_upper+1]-lifetimes[Zi][Mi_upper]))); //IN MSUN  //Highest mass of star to 'die' in current timestep from SFH bin i.

	        		if (Mi_upper_actual <= 0.0 || Mi_upper_actual > 120.0 || Mi_upper == LIFETIME_MASS_NUM-1)
	        			Mi_upper_actual = SNII_MAX_MASS; //Mi_upper_actual can be	< 0.0 for the highest Zi, as lifetimes[5][149+1] = 0.0 (i.e. doesn't exist!)

	        		if (Mi_lower_actual <= 0.85) Mi_lower_actual = AGB_MIN_MASS; //No stars below 0.85 Msun contribute to chemical enrichment.
	        		//if (Mi_upper == LIFETIME_MASS_NUM-1) Mi_upper_actual = SNII_MAX_MASS; //No stars above 120 Msun assumed to exist.
	        		if (lifetimeMasses[Mi_upper] >= SNII_MAX_MASS) Mi_upper_actual = SNII_MAX_MASS; //No stars of mass above max. SN-II progenitor assumed to exist.

	        		//if (Mi_upper_actual <= 0.0) {Mi_upper_actual = Mi_lower_actual;} //ROB: Just a condition added when artificially changing lifetimes. (08-01-13) //printf("USED!!\n");

	        		if (t_upper >= lifetimes[Zi][LIFETIME_MASS_NUM-1]) //If the longest time from SFH bin i to current timestep is shorter than the shortest possible lifetime, there's no enrichment, so skip calculations.
	        		{
	        			//*****************************************
	        			//SNe-II (Disc and Bulge):
	        			//*****************************************
	        			Zi_SNII = find_initial_metallicity_comp(Zi, i, 2); //Metallicity bin (SNe-II arrays) corresp. to metallicity Zi.

	        	    	//Check if mass range is within range for SN-II progenitor stars:
	        			Mi_lower_SNII = max_Mi_lower(Mi_lower,2); //ROB: Should we send Mi_lower_actual and Mi_upper_actual to max_Mi_lower and min_Mi_upper? (24-07-12)
	        			Mi_upper_SNII = min_Mi_upper(Mi_upper,2);

	        			if (Mi_lower_SNII <= Mi_upper_SNII)
	        			{
	        				Mi_lower_SNII = find_SNII_mass_bin(lifetimeMasses[Mi_lower_SNII]); //Mass bin (SNe-II arrays) corresp. to lowest mass of star to 'die' in current timestep from SFH bin i.
	        				Mi_upper_SNII = find_SNII_mass_bin(lifetimeMasses[Mi_upper_SNII]); //Mass bin (SNe-II arrays) corresp. to highest mass of star to 'die' in current timestep from SFH bin i.

	        				//Find true yields at the true upper and lower masses, given by 'Mi_upper_actual' and 'Mi_lower_actual':
#ifdef INDIVIDUAL_ELEMENTS
	        				find_actual_ejecta_limits(2, Mi_lower_actual, Mi_upper_actual, Mi_lower_SNII, Mi_upper_SNII, Zi_SNII,
	        						&SNIIEjectedMasses_lower_actual, &SNIIEjectedMasses_upper_actual, &SNIITotalMetals_lower_actual, &SNIITotalMetals_upper_actual,
	        						SNIIYields_lower_actual, SNIIYields_upper_actual);
#else
	        				find_actual_ejecta_limits(2, Mi_lower_actual, Mi_upper_actual, Mi_lower_SNII, Mi_upper_SNII, Zi_SNII,
	        						&SNIIEjectedMasses_lower_actual, &SNIIEjectedMasses_upper_actual, &SNIITotalMetals_lower_actual, &SNIITotalMetals_upper_actual);
#endif

	        	    		//NUMERICALLY INTEGRATE OVER THE MASS RANGE APPLICABLE FOR SNe-II:
	        	    		int j;
	        	    		for (j=Mi_lower_SNII;j<=Mi_upper_SNII;j++)
	        	    		{
	        	    			if (SNIIMasses[j] <= SNIA_MAX_MASS) //For mass range where both SN-II and SN-Ia progenitors are possible
	        	    			{
	        	    				if (j != Mi_lower_SNII && j != Mi_upper_SNII) //If mass bin j is NEITHER the lowest NOR highest mass bin to be integrated over.
	        	    				{
	        	    					L1a++;
	        	    					//SNIIRate[(STEPS*snap)+step][Zi] += (1.0-A_FACTOR) * (SNIIMasses[j+1]-SNIIMasses[j]) * ((Chabrier_IMF(SNIIMasses[j])*TheSFH[i]) + (Chabrier_IMF(SNIIMasses[j+1])*TheSFH[i]))/2.0;
	        	    					NormSNIIMassEjecRate[(STEPS*snap)+step][i][Zi] += (1.0-A_FACTOR) * ((SNIIMasses[j+1]-SNIIMasses[j]) * ((SNIIEjectedMasses[Zi_SNII][j] + SNIIEjectedMasses[Zi_SNII][j+1])/2.0)); //IN [MSun/yr] //The rate of ejection of mass from a 1 Msun stellar population.
	        	    					NormSNIIMetalEjecRate[(STEPS*snap)+step][i][Zi] += (1.0-A_FACTOR) * (SNIIMasses[j+1]-SNIIMasses[j]) * ((SNIITotalMetals[Zi_SNII][j] + SNIITotalMetals[Zi_SNII][j+1])/2.0); //IN [MSun/yr] //The rate of ejection of 'newly synthesised' metals from a 1 Msun stellar population.
                      #ifdef INDIVIDUAL_ELEMENTS
	        	    					int k;
	        	    					for (k=0;k<NUM_ELEMENTS;k++)
	        	    					{
	        	    						int kk;
											#ifndef MAINELEMENTS
	        	    						kk=k;
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (1.0-A_FACTOR) * (SNIIMasses[j+1]-SNIIMasses[j]) * ((SNIIYields[Zi_SNII][kk][j] + SNIIYields[Zi_SNII][kk][j+1])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
	        	    						#else
	        	    						switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (1.0-A_FACTOR) * (SNIIMasses[j+1]-SNIIMasses[j]) * ((SNIIYields[Zi_SNII][kk][j] + SNIIYields[Zi_SNII][kk][j+1])/2.0); //IN [Msun/yr] //The rate of ejection of individual elements from a 1 Msun stellar population.
											#endif

	        	    					}
	        	    			#endif //INDIVIDUAL_ELEMENTS

	        	    				}
	        	    				else if (j == Mi_lower_SNII && j == Mi_upper_SNII && Mi_lower_actual >= SNII_MIN_MASS) //If mass bin j is BOTH the lowest AND highest mass bin to be integrated over.
	        	    				{
	        	    					L2a++;
	        	    					//SNIIRate[(STEPS*snap)+step][Zi] += (1.0-A_FACTOR) * (Mi_upper_actual-Mi_lower_actual) * ((Chabrier_IMF(Mi_lower_actual)*TheSFH[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH[i]))/2.0;
	        	    					NormSNIIMassEjecRate[(STEPS*snap)+step][i][Zi] += (1.0-A_FACTOR) * ((Mi_upper_actual-Mi_lower_actual) * ((SNIIEjectedMasses_lower_actual + SNIIEjectedMasses_upper_actual)/2.0));
	        	    					NormSNIIMetalEjecRate[(STEPS*snap)+step][i][Zi] += (1.0-A_FACTOR) * ((Mi_upper_actual-Mi_lower_actual) * ((SNIITotalMetals_lower_actual + SNIITotalMetals_upper_actual)/2.0));
	        	    					float Sum_NormSNIIYieldRate; Sum_NormSNIIYieldRate = 0.0;
											#ifdef INDIVIDUAL_ELEMENTS
	        	    					int k;
	        	    					for (k=0;k<NUM_ELEMENTS;k++)
	        	    					{
	        	    						int kk;
											#ifndef MAINELEMENTS
	        	    						kk=k;
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (1.0-A_FACTOR) * (Mi_upper_actual-Mi_lower_actual) * ((SNIIYields_lower_actual[k] + SNIIYields_upper_actual[k])/2.0); //NB: 'SNIIYields_lower_actual' and 'SNIIYields_upper_actual' are arrays of k elements, so index with k, not kk, here.
	        	    						if (k != 0 && k != 1) {Sum_NormSNIIYieldRate += NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k];}
	        	    						#else
	        	    						switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (1.0-A_FACTOR) * (Mi_upper_actual-Mi_lower_actual) * ((SNIIYields_lower_actual[k] + SNIIYields_upper_actual[k])/2.0); //NB: 'SNIIYields_lower_actual' and 'SNIIYields_upper_actual' are arrays of k elements, so index with k, not kk, here.
											#endif
	        	    					}
											#endif //INDIVIDUAL_ELEMENTS
	        	    				}
	        	    				else if (j == Mi_lower_SNII && j == Mi_upper_SNII && Mi_lower_actual < SNII_MIN_MASS) //If mass bin j is BOTH the lowest AND highest mass bin to be integrated over, AND 'Mi_lower_actual' is below min. mass for SNe-II. (Still only counts stars from SNII_MIN_MASS to Mi_upper_actual).
	        	    				{
	        	    					L3a++;
	        	    					//printf("j = %i, Mi_lower_SNII = %i, Mi_upper_SNII = %i, Mi_lower_actual = %f, SNII_MIN_MASS = %f, Mi_upper_actual = %f\n", j, Mi_lower_SNII, Mi_upper_SNII, Mi_lower_actual, SNII_MIN_MASS, Mi_upper_actual);
	        	    					//SNIIRate[(STEPS*snap)+step][Zi] += (1.0-A_FACTOR) * (Mi_upper_actual-SNII_MIN_MASS) * ((Chabrier_IMF(SNII_MIN_MASS)*TheSFH[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH[i]))/2.0;
	        	    					NormSNIIMassEjecRate[(STEPS*snap)+step][i][Zi] += (1.0-A_FACTOR) * ((Mi_upper_actual-SNII_MIN_MASS) * ((SNIIEjectedMasses_lower_actual + SNIIEjectedMasses_upper_actual)/2.0));
	        	    					NormSNIIMetalEjecRate[(STEPS*snap)+step][i][Zi] += (1.0-A_FACTOR) * ((Mi_upper_actual-SNII_MIN_MASS) * ((SNIITotalMetals_lower_actual + SNIITotalMetals_upper_actual)/2.0));
	        	    					float Sum_NormSNIIYieldRate; Sum_NormSNIIYieldRate = 0.0;
											#ifdef INDIVIDUAL_ELEMENTS
	        	    					int k;
	        	    					for (k=0;k<NUM_ELEMENTS;k++)
	        	    					{
	        	    						int kk;
											#ifndef MAINELEMENTS
	        	    						kk=k;
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (1.0-A_FACTOR) * (Mi_upper_actual-SNII_MIN_MASS) * ((SNIIYields_lower_actual[k] + SNIIYields_upper_actual[k])/2.0);
	        	    						if (k != 0 && k != 1) {Sum_NormSNIIYieldRate += NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k];}
	        	    						#else
	        	    						switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (1.0-A_FACTOR) * (Mi_upper_actual-SNII_MIN_MASS) * ((SNIIYields_lower_actual[k] + SNIIYields_upper_actual[k])/2.0);
											#endif
	        	    					}
											#endif //INDIVIDUAL_ELEMENTS
	        	    				}
	        	    				else if (j == Mi_lower_SNII && j != Mi_upper_SNII) //If mass bin j IS the lowest but IS NOT the highest mass bin to be integrated over.
	        	    				{
	        	    					L4a++;
	        	    					//SNIIRate[(STEPS*snap)+step][Zi] += (1.0-A_FACTOR) * (SNIIMasses[j+1]-Mi_lower_actual) * ((Chabrier_IMF(Mi_lower_actual)*TheSFH[i]) + (Chabrier_IMF(SNIIMasses[j+1])*TheSFH[i]))/2.0;
	        	    					NormSNIIMassEjecRate[(STEPS*snap)+step][i][Zi] += (1.0-A_FACTOR) * ((SNIIMasses[j+1]-Mi_lower_actual) * ((SNIIEjectedMasses_lower_actual + SNIIEjectedMasses[Zi_SNII][j+1])/2.0));
	        	    					NormSNIIMetalEjecRate[(STEPS*snap)+step][i][Zi] += (1.0-A_FACTOR) * (SNIIMasses[j+1]-Mi_lower_actual) * ((SNIITotalMetals_lower_actual + SNIITotalMetals[Zi_SNII][j+1])/2.0);
											#ifdef INDIVIDUAL_ELEMENTS
	        	    					int k;
	        	    					for (k=0;k<NUM_ELEMENTS;k++)
	        	    					{
	        	    						int kk;
											#ifndef MAINELEMENTS
	        	    						kk=k;
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (1.0-A_FACTOR) * (SNIIMasses[j+1]-Mi_lower_actual) * ((SNIIYields_lower_actual[k] + SNIIYields[Zi_SNII][kk][j+1])/2.0);
	        	    						#else
	        	    						switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (1.0-A_FACTOR) * (SNIIMasses[j+1]-Mi_lower_actual) * ((SNIIYields_lower_actual[k] + SNIIYields[Zi_SNII][kk][j+1])/2.0);
											#endif
	        	    					}
											#endif //INDIVIDUAL_ELEMENTS
	        	    				}
	        	    				else if (j == Mi_upper_SNII && j != Mi_lower_SNII) //If mass bin j IS NOT the lowest but IS the highest mass bin to be integrated over.
	        	    				{
	        	    					L5a++;
	        	    					//SNIIRate[(STEPS*snap)+step][Zi] += (1.0-A_FACTOR) * (Mi_upper_actual-SNIIMasses[j]) * ((Chabrier_IMF(SNIIMasses[j])*TheSFH[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH[i]))/2.0;
	        	    					NormSNIIMassEjecRate[(STEPS*snap)+step][i][Zi] += (1.0-A_FACTOR) * ((Mi_upper_actual-SNIIMasses[j]) * ((SNIIEjectedMasses[Zi_SNII][j] + SNIIEjectedMasses_upper_actual)/2.0));
	        	    					NormSNIIMetalEjecRate[(STEPS*snap)+step][i][Zi] += (1.0-A_FACTOR) * ((Mi_upper_actual-SNIIMasses[j]) * ((SNIITotalMetals[Zi_SNII][j] + SNIITotalMetals_upper_actual)/2.0));
											#ifdef INDIVIDUAL_ELEMENTS
	        	    					int k;
	        	    					for (k=0;k<NUM_ELEMENTS;k++)
	        	    					{
	        	    						int kk;
											#ifndef MAINELEMENTS
	        	    						kk=k;
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (1.0-A_FACTOR) * (Mi_upper_actual-SNIIMasses[j]) * ((SNIIYields[Zi_SNII][kk][j] + SNIIYields_upper_actual[k])/2.0);
	        	    						#else
	        	    						switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (1.0-A_FACTOR) * (Mi_upper_actual-SNIIMasses[j]) * ((SNIIYields[Zi_SNII][kk][j] + SNIIYields_upper_actual[k])/2.0);
											#endif
	        	    					}
											#endif //INDIVIDUAL_ELEMENTS
	        	    				}
	        	    			}
	        	    			else  //For mass range where only SN-II progenitors are possible
	        	    			{
	        	    				if (j != Mi_lower_SNII && j != Mi_upper_SNII)
	        	    				{
	        	    					L1b++;
	        	    					//SNIIRate[(STEPS*snap)+step][Zi] += (SNIIMasses[j+1]-SNIIMasses[j]) * ((Chabrier_IMF(SNIIMasses[j])*TheSFH[i]) + (Chabrier_IMF(SNIIMasses[j+1])*TheSFH[i]))/2.0;
	        	    					NormSNIIMassEjecRate[(STEPS*snap)+step][i][Zi] += ((SNIIMasses[j+1]-SNIIMasses[j]) * ((SNIIEjectedMasses[Zi_SNII][j] + SNIIEjectedMasses[Zi_SNII][j+1])/2.0));
	        	    					NormSNIIMetalEjecRate[(STEPS*snap)+step][i][Zi] += (SNIIMasses[j+1]-SNIIMasses[j]) * ((SNIITotalMetals[Zi_SNII][j] + SNIITotalMetals[Zi_SNII][j+1])/2.0);
											#ifdef INDIVIDUAL_ELEMENTS
	        	    					int k;
	        	    					for (k=0;k<NUM_ELEMENTS;k++)
	        	    					{
	        	    						int kk;
											#ifndef MAINELEMENTS
	        	    						kk=k;
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (SNIIMasses[j+1]-SNIIMasses[j]) * ((SNIIYields[Zi_SNII][kk][j] + SNIIYields[Zi_SNII][kk][j+1])/2.0);
	        	    						#else
	        	    						switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (SNIIMasses[j+1]-SNIIMasses[j]) * ((SNIIYields[Zi_SNII][kk][j] + SNIIYields[Zi_SNII][kk][j+1])/2.0);
	        	    						#endif
	        	    					}
											#endif //INDIVIDUAL_ELEMENTS
	        	    				}
	        	    				else if (j == Mi_lower_SNII && j == Mi_upper_SNII)
	        	    				{
	        	    					L2b++;
	        	    					//SNIIRate[(STEPS*snap)+step][Zi] += (Mi_upper_actual-Mi_lower_actual) * ((Chabrier_IMF(Mi_lower_actual)*TheSFH[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH[i]))/2.0;
	        	    					NormSNIIMassEjecRate[(STEPS*snap)+step][i][Zi] += ((Mi_upper_actual-Mi_lower_actual) * ((SNIIEjectedMasses_lower_actual + SNIIEjectedMasses_upper_actual)/2.0));
	        	    					NormSNIIMetalEjecRate[(STEPS*snap)+step][i][Zi] += ((Mi_upper_actual-Mi_lower_actual) * ((SNIITotalMetals_lower_actual + SNIITotalMetals_upper_actual)/2.0));
											#ifdef INDIVIDUAL_ELEMENTS
	        	    					int k;
	        	    					for (k=0;k<NUM_ELEMENTS;k++)
	        	    					{
	        	    						int kk;
											#ifndef MAINELEMENTS
	        	    						kk=k;
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (Mi_upper_actual-Mi_lower_actual) * ((SNIIYields_lower_actual[k] + SNIIYields_upper_actual[k])/2.0);
	        	    						#else
	        	    						switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (Mi_upper_actual-Mi_lower_actual) * ((SNIIYields_lower_actual[k] + SNIIYields_upper_actual[k])/2.0);
											#endif
	        	    					}
											#endif //INDIVIDUAL_ELEMENTS
	        	    				}
	        	    				else if (j == Mi_lower_SNII && j != Mi_upper_SNII)
	        	    				{
	        	    					L3b++;
	        	    					//SNIIRate[(STEPS*snap)+step][Zi] += (SNIIMasses[j+1]-Mi_lower_actual) * ((Chabrier_IMF(Mi_lower_actual)*TheSFH[i]) + (Chabrier_IMF(SNIIMasses[j+1])*TheSFH[i]))/2.0;
	        	    					NormSNIIMassEjecRate[(STEPS*snap)+step][i][Zi] += ((SNIIMasses[j+1]-Mi_lower_actual) * ((SNIIEjectedMasses_lower_actual + SNIIEjectedMasses[Zi_SNII][j+1])/2.0));
	        	    					NormSNIIMetalEjecRate[(STEPS*snap)+step][i][Zi] += (SNIIMasses[j+1]-Mi_lower_actual) * ((SNIITotalMetals_lower_actual + SNIITotalMetals[Zi_SNII][j+1])/2.0);
											#ifdef INDIVIDUAL_ELEMENTS
	        	    					int k;
	        	    					for (k=0;k<NUM_ELEMENTS;k++)
	        	    					{
	        	    						int kk;
											#ifndef MAINELEMENTS
	        	    						kk=k;
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (SNIIMasses[j+1]-Mi_lower_actual) * ((SNIIYields_lower_actual[k] + SNIIYields[Zi_SNII][kk][j+1])/2.0);
	        	    						#else
	        	    						switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (SNIIMasses[j+1]-Mi_lower_actual) * ((SNIIYields_lower_actual[k] + SNIIYields[Zi_SNII][kk][j+1])/2.0);
											#endif
	        	    					}
											#endif //INDIVIDUAL_ELEMENTS
	        	    				}
	        	    				else if (j == Mi_upper_SNII && j != Mi_lower_SNII)
	        	    				{
	        	    					L4b++;
	        	    					//SNIIRate[(STEPS*snap)+step][Zi] += (Mi_upper_actual-SNIIMasses[j]) * ((Chabrier_IMF(SNIIMasses[j])*TheSFH[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH[i]))/2.0;
	        	    					NormSNIIMassEjecRate[(STEPS*snap)+step][i][Zi] += ((Mi_upper_actual-SNIIMasses[j]) * ((SNIIEjectedMasses[Zi_SNII][j] + SNIIEjectedMasses_upper_actual)/2.0));
	        	    					NormSNIIMetalEjecRate[(STEPS*snap)+step][i][Zi] += ((Mi_upper_actual-SNIIMasses[j]) * ((SNIITotalMetals[Zi_SNII][j] + SNIITotalMetals_upper_actual)/2.0));
											#ifdef INDIVIDUAL_ELEMENTS
	        	    					int k;
	        	    					for (k=0;k<NUM_ELEMENTS;k++)
	        	    					{
	        	    						int kk;
											#ifndef MAINELEMENTS
	        	    						kk=k;
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (Mi_upper_actual-SNIIMasses[j]) * ((SNIIYields[Zi_SNII][kk][j] + SNIIYields_upper_actual[k])/2.0);
	        	    						#else
	        	    						switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
	        	    						NormSNIIYieldRate[(STEPS*snap)+step][i][Zi][k] += (Mi_upper_actual-SNIIMasses[j]) * ((SNIIYields[Zi_SNII][kk][j] + SNIIYields_upper_actual[k])/2.0);
											#endif
	        	    					}
											#endif //INDIVIDUAL_ELEMENTS
	        	    				}
	        	    			} //else of "if (SNIIMasses[j] <= 16.0)"
	        	    		} //for (j=Mi_lower_SNII;j<=Mi_upper_SNII;j++)
	        	    	} //if (Mi_lower_SNII <= Mi_upper_SNII)

	        	    	//*****************************************
	        			//SNe-Ia (Disc and Bulge):
	        			//*****************************************
#ifdef DTD
	        	    	t_lower_lifetime = Mi_upper; //Lifetime bin (lifetime arrays) corresp. to longest lifetime (lowest mass) of star to 'die' in current timestep, from SFH bin i.
	        	    	t_upper_lifetime = Mi_lower; //Lifetime bin (lifetime arrays) corresp. to shortest lifetime (highest mass) of star to 'die' in current timestep, from SFH bin i.

	        	    	if ((lifetimes[Zi][Mi_lower] > 26.0*1.0e6) && (lifetimes[Zi][Mi_upper] < 21.0*1.0e9)) //P98 values: t_min = 26 Myrs (the value suggested by P98 at Z=0.02). t_max = 21 Gyrs (the value suggested by P98 at Z=0.02)
	        	    	{
	        	        	int j;
	        	        	for (j=t_lower_lifetime-1;j>=t_upper_lifetime;j--)
	        	        	{
	        	        		if (j != t_lower_lifetime && j != t_upper_lifetime) //If lifetime bin j is NEITHER the lowest NOR the highest bin to be integrated over.
        	    				{
        	    					//Calculate the normalised SN-Ia rate:
    	        	        		DTD_lower = DTDcalc(lifetimes[Zi][j+1]) * 1.0e-9; //NOTE: DTDcalc returns SNe/Gyr, not SNe/yr, hence the multiple (1.0e-9).
    	        	        		DTD_upper = DTDcalc(lifetimes[Zi][j]) * 1.0e-9;

    	        	        		//SNIaRate[(STEPS*snap)+step][Zi] += A_FACTOR*F316 * KALPHA * (lifetimes[Zi][j]-lifetimes[Zi][j+1]) * ((DTD_lower * TheSFH[i]) + (DTD_upper * TheSFH[i]))/2.0;
    	        	        		NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi] += A_FACTOR*F316 * KALPHA * (lifetimes[Zi][j]-lifetimes[Zi][j+1]) * ((DTD_lower * SNIAEJECMASS) + (DTD_upper * SNIAEJECMASS))/2.0; //IN [MSun/yr]
    	        	        		NormSNIaMetalEjecRate[(STEPS*snap)+step][i][Zi] = NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi];
											#ifdef INDIVIDUAL_ELEMENTS
    	        	        		int k;
    	        	        		for (k=0;k<NUM_ELEMENTS;k++)
    	        	        		{
    	        	        			int kk;
    									#ifndef MAINELEMENTS
    	        	        			switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=5; break; case 3: kk=6; break; case 4: kk=7; break; case 5: kk=9; break; case 6: kk=11; break; case 7: kk=13; break; case 8: kk=15; break; case 9: kk=19; break; case 10: kk=25; break;}
    	        	        			NormSNIaYieldRate[(STEPS*snap)+step][i][Zi][k] += A_FACTOR*F316 * KALPHA * (lifetimes[Zi][j]-lifetimes[Zi][j+1]) * ((DTD_lower * SNIaYields[kk]) + (DTD_upper * SNIaYields[kk]))/2.0; //IN [Msun/yr]
    									#else
    	        	        			switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=7; break; case 3: kk=11; break; case 4: kk=25; break;}
    	        	        			NormSNIaYieldRate[(STEPS*snap)+step][i][Zi][k] += A_FACTOR*F316 * KALPHA * (lifetimes[Zi][j]-lifetimes[Zi][j+1]) * ((DTD_lower * SNIaYields[kk]) + (DTD_upper * SNIaYields[kk]))/2.0; //IN [Msun/yr]
											#endif
    	        	        		}
											#endif //INDIVIDUAL_ELEMENTS
        	    				}
        	    				else if (j == t_lower_lifetime && j == t_upper_lifetime) //If lifetime bin j is BOTH the lowest AND the highest bin to be integrated over.
        	    				{
    	        	        		DTD_lower = DTDcalc(t_lower) * 1.0e-9;
    	        	        		DTD_upper = DTDcalc(t_upper) * 1.0e-9;

    	        	        		//SNIaRate[(STEPS*snap)+step][Zi] += A_FACTOR*F316 * KALPHA * (t_upper-t_lower) * ((DTD_lower * TheSFH[i]) + (DTD_upper * TheSFH[i]))/2.0;
    	        	        		NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi] += A_FACTOR*F316 * KALPHA * (t_upper-t_lower) * ((DTD_lower * SNIAEJECMASS) + (DTD_upper * SNIAEJECMASS))/2.0; //IN [MSun/yr]
    	        	        		NormSNIaMetalEjecRate[(STEPS*snap)+step][i][Zi] = NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi];
											#ifdef INDIVIDUAL_ELEMENTS
    	        	        		int k;
    	        	        		for (k=0;k<NUM_ELEMENTS;k++)
    	        	        		{
    	        	        			int kk;
    									#ifndef MAINELEMENTS
    	        	        			switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=5; break; case 3: kk=6; break; case 4: kk=7; break; case 5: kk=9; break; case 6: kk=11; break; case 7: kk=13; break; case 8: kk=15; break; case 9: kk=19; break; case 10: kk=25; break;}
    	        	        			NormSNIaYieldRate[(STEPS*snap)+step][i][Zi][k] += A_FACTOR*F316 * KALPHA * (t_upper-t_lower) * ((DTD_lower * SNIaYields[kk]) + (DTD_upper * SNIaYields[kk]))/2.0; //IN [Msun/yr]
    									#else
    	        	        			switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=7; break; case 3: kk=11; break; case 4: kk=25; break;}
    	        	        			NormSNIaYieldRate[(STEPS*snap)+step][i][Zi][k] += A_FACTOR*F316 * KALPHA * (t_upper-t_lower) * ((DTD_lower * SNIaYields[kk]) + (DTD_upper * SNIaYields[kk]))/2.0; //IN [Msun/yr]
											#endif
    	        	        		}
											#endif //INDIVIDUAL_ELEMENTS
        	    				}
        	    				else if (j == t_lower_lifetime && j != t_upper_lifetime) //If lifetime bin j IS the lowest bin, but IS NOT the highest bin to be integrated over.
        	    				{
    	        	        		DTD_lower = DTDcalc(t_lower) * 1.0e-9;
    	        	        		DTD_upper = DTDcalc(lifetimes[Zi][j]) * 1.0e-9;

    	        	        		//SNIaRate[(STEPS*snap)+step][Zi] += A_FACTOR*F316 * KALPHA * (lifetimes[Zi][j]-t_lower) * ((DTD_lower * TheSFH[i]) + (DTD_upper * TheSFH[i]))/2.0;
    	        	        		NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi] += A_FACTOR*F316 * KALPHA * (lifetimes[Zi][j]-t_lower) * ((DTD_lower * SNIAEJECMASS) + (DTD_upper * SNIAEJECMASS))/2.0; //IN [MSun/yr]
    	        	        		NormSNIaMetalEjecRate[(STEPS*snap)+step][i][Zi] = NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi];
											#ifdef INDIVIDUAL_ELEMENTS
    	        	        		int k;
    	        	        		for (k=0;k<NUM_ELEMENTS;k++)
    	        	        		{
    	        	        			int kk;
    									#ifndef MAINELEMENTS
    	        	        			switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=5; break; case 3: kk=6; break; case 4: kk=7; break; case 5: kk=9; break; case 6: kk=11; break; case 7: kk=13; break; case 8: kk=15; break; case 9: kk=19; break; case 10: kk=25; break;}
    	        	        			NormSNIaYieldRate[(STEPS*snap)+step][i][Zi][k] += A_FACTOR*F316 * KALPHA * (lifetimes[Zi][j]-t_lower) * ((DTD_lower * SNIaYields[kk]) + (DTD_upper * SNIaYields[kk]))/2.0; //IN [Msun/yr]
    									#else
    	        	        			switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=7; break; case 3: kk=11; break; case 4: kk=25; break;}
    	        	        			NormSNIaYieldRate[(STEPS*snap)+step][i][Zi][k] += A_FACTOR*F316 * KALPHA * (lifetimes[Zi][j]-t_lower) * ((DTD_lower * SNIaYields[kk]) + (DTD_upper * SNIaYields[kk]))/2.0; //IN [Msun/yr]
											#endif
    	        	        		}
											#endif //INDIVIDUAL_ELEMENTS
        	    				}
        	    				else if (j != t_lower_lifetime && j == t_upper_lifetime) //If lifetime bin j IS NOT the lowest bin, but IS the highest bin to be integrated over.
        	    				{
    	        	        		DTD_lower = DTDcalc(lifetimes[Zi][j+1]) * 1.0e-9;
    	        	        		DTD_upper = DTDcalc(t_upper) * 1.0e-9;

    	        	        		//SNIaRate[(STEPS*snap)+step][Zi] += A_FACTOR*F316 * KALPHA * (t_upper-lifetimes[Zi][j+1]) * ((DTD_lower * TheSFH[i]) + (DTD_upper * TheSFH[i]))/2.0;
    	        	        		NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi] += A_FACTOR*F316 * KALPHA * (t_upper-lifetimes[Zi][j+1]) * ((DTD_lower * SNIAEJECMASS) + (DTD_upper * SNIAEJECMASS))/2.0; //IN [MSun/yr]
    	        	        		NormSNIaMetalEjecRate[(STEPS*snap)+step][i][Zi] = NormSNIaMassEjecRate[(STEPS*snap)+step][i][Zi];
											#ifdef INDIVIDUAL_ELEMENTS
    	        	        		int k;
    	        	        		for (k=0;k<NUM_ELEMENTS;k++)
    	        	        		{
    	        	        			int kk;
    									#ifndef MAINELEMENTS
    	        	        			switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=5; break; case 3: kk=6; break; case 4: kk=7; break; case 5: kk=9; break; case 6: kk=11; break; case 7: kk=13; break; case 8: kk=15; break; case 9: kk=19; break; case 10: kk=25; break;}
    	        	        			NormSNIaYieldRate[(STEPS*snap)+step][i][Zi][k] += A_FACTOR*F316 * KALPHA * (t_upper-lifetimes[Zi][j+1]) * ((DTD_lower * SNIaYields[kk]) + (DTD_upper * SNIaYields[kk]))/2.0; //IN [Msun/yr]
    									#else
    	        	        			switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=7; break; case 3: kk=11; break; case 4: kk=25; break;}
    	        	        			NormSNIaYieldRate[(STEPS*snap)+step][i][Zi][k] += A_FACTOR*F316 * KALPHA * (t_upper-lifetimes[Zi][j+1]) * ((DTD_lower * SNIaYields[kk]) + (DTD_upper * SNIaYields[kk]))/2.0; //IN [Msun/yr]
											#endif
    	        	        		}
											#endif //INDIVIDUAL_ELEMENTS
        	    				}
	        	        	} //for (j=t_lower_lifetime;j>=t_upper_lifetime+1;j--)

	        	    	} //if ((lifetimes[Zi][Mi_lower] > 0.02*1.0e9) && (lifetimes[Zi][Mi_upper] < 11.5*1.0e9))


#endif
	        			//*****************************************
	        			//AGB Winds (Disc and Bulge):
	        			//*****************************************
	        	    	Zi_AGB = find_initial_metallicity_comp(Zi, i, 4);

	        			Mi_lower_AGB = max_Mi_lower(Mi_lower,4);
	        			Mi_upper_AGB = min_Mi_upper(Mi_upper,4);

	        	    	if (Mi_lower_AGB <= Mi_upper_AGB)
	        	    	{
	        	    		Mi_lower_AGB = find_agb_mass_bin(lifetimeMasses[Mi_lower_AGB]);
	        	    		Mi_upper_AGB = find_agb_mass_bin(lifetimeMasses[Mi_upper_AGB]);

#ifdef INDIVIDUAL_ELEMENTS
	        	    		find_actual_ejecta_limits(4, Mi_lower_actual, Mi_upper_actual, Mi_lower_AGB, Mi_upper_AGB, Zi_AGB,
	        	    				&AGBEjectedMasses_lower_actual, &AGBEjectedMasses_upper_actual, &AGBTotalMetals_lower_actual, &AGBTotalMetals_upper_actual,
	        	    				AGBYields_lower_actual, AGBYields_upper_actual);
#else
	        	    		find_actual_ejecta_limits(4, Mi_lower_actual, Mi_upper_actual, Mi_lower_AGB, Mi_upper_AGB, Zi_AGB,
	        	    				&AGBEjectedMasses_lower_actual, &AGBEjectedMasses_upper_actual, &AGBTotalMetals_lower_actual, &AGBTotalMetals_upper_actual);
#endif

	        	    		int j;
	        	    		for (j=Mi_lower_AGB;j<=Mi_upper_AGB;j++)
	        	    		{
        	    				if (j != Mi_lower_AGB && j != Mi_upper_AGB)
        	    				{
        	    					//AGBRate[(STEPS*snap)+step][Zi] += (AGBMasses[j+1]-AGBMasses[j]) * ((Chabrier_IMF(AGBMasses[j])*TheSFH[i]) + (Chabrier_IMF(AGBMasses[j+1])*TheSFH[i]))/2.0;
        	    					NormAGBMassEjecRate[(STEPS*snap)+step][i][Zi] += ((AGBMasses[j+1]-AGBMasses[j]) * ((AGBEjectedMasses[Zi_AGB][j] + AGBEjectedMasses[Zi_AGB][j+1])/2.0));
        	    					NormAGBMetalEjecRate[(STEPS*snap)+step][i][Zi] += (AGBMasses[j+1]-AGBMasses[j]) * ((AGBTotalMetals[Zi_AGB][j] + AGBTotalMetals[Zi_AGB][j+1])/2.0);
										#ifdef INDIVIDUAL_ELEMENTS
        	    					int k;
        	    					for (k=0;k<NUM_ELEMENTS;k++)
        	    					{
        	    						int kk;
										#ifndef MAINELEMENTS
        	    						kk=k;
        	    						NormAGBYieldRate[(STEPS*snap)+step][i][Zi][k] += (AGBMasses[j+1]-AGBMasses[j]) * ((AGBYields[Zi_AGB][kk][j] + AGBYields[Zi_AGB][kk][j+1])/2.0);
        	    						#else
        	    						switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
        	    						NormAGBYieldRate[(STEPS*snap)+step][i][Zi][k] += (AGBMasses[j+1]-AGBMasses[j]) * ((AGBYields[Zi_AGB][kk][j] + AGBYields[Zi_AGB][kk][j+1])/2.0);
        	    						#endif
        	    					}
										#endif //INDIVIDUAL_ELEMENTS
        	    				}
        	    				else if (j == Mi_lower_AGB && j == Mi_upper_AGB && Mi_lower_actual <= AGB_MAX_MASS)
        	    				{
        	    					//AGBRate[(STEPS*snap)+step][Zi] += (Mi_upper_actual-Mi_lower_actual) * ((Chabrier_IMF(Mi_lower_actual)*TheSFH[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH[i]))/2.0;
        	    					NormAGBMassEjecRate[(STEPS*snap)+step][i][Zi] += ((Mi_upper_actual-Mi_lower_actual) * ((AGBEjectedMasses_lower_actual + AGBEjectedMasses_upper_actual)/2.0));
        	    					NormAGBMetalEjecRate[(STEPS*snap)+step][i][Zi] += (Mi_upper_actual-Mi_lower_actual) * ((AGBTotalMetals_lower_actual + AGBTotalMetals_upper_actual)/2.0);
        	    					float Sum_NormAGBYieldRate; Sum_NormAGBYieldRate = 0.0;
										#ifdef INDIVIDUAL_ELEMENTS
        	    					int k;
        	    					for (k=0;k<NUM_ELEMENTS;k++)
        	    					{
        	    						int kk;
										#ifndef MAINELEMENTS
        	    						kk=k;
        	    						NormAGBYieldRate[(STEPS*snap)+step][i][Zi][k] += (Mi_upper_actual-Mi_lower_actual) * ((AGBYields_lower_actual[k] + AGBYields_upper_actual[k])/2.0);
        	    						if (k != 0 && k != 1) {Sum_NormAGBYieldRate += NormAGBYieldRate[(STEPS*snap)+step][i][Zi][k];}
        	    						#else
        	    						switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
        	    						NormAGBYieldRate[(STEPS*snap)+step][i][Zi][k] += (Mi_upper_actual-Mi_lower_actual) * ((AGBYields_lower_actual[k] + AGBYields_upper_actual[k])/2.0);
										#endif
        	    					}
										#endif //INDIVIDUAL_ELEMENTS
        	    				}
        	    				else if (j == Mi_lower_AGB && j != Mi_upper_AGB)
        	    				{
        	    					//AGBRate[(STEPS*snap)+step][Zi] += (AGBMasses[j+1]-Mi_lower_actual) * ((Chabrier_IMF(Mi_lower_actual)*TheSFH[i]) + (Chabrier_IMF(AGBMasses[j+1])*TheSFH[i]))/2.0;
        	    					NormAGBMassEjecRate[(STEPS*snap)+step][i][Zi] += ((AGBMasses[j+1]-Mi_lower_actual) * ((AGBEjectedMasses_lower_actual + AGBEjectedMasses[Zi_AGB][j+1])/2.0));
        	    					NormAGBMetalEjecRate[(STEPS*snap)+step][i][Zi] += (AGBMasses[j+1]-Mi_lower_actual) * ((AGBTotalMetals_lower_actual + AGBTotalMetals[Zi_AGB][j+1])/2.0);
										#ifdef INDIVIDUAL_ELEMENTS
        	    					int k;
        	    					for (k=0;k<NUM_ELEMENTS;k++)
        	    					{
        	    						int kk;
										#ifndef MAINELEMENTS
        	    						kk=k;
        	    						NormAGBYieldRate[(STEPS*snap)+step][i][Zi][k] += (AGBMasses[j+1]-Mi_lower_actual) * ((AGBYields_lower_actual[k] + AGBYields[Zi_AGB][kk][j+1])/2.0);
        	    						#else
        	    						switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
        	    						NormAGBYieldRate[(STEPS*snap)+step][i][Zi][k] += (AGBMasses[j+1]-Mi_lower_actual) * ((AGBYields_lower_actual[k] + AGBYields[Zi_AGB][kk][j+1])/2.0);
										#endif
        	    					}
										#endif //INDIVIDUAL_ELEMENTS
        	    				}
        	    				else if (j == Mi_upper_AGB && j != Mi_lower_AGB)
        	    				{
        	    					//AGBRate[(STEPS*snap)+step][Zi] += (Mi_upper_actual-AGBMasses[j]) * ((Chabrier_IMF(AGBMasses[j])*TheSFH[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH[i]))/2.0;
        	    					NormAGBMassEjecRate[(STEPS*snap)+step][i][Zi] += ((Mi_upper_actual-AGBMasses[j]) * ((AGBEjectedMasses[Zi_AGB][j] + AGBEjectedMasses_upper_actual)/2.0));
        	    					NormAGBMetalEjecRate[(STEPS*snap)+step][i][Zi] += ((Mi_upper_actual-AGBMasses[j]) * ((AGBTotalMetals[Zi_AGB][j] + AGBTotalMetals_upper_actual)/2.0));
										#ifdef INDIVIDUAL_ELEMENTS
        	    					int k;
        	    					for (k=0;k<NUM_ELEMENTS;k++)
        	    					{
        	    						int kk;
										#ifndef MAINELEMENTS
        	    						kk=k;
        	    						NormAGBYieldRate[(STEPS*snap)+step][i][Zi][k] += (Mi_upper_actual-AGBMasses[j]) * ((AGBYields[Zi_AGB][kk][j] + AGBYields_upper_actual[k])/2.0);
        	    						#else
        	    						switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
        	    						NormAGBYieldRate[(STEPS*snap)+step][i][Zi][k] += (Mi_upper_actual-AGBMasses[j]) * ((AGBYields[Zi_AGB][kk][j] + AGBYields_upper_actual[k])/2.0);
										#endif
        	    					}
										#endif //INDIVIDUAL_ELEMENTS
        	    				}
	        	    		} //for (j=Mi_lower_AGB;j<=Mi_upper_AGB;j++)
	        	    	} //if (Mi_lower_AGB <= Mi_upper_AGB)
	        	    	
	        	   } //if (t_upper >= lifetimes[Zi][Mi_lower+1])

	        	} //for (Zi=0;Zi<LIFETIME_Z_NUM;Zi++)
	          } //for (j=0;j<=width_in_timesteps;j++) //MINI_BINS	
	        } //for (i=0;i<=SFH_ibin_structure[(SFH_NBIN*snap)+step];i++)
	        	
	    } //for(step=0;step<STEPS;step++)
	} //for(snap=0;snap<(LastDarkMatterSnapShot+1);snap++)
	/*printf("L1a = %i\n", L1a);
	printf("L2a = %i\n", L2a);
	printf("L3a = %i\n", L3a);
	printf("L4a = %i\n", L4a);
	printf("L5a = %i\n", L5a);
	printf("L1b = %i\n", L1b);
	printf("L2b = %i\n", L2b);
	printf("L3b = %i\n", L3b);
	printf("L4b = %i\n", L4b);*/
	printf("Yield integrals calculated.\n");
}

int find_initial_metallicity_comp(int Zi, int sfh_bin, int table_type)
{
	int i, Zi_bin;
	double Z_in;

	Zi_bin = -1;
	i = 0;
	Z_in = lifetimeMetallicities[Zi];

	switch (table_type)
	{
		case 1: //Lifetime metallicity table
			while (Zi_bin == -1)
			{
				if (lifetimeMetallicities[i] < Z_in)
				{
					i++;
					if (i == LIFETIME_Z_NUM) Zi_bin = i-1; //If galaxy's Z is higher than max Z from table, then just take max Z from table
				}
				else Zi_bin = i;
			}
			break;
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
		//case 3 //SNIa yields are NOT metallicity dependent
		case 4: //AGB metallicity table
			while (Zi_bin == -1)
			{
				if (AGBMetallicities[i] < Z_in)
				{
					i++;
					if (i == AGB_Z_NUM) Zi_bin = i-1;
				}
				else Zi_bin = i;
			}
			break;
	}
	return Zi_bin;
}

int find_initial_mass(double lifetime, int Zi_bin)
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

int max_Mi_lower(int Mi_lower, int channel_type)
{
	switch (channel_type)
		{
			case 2: //SNII mass limits
					if (lifetimeMasses[Mi_lower] > SNII_MIN_MASS) return Mi_lower;
					else
					{
						int i;
						i = 0;
						do { i++; }
						while (lifetimeMasses[i] < SNII_MIN_MASS);
						return i;
					}
					break;
			#ifndef DTD
			case 3: //SNIa mass limits
					if (lifetimeMasses[Mi_lower] > 0.85) return Mi_lower; //NB: Lifetimes of SNe-Ia binaries depend on M2, not Mb. (ie: 0.85<=M2<=8.0)
					else
					{
						int i;
						i = 0;
						do { i++; }
						while (lifetimeMasses[i] < 0.85); //NB: Lifetimes of SNe-Ia binaries depend on M2, not Mb. (ie: 0.85<=M2<=8.0)
						return i;
					}
					break;
			#endif
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

int min_Mi_upper(int Mi_upper, int channel_type)
{
	switch (channel_type)
		{
			case 2: //SNII mass limits
					if (lifetimeMasses[Mi_upper] < SNII_MAX_MASS) return Mi_upper;
					//else return LIFETIME_MASS_NUM-1;
					else
					{
						int i;
						i = LIFETIME_MASS_NUM-1;
						do { i--; }
						while (lifetimeMasses[i] > SNII_MAX_MASS);
						return i;
					}
					break;
			#ifndef DTD
			case 3: //SNIa mass limits
					if (lifetimeMasses[Mi_upper] < 8.0) return Mi_upper; //NB: Lifetimes of SNe-Ia binaries depends on M2, not Mb. (ie: 0.85<=M2<=8.0)
					else
					{
						int i;
						i = LIFETIME_MASS_NUM-1;
						do { i--; }
						while (lifetimeMasses[i] > 8.0); //NB: Lifetimes of SNe-Ia binaries depends on M2, not Mb. (ie: 0.85<=M2<=8.0)
						return i;
					}
					break;
			#endif
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

int find_SNII_mass_bin(double masslimit)
{
    if (masslimit == SNII_MAX_MASS) return SNII_MASS_NUM-1;
    else
    {
	int Mi_bin;

	Mi_bin = -1;
	int i;
	i = 0;
	while (Mi_bin == -1)
	{
		if (SNIIMasses[i] < masslimit)
		{
			i++;
			if (i == SNII_MASS_NUM) Mi_bin = i-1; //If mass is greater than max mass for SNe-II (shouldn't be), then just return max mass (120.0 Msun)
		}
		else Mi_bin = i;
	}
	return Mi_bin;
    }
}

int find_agb_mass_bin(double masslimit)
{
	if (masslimit == AGB_MAX_MASS) return AGB_MASS_NUM-1;
	else
	{
	int Mi_bin;

	Mi_bin = -1;
	int i;
	i = 0;
	while (Mi_bin == -1)
	{
		if (AGBMasses[i] < masslimit)
		{
			i++;
			if (i == AGB_MASS_NUM) Mi_bin = i-1; //If mass is greater than max mass for AGB winds (shouldn't be), then just return max mass (5.0 Msun)
		}
		else Mi_bin = i;
	}
	return Mi_bin;
	}
}

#ifdef DTD
double DTDcalc (double timevalue)
{
#ifdef BIMODALDTD
	double timevalueM, DTDvalueM;
	timevalueM = log10(timevalue); //IN [log(YEARS)]
	if (timevalueM < 7.93) //Characteristic time == 10^7.93 yrs (NB: timevalue is in log10(yrs) here)
	{
		DTDvalueM = (1.4 - 50.0*(timevalueM-7.7)*(timevalueM-7.7)) / DTD_NORM; //When using Mannucci bi-modal DTD
		//DTDvalueM = (0.74 - 50.0*(timevalueM-7.7)*(timevalueM-7.7)) / DTD_NORM; //When using CUSTOM bi-modal DTD, with 20% in prompt component
	}
	else
	{
		DTDvalueM = (-0.8 - 0.9*(timevalueM-8.7)*(timevalueM-8.7)) / DTD_NORM; //Same eqn for Custom bi-modal DTD too
	}
	return pow(10.0,DTDvalueM); //IN [SNe/Gyr]
#endif

#ifdef CUSTOMDTD
	double timevalueM, DTDvalueM;
	timevalueM = log10(timevalue); //IN [log(YEARS)]
	if (timevalueM < 7.93) //Characteristic time == 10^7.93 yrs (NB: timevalue is in log10(yrs) here)
	{
		DTDvalueM = (0.74 - 50.0*(timevalueM-7.7)*(timevalueM-7.7)) / DTD_NORM; //When using CUSTOM bi-modal DTD, with 20% in prompt component
	}
	else
	{
		DTDvalueM = (-0.8 - 0.9*(timevalueM-8.7)*(timevalueM-8.7)) / DTD_NORM;
	}
	return pow(10.0,DTDvalueM); //IN [SNe/Gyr]
#endif

#ifdef GAUSSIANDTD
	double timevalueG, pivalue, DTDvalueG; //, tauCharac, sigmatd
	timevalueG = timevalue/1.0e9; //IN [Gyrs]
	pivalue = 3.14159;
	DTDvalueG = ((1./sqrt(2.*pivalue*SIGMA_TD*SIGMA_TD)) * exp(-(pow((timevalueG-TAUCHARAC),2.))/(2.*SIGMA_TD*SIGMA_TD)));// / DTD_NORM;
	return DTDvalueG; //IN [SNe/Gyr]
#endif

#ifdef POWERLAWDTD
	double timevalueP, DTDvalueP;
	timevalueP = timevalue/1.0e9; //IN [Gyrs]
	DTDvalueP = pow(timevalueP, DTD_SLOPE) / DTD_NORM;
	return DTDvalueP; //IN [SNe/Gyr]
#endif

#ifdef RUITERDTD
	double timevalueR, pivalue, DTDvalueR;
	timevalueR = timevalue/1.0e9; //IN [Gyrs]
	if (timevalueR <= 1.0) //Time between Gaussian and power-law components = 1.0 Gyrs
	{
		pivalue = 3.14159;
		DTDvalueR = (0.143 * (1./sqrt(2.*pivalue*SIGMA_TD*SIGMA_TD)) * exp(-(pow((timevalueR-TAUCHARAC),2.))/(2.*SIGMA_TD*SIGMA_TD))) / DTD_NORM; //0.143 factor ensures Gaussian component is 13% of total
	}
	if (timevalueR > 1.0)
	{
		DTDvalueR = pow(timevalueR, DTD_SLOPE) / DTD_NORM;
	}
	return DTDvalueR; //IN [SNe/Gyr]
#endif
}
#endif

#ifdef INDIVIDUAL_ELEMENTS
void find_actual_ejecta_limits(int channel_type, double Mi_lower_actual, double Mi_upper_actual, int Mi_lower, int Mi_upper, int Zi,
		double* EjectedMasses_lower_actual, double* EjectedMasses_upper_actual, double* TotalMetals_lower_actual, double* TotalMetals_upper_actual,
		double* Yields_lower_actual, double* Yields_upper_actual)
#else
void find_actual_ejecta_limits(int channel_type, double Mi_lower_actual, double Mi_upper_actual, int Mi_lower, int Mi_upper, int Zi,
		double* EjectedMasses_lower_actual, double* EjectedMasses_upper_actual, double* TotalMetals_lower_actual, double* TotalMetals_upper_actual)
#endif
{
	switch (channel_type)
	{
	case 2: //SNII
		if (Mi_lower == 0)
		{
			*EjectedMasses_lower_actual = SNIIEjectedMasses[Zi][0]; //If Mi_lower_actual was less than or equal to Mi_lower (e.g. 7 Msun bin for SNe-II), then EjectedMasses_lower_actual is set to ejected mass from an e.g. 7 Msun star.
		    *TotalMetals_lower_actual = SNIITotalMetals[Zi][0]; //If Mi_lower_actual was less than or equal to Mi_lower (e.g. 7 Msun bin for SNe-II), then TotalMetals_lower_actual is set to ejected mass in metals from an e.g. 7 Msun star.
#ifdef INDIVIDUAL_ELEMENTS
		    int k;
		    for (k=0;k<NUM_ELEMENTS;k++)
		    {
		    	int kk;
				#ifndef MAINELEMENTS
		    	kk=k;
		    	Yields_lower_actual[k] = SNIIYields[Zi][kk][0];
				#else
		    	switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
		    	Yields_lower_actual[k] = SNIIYields[Zi][kk][0];
				#endif
		    }
#endif //INDIVIDUAL_ELEMENTS
		}
		else
		{
			*EjectedMasses_lower_actual = SNIIEjectedMasses[Zi][Mi_lower] + ((SNIIEjectedMasses[Zi][Mi_lower+1]-SNIIEjectedMasses[Zi][Mi_lower]) * ((Mi_lower_actual-SNIIMasses[Mi_lower])/(SNIIMasses[Mi_lower+1]-SNIIMasses[Mi_lower])));
		    *TotalMetals_lower_actual = SNIITotalMetals[Zi][Mi_lower] + ((SNIITotalMetals[Zi][Mi_lower+1]-SNIITotalMetals[Zi][Mi_lower]) * ((Mi_lower_actual-SNIIMasses[Mi_lower])/(SNIIMasses[Mi_lower+1]-SNIIMasses[Mi_lower])));
#ifdef INDIVIDUAL_ELEMENTS
		    int k;
		    for (k=0;k<NUM_ELEMENTS;k++)
		    {
		    	int kk;
				#ifndef MAINELEMENTS
		    	kk=k;
		        Yields_lower_actual[k] = SNIIYields[Zi][kk][Mi_lower] + ((SNIIYields[Zi][kk][Mi_lower+1]-SNIIYields[Zi][kk][Mi_lower]) * ((Mi_lower_actual-SNIIMasses[Mi_lower])/(SNIIMasses[Mi_lower+1]-SNIIMasses[Mi_lower])));
				#else
		        switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
		        Yields_lower_actual[k] = SNIIYields[Zi][kk][Mi_lower] + ((SNIIYields[Zi][kk][Mi_lower+1]-SNIIYields[Zi][kk][Mi_lower]) * ((Mi_lower_actual-SNIIMasses[Mi_lower])/(SNIIMasses[Mi_lower+1]-SNIIMasses[Mi_lower])));
				#endif
		    }
#endif //INDIVIDUAL_ELEMENTS
		}

		if (Mi_upper == SNII_MASS_NUM-1)
		{
			*EjectedMasses_upper_actual = SNIIEjectedMasses[Zi][SNII_MASS_NUM-1]; //If Mi_upper_actual was more than or equal to 120 Msun, then EjectedMasses_upper_actual is set to ejected mass from 120 Msun star.
		    *TotalMetals_upper_actual = SNIITotalMetals[Zi][SNII_MASS_NUM-1]; //If Mi_upper_actual was more than or equal to 120 Msun, then TotalMetals_upper_actual is set to ejected mass in metals from 120 Msun star.
#ifdef INDIVIDUAL_ELEMENTS
		    int k;
		    for (k=0;k<NUM_ELEMENTS;k++)
		    {
		    	int kk;
				#ifndef MAINELEMENTS
		    	kk=k;
		    	Yields_upper_actual[k] = SNIIYields[Zi][kk][SNII_MASS_NUM-1];
				#else
		    	switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
		    	Yields_upper_actual[k] = SNIIYields[Zi][kk][SNII_MASS_NUM-1];
				#endif
		    }
#endif //INDIVIDUAL_ELEMENTS
		}
		else
		{
		    *EjectedMasses_upper_actual = SNIIEjectedMasses[Zi][Mi_upper] + ((SNIIEjectedMasses[Zi][Mi_upper+1]-SNIIEjectedMasses[Zi][Mi_upper]) * ((Mi_upper_actual-SNIIMasses[Mi_upper])/(SNIIMasses[Mi_upper+1]-SNIIMasses[Mi_upper])));
		    *TotalMetals_upper_actual = SNIITotalMetals[Zi][Mi_upper] + ((SNIITotalMetals[Zi][Mi_upper+1]-SNIITotalMetals[Zi][Mi_upper]) * ((Mi_upper_actual-SNIIMasses[Mi_upper])/(SNIIMasses[Mi_upper+1]-SNIIMasses[Mi_upper])));
#ifdef INDIVIDUAL_ELEMENTS
		    int k;
		    for (k=0;k<NUM_ELEMENTS;k++)
		    {
		    	int kk;
				#ifndef MAINELEMENTS
		    	kk=k;
		    	Yields_upper_actual[k] = SNIIYields[Zi][kk][Mi_upper] + ((SNIIYields[Zi][kk][Mi_upper+1]-SNIIYields[Zi][kk][Mi_upper]) * ((Mi_upper_actual-SNIIMasses[Mi_upper])/(SNIIMasses[Mi_upper+1]-SNIIMasses[Mi_upper])));
				#else
		        switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
		        Yields_upper_actual[k] = SNIIYields[Zi][kk][Mi_upper] + ((SNIIYields[Zi][kk][Mi_upper+1]-SNIIYields[Zi][kk][Mi_upper]) * ((Mi_upper_actual-SNIIMasses[Mi_upper])/(SNIIMasses[Mi_upper+1]-SNIIMasses[Mi_upper])));
				#endif
		     }
#endif //INDIVIDUAL_ELEMENTS
		}
	break;

	case 4: //AGB
		if (Mi_lower == 0)
		{
			*EjectedMasses_lower_actual = AGBEjectedMasses[Zi][0]; //If Mi_lower_actual was less than or equal to Mi_lower (e.g. 0.85 Msun bin for AGB), then EjectedMasses_lower_actual is set to ejected mass from an e.g. 0.85 Msun star.
		    *TotalMetals_lower_actual = AGBTotalMetals[Zi][0]; //If Mi_lower_actual was less than or equal to Mi_lower (e.g. 0.85 Msun bin for AGB), then TotalMetals_lower_actual is set to ejected mass in metals from an e.g. 0.85 Msun star.
#ifdef INDIVIDUAL_ELEMENTS
		    int k;
		    for (k=0;k<NUM_ELEMENTS;k++)
		    {
		    	int kk;
				#ifndef MAINELEMENTS
		    	kk=k;
		    	Yields_lower_actual[k] = AGBYields[Zi][kk][0];
				#else
		    	switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
		    	Yields_lower_actual[k] = AGBYields[Zi][kk][0];
				#endif
		    }
#endif //INDIVIDUAL_ELEMENTS
		}
		else
		{
			*EjectedMasses_lower_actual = AGBEjectedMasses[Zi][Mi_lower] + ((AGBEjectedMasses[Zi][Mi_lower+1]-AGBEjectedMasses[Zi][Mi_lower]) * ((Mi_lower_actual-AGBMasses[Mi_lower])/(AGBMasses[Mi_lower+1]-AGBMasses[Mi_lower])));
		    *TotalMetals_lower_actual = AGBTotalMetals[Zi][Mi_lower] + ((AGBTotalMetals[Zi][Mi_lower+1]-AGBTotalMetals[Zi][Mi_lower]) * ((Mi_lower_actual-AGBMasses[Mi_lower])/(AGBMasses[Mi_lower+1]-AGBMasses[Mi_lower])));
#ifdef INDIVIDUAL_ELEMENTS
		    int k;
		    for (k=0;k<NUM_ELEMENTS;k++)
		    {
		    	int kk;
				#ifndef MAINELEMENTS
		    	kk=k;
		        Yields_lower_actual[k] = AGBYields[Zi][kk][Mi_lower] + ((AGBYields[Zi][kk][Mi_lower+1]-AGBYields[Zi][kk][Mi_lower]) * ((Mi_lower_actual-AGBMasses[Mi_lower])/(AGBMasses[Mi_lower+1]-AGBMasses[Mi_lower])));
				#else
		        switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
		        Yields_lower_actual[k] = AGBYields[Zi][kk][Mi_lower] + ((AGBYields[Zi][kk][Mi_lower+1]-AGBYields[Zi][kk][Mi_lower]) * ((Mi_lower_actual-AGBMasses[Mi_lower])/(AGBMasses[Mi_lower+1]-AGBMasses[Mi_lower])));
				#endif
		    }
#endif //INDIVIDUAL_ELEMENTS
		}

		if (Mi_upper == AGB_MASS_NUM-1)
		{
			*EjectedMasses_upper_actual = AGBEjectedMasses[Zi][AGB_MASS_NUM-1]; //If Mi_upper_actual was more than or equal to 5 Msun, then EjectedMasses_upper_actual is set to ejected mass from 5 Msun star.
		    *TotalMetals_upper_actual = AGBTotalMetals[Zi][AGB_MASS_NUM-1]; //If Mi_upper_actual was more than or equal to 5 Msun, then TotalMetals_upper_actual is set to ejected mass in metals from 5 Msun star.
#ifdef INDIVIDUAL_ELEMENTS
		    int k;
		    for (k=0;k<NUM_ELEMENTS;k++)
		    {
		    	int kk;
				#ifndef MAINELEMENTS
		    	kk=k;
		    	Yields_upper_actual[k] = AGBYields[Zi][kk][AGB_MASS_NUM-1];
				#else
		    	switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
		    	Yields_upper_actual[k] = AGBYields[Zi][kk][AGB_MASS_NUM-1];
				#endif
		    }
#endif //INDIVIDUAL_ELEMENTS
		}
		else
		{
		    *EjectedMasses_upper_actual = AGBEjectedMasses[Zi][Mi_upper] + ((AGBEjectedMasses[Zi][Mi_upper+1]-AGBEjectedMasses[Zi][Mi_upper]) * ((Mi_upper_actual-AGBMasses[Mi_upper])/(AGBMasses[Mi_upper+1]-AGBMasses[Mi_upper])));
		    *TotalMetals_upper_actual = AGBTotalMetals[Zi][Mi_upper] + ((AGBTotalMetals[Zi][Mi_upper+1]-AGBTotalMetals[Zi][Mi_upper]) * ((Mi_upper_actual-AGBMasses[Mi_upper])/(AGBMasses[Mi_upper+1]-AGBMasses[Mi_upper])));
#ifdef INDIVIDUAL_ELEMENTS
		    int k;
		    for (k=0;k<NUM_ELEMENTS;k++)
		    {
		    	int kk;
				#ifndef MAINELEMENTS
		    	kk=k;
		    	Yields_upper_actual[k] = AGBYields[Zi][kk][Mi_upper] + ((AGBYields[Zi][kk][Mi_upper+1]-AGBYields[Zi][kk][Mi_upper]) * ((Mi_upper_actual-AGBMasses[Mi_upper])/(AGBMasses[Mi_upper+1]-AGBMasses[Mi_upper])));
				#else
		        switch(k){case 0: kk=0; break; case 1: kk=1; break; case 2: kk=4; break; case 3: kk=6; break; case 4: kk=10; break;}
		        Yields_upper_actual[k] = AGBYields[Zi][kk][Mi_upper] + ((AGBYields[Zi][kk][Mi_upper+1]-AGBYields[Zi][kk][Mi_upper]) * ((Mi_upper_actual-AGBMasses[Mi_upper])/(AGBMasses[Mi_upper+1]-AGBMasses[Mi_upper])));
				#endif
		     }
#endif //INDIVIDUAL_ELEMENTS
		}
		break;
	}
}

#endif //DETAILED_EMRICHMENT

