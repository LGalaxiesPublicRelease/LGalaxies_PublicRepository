/*
 * calc_SNe_rates.c
 *
 * Calculates the SNe rates using SFhs of one-timestep resolution. To be compared with rates using SFH bins in yield_intergrals.c
 *
 *  Created on: Aug 6, 2012
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

void SNe_rates()
{
	double previoustime, newtime, deltaT;
	int snap, step,i;
	double timet;

	int Mi_lower, Mi_upper, Mi_lower_SNII, Mi_upper_SNII, Zi_SNII, Zi_AGB, Mi_lower_AGB, Mi_upper_AGB, t_lower_lifetime, t_upper_lifetime;
	double dt, t_lower, t_upper, DTD_lower, DTD_upper;
	double Mi_lower_actual, Mi_upper_actual, SNIIEjectedMasses_lower_actual, SNIIEjectedMasses_upper_actual, SNIITotalMetals_lower_actual, SNIITotalMetals_upper_actual, SNIIYields_lower_actual[NUM_ELEMENTS], SNIIYields_upper_actual[NUM_ELEMENTS];
	double AGBEjectedMasses_lower_actual, AGBEjectedMasses_upper_actual, AGBTotalMetals_lower_actual, AGBTotalMetals_upper_actual, AGBYields_lower_actual[NUM_ELEMENTS], AGBYields_upper_actual[NUM_ELEMENTS];

	FRACCOUNTA = 0;

	int L1a,L2a,L3a,L4a,L5a,L1b,L2b,L3b,L4b;
	L1a=0;L2a=0;L3a=0;L4a=0;L5a=0;L1b=0;L2b=0;L3b=0;L4b=0;

	int counta;
	TheSFH2[0] = 1.0;
	//TheSFH2[0] = 1.0/(tau_dt[0]*UnitTime_in_years/Hubble_h);
	for(counta=1;counta<(STEPS*MAXSNAPS);counta++)
	{
		TheSFH2[counta] = 0.0;// /(tau_dt[counta]*UnitTime_in_years/Hubble_h);
	}
	//TheSFH2[0] = 0.0; TheSFH2[1] = 3.8; TheSFH2[2] = 5.0; TheSFH2[3] = 5.5; TheSFH2[4] = 5.2; TheSFH2[5] = 5.0; TheSFH2[6] = 4.8; TheSFH2[7] = 4.2; TheSFH2[8] = 4.0; TheSFH2[9] = 3.8; TheSFH2[10] = 3.2; TheSFH2[11] = 3.0; TheSFH2[12] = 2.8; TheSFH2[13] = 2.5; TheSFH2[14] = 2.2; TheSFH2[15] = 2.2; TheSFH2[16] = 2.0; TheSFH2[17] = 1.8; TheSFH2[18] = 1.8; TheSFH2[19] = 1.5;

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
	}
	else if (SNII_MAX_MASS == 70.0)
	{
		KALPHA = 1.52679;
	}
	else if (SNII_MAX_MASS == 50.0)
	{
		KALPHA = 1.56408;
		F316 = 0.0408018;
	}
	else if (SNII_MAX_MASS == 40.0)
	{
		KALPHA = 1.59203;
	}
	else
	{
		KALPHA = 1.5;
		printf("****************\n");
		printf("In calc_SNe_rates.c:\n");
		printf("\nSNII_MAX_MASS neither 40, 50, 70, 100, nor 120 Msun\n");
		printf("KALPHA set to 1.50\n");
		printf("****************\n\n");
	}

	for(snap=0;snap<MAXSNAPS-1;snap++) //LOOP OVER SNAPSHOTS
	{
	    previoustime = NumToTime(snap); //Time to z=0 from start of current snapshot
	    newtime = NumToTime(snap+1); //Time to z=0 from end of current snapshot
	    deltaT = previoustime - newtime;

	    for(step=0;step<STEPS;step++) //LOOP OVER TIMESTEPS
	    {
	    	dt = deltaT/STEPS;  //Time-width of a timestep in current snapshot
	    	timet = previoustime - (step + 0.5) * dt; //Time from middle of the current timestep to z=0

	    	for (i=0; i<(STEPS*snap)+step; i++) //LOOP OVER SFH-TIMESTEPS
	    	{
				t_lower = (tau_t[i] + (0.5*tau_dt[i]) - timet - (0.5*dt)) * UnitTime_in_years/Hubble_h; //IN YEARS //Time from middle of SFH-timestep i to high-z edge of current timestep
	        	t_upper = (tau_t[i] + (0.5*tau_dt[i]) - timet + (0.5*dt)) * UnitTime_in_years/Hubble_h; //IN YEARS //Time from middle of SFH-timestep i to low-z edge of current timestep

	        	//if (snap==0) {printf("Timestep: %i, SFH-timestep: %i, TheSFH2[%i] = %f, t_l = %f, t_u = %f [Myrs]\n", (STEPS*snap)+step, i, i, TheSFH2[i], t_lower/1.0e6, t_upper/1.0e6);}

	        	int Zi;
	        	//for (Zi=1;Zi<2;Zi++)
	        	for (Zi=0;Zi<LIFETIME_Z_NUM;Zi++) //LOOP OVER POSSIBLE INITIAL METALLICITIES (from lifetime tables)
	        	{
					Mi_lower = find_initial_mass2(t_upper, Zi); //Mass bin number (lifetime arrays) corresp. to lowest mass of star to 'die' in current timestep, from SFH bin i.
					Mi_upper = find_initial_mass2(t_lower, Zi); //Mass bin number (lifetime arrays) corresp. to highest mass of star to 'die' in current timestep, from SFH bin i.
					//if (snap==0) {printf("Zi = %i, M_l = %i, %f, M_u = %i, %f\n", Zi, Mi_lower, lifetimeMasses[Mi_lower], Mi_upper, lifetimeMasses[Mi_upper]);}
					/*if (lifetimeMasses[Mi_upper] >= SNII_MAX_MASS)
					{
						int i;
						i = 0;
						while (lifetimeMasses[i] < SNII_MAX_MASS) { i++; }
						Mi_upper = i;
					}*/

					Mi_lower_actual = lifetimeMasses[Mi_lower] + ((lifetimeMasses[Mi_lower+1]-lifetimeMasses[Mi_lower]) * ((t_upper-lifetimes[Zi][Mi_lower])/(lifetimes[Zi][Mi_lower+1]-lifetimes[Zi][Mi_lower]))); //IN MSUN  //Lowest mass of star to 'die' in current timestep from SFH bin i.
					Mi_upper_actual = lifetimeMasses[Mi_upper] + ((lifetimeMasses[Mi_upper+1]-lifetimeMasses[Mi_upper]) * ((t_lower-lifetimes[Zi][Mi_upper])/(lifetimes[Zi][Mi_upper+1]-lifetimes[Zi][Mi_upper]))); //IN MSUN  //Highest mass of star to 'die' in current timestep from SFH bin i.

					if (Mi_lower == 0) Mi_lower_actual = AGB_MIN_MASS; //No stars below 0.85 Msun contribute to chemical enrichment.
					//if (Mi_upper == LIFETIME_MASS_NUM-1) Mi_upper_actual = SNII_MAX_MASS; //No stars above 120 Msun assumed to exist.
					if (lifetimeMasses[Mi_upper] >= SNII_MAX_MASS) Mi_upper_actual = SNII_MAX_MASS; //No stars of mass above max. SN-II progenitor assumed to exist.

					if (t_upper >= lifetimes[Zi][LIFETIME_MASS_NUM-1]) //If the longest time from SFH bin i to curent timestep is shorter than the shortest possible lifetime, there's no enrichment, so skip calculations.
					{
	        			//*****************************************
	        			//SNe-II (Disc and Bulge):
	        			//*****************************************
	        	    	Zi_SNII = find_initial_metallicity_comp(Zi, i, 2); //Metallicity bin (SNe-II arrays) corresp. to metallicity Zi.

	        	    	//Check if mass range is within range for SN-II progenitor stars:
	        			Mi_lower_SNII = max_Mi_lower(Mi_lower,2); //ROB: Should we send Mi_lower_actual and Mi_upper_actual to max_Mi_lower and min_Mi_upper? (24-07-120
	        			Mi_upper_SNII = min_Mi_upper(Mi_upper,2);

	        	    	if (Mi_lower_SNII <= Mi_upper_SNII)
	        	    	{
	        	    		Mi_lower_SNII = find_SNII_mass_bin(lifetimeMasses[Mi_lower_SNII]); //Mass bin (SNe-II arrays) corresp. to lowest mass of star to 'die' in current timestep from SFH bin i.
	        	    		Mi_upper_SNII = find_SNII_mass_bin(lifetimeMasses[Mi_upper_SNII]); //Mass bin (SNe-II arrays) corresp. to highest mass of star to 'die' in current timestep from SFH bin i.

	        	    		//Find true yields at the true upper and lower masses, given by 'Mi_upper_actual' and 'Mi_lower_actual':
	        	    		find_actual_ejecta_limits(2, Mi_lower_actual, Mi_upper_actual, Mi_lower_SNII, Mi_upper_SNII, Zi_SNII,&SNIIEjectedMasses_lower_actual, &SNIIEjectedMasses_upper_actual, &SNIITotalMetals_lower_actual, &SNIITotalMetals_upper_actual, SNIIYields_lower_actual, SNIIYields_upper_actual);

	        	    		/*if(Mi_upper_actual <= SNIA_MAX_MASS)
	        	    		{
	        	    			SNIIRate2[(STEPS*snap)+step][Zi] += (1.0-A_FACTOR) * (Mi_upper_actual-Mi_lower_actual) * ((Chabrier_IMF(Mi_lower_actual)*TheSFH2[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH2[i]))/2.0;
	        	    		}
	        	    		else if (Mi_upper_actual > SNIA_MAX_MASS && Mi_lower_actual <= SNIA_MAX_MASS)
	        	    		{
	        	    			SNIIRate2[(STEPS*snap)+step][Zi] += (1.0-A_FACTOR) * (SNIA_MAX_MASS-Mi_lower_actual) * ((Chabrier_IMF(Mi_lower_actual)*TheSFH2[i]) + (Chabrier_IMF(SNIA_MAX_MASS)*TheSFH2[i]))/2.0;
	        	    			SNIIRate2[(STEPS*snap)+step][Zi] += (Mi_upper_actual-SNIA_MAX_MASS) * ((Chabrier_IMF(SNIA_MAX_MASS)*TheSFH2[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH2[i]))/2.0;
	        	    		}
	        	    		else if (Mi_lower_actual > SNIA_MAX_MASS)
	        	    		{
	        	    			SNIIRate2[(STEPS*snap)+step][Zi] += (Mi_upper_actual-Mi_lower_actual) * ((Chabrier_IMF(Mi_lower_actual)*TheSFH2[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH2[i]))/2.0;
	        	    		}*/

	        	    		//NUMERICALLY INTEGRATE OVER THE MASS RANGE APPLICABLE FOR SNe-II:
	        	    		int j;
	        	    		for (j=Mi_lower_SNII;j<=Mi_upper_SNII;j++)
	        	    		{
	        	    			if (SNIIMasses[j] <= SNIA_MAX_MASS) //For mass range where both SN-II and SN-Ia progenitors are possible
	        	    			{
	        	    				if (j != Mi_lower_SNII && j != Mi_upper_SNII) //If mass bin j is NEITHER the lowest NOR highest mass bin to be integrated over.
	        	    				{
	        	    					L1a++;
	        	    					SNIIRate2[(STEPS*snap)+step][Zi] += (1.0-A_FACTOR) * (SNIIMasses[j+1]-SNIIMasses[j]) * ((Chabrier_IMF(SNIIMasses[j])*TheSFH2[i]) + (Chabrier_IMF(SNIIMasses[j+1])*TheSFH2[i]))/2.0;
	        	    				}
	        	    				else if (j == Mi_lower_SNII && j == Mi_upper_SNII && Mi_lower_actual >= SNII_MIN_MASS) //If mass bin j is BOTH the lowest AND highest mass bin to be integrated over.
	        	    				{
	        	    					L2a++;
	        	    					SNIIRate2[(STEPS*snap)+step][Zi] += (1.0-A_FACTOR) * (Mi_upper_actual-Mi_lower_actual) * ((Chabrier_IMF(Mi_lower_actual)*TheSFH2[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH2[i]))/2.0;
	        	    				}
	        	    				else if (j == Mi_lower_SNII && j == Mi_upper_SNII && Mi_lower_actual < SNII_MIN_MASS) //If mass bin j is BOTH the lowest AND highest mass bin to be integrated over, AND 'Mi_lower_actual' is below min. mass for SNe-II. (Still only counts stars from SNII_MIN_MASS to Mi_upper_actual).
	        	    				{
	        	    					L3a++;
	        	    					SNIIRate2[(STEPS*snap)+step][Zi] += (1.0-A_FACTOR) * (Mi_upper_actual-SNII_MIN_MASS) * ((Chabrier_IMF(SNII_MIN_MASS)*TheSFH2[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH2[i]))/2.0;
	        	    				}
	        	    				else if (j == Mi_lower_SNII && j != Mi_upper_SNII) //If mass bin j IS the lowest but IS NOT the highest mass bin to be integrated over.
	        	    				{
	        	    					L4a++;
	        	    					SNIIRate2[(STEPS*snap)+step][Zi] += (1.0-A_FACTOR) * (SNIIMasses[j+1]-Mi_lower_actual) * ((Chabrier_IMF(Mi_lower_actual)*TheSFH2[i]) + (Chabrier_IMF(SNIIMasses[j+1])*TheSFH2[i]))/2.0;
	        	    				}
	        	    				else if (j == Mi_upper_SNII && j != Mi_lower_SNII) //If mass bin j IS NOT the lowest but IS the highest mass bin to be integrated over.
	        	    				{
	        	    					L5a++;
	        	    					SNIIRate2[(STEPS*snap)+step][Zi] += (1.0-A_FACTOR) * (Mi_upper_actual-SNIIMasses[j]) * ((Chabrier_IMF(SNIIMasses[j])*TheSFH2[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH2[i]))/2.0;
	        	    				}
	        	    			}
	        	    			else  //For mass range where only SN-II progenitors are possible
	        	    			{
	        	    				if (j != Mi_lower_SNII && j != Mi_upper_SNII)
	        	    				{
	        	    					L1b++;
	        	    					SNIIRate2[(STEPS*snap)+step][Zi] += (SNIIMasses[j+1]-SNIIMasses[j]) * ((Chabrier_IMF(SNIIMasses[j])*TheSFH2[i]) + (Chabrier_IMF(SNIIMasses[j+1])*TheSFH2[i]))/2.0;
	        	    				}
	        	    				else if (j == Mi_lower_SNII && j == Mi_upper_SNII)
	        	    				{
	        	    					L2b++;
	        	    					SNIIRate2[(STEPS*snap)+step][Zi] += (Mi_upper_actual-Mi_lower_actual) * ((Chabrier_IMF(Mi_lower_actual)*TheSFH2[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH2[i]))/2.0;
	        	    				}
	        	    				else if (j == Mi_lower_SNII && j != Mi_upper_SNII)
	        	    				{
	        	    					L3b++;
	        	    					SNIIRate2[(STEPS*snap)+step][Zi] += (SNIIMasses[j+1]-Mi_lower_actual) * ((Chabrier_IMF(Mi_lower_actual)*TheSFH2[i]) + (Chabrier_IMF(SNIIMasses[j+1])*TheSFH2[i]))/2.0;
	        	    				}
	        	    				else if (j == Mi_upper_SNII && j != Mi_lower_SNII)
	        	    				{
	        	    					L4b++;
	        	    					SNIIRate2[(STEPS*snap)+step][Zi] += (Mi_upper_actual-SNIIMasses[j]) * ((Chabrier_IMF(SNIIMasses[j])*TheSFH2[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH2[i]))/2.0;
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

	        	    	if ((lifetimes[Zi][Mi_lower] > 0.01*1.0e9) && (lifetimes[Zi][Mi_upper] < 20.6*1.0e9)) //Min and max lifetimes for 2ndary star in a SN-Ia-producing binary (SD scenario).
	        	    	{
	        	        	int j;
	        	        	for (j=t_lower_lifetime;j>=t_upper_lifetime;j--)
	        	        	{
        	    				if (j != t_lower_lifetime && j != t_upper_lifetime) //If lifetime bin j is NEITHER the lowest NOR the highest bin to be integrated over.
        	    				{
        	    					//Calculate the normalised SN-Ia rate:
    	        	        		DTD_lower = DTDcalc(lifetimes[Zi][j+1]) * 1.0e-9; //NOTE: DTDcalc returns SNe/Gyr, not SNe/yr, hence the multiple (1.0e-9).
    	        	        		DTD_upper = DTDcalc(lifetimes[Zi][j]) * 1.0e-9;
    	        	        		SNIaRate2[(STEPS*snap)+step][Zi] += A_FACTOR*F316 * KALPHA * (lifetimes[Zi][j]-lifetimes[Zi][j+1]) * ((DTD_lower * TheSFH2[i]) + (DTD_upper * TheSFH2[i]))/2.0;
        	    				}
        	    				else if (j == t_lower_lifetime && j == t_upper_lifetime) //If lifetime bin j is BOTH the lowest AND the highest bin to be integrated over.
        	    				{
    	        	        		DTD_lower = DTDcalc(t_lower) * 1.0e-9;
    	        	        		DTD_upper = DTDcalc(t_upper) * 1.0e-9;
    	        	        		SNIaRate2[(STEPS*snap)+step][Zi] += A_FACTOR*F316 * KALPHA * (t_upper-t_lower) * ((DTD_lower * TheSFH2[i]) + (DTD_upper * TheSFH2[i]))/2.0;
        	    				}
        	    				else if (j == t_lower_lifetime && j != t_upper_lifetime) //If lifetime bin j IS the lowest bin, but IS NOT the highest bin to be integrated over.
        	    				{
    	        	        		DTD_lower = DTDcalc(t_lower) * 1.0e-9;
    	        	        		DTD_upper = DTDcalc(lifetimes[Zi][j]) * 1.0e-9;
    	        	        		SNIaRate2[(STEPS*snap)+step][Zi] += A_FACTOR*F316 * KALPHA * (lifetimes[Zi][j]-t_lower) * ((DTD_lower * TheSFH2[i]) + (DTD_upper * TheSFH2[i]))/2.0;
        	    				}
        	    				else if (j != t_lower_lifetime && j == t_upper_lifetime) //If lifetime bin j IS NOT the lowest bin, but IS the highest bin to be integrated over.
        	    				{
    	        	        		DTD_lower = DTDcalc(lifetimes[Zi][j+1]) * 1.0e-9;
    	        	        		DTD_upper = DTDcalc(t_upper) * 1.0e-9;
    	        	        		SNIaRate2[(STEPS*snap)+step][Zi] += A_FACTOR*F316 * KALPHA * (t_upper-lifetimes[Zi][j+1]) * ((DTD_lower * TheSFH2[i]) + (DTD_upper * TheSFH2[i]))/2.0;
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

	        	    		find_actual_ejecta_limits(4, Mi_lower_actual, Mi_upper_actual, Mi_lower_AGB, Mi_upper_AGB, Zi_AGB, &AGBEjectedMasses_lower_actual, &AGBEjectedMasses_upper_actual, &AGBTotalMetals_lower_actual, &AGBTotalMetals_upper_actual, AGBYields_lower_actual, AGBYields_upper_actual);

	        	    		int j;
	        	    		for (j=Mi_lower_AGB;j<=Mi_upper_AGB;j++)
	        	    		{
        	    				if (j != Mi_lower_AGB && j != Mi_upper_AGB)
        	    				{
        	    					AGBRate2[(STEPS*snap)+step][Zi] += (AGBMasses[j+1]-AGBMasses[j]) * ((Chabrier_IMF(AGBMasses[j])*TheSFH2[i]) + (Chabrier_IMF(AGBMasses[j+1])*TheSFH2[i]))/2.0;
        	    				}
        	    				else if (j == Mi_lower_AGB && j == Mi_upper_AGB && Mi_lower_actual <= AGB_MAX_MASS)
        	    				{
        	    					AGBRate2[(STEPS*snap)+step][Zi] += (Mi_upper_actual-Mi_lower_actual) * ((Chabrier_IMF(Mi_lower_actual)*TheSFH2[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH2[i]))/2.0;
        	    				}
        	    				else if (j == Mi_lower_AGB && j != Mi_upper_AGB)
        	    				{
        	    					AGBRate2[(STEPS*snap)+step][Zi] += (AGBMasses[j+1]-Mi_lower_actual) * ((Chabrier_IMF(Mi_lower_actual)*TheSFH2[i]) + (Chabrier_IMF(AGBMasses[j+1])*TheSFH2[i]))/2.0;
        	    				}
        	    				else if (j == Mi_upper_AGB && j != Mi_lower_AGB)
        	    				{
        	    					AGBRate2[(STEPS*snap)+step][Zi] += (Mi_upper_actual-AGBMasses[j]) * ((Chabrier_IMF(AGBMasses[j])*TheSFH2[i]) + (Chabrier_IMF(Mi_upper_actual)*TheSFH2[i]))/2.0;
        	    				}
	        	    		} //for (j=Mi_lower_AGB;j<=Mi_upper_AGB;j++)
	        	    	} //if (Mi_lower_AGB <= Mi_upper_AGB)
	        	   } //if (t_upper >= lifetimes[Zi][Mi_lower+1])
	        	} //for (Zi=0;Zi<LIFETIME_Z_NUM;Zi++)
	        } //for (i=0; i<(STEPS*snap)+step; i++)
	    	//printf("%.11f, ", SNIIRate2[(STEPS*snap)+step][3]/(tau_dt[(STEPS*snap)+step]*UnitTime_in_years/Hubble_h));
	    	//printf("Timestep %i: Cosmic time %f: %.11f\n", (STEPS*snap)+step, 13.569513 - tau_t[(STEPS*snap)+step]*UnitTime_in_years/Hubble_h/1.0e9, SNIIRate2[(STEPS*snap)+step][1]);
	    	//printf("%.16f, ", SNIaRate2[(STEPS*snap)+step][3]/(tau_dt[(STEPS*snap)+step]*UnitTime_in_years/Hubble_h));
	        //printf("%.16f, ", AGBRate2[(STEPS*snap)+step][3]/(tau_dt[(STEPS*snap)+step]*UnitTime_in_years/Hubble_h));
	        //if (snap < 63) printf("%f, ", (NumToTime(snap)-(((NumToTime(snap)-NumToTime(snap+1))/STEPS)*step))*UnitTime_in_years/Hubble_h);
	        //if (snap == 63) printf("%f, ", (NumToTime(snap)-(((NumToTime(snap))/STEPS)*step))*UnitTime_in_years/Hubble_h);
	        //printf("%i, ", (STEPS*snap)+step);
	    } //for(step=0;step<STEPS;step++)
	} //for(snap=0;snap<MAXSNAPS;snap++)
	/*printf("L1a = %i\n", L1a);
	printf("L2a = %i\n", L2a);
	printf("L3a = %i\n", L3a);
	printf("L4a = %i\n", L4a);
	printf("L5a = %i\n", L5a);
	printf("L1b = %i\n", L1b);
	printf("L2b = %i\n", L2b);
	printf("L3b = %i\n", L3b);
	printf("L4b = %i\n", L4b);*/
	printf("\nSNe rates calculated.\n");
}

int find_initial_metallicity_comp2(int Zi, int sfh_bin, int table_type)
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

int find_initial_mass2(double lifetime, int Zi_bin)
{
	//printf("lifetime = %f\n", lifetime);
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

/*int max_Mi_lower2(int Mi_lower, int channel_type)
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

int min_Mi_upper2(int Mi_upper, int channel_type)
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

int find_SNII_mass_bin2(double masslimit)
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
}*/

/*int find_agb_mass_bin2(double masslimit)
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
double DTDcalc2 (double timevalue)
{
#ifdef BIMODALDTD
	double timevalueM, DTDvalueM;
	timevalueM = log10(timevalue); //IN [log(YEARS)]
	if (timevalueM < 7.93) //Characteristic time == 10^7.93 yrs (NB: timevalue is in log10(yrs) here)
	{
		DTDvalueM = (1.4 - 50.0*(timevalueM-7.7)*(timevalueM-7.7)) / DTD_NORM;
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
	//tauCharac = 1.0; //Now set in allvars.h
	//sigmatd = 0.2*TAUCHARAC; //Now set in allvars.h
	//DTDvalueG = (1./sqrt(2.*pivalue*sigmatd*sigmatd)) * exp(-(pow((timevalueG-tauCharac),2.))/(2.*sigmatd*sigmatd));
	DTDvalueG = (1./sqrt(2.*pivalue*SIGMA_TD*SIGMA_TD)) * exp(-(pow((timevalueG-TAUCHARAC),2.))/(2.*SIGMA_TD*SIGMA_TD));
	return DTDvalueG; //IN [SNe/Gyr]
#endif

#ifdef POWERLAWDTD
	double timevalueP, DTDvalueP;
	timevalueP = timevalue/1.0e9; //IN [Gyrs]
	DTDvalueP = pow(timevalueP, DTD_SLOPE) / DTD_NORM;
	return DTDvalueP; //IN [SNe/Gyr]
#endif
}
#endif*/

/*void find_actual_ejecta_limits2(int channel_type, double Mi_lower_actual, double Mi_upper_actual, int Mi_lower, int Mi_upper, int Zi, double* EjectedMasses_lower_actual, double* EjectedMasses_upper_actual, double* TotalMetals_lower_actual, double* TotalMetals_upper_actual, double* Yields_lower_actual, double* Yields_upper_actual)
{
	switch (channel_type)
	{
	case 2: //SNII
		if (Mi_lower == 0)
		{
			*EjectedMasses_lower_actual = SNIIEjectedMasses[Zi][0]; //If Mi_lower_actual was less than or equal to Mi_lower (e.g. 7 Msun bin for SNe-II), then EjectedMasses_lower_actual is set to ejected mass from an e.g. 7 Msun star.
		    *TotalMetals_lower_actual = SNIITotalMetals[Zi][0]; //If Mi_lower_actual was less than or equal to Mi_lower (e.g. 7 Msun bin for SNe-II), then TotalMetals_lower_actual is set to ejected mass in metals from an e.g. 7 Msun star.
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
		}
		else
		{
			*EjectedMasses_lower_actual = SNIIEjectedMasses[Zi][Mi_lower] + ((SNIIEjectedMasses[Zi][Mi_lower+1]-SNIIEjectedMasses[Zi][Mi_lower]) * ((Mi_lower_actual-SNIIMasses[Mi_lower])/(SNIIMasses[Mi_lower+1]-SNIIMasses[Mi_lower])));
		    *TotalMetals_lower_actual = SNIITotalMetals[Zi][Mi_lower] + ((SNIITotalMetals[Zi][Mi_lower+1]-SNIITotalMetals[Zi][Mi_lower]) * ((Mi_lower_actual-SNIIMasses[Mi_lower])/(SNIIMasses[Mi_lower+1]-SNIIMasses[Mi_lower])));
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
		}

		if (Mi_upper == SNII_MASS_NUM-1)
		{
			*EjectedMasses_upper_actual = SNIIEjectedMasses[Zi][SNII_MASS_NUM-1]; //If Mi_upper_actual was more than or equal to 120 Msun, then EjectedMasses_upper_actual is set to ejected mass from 120 Msun star.
		    *TotalMetals_upper_actual = SNIITotalMetals[Zi][SNII_MASS_NUM-1]; //If Mi_upper_actual was more than or equal to 120 Msun, then TotalMetals_upper_actual is set to ejected mass in metals from 120 Msun star.
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
		}
		else
		{
		    *EjectedMasses_upper_actual = SNIIEjectedMasses[Zi][Mi_upper] + ((SNIIEjectedMasses[Zi][Mi_upper+1]-SNIIEjectedMasses[Zi][Mi_upper]) * ((Mi_upper_actual-SNIIMasses[Mi_upper])/(SNIIMasses[Mi_upper+1]-SNIIMasses[Mi_upper])));
		    *TotalMetals_upper_actual = SNIITotalMetals[Zi][Mi_upper] + ((SNIITotalMetals[Zi][Mi_upper+1]-SNIITotalMetals[Zi][Mi_upper]) * ((Mi_upper_actual-SNIIMasses[Mi_upper])/(SNIIMasses[Mi_upper+1]-SNIIMasses[Mi_upper])));
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
		}
	break;

	case 4: //AGB
		if (Mi_lower == 0)
		{
			*EjectedMasses_lower_actual = AGBEjectedMasses[Zi][0]; //If Mi_lower_actual was less than or equal to Mi_lower (e.g. 0.85 Msun bin for AGB), then EjectedMasses_lower_actual is set to ejected mass from an e.g. 0.85 Msun star.
		    *TotalMetals_lower_actual = AGBTotalMetals[Zi][0]; //If Mi_lower_actual was less than or equal to Mi_lower (e.g. 0.85 Msun bin for AGB), then TotalMetals_lower_actual is set to ejected mass in metals from an e.g. 0.85 Msun star.
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
		}
		else
		{
			*EjectedMasses_lower_actual = AGBEjectedMasses[Zi][Mi_lower] + ((AGBEjectedMasses[Zi][Mi_lower+1]-AGBEjectedMasses[Zi][Mi_lower]) * ((Mi_lower_actual-AGBMasses[Mi_lower])/(AGBMasses[Mi_lower+1]-AGBMasses[Mi_lower])));
		    *TotalMetals_lower_actual = AGBTotalMetals[Zi][Mi_lower] + ((AGBTotalMetals[Zi][Mi_lower+1]-AGBTotalMetals[Zi][Mi_lower]) * ((Mi_lower_actual-AGBMasses[Mi_lower])/(AGBMasses[Mi_lower+1]-AGBMasses[Mi_lower])));
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
		}

		if (Mi_upper == AGB_MASS_NUM-1)
		{
			*EjectedMasses_upper_actual = AGBEjectedMasses[Zi][AGB_MASS_NUM-1]; //If Mi_upper_actual was more than or equal to 5 Msun, then EjectedMasses_upper_actual is set to ejected mass from 5 Msun star.
		    *TotalMetals_upper_actual = AGBTotalMetals[Zi][AGB_MASS_NUM-1]; //If Mi_upper_actual was more than or equal to 5 Msun, then TotalMetals_upper_actual is set to ejected mass in metals from 5 Msun star.
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
		}
		else
		{
		    *EjectedMasses_upper_actual = AGBEjectedMasses[Zi][Mi_upper] + ((AGBEjectedMasses[Zi][Mi_upper+1]-AGBEjectedMasses[Zi][Mi_upper]) * ((Mi_upper_actual-AGBMasses[Mi_upper])/(AGBMasses[Mi_upper+1]-AGBMasses[Mi_upper])));
		    *TotalMetals_upper_actual = AGBTotalMetals[Zi][Mi_upper] + ((AGBTotalMetals[Zi][Mi_upper+1]-AGBTotalMetals[Zi][Mi_upper]) * ((Mi_upper_actual-AGBMasses[Mi_upper])/(AGBMasses[Mi_upper+1]-AGBMasses[Mi_upper])));
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
		}
		break;
	}
}*/

#endif //#ifdef DETAILED_METALS_AND_MASS_RETURN


