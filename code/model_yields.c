/*
 * recipe_yields.c
 *
 *  Created on: 18.11.2011
 *      Author: robyates
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"


void update_yields_and_return_mass(int p, int centralgal, double dt, int nstep)
{
	int Zi;
	double timestep_width; //Width of current timestep in CODE UNITS
	int TimeBin; //Bin in Yield arrays corresponding to current timestep
	double Zi_disp, NormSNIIMassEjecRate_actual, NormSNIaMassEjecRate_actual, NormAGBMassEjecRate_actual, NormSNIIMetalEjecRate_actual, NormSNIaMetalEjecRate_actual, NormAGBMetalEjecRate_actual;
#ifdef INDIVIDUAL_ELEMENTS
	double NormSNIIYieldRate_actual[NUM_ELEMENTS], NormSNIaYieldRate_actual[NUM_ELEMENTS], NormAGBYieldRate_actual[NUM_ELEMENTS];
#endif
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

	TotalMassReturnedToColdDiskGas=0.0;
	TotalMassReturnedToHotGas=0.0;

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

    int i;
    for (i=0;i<=Gal[p].sfh_ibin;i++) //LOOP OVER SFH BINS
    {
    	sfh_time=Gal[p].sfh_t[i]+(0.5*Gal[p].sfh_dt[i]);
    	//time_to_ts = ((sfh_time+(0.5*Gal[p].sfh_dt[i])) - timet)*(UnitTime_in_years/Hubble_h)/1.0e6; //Time from high-z (upper) edge of SFH bin to middle of current timestep [in Myrs]
    	//tcut = 2.0*((Gal[p].Rvir/Gal[p].Vvir)/0.0001); //Maximum lifetime of stars that have their ejected put straight into the HotGas [in Myrs]



    //*****************************************
    //ENRICHMENT FROM DISK STARS INTO COLD GAS:
    //*****************************************
    if (Gal[p].sfh_DiskMass[i] > 0.0)
    {
     	//pre-calculations to speed up the code
    	DiskSFR = Gal[p].sfh_DiskMass[i]/Gal[p].sfh_dt[i];
    	step_width_times_DiskSFR = timestep_width * DiskSFR;
    	DiskSFR_physical_units = DiskSFR * (1.0e10/Hubble_h);
    	step_width_times_DiskSFR_physical_units = timestep_width * DiskSFR_physical_units;
    	inverse_DiskMass_physical_units=Hubble_h/(Gal[p].sfh_DiskMass[i]*1.0e10);
    	Disk_total_metallicity=metals_total(Gal[p].sfh_MetalsDiskMass[i])/Gal[p].sfh_DiskMass[i];


    	Zi = find_initial_metallicity(p, i, 1, 1);
    	//Interpolate the disk metallicity on the lifetimeMetallicities tables:
    	Zi_disp = (Disk_total_metallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.

    	NormSNIIMassEjecRate_actual = NormSNIIMassEjecRate[TimeBin][i][Zi] + ((NormSNIIMassEjecRate[TimeBin][i][Zi+1] - NormSNIIMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIaMassEjecRate_actual = NormSNIaMassEjecRate[TimeBin][i][Zi] + ((NormSNIaMassEjecRate[TimeBin][i][Zi+1] - NormSNIaMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormAGBMassEjecRate_actual = NormAGBMassEjecRate[TimeBin][i][Zi] + ((NormAGBMassEjecRate[TimeBin][i][Zi+1] - NormAGBMassEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIIMetalEjecRate_actual = NormSNIIMetalEjecRate[TimeBin][i][Zi] + ((NormSNIIMetalEjecRate[TimeBin][i][Zi+1] - NormSNIIMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormSNIaMetalEjecRate_actual = NormSNIaMetalEjecRate[TimeBin][i][Zi] + ((NormSNIaMetalEjecRate[TimeBin][i][Zi+1] - NormSNIaMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
    	NormAGBMetalEjecRate_actual = NormAGBMetalEjecRate[TimeBin][i][Zi] + ((NormAGBMetalEjecRate[TimeBin][i][Zi+1] - NormAGBMetalEjecRate[TimeBin][i][Zi])*Zi_disp);

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
#endif //SNIATOHOT
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

    	/*if (p==0 && Gal[p].sfh_ICM[i] > 0.0) {printf("Bulge:\n");}
    	if (p==0 && Gal[p].sfh_ICM[i] > 0.0) {printf("%.11f | %.11f %.11f %.11f\n", NormMassEjecRateSumAllTypes, NormSNIIMassEjecRate_actual, NormSNIaMassEjecRate_actual, NormAGBMassEjecRate_actual);}
    	if (p==0 && Gal[p].sfh_ICM[i] > 0.0) {printf("%.11f | %.11f %.11f %.11f | %.11f\n", max(0.0, step_width_times_BulgeSFR * NormMassEjecRateSumAllTypes), max(0.0, step_width_times_BulgeSFR * (NormSNIIMetalEjecRate_actual + (Bulge_total_metallicity * NormSNIIMassEjecRate_actual))), max(0.0, step_width_times_BulgeSFR * NormSNIaMetalEjecRate_actual), max(0.0, step_width_times_BulgeSFR * (NormAGBMetalEjecRate_actual + (Bulge_total_metallicity * NormAGBMassEjecRate_actual))), max(0.0, step_width_times_BulgeSFR * NormSNIaMetalEjecRate_actual) + max(0.0, step_width_times_BulgeSFR * (NormAGBMetalEjecRate_actual + (Bulge_total_metallicity * NormAGBMassEjecRate_actual))));}*/

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
    	//TotalMassReturnedToHotGas += 0.0;
    	//printf("BulgeToCold = %f\n\n", step_width_times_BulgeSFR * NormMassEjecRateSumAllTypes);
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

    	/*if (p==0) {printf("ICL:\n");}
    	if (p==0) {printf("%.11f | %.11f %.11f %.11f\n", NormMassEjecRateSumAllTypes, NormSNIIMassEjecRate_actual, NormSNIaMassEjecRate_actual, NormAGBMassEjecRate_actual);}
    	if (p==0) {printf("%.11f | %.11f %.11f %.11f | %.11f\n\n", max(0.0, step_width_times_ICMSFR * NormMassEjecRateSumAllTypes), max(0.0, step_width_times_ICMSFR * (NormSNIIMetalEjecRate_actual + (ICM_total_metallicity * NormSNIIMassEjecRate_actual))), max(0.0, step_width_times_ICMSFR * NormSNIaMetalEjecRate_actual), max(0.0, step_width_times_ICMSFR * (NormAGBMetalEjecRate_actual + (ICM_total_metallicity * NormAGBMassEjecRate_actual))), max(0.0, step_width_times_ICMSFR * NormSNIaMetalEjecRate_actual) + max(0.0, step_width_times_ICMSFR * (NormAGBMetalEjecRate_actual + (ICM_total_metallicity * NormAGBMassEjecRate_actual))));}*/

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

    	/*//Update ages:
        for(n=0;n<NOUT;n++)
        {
        	AgeCorrectionICM[n] += max(0.0, (sfh_time-NumToTime(ListOutputSnaps[n]))*(step_width_times_ICMSFR * NormMassEjecRateSumAllTypes));
        }*/
    } //if (Gal[p].sfh_ICM[i] > 0.0) //ICM

    } //for (i=0;i<=Gal[p].sfh_ibin;i++) //MAIN LOOP OVER SFH BINS

    /*//CALL SN-FEEDBACK RECIPE: Sending total mass returned to ColdGas to calculate FB energy:
    SN_feedback(p, centralgal, TotalMassReturnedToColdDiskGas);*/

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

#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN
    if(TotalMassReturnedToColdDiskGas>0.)
#ifndef H2_AND_RINGS
    	SN_feedback(p, centralgal, TotalMassReturnedToColdDiskGas, "ColdGas");
#else
    SN_feedback(p, centralgal, TotalMassReturnedToColdDiskGas, TotalMassReturnedToColdDiskGasr, "ColdGas", Coldmetallicityr);
#endif
    if(TotalMassReturnedToHotGas>0.)
#ifndef H2_AND_RINGS
    	SN_feedback(p, centralgal, TotalMassReturnedToHotGas, "HotGas");
#else
    SN_feedback(p, centralgal, TotalMassReturnedToHotGas, TotalMassReturnedToHotGasr, "HotGas", Hotmetallicity);
#endif
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

