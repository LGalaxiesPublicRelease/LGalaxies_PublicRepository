/*
 * model_yields.c
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

//INDIVIDUAL ELEMENTS
//All: [H][He][Cb][N][O][Ne][Mg][Si][S][Ca][Fe] or //Only [H][He][O][Mg][Fe]


void update_yields_and_return_mass(int p, int centralgal, double dt, int nstep)
{
	int Zi, igal, ii, mm;
	double timestep_width; //Width of current timestep in CODE UNITS
	int TimeBin; //Bin in Yield arrays corresponding to current timestep
	double Zi_disp;
	double sfh_time;
	//double time_to_ts; //Time from high-z (upper) edge of SFH bin to middle of current timestep (used for massive SNII to hot) [in Myrs]
	//double tcut; //Maximum lifetime of stars that have their ejected put straight into the HotGas [in Myrs]
#ifdef METALRICHWIND
#ifdef GASDENSITYFWIND
	double ColdGasSurfaceDensity;
#endif
#endif
	double fwind_SNII; //Required for SN-II metal-rich wind implementation
	double fwind_SNIa; //Required for SN-Ia metal-rich wind implementation
	double DiskSFRxStep, DiskSFRxStep_Phys, BulgeSFRxStep, BulgeSFRxStep_Phys, ICMSFRxStep, ICMSFRxStep_Phys;
	double DiskMetallicity, BulgeMetallicity, ICMMetallicity;
	double TotalMassReturnedToColdDiskGas, TotalMassReturnedToHotGas;
	double SNIIEjectaMass, SNIaEjectaMass, AGBEjectaMass;
	double SNIIUnProcessedMetals, SNIaUnProcessedMetals, AGBUnProcessedMetals, SNIIAllMetals, SNIaAllMetals, AGBAllMetals;
#ifdef INDIVIDUAL_ELEMENTS
	int ee;
	double SNIIAllElements[NUM_ELEMENTS], SNIaAllElements[NUM_ELEMENTS], AGBAllElements[NUM_ELEMENTS];
	double SNIIUnProcessedElements[NUM_ELEMENTS], SNIaUnProcessedElements[NUM_ELEMENTS], AGBUnProcessedElements[NUM_ELEMENTS];
	double DiskMetallicityElement_Phys[NUM_ELEMENTS];
	double BulgeMetallicityElement_Phys[NUM_ELEMENTS];
	double ICMMetallicityElement_Phys[NUM_ELEMENTS];
#endif
	int n; //Iterator used for loop over NOUT when updating MassWeightedAge
	double AgeCorrectionDisk[NOUT];
	double AgeCorrectionBulge[NOUT];
#ifdef H2_AND_RINGS
	int jj;
	double fractionRings[RNUM];
	double fractionRingsBulge[RNUM];
	double TotalMassReturnedToColdDiskGasr[RNUM];
#endif

	TotalMassReturnedToHotGas=0.0;
	TotalMassReturnedToColdDiskGas=0.0;
#ifdef H2_AND_RINGS
    	for(jj=0;jj<RNUM;jj++)
    	  TotalMassReturnedToColdDiskGasr[jj]=0.0;
#endif


    	for(n=0;n<NOUT;n++)
    	  {
    	    AgeCorrectionDisk[n] = 0.0;
    	    AgeCorrectionBulge[n] = 0.0;
    	  }

	timestep_width = dt; //Width of current timestep in CODE UNITS (units cancel out when dividing by SFH bin width, sfh_dt) (12-04-12)
	TimeBin = (STEPS*(Halo[Gal[p].HaloNr].SnapNum-1.0))+nstep; //TimeBin = (STEPS*Gal[p].SnapNum)+nstep; //Bin in Yield tables corresponding to current timestep //TEST!: BRUNO: Snapnum would be +1 too low for a 'jumping' galaxy (14-11-13)
	//timet = NumToTime((Halo[Gal[p].HaloNr].SnapNum-1.0)) - (nstep + 0.5) * dt; //Time from middle of the current timestep to z=0 (used here for MassWeightAge corrections)
	//NB: NumToTime(Gal[p].SnapNum) is the time to z=0 from start of current snapshot
	//    nstep is the number of the current timestep (0-19)
	//    dt is the width of one timestep within current snapshot


	//for stars dying that enrich the Hot gas directly
    if(Gal[p].Type==2)
      igal=Gal[p].CentralGal;
    else
      igal=p;

    for (ii=0;ii<=Gal[p].sfh_ibin;ii++) //LOOP OVER SFH BINS
    {
    	sfh_time=Gal[p].sfh_t[ii]+(0.5*Gal[p].sfh_dt[ii]);
    	//time_to_ts = ((sfh_time+(0.5*Gal[p].sfh_dt[i])) - timet)*(UnitTime_in_years/Hubble_h)/1.0e6; //Time from high-z (upper) edge of SFH bin to middle of current timestep [in Myrs]
    	//tcut = 2.0*((Gal[p].Rvir/Gal[p].Vvir)/0.0001); //Maximum lifetime of stars that have their ejected put straight into the HotGas [in Myrs]


    	mass_checks(p,"model_yields.c",__LINE__);


    //******************************************************
    //ENRICHMENT FROM DISK STARS (INTO COLD GAS & HOT GAS):
    //******************************************************
#ifdef H2_AND_RINGS
    	for(jj=0;jj<RNUM;jj++) 	
    		if (Gal[p].DiskMassRings[jj] > 0.0)//ROB: Discs can be destroyed (i.e. converted in to bulges). So only calculate enrichment from stars born in the disc if there is still a disc
    			if (Gal[p].sfh_DiskMassRings[jj][ii] > 0.0)    	
    			{		
#else
    if (Gal[p].DiskMass > 0.0)//ROB: Discs can be destroyed (i.e. converted in to bulges). So only calculate enrichment from stars born in the disc if there is still a disc at the current timestep. This if statement has scope until the end of the following if statement.
    if (Gal[p].sfh_DiskMass[ii] > 0.0)
    {
#endif
    

#ifdef METALRICHWIND
#ifdef GASDENSITYFWIND
#ifndef H2_AND_RINGS
    	ColdGasSurfaceDensity = (Gal[p].ColdGas*(1.0e10))/(4.0*3.14159265*Gal[p].ColdGasRadius*Gal[p].ColdGasRadius*(1.0e6*1.0e6)); //in Msun/pc^2
    	fwind_SNII = min(1.0, 1.0/(ColdGasSurfaceDensity/NORMGASDENSITY)); //Fraction of SN-II ejecta put directly into HotGas. When ColdGasSurfaceDensity <= NORMGASDENSITY (i.e. 10.0 Msun/pc), then fwind_SNII = 1.0.
    	fwind_SNIa = min(1.0, 1.0/(ColdGasSurfaceDensity/NORMGASDENSITY)); //Fraction of SN-Ia ejecta put directly into HotGas. When ColdGasSurfaceDensity <= NORMGASDENSITY (i.e. 10.0 Msun/pc), then fwind_SNII = 1.0.
#else
    	if(jj==0)
    		ColdGasSurfaceDensity = (Gal[p].ColdGasRings[jj]*(1.0e10))/(4.0*3.14159265*RingRadius[jj]*RingRadius[jj]*(1.0e6*1.0e6)); //in Msun/pc^2
    	else
    		ColdGasSurfaceDensity = (Gal[p].ColdGasRings[jj]*(1.0e10))/(4.0*3.14159265*((RingRadius[jj]*RingRadius[jj])-(RingRadius[jj-1]*RingRadius[jj-1]))*(1.0e6*1.0e6)); //in Msun/pc^2

    	fwind_SNII = min(1.0, 1.0/(ColdGasSurfaceDensity/NORMGASDENSITY)); //Fraction of SN-II ejecta put directly into HotGas.
    	fwind_SNIa = min(1.0, 1.0/(ColdGasSurfaceDensity/NORMGASDENSITY)); //Fraction of SN-Ia ejecta put directly into HotGas.

#endif	//H2_AND_RINGS
#endif	//GASDENSITYFWIND
#ifndef GASDENSITYFWIND
    fwind_SNII = FracZSNIItoHot; //FracZtoHot;
    fwind_SNIa = FracZSNIatoHot;
#endif //GASDENSITYFWIND

#else //METALRICHWIND
     	fwind_SNII = 0.0; //For all stellar ejecta (from disk) to ColdGas
     	fwind_SNIa = 0.0;
#endif //METALRICHWIND

    	//pre-calculations to speed up the code
	    //timestep_width and dt units cancel out
#ifdef H2_AND_RINGS
    	DiskSFRxStep = timestep_width * Gal[p].sfh_DiskMassRings[jj][ii]/Gal[p].sfh_dt[ii];
    	DiskSFRxStep_Phys = DiskSFRxStep * (1.0e10/Hubble_h);

    	DiskMetallicity = 0.;
    	for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
    		DiskMetallicity += Gal[p].sfh_MetalsDiskMassRings[jj][ii][mm];
    	DiskMetallicity /= Gal[p].sfh_DiskMassRings[jj][ii];

#ifdef INDIVIDUAL_ELEMENTS
    	for (ee=0;ee<NUM_ELEMENTS;ee++)
    		DiskMetallicityElement_Phys[ee] = Gal[p].sfh_DiskMass_elementsRings[jj][ii][ee] / (Gal[p].sfh_DiskMassRings[jj][ii]*1.0e10/Hubble_h);
#endif

#else
    	DiskSFRxStep = timestep_width * Gal[p].sfh_DiskMass[ii]/Gal[p].sfh_dt[ii];
    	DiskSFRxStep_Phys = DiskSFRxStep * (1.0e10/Hubble_h);

    	DiskMetallicity = 0.;
    	for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
    		DiskMetallicity += Gal[p].sfh_MetalsDiskMass[ii][mm];
    	DiskMetallicity /= Gal[p].sfh_DiskMass[ii];


#ifdef INDIVIDUAL_ELEMENTS
    	for (ee=0;ee<NUM_ELEMENTS;ee++)
    		DiskMetallicityElement_Phys[ee] = Gal[p].sfh_DiskMass_elements[ii][ee] / (Gal[p].sfh_DiskMass[ii]*1.0e10/Hubble_h);
#endif

#endif //H2_AND_RINGS

    	Zi = find_initial_metallicity(DiskMetallicity, 1, 1);
    	//Interpolate the disk metallicity on the lifetimeMetallicities tables:
    	Zi_disp = (DiskMetallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.

#ifndef INDIVIDUAL_ELEMENTS
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, DiskSFRxStep, DiskSFRxStep_Phys, DiskMetallicity,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals);
#else
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, DiskSFRxStep, DiskSFRxStep_Phys, DiskMetallicity, DiskMetallicityElement_Phys,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
				   SNIIAllElements, SNIIUnProcessedElements, SNIaAllElements, SNIaUnProcessedElements, AGBAllElements, AGBUnProcessedElements);
#endif



//If mass return or unprocessed metals larger than what is currently in the stellar population, re-adjust
    	double SNIImetals, SNIametals, AGBmetals;
#ifdef H2_AND_RINGS
		SNIImetals = Gal[p].MetalsDiskMassRings[jj][0];
		SNIametals = Gal[p].MetalsDiskMassRings[jj][1];
		AGBmetals = Gal[p].MetalsDiskMassRings[jj][2];
#else
		SNIImetals = Gal[p].MetalsDiskMass[0];
		SNIametals = Gal[p].MetalsDiskMass[1];
		AGBmetals = Gal[p].MetalsDiskMass[2];
#endif
    	if(SNIIUnProcessedMetals>SNIImetals)
    		SNIIUnProcessedMetals=SNIImetals;
    	if(SNIaUnProcessedMetals>SNIametals)
    		SNIaUnProcessedMetals=SNIametals;
    	if(AGBUnProcessedMetals>AGBmetals)
    		AGBUnProcessedMetals=AGBmetals;




    	mass_checks(p,"model_yields.c",__LINE__);


    	//************************
    	// UPDATE GAS COMPONENTS:
    	//************************

    	//Hot Gas:
     	Gal[igal].HotGas += (fwind_SNII * SNIIEjectaMass) + (fwind_SNIa * SNIaEjectaMass);
     	Gal[igal].MetalsHotGas[0] += fwind_SNII * SNIIAllMetals;
     	Gal[igal].MetalsHotGas[1] += fwind_SNIa * SNIaAllMetals;
     	mass_checks(p,"model_yields.c",__LINE__);
     	//Cold Gas:
     	Gal[p].ColdGas += (1.0-fwind_SNII)*SNIIEjectaMass + (1.0-fwind_SNIa)*SNIaEjectaMass + AGBEjectaMass;
     	Gal[p].MetalsColdGas[0] += (1.0-fwind_SNII) * SNIIAllMetals;
     	Gal[p].MetalsColdGas[1] += (1.0-fwind_SNIa) * SNIaAllMetals;
     	Gal[p].MetalsColdGas[2] += AGBAllMetals;

#ifdef H2_AND_RINGS
     	Gal[p].ColdGasRings[jj] += ((1.0-fwind_SNII) * SNIIEjectaMass + (1.0-fwind_SNIa) * SNIaEjectaMass + AGBEjectaMass);
     	Gal[p].MetalsColdGasRings[jj][0] += (1.0-fwind_SNII) * SNIIAllMetals;
     	Gal[p].MetalsColdGasRings[jj][1] += (1.0-fwind_SNIa) * SNIaAllMetals;
     	Gal[p].MetalsColdGasRings[jj][2] +=  AGBAllMetals;
#endif
    	//Total:
     	TotalMassReturnedToHotGas += (fwind_SNII * SNIIEjectaMass) + (fwind_SNIa * SNIaEjectaMass);
     	TotalMassReturnedToColdDiskGas += (1.0-fwind_SNII) * SNIIEjectaMass + (1.0-fwind_SNIa) * SNIaEjectaMass + AGBEjectaMass; //Only use energy from SNe that eject into ColdGas to reheat
#ifdef H2_AND_RINGS
     	TotalMassReturnedToColdDiskGasr[jj] += ((1.0-fwind_SNII) * SNIIEjectaMass + (1.0-fwind_SNIa) * SNIaEjectaMass + AGBEjectaMass);
#endif //H2_AND_RINGS



#ifdef INDIVIDUAL_ELEMENTS
     	for(ee=0;ee<NUM_ELEMENTS;ee++)
     	  {
     	    Gal[igal].HotGas_elements[ee] += fwind_SNII * SNIIAllElements[ee] + fwind_SNIa * SNIaAllElements[ee];
     	    Gal[p].ColdGas_elements[ee] += (1.0-fwind_SNII) * SNIIAllElements[ee] + (1.0-fwind_SNIa) * SNIaAllElements[ee] + AGBAllElements[ee];
#ifdef H2_AND_RINGS
     	    Gal[p].ColdGasRings_elements[jj][ee] += (1.0-fwind_SNII[jj]) * SNIIAllElements[ee] + (1.0-fwind_SNIa[jj]) * SNIaAllElements[ee] + AGBAllElements[ee];
#endif // H2_AND_RINGS
     	  }
#endif //INDIVIDUAL_ELEMENTS



     	mass_checks(p,"model_yields.c",__LINE__);


     	//*****************************
     	// UPDATE DISK MASS COMPONENTS:
     	//*****************************
     	/*ROB (13-02-13): All the mass/metals/elements in the stars that die in this timestep are lost from the stellar component.
     	 * i.e. All the mass/metals/elements in the stars at birth are removed...
     	 *...Some goes to the gas (+ newly synthesised component), the rest goes into the 'stellar remnants' which are not tracked and do not contribute to the stellar component's mass/metals/elements budget.*/
     	Gal[p].DiskMass -= SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
     	Gal[p].MetalsDiskMass[0] -= SNIIUnProcessedMetals;
     	Gal[p].MetalsDiskMass[1] -= SNIaUnProcessedMetals;
     	Gal[p].MetalsDiskMass[2] -= AGBUnProcessedMetals;

#ifdef H2_AND_RINGS
     	Gal[p].DiskMassRings[jj]-= (SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
     	Gal[p].MetalsDiskMassRings[jj][0]-= SNIIUnProcessedMetals;
     	Gal[p].MetalsDiskMassRings[jj][1]-= SNIaUnProcessedMetals;
     	Gal[p].MetalsDiskMassRings[jj][2]-= AGBUnProcessedMetals;
#endif
	    mass_checks(p,"model_yields.c",__LINE__);
#ifdef INDIVIDUAL_ELEMENTS

     	for(ee=0;ee<NUM_ELEMENTS;ee++)
     	  {
     	    Gal[p].DiskMass_elements[ee] -= SNIIUnProcessedElements[ee]+SNIaUnProcessedElements[ee]+AGBUnProcessedElements[ee];
#ifdef H2_AND_RINGS
     	    Gal[p].DiskMassRings_elements[jj][ee] -= (SNIIUnProcessedElements[ee]+SNIaUnProcessedElements[ee]+AGBUnProcessedElements[ee]);
#endif
     	  }
#endif //INDIVIDUAL_ELEMENTS

     	//Update ages:
     	for(n=0;n<NOUT;n++)
     	  {
    	    AgeCorrectionDisk[n] += (sfh_time-NumToTime(ListOutputSnaps[n]))*(SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
    	    if (AgeCorrectionDisk[n] < 0.0) AgeCorrectionDisk[n] = 0.0;
    	  }

    	mass_checks(p,"model_yields.c",__LINE__);
    } //if (Gal[p].sfh_DiskMass[i] > 0.0) -> all disk properties updated













    //******************************************************
    //ENRICHMENT FROM BULGE STARS (INTO COLD GAS & HOT GAS):
    //******************************************************
#ifdef H2_AND_RINGS
    	for(jj=0;jj<RNUM;jj++)
    		if (Gal[p].sfh_BulgeMassRings[jj][ii] > 0.0)
    		{
    			BulgeSFRxStep = timestep_width * Gal[p].sfh_BulgeMassRings[jj][ii]/Gal[p].sfh_dt[ii];
    			BulgeSFRxStep_Phys = BulgeSFRxStep * (1.0e10/Hubble_h);
    			//printf("sfr=%0.5e\n",DiskSFRxStep_Phys);
    			BulgeMetallicity = 0.;
    			for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
    				BulgeMetallicity += Gal[p].sfh_MetalsBulgeMassRings[jj][ii][mm];
    			BulgeMetallicity /= Gal[p].sfh_BulgeMassRings[jj][ii];

#ifdef INDIVIDUAL_ELEMENTS
    			for (ee=0;ee<NUM_ELEMENTS;ee++)
    				BulgeMetallicityElement_Phys[ee] = Gal[p].sfh_BulgeMass_elementsRings[jj][ii][ee] / (Gal[p].sfh_BulgeMassRings[jj][ii]*1.0e10/Hubble_h);
#endif

#else //H2_AND_RINGS

    	if (Gal[p].sfh_BulgeMass[ii] > 0.0)
    	  {
    		BulgeSFRxStep = timestep_width * Gal[p].sfh_BulgeMass[ii]/Gal[p].sfh_dt[ii];
    		BulgeSFRxStep_Phys = BulgeSFRxStep * (1.0e10/Hubble_h);

    		BulgeMetallicity = 0.;
    		for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
    			BulgeMetallicity += Gal[p].sfh_MetalsBulgeMass[ii][mm];
    		BulgeMetallicity /= Gal[p].sfh_BulgeMass[ii];


#ifdef INDIVIDUAL_ELEMENTS
    		for (ee=0;ee<NUM_ELEMENTS;ee++)
    			BulgeMetallicityElement_Phys[ee] = Gal[p].sfh_BulgeMass_elements[ii][ee] / (Gal[p].sfh_BulgeMass[ii]*1.0e10/Hubble_h);
#endif

#endif //H2_AND_RINGS

    	Zi = find_initial_metallicity(BulgeMetallicity, 1, 2);
    	//Interpolate the bulge luminosity on the lifetimeMetallicities tables:
    	Zi_disp = (BulgeMetallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.


#ifndef INDIVIDUAL_ELEMENTS
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, BulgeSFRxStep, BulgeSFRxStep_Phys, BulgeMetallicity,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals);

#else
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, BulgeSFRxStep, BulgeSFRxStep_Phys, BulgeMetallicity, BulgeMetallicityElement_Phys,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
				   SNIIAllElements, SNIIUnProcessedElements, SNIaAllElements, SNIaUnProcessedElements, AGBAllElements, AGBUnProcessedElements);
#endif

    	//If mass return or unprocessed metals larger than what is currently in the stellar population, re-adjust
    	double SNIImetals, SNIametals, AGBmetals;

#ifdef H2_AND_RINGS
    	SNIImetals = Gal[p].MetalsBulgeMassRings[jj][0];
    	SNIametals = Gal[p].MetalsBulgeMassRings[jj][1];
    	AGBmetals = Gal[p].MetalsBulgeMassRings[jj][2];
#else
    	SNIImetals = Gal[p].MetalsBulgeMass[0];
    	SNIametals = Gal[p].MetalsBulgeMass[1];
    	AGBmetals = Gal[p].MetalsBulgeMass[2];
#endif
    	if(SNIIUnProcessedMetals>SNIImetals)
    		SNIIUnProcessedMetals=SNIImetals;
    	if(SNIaUnProcessedMetals>SNIametals)
    		SNIaUnProcessedMetals=SNIametals;
    	if(AGBUnProcessedMetals>AGBmetals)
    		AGBUnProcessedMetals=AGBmetals;



    	//*****************************
    	// UPDATE HOT GAS COMPONENTS:
    	//*****************************
#ifndef BULGE_TO_COLD
    	Gal[igal].HotGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	//If there was no hotgas left in the galaxy it will probably be stripped next step.
    	//Give it a fake HotRadius for now to avoid crash at mass checks
    	if(Gal[igal].HotGas>0. && Gal[igal].HotRadius==0)
    	  Gal[igal].HotRadius=1.e-10;
    	TotalMassReturnedToHotGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	Gal[igal].MetalsHotGas[0] += SNIIAllMetals;
    	Gal[igal].MetalsHotGas[1] += SNIaAllMetals;
    	Gal[igal].MetalsHotGas[2] += AGBAllMetals;

#ifdef INDIVIDUAL_ELEMENTS
    	for(ee=0;ee<NUM_ELEMENTS;ee++)
    	  Gal[igal].HotGas_elements[ee] += SNIIAllElements[ee]+SNIaAllElements[ee]+AGBAllElements[ee];
#endif //INDIVIDUAL_ELEMENTS



#else //BULGE_TO_COLD
    	Gal[p].ColdGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	Gal[p].MetalsColdGas[0] += SNIIAllMetals;
    	Gal[p].MetalsColdGas[1] += SNIaAllMetals;
    	Gal[p].MetalsColdGas[2] += AGBAllMetals;

#ifdef H2_AND_RINGS
    	Gal[p].ColdGasRings[jj] += (SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
    	Gal[p].MetalsColdGasRings[jj][0] += SNIIAllMetals;
    	Gal[p].MetalsColdGasRings[jj][1] += SNIaAllMetals;
    	Gal[p].MetalsColdGasRings[jj][2] += AGBAllMetals;
#endif

     	TotalMassReturnedToColdDiskGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
#ifdef H2_AND_RINGS
     	TotalMassReturnedToColdDiskGasr[jj]+= (SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
#endif

#ifdef INDIVIDUAL_ELEMENTS
     	for(ee=0;ee<NUM_ELEMENTS;ee++)
     	  {
     	    Gal[p].ColdGas_elements[ee] += SNIIAllElements[ee]+SNIaAllElements[ee]+AGBAllElements[ee];
#ifdef H2_AND_RINGS
     	    Gal[p].ColdGasRings_elements[jj][ee]  += (SNIIAllElements[ee] + SNIaAllElements[ee] + AGBAllElements[ee]);
#endif//H2_AND_RINGS
     	  }
#endif //INDIVIDUAL_ELEMENTS

#endif //BULGE_TO_COLD







    	//*****************************
    	// UPDATE BULGE MASS COMPONENTS:
    	//*****************************
    	Gal[p].BulgeMass -= SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	Gal[p].MetalsBulgeMass[0] -= SNIIUnProcessedMetals;
    	Gal[p].MetalsBulgeMass[1] -= SNIaUnProcessedMetals;
    	Gal[p].MetalsBulgeMass[2] -= AGBUnProcessedMetals;

#ifdef H2_AND_RINGS
    	Gal[p].BulgeMassRings[jj]-= (SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
    	Gal[p].MetalsBulgeMassRings[jj][0]-= SNIIUnProcessedMetals;
    	Gal[p].MetalsBulgeMassRings[jj][1]-= SNIaUnProcessedMetals;
    	Gal[p].MetalsBulgeMassRings[jj][2]-= AGBUnProcessedMetals;
#endif

#ifdef INDIVIDUAL_ELEMENTS
    	for(ee=0;ee<NUM_ELEMENTS;ee++)
    	  {
    	    Gal[p].BulgeMass_elements[ee] -= SNIIUnProcessedElements[ee]+SNIaUnProcessedElements[ee]+AGBUnProcessedElements[ee];
#ifdef H2_AND_RINGS
    	    Gal[p].BulgeMassRings_elements[jj][ee] -= (SNIIUnProcessedElements[ee]+SNIaUnProcessedElements[ee]+AGBUnProcessedElements[ee]);
#endif
    	  }
#endif //INDIVIDUAL_ELEMENTS

    	//Update ages:
        for(n=0;n<NOUT;n++)
        {
        	AgeCorrectionBulge[n] += (sfh_time-NumToTime(ListOutputSnaps[n]))*(SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
        	if (AgeCorrectionBulge[n] < 0.0) AgeCorrectionBulge[n] = 0.0;
        }

        mass_checks(p,"model_yields.c",__LINE__);
    } //if (Gal[p].sfh_BulgeMass[i] > 0.0) //all BULGE properties updated









    //*****************************************
    //ENRICHMENT FROM ICL STARS INTO HOT GAS:
    //*****************************************

    if (Gal[p].sfh_ICM[ii] > 0.0)
      {
    	//pre-calculations to speed up the code
    	//Note: This is NOT really an SFR, as no stars are formed in the ICM. Rather, this is a star-transfer rate from satellite disruption to the stellar halo.
    	ICMSFRxStep = timestep_width * Gal[p].sfh_ICM[ii]/Gal[p].sfh_dt[ii];
    	ICMSFRxStep_Phys = ICMSFRxStep * (1.0e10/Hubble_h) ;
    	ICMMetallicity=0.;
    	for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
    	  ICMMetallicity += Gal[p].sfh_MetalsICM[ii][mm];
    	ICMMetallicity /= Gal[p].sfh_ICM[ii];
#ifdef INDIVIDUAL_ELEMENTS
    	for (ee=0;ee<NUM_ELEMENTS;ee++)
    	  ICMMetallicityElement_Phys[ee] = Gal[p].sfh_ICM_elements[ii][ee] / (Gal[p].sfh_ICM[ii]*1.0e10/Hubble_h);
#endif

    	Zi = find_initial_metallicity(ICMMetallicity, 1, 3);
    	//Interpolate the ICM metallicity on the lifetimeMetallicities tables:
    	Zi_disp = (ICMMetallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.


#ifndef INDIVIDUAL_ELEMENTS
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, ICMSFRxStep, ICMSFRxStep_Phys, ICMMetallicity,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals);

#else
    	compute_actual_eject_rates(TimeBin, ii, Zi, Zi_disp, Gal[p].sfh_ibin, ICMSFRxStep, ICMSFRxStep_Phys, ICMMetallicity, ICMMetallicityElement_Phys,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
				   SNIIAllElements, SNIIUnProcessedElements, SNIaAllElements, SNIaUnProcessedElements, AGBAllElements, AGBUnProcessedElements);
#endif

    	//If mass return or unprocessed metals larger than what is currently in the stellar population, re-adjust
    	double SNIImetals, SNIametals, AGBmetals;

    	SNIImetals = Gal[p].MetalsICM[0];
    	SNIametals = Gal[p].MetalsICM[1];
    	AGBmetals = Gal[p].MetalsICM[2];

    	if(SNIIUnProcessedMetals>SNIImetals)
    		SNIIUnProcessedMetals=SNIImetals;
    	if(SNIaUnProcessedMetals>SNIametals)
    		SNIaUnProcessedMetals=SNIametals;
    	if(AGBUnProcessedMetals>AGBmetals)
    		AGBUnProcessedMetals=AGBmetals;

    	//*****************************
    	// UPDATE HOT GAS COMPONENTS:
    	//*****************************
    	Gal[igal].HotGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	//If there was no hotgas left in the galaxy it will probably be stripped next step.
    	//Give it a fake HotRadius for now to avoid crash at mass checks
    	if(Gal[igal].HotGas>0. && Gal[igal].HotRadius==0)
    	  Gal[igal].HotRadius=1.e-10;

    	TotalMassReturnedToHotGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	Gal[igal].MetalsHotGas[0] += SNIIAllMetals;
    	Gal[igal].MetalsHotGas[1] += SNIaAllMetals;
    	Gal[igal].MetalsHotGas[2] += AGBAllMetals;

#ifdef INDIVIDUAL_ELEMENTS
    	for(ee=0;ee<NUM_ELEMENTS;ee++)
    	  Gal[igal].HotGas_elements[ee] += SNIIAllElements[ee] + SNIaAllElements[ee] + AGBAllElements[ee];
#endif //INDIVIDUAL_ELEMENTS


    	//*****************************
    	// UPDATE ICM COMPONENTS:
    	//*****************************
    	Gal[p].ICM -= SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	Gal[p].MetalsICM[0] -= SNIIUnProcessedMetals;
    	Gal[p].MetalsICM[1] -= SNIaUnProcessedMetals;
    	Gal[p].MetalsICM[2] -= AGBUnProcessedMetals;

#ifdef INDIVIDUAL_ELEMENTS
    	for(ee=0;ee<NUM_ELEMENTS;ee++)
    	  Gal[p].ICM_elements[ee] -= SNIIUnProcessedElements[ee] + SNIaUnProcessedElements[ee] + AGBUnProcessedElements[ee];
#endif //INDIVIDUAL_ELEMENTS

    	//Update ages:
        //for(n=0;n<NOUT;n++)
        //{
        //	AgeCorrectionICM[n] += (sfh_time-NumToTime(ListOutputSnaps[n]))*(ICMSFRxStep * NormMassEjecRateAllTypes);
        //}

    	mass_checks(p,"model_yields.c",__LINE__);
      } //if (Gal[p].sfh_ICM[i] > 0.0) //all ICM properties updated



    } //for (i=0;i<=Gal[p].sfh_ibin;i++) //MAIN LOOP OVER SFH BINS








    //Update Mass-weighted ages:
    for(n=0;n<NOUT;n++)
      {
    	Gal[p].MassWeightAge[n] -= (AgeCorrectionDisk[n]+AgeCorrectionBulge[n]);
      }



#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN

    if(TotalMassReturnedToColdDiskGas>0.)
      {
	if(TotalMassReturnedToColdDiskGas>Gal[p].ColdGas)
	  TotalMassReturnedToColdDiskGas=Gal[p].ColdGas;
#ifndef H2_AND_RINGS
	SN_feedback(p, centralgal, TotalMassReturnedToColdDiskGas, "ColdGas");
#else
	SN_feedback(p, centralgal, TotalMassReturnedToColdDiskGas, TotalMassReturnedToColdDiskGasr, "ColdGas");
#endif
      }

    //this mass will only result in ejection, no reheating
    if(TotalMassReturnedToHotGas>0.)
      {
    	if(TotalMassReturnedToHotGas>Gal[p].HotGas)
    		TotalMassReturnedToHotGas=Gal[p].HotGas;
#ifndef H2_AND_RINGS
    	SN_feedback(p, centralgal, TotalMassReturnedToHotGas, "HotGas");
#else
    	double HotGasRings[RNUM];
		for(jj=0;jj<RNUM;jj++)
			//HotGasRings[jj]=TotalMassReturnedToHotGas*fractionRingsBulge[jj];
			HotGasRings[jj]=0.;
		SN_feedback(p, centralgal, TotalMassReturnedToHotGas, HotGasRings, "HotGas");
#endif
      }

#endif //FEEDBACK_COUPLED_WITH_MASS_RETURN

      mass_checks(p,"model_yields.c",__LINE__);

      /*double diskspinpar=sqrt(Gal[p].DiskSpin[0] * Gal[p].DiskSpin[0] +
     			     Gal[p].DiskSpin[1] * Gal[p].DiskSpin[1] +
     			     Gal[p].DiskSpin[2] * Gal[p].DiskSpin[2] );
      int ii;
      if (Gal[p].ColdGas+TotalMassReturnedToColdDiskGas > 1.e-8)
        for(ii=0;ii<3;ii++)
          Gal[p].ColdGasSpin[ii]=((Gal[p].ColdGasSpin[ii])*(Gal[p].ColdGas)+1./sqrt(3)*diskspinpar*TotalMassReturnedToColdDiskGas)/(Gal[p].ColdGas+TotalMassReturnedToColdDiskGas);*/

      if (Gal[p].ColdGas+TotalMassReturnedToColdDiskGas > 1.e-8)
	for (ii = 0; ii < 3; ii++)
	  Gal[p].ColdGasSpin[ii]=((Gal[p].ColdGasSpin[ii])*(Gal[p].ColdGas)+TotalMassReturnedToColdDiskGas*Gal[p].DiskSpin[ii])/(Gal[p].ColdGas+TotalMassReturnedToColdDiskGas);

      if (DiskRadiusModel == 0)
	{
	  Gal[p].ColdGasRadius = get_gas_disk_radius(p);
	  Gal[p].DiskRadius = get_stellar_disk_radius(p);
	}
}


#ifndef INDIVIDUAL_ELEMENTS
void compute_actual_eject_rates(int TimeBin, int ii, int Zi, double Zi_disp, int sfh_ibin, double SFRxStep, double SFRxStep_Phys, double Metallicity,
				 double *SNIIEjectaMass, double *SNIIAllMetals, double *SNIIUnProcessedMetals,
				 double *SNIaEjectaMass, double *SNIaAllMetals, double *SNIaUnProcessedMetals,
				 double *AGBEjectaMass, double *AGBAllMetals, double *AGBUnProcessedMetals)
#else
void compute_actual_eject_rates(int TimeBin, int ii, int Zi, double Zi_disp, int sfh_ibin, double SFRxStep, double SFRxStep_Phys, double Metallicity, double *MetallicityElement_Phys,
				 double *SNIIEjectaMass, double *SNIIAllMetals, double *SNIIUnProcessedMetals,
				 double *SNIaEjectaMass, double *SNIaAllMetals, double *SNIaUnProcessedMetals,
				 double *AGBEjectaMass, double *AGBAllMetals, double *AGBUnProcessedMetals,
				 double *SNIIAllElements, double *SNIIUnProcessedElements,
				 double *SNIaAllElements, double *SNIaUnProcessedElements,
				 double *AGBAllElements, double *AGBUnProcessedElements)
#endif
{
  double NormSNIIMassEjecRate_actual, NormSNIaMassEjecRate_actual, NormAGBMassEjecRate_actual;
  double NormSNIIMetalEjecRate_actual, NormSNIaMetalEjecRate_actual, NormAGBMetalEjecRate_actual;
#ifdef INDIVIDUAL_ELEMENTS
  int ee;
  double NormSNIIYieldRate_actual[NUM_ELEMENTS], NormSNIaYieldRate_actual[NUM_ELEMENTS], NormAGBYieldRate_actual[NUM_ELEMENTS];
#endif

//pre-calculations to speed up the code
  NormSNIIMassEjecRate_actual = NormSNIIMassEjecRate[TimeBin][ii][Zi] + ((NormSNIIMassEjecRate[TimeBin][ii][Zi+1] - NormSNIIMassEjecRate[TimeBin][ii][Zi])*Zi_disp);
  NormSNIaMassEjecRate_actual = NormSNIaMassEjecRate[TimeBin][ii][Zi] + ((NormSNIaMassEjecRate[TimeBin][ii][Zi+1] - NormSNIaMassEjecRate[TimeBin][ii][Zi])*Zi_disp);
  NormAGBMassEjecRate_actual = NormAGBMassEjecRate[TimeBin][ii][Zi] + ((NormAGBMassEjecRate[TimeBin][ii][Zi+1] - NormAGBMassEjecRate[TimeBin][ii][Zi])*Zi_disp);
  NormSNIIMetalEjecRate_actual = NormSNIIMetalEjecRate[TimeBin][ii][Zi] + ((NormSNIIMetalEjecRate[TimeBin][ii][Zi+1] - NormSNIIMetalEjecRate[TimeBin][ii][Zi])*Zi_disp);
  NormSNIaMetalEjecRate_actual = NormSNIaMetalEjecRate[TimeBin][ii][Zi] + ((NormSNIaMetalEjecRate[TimeBin][ii][Zi+1] - NormSNIaMetalEjecRate[TimeBin][ii][Zi])*Zi_disp);
  NormAGBMetalEjecRate_actual = NormAGBMetalEjecRate[TimeBin][ii][Zi] + ((NormAGBMetalEjecRate[TimeBin][ii][Zi+1] - NormAGBMetalEjecRate[TimeBin][ii][Zi])*Zi_disp);

#ifdef INDIVIDUAL_ELEMENTS //Work out the actual yield of element k, by interpolating between the yields in the look-up table created by yield_integrals.c.
  for (ee=0;ee<NUM_ELEMENTS;ee++)
    {
      NormSNIIYieldRate_actual[ee] = NormSNIIYieldRate[TimeBin][ii][Zi][ee] + ((NormSNIIYieldRate[TimeBin][ii][Zi+1][ee] - NormSNIIYieldRate[TimeBin][ii][Zi][ee])*Zi_disp);
      NormSNIaYieldRate_actual[ee] = NormSNIaYieldRate[TimeBin][ii][Zi][ee] + ((NormSNIaYieldRate[TimeBin][ii][Zi+1][ee] - NormSNIaYieldRate[TimeBin][ii][Zi][ee])*Zi_disp);
      NormAGBYieldRate_actual[ee] = NormAGBYieldRate[TimeBin][ii][Zi][ee] + ((NormAGBYieldRate[TimeBin][ii][Zi+1][ee] - NormAGBYieldRate[TimeBin][ii][Zi][ee])*Zi_disp);
    }
#endif

#ifdef FAST_TESTING_MODE
    	/* if(NormAGBMassEjecRate_actual>0.1)
    	          NormAGBMassEjecRate_actual=0.01;

    	        if(NormSNIIMassEjecRate_actual+NormSNIaMassEjecRate_actual+NormAGBMassEjecRate_actual>1.)
    	          {
    	            NormSNIIMassEjecRate_actual=0.43;
    	            NormSNIaMassEjecRate_actual=0.3;
    	            NormAGBMassEjecRate_actual=0.1;
    	          }*/
#endif

#ifdef INSTANTANEOUS_RECYCLE //to recover results from instantaneous recycling approximation
  reset_ejection_rates(ii, sfh_ibin,
		       &NormSNIIMassEjecRate_actual, &NormSNIIMetalEjecRate_actual,
		       &NormSNIaMassEjecRate_actual, &NormAGBMassEjecRate_actual,
		       &NormSNIaMetalEjecRate_actual, &NormAGBMetalEjecRate_actual);
#endif //INSTANTANEOUS_RECYCLE


  //SNII
  *SNIIEjectaMass = SFRxStep * NormSNIIMassEjecRate_actual;
#ifdef INSTANTANEOUS_RECYCLE
  *SNIIAllMetals = SFRxStep * (NormSNIIMetalEjecRate_actual + Metallicity * NormSNIIMassEjecRate_actual);
#else
#ifdef PORTINARI
  *SNIIAllMetals = SFRxStep * (NormSNIIMetalEjecRate_actual + Metallicity * NormSNIIMassEjecRate_actual);
#endif
#ifdef CHIEFFI
  //ROB: No unsynth component required for SN-II ejecta, when using the Chieffi & Limongi 92007) yield tables/
  *SNIIAllMetals = SFRxStep * NormSNIIMetalEjecRate_actual;
#endif
#endif //INSTANTANEOUS_RECYCLE
  *SNIIUnProcessedMetals = SFRxStep * Metallicity * NormSNIIMassEjecRate_actual;

  //SNIa
  *SNIaEjectaMass = SFRxStep * NormSNIaMassEjecRate_actual;
  //SNIa yields are written in a different way to SNIIa and AGB so the following line is correct
  *SNIaAllMetals = SFRxStep * NormSNIaMetalEjecRate_actual;
  *SNIaUnProcessedMetals = SFRxStep * Metallicity * NormSNIaMassEjecRate_actual;

  //AGB
  *AGBEjectaMass = SFRxStep * NormAGBMassEjecRate_actual;
  *AGBAllMetals = SFRxStep * (NormAGBMetalEjecRate_actual + (Metallicity * NormAGBMassEjecRate_actual));
  *AGBUnProcessedMetals = SFRxStep * Metallicity * NormAGBMassEjecRate_actual;

#ifdef INDIVIDUAL_ELEMENTS
  for(ee=0;ee<NUM_ELEMENTS;ee++)
    {
#ifdef PORTINARI
      SNIIAllElements[ee] = SFRxStep_Phys * (NormSNIIYieldRate_actual[ee] + MetallicityElement_Phys[ee] * NormSNIIMassEjecRate_actual);
#endif
#ifdef CHIEFFI
      SNIIAllElements[ee] = SFRxStep_Phys * NormSNIIYieldRate_actual[ee];
#endif
      SNIIUnProcessedElements[ee] = SFRxStep_Phys * MetallicityElement_Phys[ee] * NormSNIIMassEjecRate_actual;

      SNIaAllElements[ee] = SFRxStep_Phys * NormSNIaYieldRate_actual[ee];
      SNIaUnProcessedElements[ee] = SFRxStep_Phys * MetallicityElement_Phys[ee] * NormSNIaMassEjecRate_actual;

      AGBAllElements[ee] = SFRxStep_Phys * (NormAGBYieldRate_actual[ee] + MetallicityElement_Phys[ee] * NormAGBMassEjecRate_actual);
      AGBUnProcessedElements[ee] = SFRxStep_Phys * MetallicityElement_Phys[ee] * NormAGBMassEjecRate_actual;
    }
#endif

}








int find_initial_metallicity(double metallicity, int table_type, int component_type)
{
	int i, Zi_bin;
	Zi_bin = -1;
	i = 0;

	switch (table_type)
	{
	case 1: //Lifetime metallicity table
	  while (Zi_bin == -1)
		{
		  if (lifetimeMetallicities[i] < metallicity)
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
		  if (SNIIMetallicities[i] < metallicity)
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
		  if (AGBMetallicities[i] < metallicity)
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


#ifdef INSTANTANEOUS_RECYCLE //to recover results from instantaneous recycling approximation
 void reset_ejection_rates(int ii, int sfh_ibin,
		 double *NormSNIIMassEjecRate_actual, double *NormSNIIMetalEjecRate_actual,
		 double *NormSNIaMassEjecRate_actual, double *NormAGBMassEjecRate_actual,
		 double *NormSNIaMetalEjecRate_actual, double *NormAGBMetalEjecRate_actual)
 {
    	if(ii==sfh_ibin)
    	{
    	    *NormSNIIMassEjecRate_actual = RecycleFraction;
    	    *NormSNIIMetalEjecRate_actual = Yield;
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

