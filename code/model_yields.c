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
	int Zi, igal, ii;
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
	double fwind; //Required for metal-rich wind implementation
	double DiskSFRxStep, DiskSFRxStep_Phys, BulgeSFRxStep, BulgeSFRxStep_Phys, ICMSFRxStep, ICMSFRxStep_Phys;
	double DiskMetallicity, BulgeMetallicity, ICMMetallicity;
	double TotalMassReturnedToColdDiskGas, TotalMassReturnedToHotGas;
#ifdef H2_AND_RINGS
	double TotalMassReturnedToColdDiskGasr[RNUM];
#endif
	double SNIIEjectaMass, SNIaEjectaMass, AGBEjectaMass;
	double SNIIUnProcessedMetals, SNIaUnProcessedMetals, AGBUnProcessedMetals, SNIIAllMetals, SNIaAllMetals, AGBAllMetals;
#ifdef INDIVIDUAL_ELEMENTS
	double SNIIAllElements[NUM_ELEMENTS], SNIaAllElements[NUM_ELEMENTS], AGBAllElements[NUM_ELEMENTS];
	double SNIIUnProcessedElements[NUM_ELEMENTS], SNIaUnProcessedElements[NUM_ELEMENTS], AGBUnProcessedElements[NUM_ELEMENTS];
	int kk;
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
#ifdef METALRICHWIND
#ifdef GASDENSITYFWIND
	ColdGasSurfaceDensity = (Gal[p].ColdGas*(1.0e10/Hubble_h))/(4.0*3.14159265*Gal[p].ColdGasRadius*Gal[p].ColdGasRadius/Hubble_h);
	fwind = min(1.0, 1.0/(ColdGasSurfaceDensity/5.0e12)); //Fraction of SN-II ejecta put directly into HotGas
	if (Gal[p].ColdGas != (float)Gal[p].ColdGas) {fwind = 1.0;}
#endif
#ifndef GASDENSITYFWIND
	fwind = FracZtoHot;
#endif
#endif
#ifndef METALRICHWIND
	fwind = 0.0; //For all stellar ejecta (from disk) to ColdGas
#endif

	//for stars dying that enrich the Hot gas directly
    if(Gal[p].Type==2)
      igal=Gal[p].CentralGal;
    else
      igal=p;

    int i;
    for (i=0;i<=Gal[p].sfh_ibin;i++) //LOOP OVER SFH BINS
    {
    	sfh_time=Gal[p].sfh_t[i]+(0.5*Gal[p].sfh_dt[i]);
    	//time_to_ts = ((sfh_time+(0.5*Gal[p].sfh_dt[i])) - timet)*(UnitTime_in_years/Hubble_h)/1.0e6; //Time from high-z (upper) edge of SFH bin to middle of current timestep [in Myrs]
    	//tcut = 2.0*((Gal[p].Rvir/Gal[p].Vvir)/0.0001); //Maximum lifetime of stars that have their ejected put straight into the HotGas [in Myrs]


    	mass_checks(p,"model_yields.c",__LINE__);
    	

    //******************************************************
    //ENRICHMENT FROM DISK STARS (INTO COLD GAS & HOT GAS):
    //******************************************************
    if (Gal[p].DiskMass > 0.0)
    if (Gal[p].sfh_DiskMass[i] > 0.0)
      {
	//pre-calculations to speed up the code
	//timestep_width and dt units cancel out
    	DiskSFRxStep = timestep_width * Gal[p].sfh_DiskMass[i]/Gal[p].sfh_dt[i];
    	DiskSFRxStep_Phys = DiskSFRxStep * (1.0e10/Hubble_h);
    	DiskMetallicity = 0.;
    	for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
    	  DiskMetallicity += Gal[p].sfh_MetalsDiskMass[i][ii];
    	DiskMetallicity /= Gal[p].sfh_DiskMass[i];

#ifdef INDIVIDUAL_ELEMENTS
    	for (kk=0;kk<NUM_ELEMENTS;kk++)
    	  DiskMetallicityElement_Phys[kk] = Gal[p].sfh_DiskMass_elements[i][kk] / (Gal[p].sfh_DiskMass[i]*1.0e10/Hubble_h);
#endif

    	Zi = find_initial_metallicity(p, i, 1, 1);
    	//Interpolate the disk metallicity on the lifetimeMetallicities tables:
    	Zi_disp = (DiskMetallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.


#ifndef INDIVIDUAL_ELEMENTS
    	compute_actual_eject_rates(TimeBin, i, Zi, Zi_disp, Gal[p].sfh_ibin, DiskSFRxStep, DiskSFRxStep_Phys, DiskMetallicity,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals);
#else
    	compute_actual_eject_rates(TimeBin, i, Zi, Zi_disp, Gal[p].sfh_ibin, DiskSFRxStep, DiskSFRxStep_Phys, DiskMetallicity, DiskMetallicityElement_Phys,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
				   SNIIAllElements, SNIIUnProcessedElements, SNIaAllElements, SNIaUnProcessedElements, AGBAllElements, AGBUnProcessedElements);
#endif


#ifdef H2_AND_RINGS
    	for(jj=0;jj<RNUM;jj++)
    	  fractionRings[jj]=Gal[p].sfh_DiskMassRings[jj][i]/Gal[p].sfh_DiskMass[i];
    	//fractionRings[jj]=Gal[p].DiskMassRings[jj]/Gal[p].DiskMass;
#endif






    	//************************
    	// UPDATE GAS COMPONENTS:
    	//************************

#ifndef SNIATOHOT
     	//Hot Gas
     	Gal[igal].HotGas += fwind * SNIIEjectaMass; //i.e. *only* SN-II ejecta could make it to the HotGas (if METALRICHWINDS is on). SN-Ia (and AGB) ejecta go to the ColdGas.
     	//If there was no hotgas left in the galaxy it will probably be stripped next step. Give it a fake HotRadius for now to avoid crash at mass checks
     	if(Gal[igal].HotGas>0. && Gal[igal].HotRadius==0)
     	  Gal[igal].HotRadius=1.e-10;
     	Gal[igal].MetalsHotGas[0] += fwind * SNIIAllMetals;

     	//Cold Gas
     	//Gal[p].ColdGas += (DiskSFRxStep * NormMassEjecRateAllTypes)-fwind * SNIIEjecta;
     	Gal[p].ColdGas += (1.0-fwind) * SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
     	Gal[p].MetalsColdGas[0] += (1.0-fwind) * SNIIAllMetals;
     	Gal[p].MetalsColdGas[1] += SNIaAllMetals;
     	Gal[p].MetalsColdGas[2] += AGBAllMetals;
#ifdef H2_AND_RINGS
     	for(jj=0;jj<RNUM;jj++)
     	  {
     	    //Gal[p].ColdGasRings[jj] += fractionRings[jj] * ((DiskSFRxStep * NormMassEjecRateAllTypes)-fwind * SNIIEjecta);
     	    Gal[p].ColdGasRings[jj] += fractionRings[jj] * ((1.0-fwind) * SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
     	    Gal[p].MetalsColdGasRings[jj][0] += fractionRings[jj] * (1.0-fwind) * SNIIAllMetals;
     	    Gal[p].MetalsColdGasRings[jj][1] += fractionRings[jj] * SNIaAllMetals;
     	    Gal[p].MetalsColdGasRings[jj][2] += fractionRings[jj] * AGBAllMetals;
     	  }
#endif

     	//Total
     	TotalMassReturnedToHotGas += fwind * SNIIEjectaMass;
     	//TotalMassReturnedToColdDiskGas += (DiskSFRxStep * NormMassEjecRateAllTypes)-fwind * SNIIEjecta; //Only use energy from SNe that eject into ColdGas to reheat
     	TotalMassReturnedToColdDiskGas += (1.0-fwind) * SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
#ifdef H2_AND_RINGS
     	for(jj=0;jj<RNUM;jj++)
     	  TotalMassReturnedToColdDiskGasr[jj]+=fractionRings[jj] * ((1.0-fwind) * SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
#endif


#else //ifdef SNIATOHOT

     	//Hot Gas
     	Gal[igal].HotGas += fwind * SNIIEjectaMass + SNIaEjectaMass;
     	//If there was no hotgas left in the galaxy it will probably be stripped next step. Give it a fake HotRadius for now to avoid crash at mass checks
     	//if(Gal[igal].HotGas>0. && Gal[igal].HotRadius==0)
     	//  Gal[igal].HotRadius=1.e-10;
     	Gal[igal].MetalsHotGas[0] += fwind * SNIIAllMetals;
     	Gal[igal].MetalsHotGas[1] += SNIaAllMetals;

     	//Cold Gas
     	Gal[p].ColdGas += (1.0-fwind) * SNIIEjectaMass + AGBEjectaMass;
     	Gal[p].MetalsColdGas[0] += (1.0-fwind) * SNIIAllMetals;
     	Gal[p].MetalsColdGas[2] += AGBAllMetals;
#ifdef H2_AND_RINGS
     	for(jj=0;jj<RNUM;jj++)
     	  {
     	    Gal[p].ColdGasRings[jj] += fractionRings[jj] * ((1.0-fwind) * SNIIEjectaMass + AGBEjectaMass);
     	    Gal[p].MetalsColdGasRings[jj][0] += fractionRings[jj] * (1.0-fwind) * SNIIAllMetals;
     	    Gal[p].MetalsColdGasRings[jj][2] += fractionRings[jj] * AGBAllMetals;
     	  }
#endif

     	//Total
     	TotalMassReturnedToHotGas += fwind * SNIIEjectaMass + SNIaEjectaMass;
     	TotalMassReturnedToColdDiskGas += (1.0-fwind) * SNIIEjectaMass + AGBEjectaMass; //Only use energy from SNe that eject into ColdGas to reheat
#ifdef H2_AND_RINGS
     	for(jj=0;jj<RNUM;jj++)
     	  TotalMassReturnedToColdDiskGasr[jj]+=fractionRings[jj] * ((1.0-fwind) * SNIIEjectaMass + AGBEjectaMass);
#endif


#endif
     	


#ifdef INDIVIDUAL_ELEMENTS



#ifndef SNIATOHOT
     	for(kk=0;kk<NUM_ELEMENTS;kk++)
     	  {
     	    Gal[igal].HotGas_elements[kk] += fwind * SNIIAllElements[kk];
     	    Gal[p].ColdGas_elements[kk] += (1.0-fwind) * SNIIAllElements[kk] + SNIaAllElements[kk] + AGBAllElements[kk];
#ifdef H2_AND_RINGS
     	    for(jj=0;jj<RNUM;jj++)
     	      Gal[p].ColdGasRings_elements[jj][kk] += fractionRings[jj] * ((1.0-fwind) * SNIIAllElements[kk] + SNIaAllElements[kk] + AGBAllElements[kk]);
#endif
     	  }
#else //SNIATOHOT
     	for(kk=0;kk<NUM_ELEMENTS;kk++)
     	  {
     	    Gal[igal].HotGas_elements[kk] += fwind * SNIIAllElements[kk] + SNIaAllElements[kk];
     	    Gal[p].ColdGas_elements[kk] += (1.0-fwind) * SNIIAllElements[kk] + AGBAllElements[kk];
#ifdef H2_AND_RINGS
     	    for(jj=0;jj<RNUM;jj++)
     	      Gal[p].ColdGasRings_elements[jj][kk] += fractionRings[jj] * ((1.0-fwind) * SNIIAllElements[kk] + AGBAllElements[kk]);
#endif
     	  }
#endif //SNIATOHOT

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
     	for(jj=0;jj<RNUM;jj++)
     	  {
     	    Gal[p].DiskMassRings[jj]-= fractionRings[jj] * (SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
     	    Gal[p].MetalsDiskMassRings[jj][0]-= fractionRings[jj] * SNIIUnProcessedMetals;
     	    Gal[p].MetalsDiskMassRings[jj][1]-= fractionRings[jj] * SNIaUnProcessedMetals;
     	    Gal[p].MetalsDiskMassRings[jj][2]-= fractionRings[jj] * AGBUnProcessedMetals;
     	  }
#endif
	//mass_checks(p,"model_yields.c",__LINE__);
#ifdef INDIVIDUAL_ELEMENTS

     	for(kk=0;kk<NUM_ELEMENTS;kk++)
     	  {
     	    Gal[p].DiskMass_elements[kk] -= SNIIUnProcessedElements[kk]+SNIaUnProcessedElements[kk]+AGBUnProcessedElements[kk];
#ifdef H2_AND_RINGS
     	    for(jj=0;jj<RNUM;jj++)
     	      Gal[p].DiskMassRings_elements[jj][kk] -= fractionRings[jj] * (SNIIUnProcessedElements[kk]+SNIaUnProcessedElements[kk]+AGBUnProcessedElements[kk]);
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

   if (Gal[p].sfh_BulgeMass[i] > 0.0)
    {
    	//pre-calculations to speed up the code
    	//Note: This is NOT really a SFR, as no stars are formed in the bulge. Rather, this is a star-transfer rate from the disc (or mergers) to the bulge.
    	BulgeSFRxStep = timestep_width * Gal[p].sfh_BulgeMass[i]/Gal[p].sfh_dt[i];
    	BulgeSFRxStep_Phys = BulgeSFRxStep * (1.0e10/Hubble_h);
    	BulgeMetallicity = 0.;
    	for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
    	  BulgeMetallicity += Gal[p].sfh_MetalsBulgeMass[i][ii];
    	BulgeMetallicity /= Gal[p].sfh_BulgeMass[i];
#ifdef INDIVIDUAL_ELEMENTS
    	for (kk=0;kk<NUM_ELEMENTS;kk++)
    	  BulgeMetallicityElement_Phys[kk] = Gal[p].sfh_BulgeMass_elements[i][kk] / (Gal[p].sfh_BulgeMass[i]*1.0e10/Hubble_h);
#endif

    	Zi = find_initial_metallicity(p, i, 1, 2);
    	//Interpolate the bulge luminosity on the lifetimeMetallicities tables:
    	Zi_disp = (BulgeMetallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.


#ifndef INDIVIDUAL_ELEMENTS
    	compute_actual_eject_rates(TimeBin, i, Zi, Zi_disp, Gal[p].sfh_ibin, BulgeSFRxStep, BulgeSFRxStep_Phys, BulgeMetallicity,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals);
#else
    	compute_actual_eject_rates(TimeBin, i, Zi, Zi_disp, Gal[p].sfh_ibin, BulgeSFRxStep, BulgeSFRxStep_Phys, BulgeMetallicity, BulgeMetallicityElement_Phys,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
				   SNIIAllElements, SNIIUnProcessedElements, SNIaAllElements, SNIaUnProcessedElements, AGBAllElements, AGBUnProcessedElements);
#endif




//#ifdef BULGE_TO_COLD
#ifdef H2_AND_RINGS
    	//distribution of mass being deposited from the bulge into the cold gas - Jaffe profile
#ifndef RINGS_IN_BULGES
	double rb=Gal[p].BulgeSize, TotMassInsideRings=0.;
	
	if(rb>0.)
	  TotMassInsideRings=(RingRadius[RNUM-1]/rb)/(1+RingRadius[RNUM-1]/rb);
	
	if(TotMassInsideRings>0.)
	  {
	    fractionRingsBulge[0]=(RingRadius[0]/rb)/(1+RingRadius[0]/rb)/TotMassInsideRings;
	    for(jj=1; jj<RNUM; jj++)
		fractionRingsBulge[jj]= ((RingRadius[jj]/rb)/(1+RingRadius[jj]/rb)-(RingRadius[jj-1]/rb)/(1+RingRadius[jj-1]/rb))/TotMassInsideRings;
	  }
	else
	  for(jj=0; jj<RNUM; jj++)
	    //    fractionRingsBulge[jj]=1./RNUM;
		fractionRingsBulge[jj]=0.;
#else
	for(jj=0;jj<RNUM;jj++)
	  fractionRingsBulge[jj]=Gal[p].sfh_BulgeMassRings[jj][i]/Gal[p].sfh_BulgeMass[i];
	  //fractionRingsBulge[jj]=Gal[p].BulgeMassRings[jj]/Gal[p].BulgeMass;
#endif
#endif
//#endif



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
    	for(kk=0;kk<NUM_ELEMENTS;kk++)
    	  Gal[igal].HotGas_elements[kk] += SNIIAllElements[kk]+SNIaAllElements[kk]+AGBAllElements[kk];
#endif //INDIVIDUAL_ELEMENTS



#else //BULGE_TO_COLD
    	Gal[p].ColdGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	Gal[p].MetalsColdGas[0] += SNIIAllMetals;
    	Gal[p].MetalsColdGas[1] += SNIaAllMetals;
    	Gal[p].MetalsColdGas[2] += AGBAllMetals;

#ifdef H2_AND_RINGS
     	for(jj=0;jj<RNUM;jj++)
     	  {
     	    Gal[p].ColdGasRings[jj] += fractionRingsBulge[jj] *  (SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
     	    Gal[p].MetalsColdGasRings[jj][0] += fractionRingsBulge[jj] *  SNIIAllMetals;
     	    Gal[p].MetalsColdGasRings[jj][1] += fractionRingsBulge[jj] *  SNIaAllMetals;
     	    Gal[p].MetalsColdGasRings[jj][2] += fractionRingsBulge[jj] *  AGBAllMetals;
     	  }
#endif

     	TotalMassReturnedToColdDiskGas += SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
#ifdef H2_AND_RINGS
     	for(jj=0;jj<RNUM;jj++)
     	  TotalMassReturnedToColdDiskGasr[jj]+=fractionRingsBulge[jj] *  (SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
#endif

#ifdef INDIVIDUAL_ELEMENTS
     	for(kk=0;kk<NUM_ELEMENTS;kk++)
     	  {
     	    Gal[p].ColdGas_elements[kk] += SNIIAllElements[kk]+SNIaAllElements[kk]+AGBAllElements[kk];
#ifdef H2_AND_RINGS
     	    for(jj=0;jj<RNUM;jj++)
     	      Gal[p].ColdGasRings_elements[jj][kk]  += fractionRingsBulge[jj] * (SNIIAllElements[kk] + SNIaAllElements[kk] + AGBAllElements[kk]);
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

#ifdef RINGS_IN_BULGES
     	for(jj=0;jj<RNUM;jj++)
     	  {
     	    Gal[p].BulgeMassRings[jj]-= fractionRingsBulge[jj] * (SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass);
     	    Gal[p].MetalsBulgeMassRings[jj][0]-= fractionRingsBulge[jj] * SNIIUnProcessedMetals;
     	    Gal[p].MetalsBulgeMassRings[jj][1]-= fractionRingsBulge[jj] * SNIaUnProcessedMetals;
     	    Gal[p].MetalsBulgeMassRings[jj][2]-= fractionRingsBulge[jj] * AGBUnProcessedMetals;
     	  }
#endif

#ifdef INDIVIDUAL_ELEMENTS
    	for(kk=0;kk<NUM_ELEMENTS;kk++)
    	  {
    	    Gal[p].BulgeMass_elements[kk] -= SNIIUnProcessedElements[kk]+SNIaUnProcessedElements[kk]+AGBUnProcessedElements[kk];
#ifdef RINGS_IN_BULGES
     	    for(jj=0;jj<RNUM;jj++)
     	      Gal[p].BulgeMassRings_elements[jj][kk] -= fractionRingsBulge[jj] * (SNIIUnProcessedElements[kk]+SNIaUnProcessedElements[kk]+AGBUnProcessedElements[kk]);
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








#ifdef CHECK_SFH_AND_RINGS

/*********************************************
 *
 * Ensure that no values became negative due
 * to precision
 *
 *********************************************/
int mm;

#ifdef CHECK_NO_NEGATIVE_VALUES
  //cold gas galaxy p
  if(Gal[p].ColdGas < 0.0)
    if(Gal[p].ColdGas > -1e-3)
      {
	Gal[p].ColdGas = 0.;
	for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
	  Gal[p].MetalsColdGas[mm] = 0.;
#ifdef INDIVIDUAL_ELEMENTS
	for(kk=0;kk<NUM_ELEMENTS;kk++)
	  Gal[q].ColdGas_elements[kk] = 0.;
#endif
#ifdef H2_AND_RINGS
	for (jj=0; jj<RNUM; jj++)
   	  {
   	    Gal[p].ColdGasRings[jj]=0.;
   	    for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
   	      Gal[p].MetalsColdGasRings[jj][mm] = 0.;
#ifdef INDIVIDUAL_ELEMENTS
   	 for(kk=0;kk<NUM_ELEMENTS;kk++)
   	    Gal[p].ColdGasRings_elements[jj][kk] = 0.;
#endif
   	  }
#endif
      }



  //disk mass galaxy p
  if(Gal[p].DiskMass < 0.0)
    if(Gal[p].DiskMass > -1e-3)
      {
  	Gal[p].DiskMass = 0.;
  	for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
  	  Gal[p].MetalsDiskMass[mm] = 0.;
#ifdef INDIVIDUAL_ELEMENTS
  	for(kk=0;kk<NUM_ELEMENTS;kk++)
  	Gal[q].DiskMass_elements[kk] = 0.;
#endif
#ifdef H2_AND_RINGS
  	for (jj=0; jj<RNUM; jj++)
  	  {
  	    Gal[p].DiskMassRings[jj]=0.;
  	    for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
  	      Gal[p].MetalsDiskMassRings[jj][mm] = 0.;
#ifdef INDIVIDUAL_ELEMENTS
  	  for(kk=0;kk<NUM_ELEMENTS;kk++)
  	    Gal[p].DiskMassRings_elements[jj][kk] = 0.;
#endif
  	  }
#endif

#ifdef STAR_FORMATION_HISTORY
       for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
	 {
	   Gal[p].sfh_DiskMass[ii] = 0.;
	   for(jj=0;jj<RNUM;jj++)
	     Gal[p].sfh_DiskMassRings[jj][ii] =0.;
	   for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
	     Gal[p].sfh_MetalsDiskMass[ii][mm] =0.;
#ifdef INDIVIDUAL_ELEMENTS
	   for(kk=0;kk<NUM_ELEMENTS;kk++)
	     Gal[p].sfh_DiskMass_elements[ii][kk] =0.;
#endif
	 }
#endif
       }


  //bulge mass galaxy p
  if(Gal[p].BulgeMass < 0.0)
    if(Gal[p].BulgeMass > -1e-3)
      {
	Gal[p].BulgeMass = 0.;
	for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
	  Gal[p].MetalsBulgeMass[mm] = 0.;
#ifdef INDIVIDUAL_ELEMENTS
    	for(kk=0;kk<NUM_ELEMENTS;kk++)
    	  Gal[q].BulgeMass_elements[kk] = 0.;
#endif
#ifdef H2_AND_RINGS
#ifdef RINGS_IN_BULGES
    	for (jj=0; jj<RNUM; jj++)
    	  {
    	    Gal[p].BulgeMassRings[jj]=0.;
    	    for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
    	      Gal[p].MetalsBulgeMassRings[jj][mm] = 0.;
#ifdef INDIVIDUAL_ELEMENTS
    	    for(kk=0;kk<NUM_ELEMENTS;kk++)
    	      Gal[p].BulgeMassRings_elements[jj][kk] = 0.;
#endif
    	  }
#endif
#endif

#ifdef STAR_FORMATION_HISTORY
       for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
	 {
	   Gal[p].sfh_BulgeMass[ii] = 0.;
	   for(jj=0;jj<RNUM;jj++)
	     Gal[p].sfh_BulgeMassRings[jj][ii] =0.;
	   for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
	     Gal[p].sfh_MetalsBulgeMass[ii][mm] =0.;
#ifdef INDIVIDUAL_ELEMENTS
	   for(kk=0;kk<NUM_ELEMENTS;kk++)
	     Gal[p].sfh_BulgeMass_elements[ii][kk] =0.;
#endif
	 }
#endif
      }


  //set all negative rings to zero
  for (jj=0; jj<RNUM; jj++)
    {
      //gasmass
      if(Gal[p].ColdGasRings[jj]<0.)
     	{
     	  Gal[p].ColdGasRings[jj]=0.;
          for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
            Gal[p].MetalsColdGasRings[jj][mm] = 0.;
     #ifdef INDIVIDUAL_ELEMENTS
       	  for(kk=0;kk<NUM_ELEMENTS;kk++)
       	    Gal[p].ColdGasRings_elements[jj][kk] = 0.;
     #endif
     	}

      for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
     	if(Gal[p].MetalsColdGasRings[jj][mm]<0.)
                Gal[p].MetalsColdGasRings[jj][mm] = 0.;
     #ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
     	if(Gal[p].ColdGasRings_elements[jj][kk] < 0.)
     	  Gal[p].ColdGasRings_elements[jj][kk] = 0.;
     #endif

      //diskmass
      if(Gal[p].DiskMassRings[jj]<0.)
	{
	  Gal[p].DiskMassRings[jj]=0.;
#ifdef STAR_FORMATION_HISTORY
	  for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
	     Gal[p].sfh_DiskMassRings[jj][ii] =0.;
#endif
       	  for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
       	      Gal[p].MetalsDiskMassRings[jj][mm] = 0.;
#ifdef INDIVIDUAL_ELEMENTS
  	  for(kk=0;kk<NUM_ELEMENTS;kk++)
  	    Gal[p].DiskMassRings_elements[jj][kk] = 0.;
#endif
	}

      for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
	if(Gal[p].MetalsDiskMassRings[jj][mm]<0.)
           Gal[p].MetalsDiskMassRings[jj][mm] = 0.;
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	if(Gal[p].DiskMassRings_elements[jj][kk] < 0.)
	  Gal[p].DiskMassRings_elements[jj][kk] = 0.;
#endif

#ifdef RINGS_IN_BULGES
      //bulgemass
      if(Gal[p].BulgeMassRings[jj]<0.)
      	{
      	  Gal[p].BulgeMassRings[jj]=0.;
#ifdef STAR_FORMATION_HISTORY
	  for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
	     Gal[p].sfh_BulgeMassRings[jj][ii] =0.;
#endif
          for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
            Gal[p].MetalsBulgeMassRings[jj][mm] = 0.;
#ifdef INDIVIDUAL_ELEMENTS
        	  for(kk=0;kk<NUM_ELEMENTS;kk++)
        	    Gal[p].BulgeMassRings_elements[jj][kk] = 0.;
#endif
      	}

      for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
      	if(Gal[p].MetalsBulgeMassRings[jj][mm]<0.)
                 Gal[p].MetalsBulgeMassRings[jj][mm] = 0.;
      #ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
      	if(Gal[p].BulgeMassRings_elements[jj][kk] < 0.)
      	  Gal[p].BulgeMassRings_elements[jj][kk] = 0.;
      #endif
#endif

    }

#endif //CHECK_NO_NEGATIVE_VALUES







  /*************************************************
   *
   * set rings to 0 if total is 0
   *
   */
    //set metals and elements to 0
    //galaxy p cold gas
  //set rings to 0
 #ifdef H2_AND_RINGS
   if(Gal[p].ColdGas == 0. || Gal[p].DiskMass == 0. || Gal[p].BulgeMass == 0.)
     for (jj=0; jj<RNUM; jj++)
       {
 	//galaxy p cold gas
 	if(Gal[p].ColdGas == 0.)
 	  {
 	    Gal[p].ColdGasRings[jj]=0.;
 	    for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
 	      Gal[p].MetalsColdGasRings[jj][mm] = 0.;
 #ifdef INDIVIDUAL_ELEMENTS
 	   for(kk=0;kk<NUM_ELEMENTS;kk++)
 	    Gal[p].ColdGasRings_elements[jj][kk] = 0.;
 #endif
 	  }

 	//galaxy p disk mass
 	if(Gal[p].DiskMass == 0.)
 	  {
 	    Gal[p].DiskMassRings[jj]=0.;
 	    for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
 	      Gal[p].MetalsDiskMassRings[jj][mm] = 0.;
 #ifdef INDIVIDUAL_ELEMENTS
 	   for(kk=0;kk<NUM_ELEMENTS;kk++)
 	    Gal[p].DiskMassRings_elements[jj][kk] = 0.;
 #endif
 	  }


#ifdef RINGS_IN_BULGES
 	//galaxy p Bulge mass
 	if(Gal[p].BulgeMass == 0.)
 	  {
 	    Gal[p].BulgeMassRings[jj]=0.;
 	    for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
 	      Gal[p].MetalsBulgeMassRings[jj][mm] = 0.;
#ifdef INDIVIDUAL_ELEMENTS
 	    for(kk=0;kk<NUM_ELEMENTS;kk++)
 	      Gal[p].BulgeMassRings_elements[jj][kk] = 0.;
#endif
 	  }
#endif //RINGS_IN_BULGES

       }
 #endif



/*************************************************
 *
 * check to make sure sum of rings = total mass
 * if not update the rings
 *
 * // Mass += fractionRings[jj]/RNUM*Gal[q].HotGas;
 * The correction is never larger than 10% except when
 * total is non-zero and sum of rings is zero or vice
 * versa. This only happens for very small masses due
 * to precision.
 *
 * If the correction is not always made, over enough
 * calculations large difference can appear between
 * total and rings sum.
 */

  double sum_rings, sum_rings_metals, ring_sum_minus_tot=0.;
  double fract, fract_metals;

  //galaxy p cold gas
  if(Gal[p].ColdGas>0.)
    {
      sum_rings=0.;
      sum_rings_metals=0.;
      for (jj=0; jj<RNUM; jj++)
	{
	  sum_rings+=Gal[p].ColdGasRings[jj];
	  for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
	    sum_rings_metals+= Gal[p].MetalsColdGasRings[jj][mm];
	} // Mass += fractionRings[jj]/RNUM*Gal[q].HotGas;

      if(sum_rings>0.)
	fract=Gal[p].ColdGas/sum_rings;
      else
	{
	  //if there is no mass in rings, reset total to 0
	  fract=0.;
	  Gal[p].ColdGas=0.;
	}

      if(sum_rings_metals>0.)
	{
	  fract_metals=0.;
	  for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
	    fract_metals+= Gal[p].MetalsColdGas[mm];
	  fract_metals /= sum_rings_metals;
	}
      else
	{
	  //if there is no mass in rings, reset total to 0
	  fract_metals=0.;
	  for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
		Gal[p].MetalsColdGas[mm]=0.;
	}

      if(sum_rings/Gal[p].ColdGas !=1.)
	{
	  for (jj=0;jj<RNUM;jj++)
	    {
	      Gal[p].ColdGasRings[jj] *= fract;
	      for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
		Gal[p].MetalsColdGasRings[jj][mm] = (Gal[p].MetalsColdGasRings[jj][mm]*fract_metals);

	    }

	  ring_sum_minus_tot=-Gal[p].ColdGas;
	  for (jj=0; jj<RNUM; jj++) // Mass += fractionRings[jj]/RNUM*Gal[q].HotGas;
	    ring_sum_minus_tot+=Gal[p].ColdGasRings[jj];
	  if((ring_sum_minus_tot < -1e-4 && ring_sum_minus_tot/Gal[p].ColdGas < -1e-4) ||
	      (ring_sum_minus_tot >  1e-4 && ring_sum_minus_tot/Gal[p].ColdGas >  1e-4))
	    {
	      printf("ring_sum_minus_tot=%0.10e 1e-5*ColdGas=%0.10e\n", ring_sum_minus_tot, 1e-5*Gal[p].ColdGas);
	      terminate("");
	    }
	}
    }


  //printf("DiskMass= %e\n",Gal[p].DiskMass);
   //galaxy p disk mass
  if(Gal[p].DiskMass>0.)
    {
      sum_rings=0.;
      sum_rings_metals=0;
      for (jj=0; jj<RNUM; jj++)
	{
	  sum_rings+=Gal[p].DiskMassRings[jj];
	  for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
	    sum_rings_metals+= Gal[p].MetalsDiskMassRings[jj][mm];
     	}

      //printf("sum_rings= %e\n",sum_rings);

      if(sum_rings>0.)
	fract=Gal[p].DiskMass/sum_rings;
      else
	{
	 //if there is no mass in rings, reset total to 0
	 fract=0.;
	Gal[p].DiskMass=0.;
	}
      /*if(fract>0. && Gal[p].DiskMass>1e-6 && sum_rings>1e-6  && (fract>1.3 || fract<0.7))
      	{
      	  printf("fract=%0.5f %s %d\n",fract,call_function, call_line);
      	  terminate("Total and Ring Sum for Metals too different Disk Mass");
      	}*/
      if(sum_rings_metals>0.)
      	{
       	 fract_metals=0.;
         for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
           fract_metals+= Gal[p].MetalsDiskMass[mm];
	 fract_metals /= sum_rings_metals;
        }
      else
	{
	//if there is no mass in rings, reset total to 0
	fract_metals=0.;
	for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
	  Gal[p].MetalsDiskMass[mm]=0.;
	}

      //printf("fract= %e\n",fract);

      if(sum_rings/Gal[p].DiskMass !=1.)
	{
	  for (jj=0;jj<RNUM;jj++)
	    {
	      Gal[p].DiskMassRings[jj] *= fract;
	      for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
		Gal[p].MetalsDiskMassRings[jj][mm] = (Gal[p].MetalsDiskMassRings[jj][mm]*fract_metals);

	    }

	  ring_sum_minus_tot=-Gal[p].DiskMass;
	  for (jj=0; jj<RNUM; jj++)
	    ring_sum_minus_tot+=Gal[p].DiskMassRings[jj];
	  if((ring_sum_minus_tot < -1e-6 && ring_sum_minus_tot < -1e-6*Gal[p].DiskMass) ||
	      (ring_sum_minus_tot >  1e-6 && ring_sum_minus_tot >  1e-6*Gal[p].DiskMass))
	    {
	      printf("ring_sum_minus_tot=%0.10e 1e-5*ColdGas=%0.10e\n", ring_sum_minus_tot, 1e-5*Gal[p].DiskMass);
	      terminate("");
	    }
	}
    }


#endif //CHECK_SFH_AND_RINGS







    //*****************************************
    //ENRICHMENT FROM ICL STARS INTO HOT GAS:
    //*****************************************

    if (Gal[p].sfh_ICM[i] > 0.0)
      {
    	//pre-calculations to speed up the code
    	//Note: This is NOT really an SFR, as no stars are formed in the ICM. Rather, this is a star-transfer rate from satellite disruption to the stellar halo.
    	ICMSFRxStep = timestep_width * Gal[p].sfh_ICM[i]/Gal[p].sfh_dt[i];
    	ICMSFRxStep_Phys = ICMSFRxStep * (1.0e10/Hubble_h) ;
    	ICMMetallicity=0.;
    	for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
    	  ICMMetallicity += Gal[p].sfh_MetalsICM[i][ii];
    	ICMMetallicity /= Gal[p].sfh_ICM[i];
#ifdef INDIVIDUAL_ELEMENTS
    	for (kk=0;kk<NUM_ELEMENTS;kk++)
    	  ICMMetallicityElement_Phys[kk] = Gal[p].sfh_ICM_elements[i][kk] / (Gal[p].sfh_ICM[i]*1.0e10/Hubble_h);
#endif

    	Zi = find_initial_metallicity(p, i, 1, 3);
    	//Interpolate the ICM metallicity on the lifetimeMetallicities tables:
    	Zi_disp = (ICMMetallicity - lifetimeMetallicities[Zi])/(lifetimeMetallicities[Zi+1] - lifetimeMetallicities[Zi]);
    	if (Zi_disp < 0.0) Zi_disp = 0.0; //Don't want to extrapolate yields down below lifetimeMetallicities[0]=0.0004. Instead, assume constant yield below this metallicity.


#ifndef INDIVIDUAL_ELEMENTS
    	compute_actual_eject_rates(TimeBin, i, Zi, Zi_disp, Gal[p].sfh_ibin, ICMSFRxStep, ICMSFRxStep_Phys, ICMMetallicity,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals);
#else
    	compute_actual_eject_rates(TimeBin, i, Zi, Zi_disp, Gal[p].sfh_ibin, ICMSFRxStep, ICMSFRxStep_Phys, ICMMetallicity, ICMMetallicityElement_Phys,
				   &SNIIEjectaMass, &SNIIAllMetals, &SNIIUnProcessedMetals,
				   &SNIaEjectaMass, &SNIaAllMetals, &SNIaUnProcessedMetals,
				   &AGBEjectaMass, &AGBAllMetals, &AGBUnProcessedMetals,
				   SNIIAllElements, SNIIUnProcessedElements, SNIaAllElements, SNIaUnProcessedElements, AGBAllElements, AGBUnProcessedElements);
#endif




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
    	for(kk=0;kk<NUM_ELEMENTS;kk++)
    	  Gal[igal].HotGas_elements[kk] += SNIIAllElements[kk] + SNIaAllElements[kk] + AGBAllElements[kk];
#endif //INDIVIDUAL_ELEMENTS


    	//*****************************
    	// UPDATE ICM COMPONENTS:
    	//*****************************
    	Gal[p].ICM -= SNIIEjectaMass + SNIaEjectaMass + AGBEjectaMass;
    	Gal[p].MetalsICM[0] -= SNIIUnProcessedMetals;
    	Gal[p].MetalsICM[1] -= SNIaUnProcessedMetals;
    	Gal[p].MetalsICM[2] -= AGBUnProcessedMetals;

#ifdef INDIVIDUAL_ELEMENTS
    	for(kk=0;kk<NUM_ELEMENTS;kk++)
    	  Gal[p].ICM_elements[kk] -= SNIIUnProcessedElements[kk] + SNIaUnProcessedElements[kk] + AGBUnProcessedElements[kk];
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
void compute_actual_eject_rates(int TimeBin, int i, int Zi, double Zi_disp, int sfh_ibin, double SFRxStep, double SFRxStep_Phys, double Metallicity,
				 double *SNIIEjectaMass, double *SNIIAllMetals, double *SNIIUnProcessedMetals,
				 double *SNIaEjectaMass, double *SNIaAllMetals, double *SNIaUnProcessedMetals,
				 double *AGBEjectaMass, double *AGBAllMetals, double *AGBUnProcessedMetals)
#else
void compute_actual_eject_rates(int TimeBin, int i, int Zi, double Zi_disp, int sfh_ibin, double SFRxStep, double SFRxStep_Phys, double Metallicity, double *MetallicityElement_Phys,
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
  int kk;
  double NormSNIIYieldRate_actual[NUM_ELEMENTS], NormSNIaYieldRate_actual[NUM_ELEMENTS], NormAGBYieldRate_actual[NUM_ELEMENTS];
#endif

//pre-calculations to speed up the code
  NormSNIIMassEjecRate_actual = NormSNIIMassEjecRate[TimeBin][i][Zi] + ((NormSNIIMassEjecRate[TimeBin][i][Zi+1] - NormSNIIMassEjecRate[TimeBin][i][Zi])*Zi_disp);
  NormSNIaMassEjecRate_actual = NormSNIaMassEjecRate[TimeBin][i][Zi] + ((NormSNIaMassEjecRate[TimeBin][i][Zi+1] - NormSNIaMassEjecRate[TimeBin][i][Zi])*Zi_disp);
  NormAGBMassEjecRate_actual = NormAGBMassEjecRate[TimeBin][i][Zi] + ((NormAGBMassEjecRate[TimeBin][i][Zi+1] - NormAGBMassEjecRate[TimeBin][i][Zi])*Zi_disp);
  NormSNIIMetalEjecRate_actual = NormSNIIMetalEjecRate[TimeBin][i][Zi] + ((NormSNIIMetalEjecRate[TimeBin][i][Zi+1] - NormSNIIMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
  NormSNIaMetalEjecRate_actual = NormSNIaMetalEjecRate[TimeBin][i][Zi] + ((NormSNIaMetalEjecRate[TimeBin][i][Zi+1] - NormSNIaMetalEjecRate[TimeBin][i][Zi])*Zi_disp);
  NormAGBMetalEjecRate_actual = NormAGBMetalEjecRate[TimeBin][i][Zi] + ((NormAGBMetalEjecRate[TimeBin][i][Zi+1] - NormAGBMetalEjecRate[TimeBin][i][Zi])*Zi_disp);

#ifdef INDIVIDUAL_ELEMENTS //Work out the actual yield of element k, by interpolating between the yields in the look-up table created by yield_integrals.c.
  for (kk=0;kk<NUM_ELEMENTS;kk++)
    {
      NormSNIIYieldRate_actual[kk] = NormSNIIYieldRate[TimeBin][i][Zi][kk] + ((NormSNIIYieldRate[TimeBin][i][Zi+1][kk] - NormSNIIYieldRate[TimeBin][i][Zi][kk])*Zi_disp);
      NormSNIaYieldRate_actual[kk] = NormSNIaYieldRate[TimeBin][i][Zi][kk] + ((NormSNIaYieldRate[TimeBin][i][Zi+1][kk] - NormSNIaYieldRate[TimeBin][i][Zi][kk])*Zi_disp);
      NormAGBYieldRate_actual[kk] = NormAGBYieldRate[TimeBin][i][Zi][kk] + ((NormAGBYieldRate[TimeBin][i][Zi+1][kk] - NormAGBYieldRate[TimeBin][i][Zi][kk])*Zi_disp);
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
  reset_ejection_rates(i, sfh_ibin,
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
  for(kk=0;kk<NUM_ELEMENTS;kk++)
    {
#ifdef PORTINARI
      SNIIAllElements[kk] = SFRxStep_Phys * (NormSNIIYieldRate_actual[kk] + MetallicityElement_Phys[kk] * NormSNIIMassEjecRate_actual);
#endif
#ifdef CHIEFFI
      SNIIAllElements[kk] = SFRxStep_Phys * NormSNIIYieldRate_actual[kk];
#endif
      SNIIUnProcessedElements[kk] = SFRxStep_Phys * MetallicityElement_Phys[kk] * NormSNIIMassEjecRate_actual;

      SNIaAllElements[kk] = SFRxStep_Phys * NormSNIaYieldRate_actual[kk];
      SNIaUnProcessedElements[kk] = SFRxStep_Phys * MetallicityElement_Phys[kk] * NormSNIaMassEjecRate_actual;

      AGBAllElements[kk] = SFRxStep_Phys * (NormAGBYieldRate_actual[kk] + MetallicityElement_Phys[kk] * NormAGBMassEjecRate_actual);
      AGBUnProcessedElements[kk] = SFRxStep_Phys * MetallicityElement_Phys[kk] * NormAGBMassEjecRate_actual;
    }
#endif

}








int find_initial_metallicity(int p, int sfh_bin, int table_type, int component_type)
{
  int ii;
  if (component_type == 1) //Disk stars
	{
	int i, Zi_bin;
	double initMetals=0., Z_disk;

	for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	  initMetals += Gal[p].sfh_MetalsDiskMass[sfh_bin][ii];
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
		double initMetals=0., Z_bulge;

		for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
			  initMetals += Gal[p].sfh_MetalsBulgeMass[sfh_bin][ii];
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
					if (lifetimeMetallicities[i] < Z_bulge) //Gal[p].sfh_MetalsDiskMass[sfh_bin][0]/Gal[p].sfh_DiskMass[sfh_bin])
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
			double initMetals=0., Z_ICM;

			for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
				  initMetals += Gal[p].sfh_MetalsICM[sfh_bin][ii];
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

