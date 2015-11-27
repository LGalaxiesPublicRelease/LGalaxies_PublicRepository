#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"

#include "mcmc_vars.h"
#include "mcmc_proto.h"


#ifdef MCMC
/** @brief Writes galaxies into a structure to be used by the MCMC */

void save_galaxy_for_mcmc(int gal_index)
{
  int snap, fof, j, aa;
  float low_mass_limit, high_mass_limit, StellarMass;
  double log10_Hubble_h, PMass;

#ifdef MR_PLUS_MRII
  if(Switch_MR_MRII==1)
    {
      low_mass_limit=9.5;
      high_mass_limit=14.0;
    }
  else
    if(Switch_MR_MRII==2)
      {
	low_mass_limit=6.0;
	high_mass_limit=9.5;
      }
#else
  low_mass_limit=7.27;
  high_mass_limit=13.0;
#ifdef MRII
  low_mass_limit=6.0;
  high_mass_limit=11.27;
#endif
#endif


  log10_Hubble_h=log10(Hubble_h);

  for(snap=0;snap<NOUT;snap++)
    {
#ifndef HALOMODEL
      StellarMass=log10(1E10 * (HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)*Hubble_h);
#else
      //no h factor in masses for OPT+=DHALOMODEL
      StellarMass=log10(1E10 * (HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)/Hubble_h);
#endif
      //THE ERROR IS NOW INCLUDED IN mcmc_likelihood.c
      //StellarMass+=gassdev(&MCMCseed)*0.08*(1+MCMCConstraintsZZ[snap]);

      for(fof=0;fof<NFofsInSample[snap]; fof++)
	if( StellarMass > low_mass_limit &&  StellarMass < high_mass_limit &&
	    HaloIDs[HaloGal[gal_index].HaloNr].FirstHaloInFOFgroup == MCMC_FOF[fof].FoFID[snap])
	  {
	    MCMC_GAL[TotMCMCGals[snap]].StellarMass[snap] = StellarMass;
	    MCMC_GAL[TotMCMCGals[snap]].ColdGas[snap] = log10(1E10 * (HaloGal[gal_index].ColdGas*Hubble_h));
	    MCMC_GAL[TotMCMCGals[snap]].BulgeMass[snap] = log10(1E10 * HaloGal[gal_index].BulgeMass*Hubble_h);
	    MCMC_GAL[TotMCMCGals[snap]].BlackHoleMass[snap] = log10(1E10 * HaloGal[gal_index].BlackHoleMass); //black hole in units of h^-1
	    //in units of Solar Masses per yr
	    MCMC_GAL[TotMCMCGals[snap]].Sfr[snap]= log10(HaloGal[gal_index].Sfr * UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS*Hubble_h);

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef POST_PROCESS_MAGS
	    struct GALAXY_OUTPUT galaxy_output;

	    //in case of postprocess magnitudes they are only calculates here, inside prepare
	    prepare_galaxy_for_output(snap, &HaloGal[gal_index], &galaxy_output);

	    MCMC_GAL[TotMCMCGals[snap]].MagU[snap] = galaxy_output.MagDust[0]-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].MagB[snap] = galaxy_output.MagDust[1]-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].MagV[snap] = galaxy_output.MagDust[2]-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].MagJ[snap] = galaxy_output.MagDust[3]-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].MagK[snap] = galaxy_output.MagDust[4]-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].Magu[snap] = galaxy_output.MagDust[5]-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].Magr[snap] = galaxy_output.MagDust[6]-5.*log10_Hubble_h;


#else
	    MCMC_GAL[TotMCMCGals[snap]].MagU[snap] = lum_to_mag(HaloGal[gal_index].LumDust[0][snap])-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].MagB[snap] = lum_to_mag(HaloGal[gal_index].LumDust[1][snap])-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].MagV[snap] = lum_to_mag(HaloGal[gal_index].LumDust[2][snap])-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].MagJ[snap] = lum_to_mag(HaloGal[gal_index].LumDust[3][snap])-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].MagK[snap] = lum_to_mag(HaloGal[gal_index].LumDust[4][snap])-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].Magu[snap] = lum_to_mag(HaloGal[gal_index].LumDust[5][snap])-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].Magr[snap] = lum_to_mag(HaloGal[gal_index].LumDust[6][snap])-5.*log10_Hubble_h;
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

#ifdef HALOMODEL
	    if(snap==0)
	      {
		MCMC_GAL[TotMCMCGals[snap]].fofid[snap] = fof;
		//MCMC_GAL[TotMCMCGals[snap]].M_Crit200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].Len*PartMass*1.e10);
		MCMC_GAL[TotMCMCGals[snap]].M_Crit200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].M_Crit200*1.e10);
		MCMC_GAL[TotMCMCGals[snap]].M_Mean200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].M_Mean200*1.e10);
#ifdef MCRIT
		MCMC_GAL[TotMCMCGals[snap]].M_Mean200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].M_Crit200*1.e10);
#endif
		MCMC_GAL[TotMCMCGals[snap]].x[snap] = HaloGal[gal_index].Pos[0];
		MCMC_GAL[TotMCMCGals[snap]].y[snap] = HaloGal[gal_index].Pos[1];
		MCMC_GAL[TotMCMCGals[snap]].z[snap] = HaloGal[gal_index].Pos[2];
		MCMC_GAL[TotMCMCGals[snap]].Type[snap] = HaloGal[gal_index].Type;
		MCMC_GAL[TotMCMCGals[snap]].ngal[snap] = 0;
	      }
#endif

	    MCMC_GAL[TotMCMCGals[snap]].Weight[snap] = MCMC_FOF[fof].Weight[snap];

#ifdef HALOMODEL
	    //NOW GET PROPERTIES FOR FOF GROUPS - done for the particular fof that current galaxy resides in
	    ++MCMC_FOF[fof].NGalsInFoF[snap];

	    if(HaloGal[gal_index].Type==0)
	      {
		MCMC_FOF[fof].IndexOfCentralGal[snap]=TotMCMCGals[snap];
		//MCMC_FOF[fof].M_Crit200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].Len*PartMass*1.e10);
		MCMC_FOF[fof].M_Crit200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].M_Crit200*1.e10);
		MCMC_FOF[fof].M_Mean200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].M_Mean200*1.e10);
#ifdef MCRIT
		MCMC_FOF[fof].M_Mean200[snap] = log10(Halo[HaloGal[gal_index].HaloNr].M_Crit200*1.e10);
#endif
    		  }
#endif //HALOMODEL



	    ++TotMCMCGals[snap];

	    if(TotMCMCGals[snap] > MCMCAllocFactor)
	      terminate("Maximum number of galaxies in MCMC structure reached. Increase MCMCSmartFactor\n");
    	  }
      }

    return;
}
#endif


