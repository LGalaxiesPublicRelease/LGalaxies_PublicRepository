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
      StellarMass=log10(1E10 * (HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)*Hubble_h);
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
	    MCMC_GAL[TotMCMCGals[snap]].MagV[snap] = galaxy_output.MagDust[1]-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].MagJ[snap] = galaxy_output.MagDust[2]-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].Magu[snap] = galaxy_output.MagDust[3]-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].Magr[snap] = galaxy_output.MagDust[4]-5.*log10_Hubble_h;
#else
	    MCMC_GAL[TotMCMCGals[snap]].MagU[snap] = lum_to_mag(HaloGal[gal_index].LumDust[0][snap])-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].MagV[snap] = lum_to_mag(HaloGal[gal_index].LumDust[1][snap])-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].MagJ[snap] = lum_to_mag(HaloGal[gal_index].LumDust[2][snap])-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].Magu[snap] = lum_to_mag(HaloGal[gal_index].LumDust[3][snap])-5.*log10_Hubble_h;
	    MCMC_GAL[TotMCMCGals[snap]].Magr[snap] = lum_to_mag(HaloGal[gal_index].LumDust[4][snap])-5.*log10_Hubble_h;
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

	    MCMC_GAL[TotMCMCGals[snap]].Weight[snap] = MCMC_FOF[fof].Weight[snap];

	    ++TotMCMCGals[snap];

	    if(TotMCMCGals[snap] > MCMCAllocFactor)
	      terminate("Maximum number of galaxies in MCMC structure reached. Increase MCMCSmartFactor\n");
    	  }
      }

    return;
}
#endif


