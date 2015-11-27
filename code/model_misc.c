#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <sys/stat.h>

#include "allvars.h"
#include "proto.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

/**@file model_misc.c
 * @brief model_misc.c contains a mix of recipes used to: calculate disk
 *        sizes, initiate a galaxy structure, get the metallicity, add
 *        luminosities, convert snap to age, convert snap to z, calculate
 *        max of two numbers, get virial mass, get virial velocity, get virial
 *        radius, luminosity to mass, update central galaxy,
 *        update type 1 and type2, transfer stars and gas between galaxies. */

/** This routine is no longer called (after Guo2010)*/
double get_disk_radius(int halonr, int p)
{
  double SpinParameter;


  if(DiskRadiusModel == 1)
    {
      /*  See Mo, Mao & White (1998) eq12, and using a Bullock style lambda.  Since this is the scale length
         we take the typical star forming region as 3 times this using the Milky Way as an approximate guide */
      SpinParameter =
	sqrt(Halo[halonr].Spin[0] * Halo[halonr].Spin[0] + Halo[halonr].Spin[1] * Halo[halonr].Spin[1] +
	     Halo[halonr].Spin[2] * Halo[halonr].Spin[2]) / (1.414 * Gal[p].Vvir * Gal[p].Rvir);
      return 3.0 * (SpinParameter / 1.414) * Gal[p].Rvir;
    }
  else
    /*  simple prescription */
    return Gal[p].Rvir / 10.0;


}

/** @brief Updates the gas disk radius.
 *
 *  The gas disk is assumed to be thin, in centrifugal equilibrium and to have
 *  an exponential density profile:
 *
 *  \f$ \Sigma(R_{\rm{gas}})=
 *      \Sigma_{\rm{gas0}}e^{-\frac{R_{\rm{gas}}}{R_{\rm{gas,d}}}}, \f$
 *
 *  where \f$R_{\rm{gas,d}}\f$ is the scale length of the gas disk and
 *  \f$\Sigma_{\rm{gas0}}\f$ is the corresponding central surface density.
 *
 *  Assuming a flat circular velocity curve (galaxy with a negligible self-gravity)
 *  in an isothermal dark matter halo, the gas disk scale length is given by:
 *
 *  \f$ R_{\rm{gas,d}}=\frac{J_{\rm{gas}}/M_{\rm{gas}}}{2V_{\rm{max}}}\f$,
 *
 *  assuming conservation of the angular momentum of the cooling gas and that the
 *  maximum circular velocity of satellite galaxies does not change after infall
 *  (inner dark matter regions are compact and don't change). */

void get_gas_disk_radius(int p)
{
  double dgas;
  if(Gal[p].Type == 0)
    dgas =
      3.0 * sqrt(Gal[p].GasSpin[0] * Gal[p].GasSpin[0] + Gal[p].GasSpin[1] * Gal[p].GasSpin[1] +
		 Gal[p].GasSpin[2] * Gal[p].GasSpin[2]) / 2.0 / Gal[p].Vmax;
  else
    dgas =
      3.0 * sqrt(Gal[p].GasSpin[0] * Gal[p].GasSpin[0] + Gal[p].GasSpin[1] * Gal[p].GasSpin[1] +
		 Gal[p].GasSpin[2] * Gal[p].GasSpin[2]) / 2.0 / Gal[p].InfallVmax;

  Gal[p].GasDiskRadius = dgas;

  //if the spin is not available for the halo
  if(Gal[p].GasSpin[0]==0 && Gal[p].GasSpin[1]==0 && Gal[p].GasSpin[2]==0)
  	Gal[p].GasDiskRadius = Gal[p].Rvir / 10.0;
}


/** @brief Updates the stellar disk radius.
 *
 *  The stellar disk is assumed to be thin, in centrifugal equilibrium and to have
 *  an exponential density profile:
 *
 *  \f$ \Sigma(R_{\star})=
 *      \Sigma_{\star\rm{0}}e^{-\frac{R_{\star}}{R_{\rm{\star,d}}}}, \f$
 *
 *  where \f$R_{\rm{\star,d}}\f$ is the scale length of the gas disk and
 *  \f$\Sigma_{\star\rm{0}}\f$ is the corresponding central surface density.
 *
 *  Assuming a flat circular velocity curve (galaxy with a negligible self-gravity)
 *  in an isothermal dark matter halo, the gas disk scale length is given by:
 *
 *  \f$ R_{\rm{\star,d}}=\frac{J_{\star}/M_{\star}}{2V_{\rm{max}}}\f$,
 *
 *  assuming that the maximum circular velocity of satellite galaxies does not
 *  change after infall (inner dark matter regions are compact and don't change). */

void get_stellar_disk_radius(int p)
{
  double dstar;

  if(Gal[p].Type == 0)
    dstar = 3.0 * sqrt(Gal[p].StellarSpin[0] * Gal[p].StellarSpin[0] +
		       Gal[p].StellarSpin[1] * Gal[p].StellarSpin[1] +
		       Gal[p].StellarSpin[2] * Gal[p].StellarSpin[2]) / 2.0 / Gal[p].Vmax;
  else
    dstar = 3.0 * sqrt(Gal[p].StellarSpin[0] * Gal[p].StellarSpin[0] +
		       Gal[p].StellarSpin[1] * Gal[p].StellarSpin[1] +
		       Gal[p].StellarSpin[2] * Gal[p].StellarSpin[2]) / 2.0 / Gal[p].InfallVmax;

  Gal[p].StellarDiskRadius = dstar;

  //if the spin is not available for the halo
   if(dstar==0)
     Gal[p].StellarDiskRadius = Gal[p].Rvir / 10.0;
}

/** @brief Initiates the value of the disk radius.
 *
 *  First determination of radius in Guo2010 (same as in Delucia2007), after this,
 *  the disks are updated using get_gas_disk_radius and get_stellar_disk_radius.
 *  Two options are available:
 *
 *    If DiskRadiusModel = 2 then \f$R_{\rm{disk}}=\frac{R_{\rm{vir}}}{10}\f$
 *
 *    If DiskRadiusModel = 0 or 1 then the Mo, Mao & White (1998) formalism is
 *    used with a Bullock style \f$\lambda\f$:
 *
 *    \f$ R_d=\frac{1}{\sqrt{2}}\frac{j_d}{m_d}\lambda r_{200}\f$
 *
 *    and using the Milky Way as an approximate guide \f$R_{\rm{disk}}=3R_d\f$. */

double get_initial_disk_radius(int halonr, int p)
{
  double SpinParameter, dgas;

  if(DiskRadiusModel == 0 || DiskRadiusModel == 1)
    {
      /*spin parameter */
      SpinParameter = sqrt(Halo[halonr].Spin[0] * Halo[halonr].Spin[0] +
			   Halo[halonr].Spin[1] * Halo[halonr].Spin[1] +
  			   Halo[halonr].Spin[2] * Halo[halonr].Spin[2]) / (1.414 * Gal[p].Vvir * Gal[p].Rvir);

      if(Gal[p].GasSpin[0]==0 && Gal[p].GasSpin[1]==0 && Gal[p].GasSpin[2]==0)
	dgas = Gal[p].Rvir / 10.0;
      else
	dgas = 3.0 * (SpinParameter / 1.414) * Gal[p].Rvir;

      return dgas;
    }
  else
    /*  simpler prescription */
    return Gal[p].Rvir / 10.0;
}



/** @brief Initializes the Galaxy Structure by setting all its
 *         elements to zero. */
void init_galaxy(int p, int halonr)
{
  int j, outputbin;

  /* make explicitly sure that the whole galaxy structure has defined 0 values */
  memset(&Gal[p], 0, sizeof(struct GALAXY));

  Gal[p].NextGalaxy = -1;
#ifdef GALAXYTREE
  Gal[p].FirstProgGal = -1;
#endif

  if(halonr != Halo[halonr].FirstHaloInFOFgroup)
    {
      terminate("Hah?\n");
    }

  Gal[p].Type = 0;

  Gal[p].HaloNr = halonr;
  Gal[p].MostBoundID = Halo[halonr].MostBoundID;
  Gal[p].SnapNum = Halo[halonr].SnapNum - 1;
#ifdef HALOPROPERTIES
  Gal[p].HaloM_Mean200 = Halo[halonr].M_Mean200;
  Gal[p].HaloM_Crit200 = Halo[halonr].M_Crit200;
  Gal[p].HaloM_TopHat = Halo[halonr].M_TopHat;
  Gal[p].HaloVelDisp = Halo[halonr].VelDisp;
  Gal[p].HaloVmax = Halo[halonr].Vmax;
#endif

  for(j = 0; j < 3; j++)
    {
      Gal[p].Pos[j] = Halo[halonr].Pos[j];
      Gal[p].Vel[j] = Halo[halonr].Vel[j];
      Gal[p].GasSpin[j] = Halo[halonr].Spin[j];
      Gal[p].StellarSpin[j] = Halo[halonr].Spin[j];
      Gal[p].HaloSpin[j] = Halo[halonr].Spin[j];
      Gal[p].MergCentralPos[j] = Gal[p].Pos[j];
      Gal[p].DistanceToCentralGal[j]=0.0;
#ifdef HALOPROPERTIES
      Gal[p].HaloPos[j] = Halo[halonr].Pos[j];
      Gal[p].HaloVel[j] = Halo[halonr].Vel[j];
#endif
    }

  Gal[p].Len = Halo[halonr].Len;
  Gal[p].Vmax = Halo[halonr].Vmax;
  Gal[p].InfallVmax = Halo[halonr].Vmax;
  Gal[p].InfallVmaxPeak = Gal[p].InfallVmax;
  Gal[p].Vvir = get_virial_velocity(halonr);
  Gal[p].Mvir = get_virial_mass(halonr);
  Gal[p].Rvir = get_virial_radius(halonr);
  Gal[p].MergeSat = 0.0;
  Gal[p].InfallSnap = Halo[halonr].SnapNum;

  Gal[p].ColdGas = 0.0;
  Gal[p].DiskMass = 0.0;
  Gal[p].BulgeMass = 0.0;
  Gal[p].HotGas = 0.0;
  Gal[p].EjectedMass = 0.0;
  Gal[p].ICM = 0.0;
#ifdef TRACK_BURST
  Gal[p].BurstMass = 0.0;
#endif

  Gal[p].BlackHoleMass = 0.0;
  /*ram pressure*/
  Gal[p].HotRadius=Gal[p].Rvir;
#ifdef GALAXYTREE
  Gal[p].DisruptOn = 0;
#endif

#ifdef DETAILED_METALS_AND_MASS_RETURN
  Gal[p].MetalsColdGas = metals_init();
  Gal[p].MetalsDiskMass = metals_init();
  Gal[p].MetalsBulgeMass = metals_init();
  Gal[p].MetalsHotGas = metals_init();
  Gal[p].MetalsEjectedMass = metals_init();
  Gal[p].MetalsICM = metals_init();
#ifdef METALS_SELF
  Gal[p].MetalsHotGasSelf = metals_init();
#endif
#else
  Gal[p].MetalsColdGas = 0.0;
  Gal[p].MetalsDiskMass = 0.0;
  Gal[p].MetalsBulgeMass = 0.0;
  Gal[p].MetalsHotGas = 0.0;
  Gal[p].MetalsEjectedMass = 0.0;
  Gal[p].MetalsICM = 0.0;
#ifdef METALS_SELF
  Gal[p].MetalsHotGasSelf = 0.0;
#endif
#endif //DETAILED_METALS_AND_MASS_RETURN

  //inclination defined as the angle between galaxy spin and the z-axis
  Gal[p].CosInclination = 0.0;

  Gal[p].PrimordialAccretionRate = 0.0;
  Gal[p].CoolingRate = 0.0;
  Gal[p].CoolingRate_beforeAGN = 0.0;
  Gal[p].CoolingRadius = 0.0;
  Gal[p].CoolingGas = 0.0;
  Gal[p].QuasarAccretionRate=0.0;
  Gal[p].RadioAccretionRate=0.0;
  Gal[p].AGNheatingFromCentral = 0.0;

  Gal[p].Sfr = 0.0;
  Gal[p].SfrBulge = 0.0;

  Gal[p].StarMerge=0.0;

  Gal[p].XrayLum = 0.0;
  Gal[p].GasDiskRadius = get_initial_disk_radius(halonr, p);
  Gal[p].StellarDiskRadius = Gal[p].GasDiskRadius;
  Gal[p].BulgeSize = 0.0;

  Gal[p].OriMergTime = 0.0;
  Gal[p].MergTime = 0.0;

  for(outputbin = 0; outputbin < NOUT; outputbin++)
    Gal[p].MassWeightAge[outputbin] = 0.0;

#ifndef  POST_PROCESS_MAGS
#ifdef OUTPUT_REST_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++) {
    for(j = 0; j < NMAG; j++) {
      Gal[p].Lum[j][outputbin]         = 0.0;
      Gal[p].YLum[j][outputbin]        = 0.0;
      Gal[p].LumBulge[j][outputbin]    = 0.0;
      Gal[p].YLumBulge[j][outputbin]   = 0.0;
      Gal[p].LumDust[j][outputbin]     = 0.0;
#ifdef ICL
      Gal[p].ICLLum[j][outputbin]      = 0.0;
#endif
    }
  }
#endif
#ifdef COMPUTE_OBS_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++) {
    for(j = 0; j < NMAG; j++) {
      Gal[p].ObsLum[j][outputbin]        = 0.0;
      Gal[p].ObsYLum[j][outputbin]       = 0.0;
      Gal[p].ObsLumBulge[j][outputbin]   = 0.0;
      Gal[p].ObsYLumBulge[j][outputbin]  = 0.0;
      Gal[p].ObsLumDust[j][outputbin]    = 0.0;
#ifdef ICL
      Gal[p].ObsICL[j][outputbin]        = 0.0;
#endif
	  
#ifdef OUTPUT_MOMAF_INPUTS
      Gal[p].dObsLum[j][outputbin]       = 0.0;
      Gal[p].dObsYLum[j][outputbin]      = 0.0;
      Gal[p].dObsLumBulge[j][outputbin]  = 0.0;
      Gal[p].dObsYLumBulge[j][outputbin] = 0.0;
      Gal[p].dObsLumDust[j][outputbin]   = 0.0;
#ifdef ICL
      Gal[p].dObsICL[j][outputbin]        = 0.0;
#endif
#endif
    }
  }
#endif
#endif //POST_PROCESS_MAGS

#ifdef GALAXYTREE
  Gal[p].FirstProgGal = -1;
#endif

#ifdef STAR_FORMATION_HISTORY
  sfh_initialise(p);
#endif //STAR_FORMATION_HISTORY

#ifdef INDIVIDUAL_ELEMENTS
  Gal[p].ColdGas_elements = elements_init();
  Gal[p].DiskMass_elements = elements_init();
  Gal[p].BulgeMass_elements = elements_init();
  Gal[p].HotGas_elements = elements_init();
  Gal[p].EjectedMass_elements = elements_init();
  Gal[p].ICM_elements = elements_init();
#endif
}

/**@brief Whenever star formation occurs, calculates the luminosity corresponding
  *        to the mass of stars formed, considering the metallicity and age of the
  *        material.
  *
  * The semi-analytic code uses look up tables produced by Evolutionary Population
  * Synthesis Models to convert the mass formed on every star formation episode
  * into a luminosity. Each of These tables corresponds to a simple stellar
  * population i.e, a population with a single metallicity. For a given IMF,
  * metatillicty and age, the tables give the luminosity for a
  * \f$ 10^{11}M_\odot\f$ burst. The default model uses a Chabrier IMF and
  * stellar populations from Bruzual & Charlot 2003 with 6 different metallicites.
  *
  * The magnitudes are immediately calculated for each output bin, so that we know
  * the age of each population that contributed to a galaxy total population: the
  * age between creation and output. Apart from the different ages of the populations
  * at a given output bin, if the option COMPUTE_OBS_MAGS is turned on, then we also
  * need to know the K-corrections (going the opposite directions as in observations)
  * that will affect each population.
  *
  * For each metallicity there is a look up table which has the different magnitudes
  * for each age and then this is k-corrected to all the snapshots.
  *
  * If MetallicityOption = 0 -> only solar metallicity.
  * If MetallicityOption = 1 -> 6 metallicities.
  * */
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef  POST_PROCESS_MAGS
void add_to_luminosities(int p, double mstars, double time, double dt, double metallicity)
{
  int outputbin, metindex, tabindex, j;
  double f1, f2, fmet1, fmet2, LuminosityToAdd, dLuminosityToAdd;
  double X1, age, tbc;
 	int N_AgeBins=1, ii;
  double upper_time;

  /* Time bellow which the luminosities are corrected for extinction due to
   * molecular birth clouds.  */
   tbc = 10.0 / UnitTime_in_Megayears * Hubble_h;


  /* mstars converted from 1.e10Msun/h to 1.e11 Msun */
  X1 = mstars/N_AgeBins * 0.1 / Hubble_h;

  /* now we have to change the luminosities accordingly. */
  /* note: we already know at which place we have to look up the tables,
   * since we know the output times, the current time and the metallicity.
   * find_interpolated_lum() finds the 2 closest points in the SPS table
   * in terms of age and metallicity. Time gives the time_to_present for
   * the current step while NumToTime(ListOutputSnaps[outputbin]) gives
   * the time of the output snap - units Mpc/Km/s/h */
  upper_time=time+dt/2.;

  //if one wants to have finner bins for the star formation then the STEPS
  //of the calculation. if N_AgeBins=1 it doesn't do anything
  for(ii=0;ii<N_AgeBins;ii++)
  {
  	time=upper_time-ii*dt/((float)N_AgeBins)-dt/((float)N_AgeBins)/2.;

#ifdef OUTPUT_REST_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++)
    {
      find_interpolated_lum(time, NumToTime(ListOutputSnaps[outputbin]), metallicity,
			    &metindex, &tabindex, &f1, &f2, &fmet1, &fmet2);

      if(MetallicityOption == 0)
	metindex = 4;		// reset met index to use only solar metallicity

      age = time - NumToTime(ListOutputSnaps[outputbin]);
      /* For rest-frame, there is no K-correction on magnitudes,
       * hence the 0 in LumTables[j][metindex][0][tabindex] */
      for(j = 0; j < NMAG; j++)
        {
    	  //interpolation between the points found by find_interpolated_lum
    	  LuminosityToAdd = X1 * (fmet1 * (f1 * LumTables[j][metindex][0][tabindex] +
    			                   f2 * LumTables[j][metindex][0][tabindex + 1]) +
    			          fmet2 * (f1 * LumTables[j][metindex + 1][0][tabindex] +
					   f2 * LumTables[j][metindex + 1][0][tabindex + 1]));
    	  Gal[p].Lum[j][outputbin] += LuminosityToAdd;

    	  /*luminosity used for extinction due to young birth clouds */
    	  if(age <= tbc)
	        Gal[p].YLum[j][outputbin] += LuminosityToAdd;
        }

    }
#endif //OUTPUT_REST_MAGS

#ifdef COMPUTE_OBS_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++)
    {
      find_interpolated_lum(time, NumToTime(ListOutputSnaps[outputbin]), metallicity,
			    &metindex, &tabindex, &f1, &f2, &fmet1, &fmet2);

      if(MetallicityOption == 0)
	    metindex = 4;		// reset met index to use only solar metallicity

      int zindex = ((LastDarkMatterSnapShot+1) - 1) - ListOutputSnaps[outputbin];

      age = time - NumToTime(ListOutputSnaps[outputbin]);

      /* Note the zindex in LumTables[][][][] meaning the magnitudes are now
       * "inversely k-corrected to get observed frame at output bins" */
      for(j = 0; j < NMAG; j++)
        {
    	  //interpolation between the points found by find_interpolated_lum
    	  LuminosityToAdd = X1 * (fmet1 * (f1 * LumTables[j][metindex][zindex][tabindex] +
    			                   f2 * LumTables[j][metindex][zindex][tabindex + 1]) +
			          fmet2 * (f1 * LumTables[j][metindex + 1][zindex][tabindex] +
			                   f2 * LumTables[j][metindex + 1][zindex][tabindex + 1]));
    	  Gal[p].ObsLum[j][outputbin] += LuminosityToAdd;

#ifdef OUTPUT_MOMAF_INPUTS
    	  dLuminosityToAdd = X1 * (fmet1 * (f1 * LumTables[j][metindex][zindex + 1][tabindex] +
				            f2 * LumTables[j][metindex][zindex + 1][tabindex + 1]) +
			           fmet2 * (f1 * LumTables[j][metindex + 1][zindex + 1][tabindex] +
				            f2 * LumTables[j][metindex + 1][zindex + 1][tabindex + 1]));
    	  Gal[p].dObsLum[j][outputbin] += dLuminosityToAdd;
#endif

    	  /*luminosity used for extinction due to young birth clouds */
    	  if(age <= tbc)
    	    {
    	      Gal[p].ObsYLum[j][outputbin] += LuminosityToAdd;
#ifdef OUTPUT_MOMAF_INPUTS
    	      Gal[p].dObsYLum[j][outputbin] += dLuminosityToAdd;
#endif
    	    }

        }
    }
#endif //COMPUTE_OBS_MAGS

  }//end loop on small age bins

}
#endif  //POST_PROCESS_MAGS
#endif  //COMPUTE_SPECPHOT_PROPERTIES



/**@brief gives the time from a given snapshot to z=0 (time in code_units/h).*/
double NumToTime(int num)
{
  return Age[num];
}


/**@brief Calculates the virial mass: \f$M_{\rm{crit200}}\f$ for central halos
 *        with \f$M_{\rm{crit200}}\f$ or Len*PartMass for central halos without. */

double get_virial_mass(int halonr)
{
  if(halonr == Halo[halonr].FirstHaloInFOFgroup && Halo[halonr].M_Crit200)
    return Halo[halonr].M_Crit200;	/* take spherical overdensity mass estimate */
  else
    return Halo[halonr].Len * PartMass;
}


/**@brief Calculates the virial velocity from the virial mass and radius.
 *
 * Calculates virial velocity:
 *    \f$ V_{\rm{vir}}=\frac{GM_{\rm{vir}}}{R_{\rm{vir}}} \f$*/

double get_virial_velocity(int halonr)
{
  return sqrt(G * get_virial_mass(halonr) / get_virial_radius(halonr));
}


double hubble_of_z(int halonr)
{
  double zplus1;

  zplus1 = 1 + ZZ[Halo[halonr].SnapNum];

  /*get H for current z*/
  return Hubble * sqrt(Omega * zplus1 * zplus1 * zplus1 + (1 - Omega - OmegaLambda) * zplus1 * zplus1 +
		       OmegaLambda);
}

/**@brief Calculates virial radius from a critical overdensity
 *
 * Calculates virial radius using:
 * \f$ M_{\rm{vir}}=\frac{4}{3}\pi R_{\rm{vir}}^3 \rho_c \Delta_c\f$.
 *
 * From which, assuming \f$ \Delta_c=200\f$, *
 * \f$ R_{\rm{vir}}=\left( \frac{3M_{\rm{vir}}}{4\pi 200 \rho_c}\right)^{1/3}\f$
 */
double get_virial_radius(int halonr)
{
  double hubble_z, rhocrit, fac;

  /*get H for current z*/
  hubble_z = hubble_of_z(halonr);

  rhocrit = 3 * hubble_z * hubble_z / (8 * M_PI * G);
  fac = 1 / (200 * 4 * M_PI / 3.0 * rhocrit); 
  return pow(get_virial_mass(halonr) * fac, 1.0 / 3);
}


/**@brief Converts luminosities into magnitudes
 *
 * Converts luminosities into magnitudes:
 * \f$ M=-2.5\mathrm{log}_{10}(L) \f$ */
double lum_to_mag(double lum)
{
  if(lum > 0)
    return -2.5 * log10(lum);
  else
    return 99.0;
}



/**@brief Updates properties of central galaxies.
 *
 *   \f$M_{\rm{vir}}\f$, \f$R_{\rm{vir}}\f$ and \f$V_{\rm{vir}}\f$ are only
 *   updated for type 0's. Once galaxies become satellites these quantities
 *   stay unchanged, so will be the values at infall.
 *
 *   If type = 0 then the HotRadius is the Viral Radius, which will be used in
 *   the cooling recipe.
 *
 *   Other infall information will not be used for central galaxies so we do not
 *   care whether they carry the correct values. */
void update_centralgal(int ngal,int halonr)
{
  int j;
  Gal[ngal].Type = 0;
 
  Gal[ngal].InfallVmax = Halo[halonr].Vmax;
  if(Gal[ngal].InfallVmaxPeak < Gal[ngal].InfallVmax)
  	Gal[ngal].InfallVmaxPeak = Gal[ngal].InfallVmax;
  Gal[ngal].Rvir = get_virial_radius(halonr);
  Gal[ngal].Vvir = get_virial_velocity(halonr);
  Gal[ngal].Mvir = get_virial_mass(halonr);
  Gal[ngal].InfallSnap = Halo[halonr].SnapNum;
  Gal[ngal].InfallHotGas=Gal[ngal].HotGas;
  //Gal[ngal].InfallHotGasRadius=Gal[ngal].Rvir;
  
  /* if type =0 then hotradius =viral radius, this will be used in the cooling recipe; */
  Gal[ngal].HotRadius = Gal[ngal].Rvir;
  Gal[ngal].MergeOn= 0;
  for (j=0;j<3;j++)
  Gal[ngal].HaloSpin[j] = Halo[halonr].Spin[j];
}


/**@brief Updates properties of type 1 galaxies.
 *
 * A dynamical friction decay time scale is calculated
 * for type 1's (as is done for type 2 - introduced for millennium II where the
 * increased resolution means type 1 always retain some dark matter and orbit
 * around for a long time). This is only calculated when the baryonic mass of
 * the type 1 becomes larger than its dark matter mass. The code finds the type
 * 0 to which this galaxy should merge and then sets up the merging clock.
 * */
void update_type_1(int ngal, int halonr, int prog)
{
  int current,descendant,firstdes;

  Gal[ngal].Type = 1;

  if(Gal[ngal].MergeOn == 0)
  {
    /*If baryonic mass > dark matter mass*/
    if (Gal[ngal].ColdGas+Gal[ngal].DiskMass+Gal[ngal].BulgeMass > Halo[halonr].Len*PartMass)
    {

    	current = halonr;
      descendant = Halo[halonr].Descendant;
      firstdes = Halo[Halo[halonr].FirstHaloInFOFgroup].Descendant;

      /* In case this is the last snapnum (firstdes == -1), it means that we tracked all
       * the way down to redshift =0 and mergeon should be trun on. Otherwise, it is the
       * case that the current halo and the corresponding fof central subhalo are
       * "mysteriously" lost in the dark matter simulation at an intermediate redshift
       * and this galaxy would not be treated further anyway further. Thus the mergeon
       * value is irrelevant. Here mergeon is set to 1. */
      if(descendant == -1)
	    Gal[ngal].MergeOn = 1;

	  /* checks when the galaxy "disappears" (when it merges with the type 0) in order to get
	   * the type 0 ID into which this type 1 will be merged. */
      while(descendant >= 0 && firstdes >= 0)
      {
      	if(firstdes != Halo[firstdes].FirstHaloInFOFgroup)
      		break;

      	if(Halo[descendant].FirstHaloInFOFgroup != Halo[firstdes].FirstHaloInFOFgroup)
      		break;

      	if(descendant != Halo[descendant].FirstHaloInFOFgroup && current == Halo[descendant].FirstProgenitor)
      		if(Gal[ngal].ColdGas + Gal[ngal].DiskMass + Gal[ngal].BulgeMass < Halo[descendant].Len * PartMass)
      			break;


      	if(descendant == Halo[descendant].FirstHaloInFOFgroup && current == Halo[descendant].FirstProgenitor)
      		break;

      	if(descendant == Halo[descendant].FirstHaloInFOFgroup && current != Halo[descendant].FirstProgenitor)
      	{
      		Gal[ngal].MergeOn = 1;
      		break;
      	}

      	if(descendant != Halo[descendant].FirstHaloInFOFgroup && current != Halo[descendant].FirstProgenitor)
      	  break;
	  
      	current=descendant;
      	firstdes = Halo[firstdes].Descendant;
      	descendant=Halo[descendant].Descendant;
	      
      	/* In case this is the last snapnum (firstdes == -1), it means that we tracked all
      	 * the way down to redshift =0 and mergeon should be trun on. Otherwise, it is the
      	 * case that the current halo and the corresponding fof central subhalo are
      	 * "mysteriously" lost in the dark matter simulation at an intermediate redshift
      	 * and this galaxy would not be treated further anyway further. Thus the mergeon
      	 * value is irrelevant. Here mergeon is set to 1. */
      	if(firstdes == -1)
      	  {
      	    if (descendant == -1)
      	      Gal[ngal].MergeOn = 1;
      	    break;
      	  }
      }
	  
   
   
      /*Sets up the dynamical friction decay merging clock as for type 2 galaxies. */
      if (descendant < 0 || Gal[ngal].MergeOn == 1)
      {
      	Gal[ngal].MergeOn = 1;
      	//In case central galaxy has no progenitor
      	if (Halo[Halo[halonr].FirstHaloInFOFgroup].FirstProgenitor == -1 )
      		Gal[ngal].MergTime = estimate_merging_time(prog,Halo[halonr].FirstHaloInFOFgroup,ngal);
      	else
      		Gal[ngal].MergTime = estimate_merging_time(prog,Halo[Halo[halonr].FirstHaloInFOFgroup].FirstProgenitor,ngal);
      	Gal[ngal].MergTime -= NumToTime(Halo[halonr].SnapNum) - NumToTime(Halo[prog].SnapNum);
	  	 //to calculate the position of type 2
      	Gal[ngal].OriMergTime=Gal[ngal].MergTime;
      	Gal[ngal].OriMvir = get_virial_mass(prog);
      	Gal[ngal].OriRvir = get_virial_radius(prog);
      }
    }
  }
  
  /*Mvir, Rvir and Vvir keep their value fixed after infall*/
}	     

  
/**@brief Updates properties of type 2 galaxies.
 *
 *  Sets Hot Radius to 0, since all the hot gas has been stripped.
 *  Calls estimate_merging_time to get the merging time scale, calculated for
 *  the orbital decay due to dynamical friction, since this galaxy has lost its
 *  dark matter halo and its position cannot be tracked. */
void update_type_2(int ngal,int halonr, int prog,int mostmassive)
{

 if(Gal[ngal].Type != 2)
   {
     int j;
     for(j=0; j<3; j++)
       {
	 Gal[ngal].Pos_notupdated[j] = Gal[ngal].Pos[j];
	 Gal[ngal].Vel_notupdated[j] = Gal[ngal].Vel[j];
       }
   }

  Gal[ngal].Type = 2;

  Gal[ngal].HotRadius = 0.0;

  /* Estimate remaining merging timescale. */
  if (Gal[ngal].MergeOn == 0)
    {
      //if central galaxy has no progenitor
      if (mostmassive == -1)
	mostmassive = halonr;

      Gal[ngal].MergTime = estimate_merging_time(prog,mostmassive,ngal);
      Gal[ngal].MergTime -= NumToTime(Halo[halonr].SnapNum) - NumToTime(Halo[prog].SnapNum);
      //to calculate the position of type 2
      Gal[ngal].OriMergTime=Gal[ngal].MergTime;
      Gal[ngal].OriMvir = get_virial_mass(prog);
      Gal[ngal].OriRvir = get_virial_radius(prog);
    }
}

void transfer_stars(int p, char cp[], int q, char cq[], double fraction) {

  /* Transfers a fraction of component cq of galaxy q onto component cp of galaxy p.
   * cp and cq must each be one of:
   *   Disk
   *   Bulge
   *   ICM
   * If -DTRACK_BURST is set then can also specify BurstMass as an option.  This is
   * a little different in that it is not a separate component, so it should only
   * be transferred if both cq and cp are BurstMass
   *
   */

  float Mass;
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals Metals;
#ifdef INDIVIDUAL_ELEMENTS
  struct elements Yield;
  struct elements sfh_Elements[SFH_NBIN];
#endif
#else
  float Metals;
#endif //DETAILED_METALS_AND_MASS_RETURN

#ifdef STAR_FORMATION_HISTORY
  int i;
  float sfh_Mass[SFH_NBIN];
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals sfh_Metals[SFH_NBIN];
#else
  float sfh_Metals[SFH_NBIN];
#endif
#endif

  /* Sanity checks */
  if (fraction > 1.) {
    printf("\n*** transfer_stars: fraction>1 ***\n");
    exit(1);
  }
#ifdef STAR_FORMATION_HISTORY
  if (Gal[p].sfh_ibin != Gal[q].sfh_ibin) {
    printf("\n*** transfer_stars: inconsistent itimes ***\n");
    for(i=0;i<SFH_NBIN;i++)
      printf("Bin[%d] time_1=%e dt_1=%e Nbins_1=%d time_2=%e dt_2=%e Nbins_2=%d\n",i,
      		  Gal[p].sfh_t[i],Gal[p].sfh_dt[i],Gal[p].sfh_Nbins[i],
      		  Gal[q].sfh_t[i],Gal[q].sfh_dt[i],Gal[q].sfh_Nbins[i]);
    exit(1);
  }
#endif
#ifdef TRACK_BURST
  if ((strcmp(cp,"BurstMass")==0 && !strcmp(cq,"BurstMass")==0) ||
      (strcmp(cq,"BurstMass")==0 && !strcmp(cp,"BurstMass")==0)) {
    printf("\n*** transfer_stars: used incorrectly with BurstMass ***\n");
    exit(1);
  }
#endif

  //Mass and metals to be transfered
  if (strcmp(cq,"Disk")==0)
  {
    Mass = fraction*Gal[q].DiskMass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++)
    	sfh_Mass[i]=fraction*Gal[q].sfh_DiskMass[i];
#endif
    Metals=metals_add(metals_init(),Gal[q].MetalsDiskMass,fraction);
#ifdef INDIVIDUAL_ELEMENTS
    Yield = elements_add(elements_init(),Gal[q].DiskMass_elements,fraction);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++)
      {
	sfh_Metals[i]=metals_add(metals_init(),Gal[q].sfh_MetalsDiskMass[i],fraction);
#ifdef INDIVIDUAL_ELEMENTS
	sfh_Elements[i]=elements_add(elements_init(),Gal[q].sfh_ElementsDiskMass[i],fraction);
#endif
      }
#endif  
  }
  else if (strcmp(cq,"Bulge")==0)
  {
    Mass=fraction*Gal[q].BulgeMass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++) sfh_Mass[i]=fraction*Gal[q].sfh_BulgeMass[i];
#endif
    Metals=metals_add(metals_init(),Gal[q].MetalsBulgeMass,fraction);
#ifdef INDIVIDUAL_ELEMENTS
    Yield = elements_add(elements_init(),Gal[q].BulgeMass_elements,fraction);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++)
      {
	sfh_Metals[i]=metals_add(metals_init(),Gal[q].sfh_MetalsBulgeMass[i],fraction);
#ifdef INDIVIDUAL_ELEMENTS
	sfh_Elements[i]=elements_add(elements_init(),Gal[q].sfh_ElementsBulgeMass[i],fraction);
#endif
      }
#endif  
  }
  else if (strcmp(cq,"ICM")==0)
  {
    Mass=fraction*Gal[q].ICM;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++) sfh_Mass[i]=fraction*Gal[q].sfh_ICM[i];
#endif
    Metals=metals_add(metals_init(),Gal[q].MetalsICM,fraction);
#ifdef INDIVIDUAL_ELEMENTS
    Yield = elements_add(elements_init(),Gal[q].ICM_elements,fraction);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++)
      {
	sfh_Metals[i]=metals_add(metals_init(),Gal[q].sfh_MetalsICM[i],fraction);
#ifdef INDIVIDUAL_ELEMENTS
	sfh_Elements[i]=elements_add(elements_init(),Gal[q].sfh_ElementsICM[i],fraction);
#endif
      }
#endif  
#ifdef TRACK_BURST
  }
  else if (strcmp(cq,"Burst")==0)
  {
    Mass = fraction*Gal[q].BurstMass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++) sfh_Mass[i]=fraction*Gal[q].sfh_BurstMass[i];
#endif
#endif
  }
  else
  {
    printf("Unknown component type %s in call to transfer_stars\n",cq);
    exit(1);
  }

  //Add to galaxy p
  if (strcmp(cp,"Disk")==0)
  {
    Gal[p].DiskMass += Mass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[p].sfh_ibin; i++) Gal[p].sfh_DiskMass[i] += sfh_Mass[i];
#endif
    Gal[p].MetalsDiskMass=metals_add(Gal[p].MetalsDiskMass,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
    Gal[p].DiskMass_elements = elements_add(Gal[p].DiskMass_elements,Yield,1.);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[p].sfh_ibin; i++)
      {
	Gal[p].sfh_MetalsDiskMass[i]=metals_add(Gal[p].sfh_MetalsDiskMass[i],sfh_Metals[i],1.);
#ifdef INDIVIDUAL_ELEMENTS
	Gal[p].sfh_ElementsDiskMass[i]=elements_add(Gal[p].sfh_ElementsDiskMass[i],sfh_Elements[i],1.);
#endif
      }
#endif  
  }
  else if (strcmp(cp,"Bulge")==0)
  {
    Gal[p].BulgeMass += Mass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[p].sfh_ibin; i++) Gal[p].sfh_BulgeMass[i] += sfh_Mass[i];
#endif
    Gal[p].MetalsBulgeMass=metals_add(Gal[p].MetalsBulgeMass,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
    Gal[p].BulgeMass_elements = elements_add(Gal[p].BulgeMass_elements,Yield,1.);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[p].sfh_ibin; i++)
      {
	Gal[p].sfh_MetalsBulgeMass[i]=metals_add(Gal[p].sfh_MetalsBulgeMass[i],sfh_Metals[i],1.);
#ifdef INDIVIDUAL_ELEMENTS
	Gal[p].sfh_ElementsBulgeMass[i]=elements_add(Gal[p].sfh_ElementsBulgeMass[i],sfh_Elements[i],1.);
#endif
      }
#endif  
  }
  else if (strcmp(cp,"ICM")==0)
  {
    Gal[p].ICM += Mass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[p].sfh_ibin; i++) Gal[p].sfh_ICM[i] += sfh_Mass[i];
#endif
    Gal[p].MetalsICM=metals_add(Gal[p].MetalsICM,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
    Gal[p].ICM_elements = elements_add(Gal[p].ICM_elements,Yield,1.);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[p].sfh_ibin; i++)
      {
	Gal[p].sfh_MetalsICM[i]=metals_add(Gal[p].sfh_MetalsICM[i],sfh_Metals[i],1.);
#ifdef INDIVIDUAL_ELEMENTS
	Gal[p].sfh_ElementsICM[i]=elements_add(Gal[p].sfh_ElementsICM[i],sfh_Elements[i],1.);
#endif
      }
#endif  
#ifdef TRACK_BURST
  }
  else if (strcmp(cp,"Burst")==0)
  {
    Gal[p].BurstMass += Mass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[p].sfh_ibin; i++) Gal[p].sfh_BurstMass[i] += sfh_Mass[i];
#endif
#endif
  }
  else
  {
    printf("Unknown component type %s in call to transfer_stars\n",cp);
    exit(1);
  }

  //Subtract from galaxy q; 
  if (strcmp(cq,"Disk")==0)
  {
    Gal[q].DiskMass -= Mass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++) Gal[q].sfh_DiskMass[i] -= sfh_Mass[i];
#endif
    Gal[q].MetalsDiskMass=metals_add(Gal[q].MetalsDiskMass,Metals,-1.);
#ifdef INDIVIDUAL_ELEMENTS
    Gal[q].DiskMass_elements = elements_add(Gal[q].DiskMass_elements,Yield,-1.);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++) 
      {
    	Gal[q].sfh_MetalsDiskMass[i]=metals_add(Gal[q].sfh_MetalsDiskMass[i],sfh_Metals[i],-1.);
#ifdef INDIVIDUAL_ELEMENTS
    	Gal[q].sfh_ElementsDiskMass[i]=elements_add(Gal[q].sfh_ElementsDiskMass[i],sfh_Elements[i],-1.);
#endif
      }
#endif  
  }
  else if (strcmp(cq,"Bulge")==0)
  {
    Gal[q].BulgeMass -= Mass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++) Gal[q].sfh_BulgeMass[i] -= sfh_Mass[i];
#endif
    Gal[q].MetalsBulgeMass=metals_add(Gal[q].MetalsBulgeMass,Metals,-1.);
#ifdef INDIVIDUAL_ELEMENTS
    Gal[q].BulgeMass_elements = elements_add(Gal[q].BulgeMass_elements,Yield,-1.);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++)
      {
	Gal[q].sfh_MetalsBulgeMass[i]=metals_add(Gal[q].sfh_MetalsBulgeMass[i],sfh_Metals[i],-1.);
#ifdef INDIVIDUAL_ELEMENTS
	Gal[q].sfh_ElementsBulgeMass[i]=elements_add(Gal[q].sfh_ElementsBulgeMass[i],sfh_Elements[i],-1.);
#endif
      }
#endif  
  }
  else if (strcmp(cq,"ICM")==0)
  {
    Gal[q].ICM -= Mass;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++) Gal[q].sfh_ICM[i] -= sfh_Mass[i];
#endif
    Gal[q].MetalsICM=metals_add(Gal[q].MetalsICM,Metals,-1.);
#ifdef INDIVIDUAL_ELEMENTS
    Gal[q].ICM_elements = elements_add(Gal[q].ICM_elements,Yield,-1.);
#endif
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++)
      {
	Gal[q].sfh_MetalsICM[i]=metals_add(Gal[q].sfh_MetalsICM[i],sfh_Metals[i],-1.);
#ifdef INDIVIDUAL_ELEMENTS
	Gal[q].sfh_ElementsICM[i]=elements_add(Gal[q].sfh_ElementsICM[i],sfh_Elements[i],-1.);
#endif
      }
#endif  
#ifdef TRACK_BURST
  }
  else if (strcmp(cq,"Burst")==0)
  {
    Gal[q].BurstMass -=0.;
#ifdef STAR_FORMATION_HISTORY
    for (i=0; i<=Gal[q].sfh_ibin; i++) Gal[q].sfh_BurstMass[i] -= sfh_Mass[i];
#endif
#endif
  } else
  {
    printf("Unknown component type %s in call to transfer_stars\n",cq);
    exit(1);
  }

  return;
}


void transfer_gas(int p, char cp[], int q, char cq[], double fraction, char call_function[], int call_line) {

  /* Transfers a fraction of component cq of galaxy q onto component cp of galaxy p.
   * cp and cq must each be one of:
   *   Cold
   *   Hot
   *   Ejected
   */
  float Mass;
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals Metals;
#ifdef INDIVIDUAL_ELEMENTS
  struct elements Yield;
#endif
#else
  float Metals;
#endif

  /* Sanity check */

  if (fraction > 1.) //1.000001
  {
  	char sbuf[1000];
  	sprintf(sbuf, "\nparent call from: %s, line %d \ntransfer_gas: fraction>1\nfraction = %.11f\nFrom '%s' to '%s\n",call_function, call_line, fraction,cq, cp);
  	terminate(sbuf);
  }
  //if (fraction > 1. && fraction < 1.000001) { //1.0000001
	//  fraction = 1.0;
    /*printf("\n*** fraction forced to 1.0 ***\n");
    printf("*** fraction was = %.11f ***\n", fraction);
    printf("*** From '%s' to '%s' ***\n\n", cq, cp);*/
 // }

  //Mass and Metals to be transfered
  if (strcmp(cq,"Cold")==0)
    {
      Mass = fraction*Gal[q].ColdGas;
      Metals = metals_add(metals_init(),Gal[q].MetalsColdGas,fraction);
#ifdef INDIVIDUAL_ELEMENTS
      Yield = elements_add(elements_init(),Gal[q].ColdGas_elements,fraction);
#endif
    }
  else if (strcmp(cq,"Hot")==0)
    {
      Mass=fraction*Gal[q].HotGas;
      Metals = metals_add(metals_init(),Gal[q].MetalsHotGas,fraction);
#ifdef INDIVIDUAL_ELEMENTS
    Yield = elements_add(elements_init(),Gal[q].HotGas_elements,fraction);
#endif
    }
  else if (strcmp(cq,"Ejected")==0)
    {
      Mass=fraction*Gal[q].EjectedMass;
      Metals = metals_add(metals_init(),Gal[q].MetalsEjectedMass,fraction);
#ifdef INDIVIDUAL_ELEMENTS
    Yield = elements_add(elements_init(),Gal[q].EjectedMass_elements,fraction);
#endif
    }
  else
    {
      printf("Unknown component type %s in call to transfer_gas\n",cq);
      exit(1);
    }

  //Add to galaxy p
  if (strcmp(cp,"Cold")==0)
    {
      Gal[p].ColdGas += Mass;
      Gal[p].MetalsColdGas = metals_add(Gal[p].MetalsColdGas,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
      Gal[p].ColdGas_elements = elements_add(Gal[p].ColdGas_elements,Yield,1.);
#endif
    }
  else if (strcmp(cp,"Hot")==0)
    {
      Gal[p].HotGas += Mass;
      Gal[p].MetalsHotGas = metals_add(Gal[p].MetalsHotGas,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
      Gal[p].HotGas_elements = elements_add(Gal[p].HotGas_elements,Yield,1.);
#endif
#ifdef METALS_SELF
      if (p==q) Gal[p].MetalsHotGasSelf = metals_add(Gal[p].MetalsHotGasSelf,Metals,1.);
#endif
    }
  else if (strcmp(cp,"Ejected")==0)
    {
      Gal[p].EjectedMass += Mass;
      Gal[p].MetalsEjectedMass = metals_add(Gal[p].MetalsEjectedMass,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
      Gal[p].EjectedMass_elements = elements_add(Gal[p].EjectedMass_elements,Yield,1.);
#endif
    }
  else
    {
      printf("Unknown component type %s in call to transfer_gas\n",cp);
      exit(1);
    }

  //Subtract from galaxy q;
  if (strcmp(cq,"Cold")==0)   {
    Gal[q].ColdGas -= Mass;
    Gal[q].MetalsColdGas = metals_add(Gal[q].MetalsColdGas,Metals,-1.);
#ifdef INDIVIDUAL_ELEMENTS
    Gal[q].ColdGas_elements = elements_add(Gal[q].ColdGas_elements,Yield,-1.);
#endif
  } else if (strcmp(cq,"Hot")==0) {
    Gal[q].HotGas -= Mass;
    Gal[q].MetalsHotGas = metals_add(Gal[q].MetalsHotGas,Metals,-1.);
#ifdef INDIVIDUAL_ELEMENTS
    Gal[q].HotGas_elements = elements_add(Gal[q].HotGas_elements,Yield,-1.);
#endif
#ifdef METALS_SELF
    Gal[q].MetalsHotGasSelf = metals_add(Gal[q].MetalsHotGasSelf,Metals,-1.);
#endif
  } else if (strcmp(cq,"Ejected")==0) {
    Gal[q].EjectedMass -= Mass;
    Gal[q].MetalsEjectedMass = metals_add(Gal[q].MetalsEjectedMass,Metals,-1.);
#ifdef INDIVIDUAL_ELEMENTS
    Gal[q].EjectedMass_elements = elements_add(Gal[q].EjectedMass_elements,Yield,-1.);
#endif
  } else {
    printf("Unknown component type %s in call to transfer_gas\n",cq);
    exit(1);
  }

  return;
}



void mass_checks(char string[], int igal) {

  /* Some sanity checks on the masses of different components. 
   * If due to rounding error, then apply a correction;
   * otherwise print error message and exit
   */
	//ROB: Should probably make some of these for the elements

#ifndef MASS_CHECKS
  return;
#endif
  
#ifdef STAR_FORMATION_HISTORY
  int i;
  double sfh_sum;
#endif

  if(Gal[igal].ColdGas < 1.e-8)
    {
      Gal[igal].ColdGas = 0.;
      Gal[igal].MetalsColdGas = metals_init();
    }

  //check if the gas mass is less than 0
  if(Gal[igal].ColdGas < 0.0) {
    if (Gal[igal].ColdGas > -1e-7)
      Gal[igal].ColdGas = 0.;
    else {
      printf("\n*** Mass check error, called from: %s, ColdGas < 0. ***\n",string);
      printf("                ColdGas[%d] = %g\n",igal,Gal[igal].ColdGas);
      terminate("");
    }
  }

  //check if the mass in metals is less than 0
  if(metals_total(Gal[igal].MetalsColdGas) < 0.0) {
    if (metals_total(Gal[igal].MetalsColdGas) > -1e-7)
      Gal[igal].MetalsColdGas = metals_init();
    else {
      printf("\n*** Mass check error, called from: %s, MetalsColdGas < 0. ***\n",string);
      printf("                MetalsColdGas[%d] = %g\n",igal,metals_total(Gal[igal].MetalsColdGas));
      terminate("");
    }
  }

  //check if the mass in metals is greater than the gas mass
  if(metals_total(Gal[igal].MetalsColdGas) > Gal[igal].ColdGas) {
    if (metals_total(Gal[igal].MetalsColdGas) < 1e-7)
      Gal[igal].MetalsColdGas = metals_add(metals_init(),Gal[igal].MetalsColdGas,
					   Gal[igal].ColdGas/metals_total(Gal[igal].MetalsColdGas));
    else {
      printf("\n*** Mass check error, called from: %s, MetalsColdGas > ColdGas ***\n",string);
      printf("          MetalsColdGas[%d] = %g\n",igal,metals_total(Gal[igal].MetalsColdGas));
      printf("                ColdGas[%d] = %g\n",igal,Gal[igal].ColdGas);
      terminate("");
    }
  }

  if(Gal[igal].HotGas < 0.0) {
    if (Gal[igal].HotGas > -1e-7)
      Gal[igal].HotGas = 0.;
    else {
      printf("\n*** Mass check error, called from: %s, HotGas < 0. ***\n",string);
      printf("                HotGas[%d] = %g\n",igal,Gal[igal].HotGas);
      terminate("");
    }
  }

  if(metals_total(Gal[igal].MetalsHotGas) < 0.0) {
    if (metals_total(Gal[igal].MetalsHotGas) > -1e-7)
      Gal[igal].MetalsHotGas = metals_init();
    else {
      printf("\n*** Mass check error, called from: %s, MetalsHotGas < 0. ***\n",string);
      printf("                MetalsHotGas[%d] = %g\n",igal,metals_total(Gal[igal].MetalsHotGas));
      terminate("");
    }
  }

  if(metals_total(Gal[igal].MetalsHotGas) > Gal[igal].HotGas) {
    if (metals_total(Gal[igal].MetalsHotGas) < 1e-7)
      Gal[igal].MetalsHotGas = metals_add(metals_init(),Gal[igal].MetalsHotGas,
					   Gal[igal].HotGas/metals_total(Gal[igal].MetalsHotGas));
   else {
      printf("\n***  Mass check error, called from: %s, MetalsHotGas > HotGas ***\n",string);
      printf("          MetalsHotGas[%d] = %g\n",igal,metals_total(Gal[igal].MetalsHotGas));
      printf("                HotGas[%d] = %g\n",igal,Gal[igal].HotGas);
      printf("          MetalsHotGas[%d] = %.11f\n",igal,metals_total(Gal[igal].MetalsHotGas));
      printf("                HotGas[%d] = %.11f\n",igal,Gal[igal].HotGas);
      printf("             BulgeMass[%d] = %g\n",igal,Gal[igal].BulgeMass);
      printf("           EjectedMass[%d] = %g\n",igal,Gal[igal].EjectedMass);
      printf("                  Snapnum = %i\n",Gal[igal].SnapNum);
      terminate("");
    }
  }

  if(Gal[igal].EjectedMass < 0.0) {
    if (Gal[igal].EjectedMass > -1e-7)
      Gal[igal].EjectedMass = 0.;
    else {
      printf("\n*** Mass check error, called from: %s, EjectedMass < 0. ***\n",string);
      printf("                EjectedMass[%d] = %g\n",igal,Gal[igal].EjectedMass);
      terminate("");
    }
  }

  if(metals_total(Gal[igal].MetalsEjectedMass) < 0.0) {
    if (metals_total(Gal[igal].MetalsEjectedMass) > -1e-7)
      Gal[igal].MetalsEjectedMass = metals_init();
    else {
      printf("\n*** Mass check error, called from: %s, MetalsEjectedMass < 0. ***\n",string);
      printf("                MetalsEjectedMass[%d] = %g\n",igal,metals_total(Gal[igal].MetalsEjectedMass));
      terminate("");
    }
  }

  if(metals_total(Gal[igal].MetalsEjectedMass) > Gal[igal].EjectedMass) {
    if (metals_total(Gal[igal].MetalsEjectedMass) < 1e-7)
      Gal[igal].MetalsEjectedMass = metals_add(metals_init(),Gal[igal].MetalsEjectedMass,
					   Gal[igal].EjectedMass/metals_total(Gal[igal].MetalsEjectedMass));
    else {
      printf("\n*** Mass check error, called from: %s, MetalsEjectedMass > EjectedMass ***\n",string);
      printf("          MetalsEjectedMass[%d] = %g\n",igal,metals_total(Gal[igal].MetalsEjectedMass));
      printf("                EjectedMass[%d] = %g\n",igal,Gal[igal].EjectedMass);
      terminate("");
    }
  }

  if(Gal[igal].DiskMass < 0.0) {
    if (Gal[igal].DiskMass > -1e-7)
      Gal[igal].DiskMass = 0.;
    else {
      printf("\n*** Mass check error, called from: %s, DiskMass < 0. ***\n",string);
      printf("                DiskMass[%d] = %g\n",igal,Gal[igal].DiskMass);
      terminate("");
    }
  }

  if(Gal[igal].BulgeMass < 0.0) {
    if (Gal[igal].BulgeMass > -1e-7)
      Gal[igal].BulgeMass = 0.;
    else {
      printf("\n*** Mass check error, called from: %s, BulgeMass < 0. ***\n",string);
      printf("                BulgeMass[%d] = %g\n",igal,Gal[igal].BulgeMass);
      terminate("");
    }
  }

  if(Gal[igal].ICM < 0.0) {
    if (Gal[igal].ICM > -1e-7)
      Gal[igal].ICM = 0.;
    else {
      printf("\n*** Mass check error, called from: %s, ICM < 0. ***\n",string);
      printf("                ICM[%d] = %g\n",igal,Gal[igal].ICM);
      terminate("");
    }
  }

#ifdef TRACK_BURST
  if(Gal[igal].BurstMass < 0.0) {
    if (Gal[igal].BurstMass > -1e-7)
      Gal[igal].BurstMass = 0.;
    else {
      printf("\n*** Mass check error, called from: %s, BurstMass < 0. ***\n",string);
      printf("                BurstMass[%d] = %g\n",igal,Gal[igal].BurstMass);
      terminate("");
    }
  }
#endif

  /* If DETAILED_METALS_AND_MASS_RETURN, sfh stores accumulation of 'stars', not 'stars-recycFrac'.
   * Therefore, it's sum doesn't equal DiskMass any more.*/
#ifndef DETAILED_METALS_AND_MASS_RETURN
#ifdef STAR_FORMATION_HISTORY
  sfh_sum=-Gal[igal].DiskMass;
  for (i=0; i<=Gal[igal].sfh_ibin; i++) sfh_sum+=Gal[igal].sfh_DiskMass[i];
  if((sfh_sum < -1e-4 && sfh_sum < -1e-4*Gal[igal].DiskMass) ||
      (sfh_sum >  1e-4 && sfh_sum >  1e-4*Gal[igal].DiskMass))
    {
      printf("                     sfh_sum = %g\n",sfh_sum);
      printf("                DiskMass[%d] = %g\n",igal,Gal[igal].DiskMass);
      printf("            sfh_DiskMass[%d] = %g\n",igal,sfh_sum+Gal[igal].DiskMass);
      char sbuf[1000];
      sprintf(sbuf, "\n*** Mass check error, called from: %s, Inconsistent sfh for DiskMass.*** \n",string);
      terminate(sbuf);
    }

  sfh_sum=-Gal[igal].BulgeMass;
  for (i=0; i<=Gal[igal].sfh_ibin; i++) sfh_sum+=Gal[igal].sfh_BulgeMass[i];
  if((sfh_sum < -1e-4 && sfh_sum < -1e-4*Gal[igal].BulgeMass) ||
      (sfh_sum >  1e-4 && sfh_sum >  1e-4*Gal[igal].BulgeMass))
    {
      printf("                     sfh_sum = %g\n",sfh_sum);
      printf("                BulgeMass[%d] = %g\n",igal,Gal[igal].BulgeMass);
      printf("            sfh_BulgeMass[%d] = %g\n",igal,sfh_sum+Gal[igal].BulgeMass);
      char sbuf[1000];
      sprintf(sbuf, "\n*** Mass check error, called from: %s, Inconsistent sfh for BulgeMass. ***\n",string);
      terminate(sbuf);
    }

  sfh_sum=-Gal[igal].ICM;
  for (i=0; i<=Gal[igal].sfh_ibin; i++) sfh_sum+=Gal[igal].sfh_ICM[i];
  if(sfh_sum < -1e-4 || sfh_sum > 1e-4)
    {
      printf("                     sfh_sum = %g\n",sfh_sum);
      printf("                ICM[%d] = %g\n",igal,Gal[igal].ICM);
      printf("            sfh_ICM[%d] = %g\n",igal,sfh_sum+Gal[igal].ICM);
      for (i=0; i<=Gal[igal].sfh_ibin; i++)
        printf("%d %f\n",i,Gal[igal].sfh_ICM[i]);
      char sbuf[1000];
      sprintf(sbuf, "\n*** Mass check error, called from: %s, Inconsistent sfh for ICM. ***\n",string);
      terminate(sbuf);
    }
#endif //STAR_FORMATION_HISTORY
#endif //DETAILED_ENRICHEMENT

  return;
}


double separation_gal(int p, int q) {

  /* Calculates the separation of galaxies p and q, allowing for wrapping */

  int i;
  double sep1,sep2;

  sep2=0.;
  for (i=0; i<3; i++)
    {
      sep1 =  wrap(Gal[p].Pos[i] - Gal[q].Pos[i],BoxSize);
      sep2+=sep1*sep1;
    }
  return sqrt(sep2);
}

double separation_halo(int p, int q) {

  /* Calculates the separation of galaxies p and q, allowing for wrapping */

  int i;
  double sep1,sep2;

  sep2=0.;
  for (i=0; i<3; i++)
    {
	  sep1 =  wrap(Halo[p].Pos[i] - Halo[q].Pos[i],BoxSize);
	  sep2+=sep1*sep1;
    }
  return sqrt(sep2);
}

float get_nr_files_to_process(int ThisTask)
{
  int nfiles, filenr, file;
  time_t start;


  nfiles=0;
  time(&start);


#ifndef MCMC
#ifndef OVERWRITE_OUTPUT
  /* a small delay so that processors dont use the same file */
#ifdef PARALLEL
  time_t current;

  if(ThisTask!=0)
    {
	  do
		  time(&current);
	  while(difftime(current, start) < 10.0);
    }
#endif
#endif
#endif

  if(ThisTask==0)
    {
      for(filenr = FirstFile; filenr <= LastFile; filenr++)
	{
#ifdef SPECIFYFILENR
	  file = ListInputFilrNr[filenr];
#else
	  file=filenr;
#endif

#ifndef OVERWRITE_OUTPUT
	  char buf[1000];
#ifdef GALAXYTREE
	  sprintf(buf, "%s/%s_galtree_%d", FinalOutputDir, FileNameGalaxies, file);
#else
	  sprintf(buf, "%s/%s_z%1.2f_%d", FinalOutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[0]], file);
#endif
	  struct stat filestatus;
	  if(stat(buf, &filestatus) != 0)	// seems to exist
#endif
	    nfiles+=1;
	}
    }
#ifdef PARALLEL
  MPI_Bcast(&nfiles,1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  return nfiles;
}

void assign_files_to_tasks(int *FileToProcess, int *TaskToProcess, int ThisTask, int NTask, int nfiles)
{
  int i,j, filenr, file;

  if(ThisTask==0)
    {
      i=0;
      j=0;
      for(filenr = FirstFile; filenr <= LastFile; filenr++)
	{
#ifdef SPECIFYFILENR
	  file = ListInputFilrNr[filenr];
#else
	  file=filenr;
#endif
#ifndef OVERWRITE_OUTPUT
	  char buf[1000];
#ifdef GALAXYTREE
	  sprintf(buf, "%s/%s_galtree_%d", FinalOutputDir, FileNameGalaxies, file);
#else
	  sprintf(buf, "%s/%s_z%1.2f_%d", FinalOutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[0]], file);
#endif
	  struct stat filestatus;
	  if(stat(buf, &filestatus) != 0)	// doesn't exist
	    {
#endif
	      FileToProcess[i]=file;
#ifdef PARALLEL
	      TaskToProcess[i]=j;
#else
	      TaskToProcess[i]=0;
#endif
	      i+=1;
	      j+=1;
	      if(j==NTask)
		j=0;
#ifndef OVERWRITE_OUTPUT
	    }
#endif
	}
    }
#ifdef PARALLEL
  MPI_Bcast(FileToProcess,sizeof(int) * nfiles, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(TaskToProcess,sizeof(int) * nfiles, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
}


//MATH MISC - PROBABLY SHOULD GO INTO SEPARATE FILE
//Finds interpolation point
//the value j so that xx[j]<x<xx[jj+1]
void locate(double *xx, int n, double x, int *j)
{
  unsigned long ju,jm,jl;
  int ascnd;

  jl=0;
  ju=n+1;
  ascnd=(xx[n] >= xx[1]);

  while (ju-jl > 1)
    {
      jm=(ju+jl) >> 1;
      if (x >= xx[jm] == ascnd)
	jl=jm;
      else
	ju=jm;
    }

  if (x == xx[1]) *j=1;
  else if(x == xx[n]) *j=n-1;
  else *j=jl;

}



//!*********************************************************
//!Simpsons quadratures for the signal
//!********************************************************
//erg/s/A --> erg/s/Hz --> erg/s


double integrate(double *flux, int Grid_Length)
{
  double sum[3], I[3], f[4];
  double integral=0.0;
  int i,k;

  for(i=0;i<3;i++)sum[i]=0.0;
  for(i=0;i<3;i++)I[i]=0.0;
  for(i=0;i<4;i++)f[i]=0.0;

  for(i=0;i<Grid_Length/2-2;i++)
    {
      k=2*i+1;                    //odd indexes
      f[2]=flux[k];
      sum[1]=sum[1]+f[2];
    }
  I[1]=sum[1]*2./3.;

  for(i=0;i<Grid_Length/2-1;i++)
    {
      k=2*i;
      f[3]=flux[k] ;     //even indexes
      sum[2]=sum[2]+f[3];
    }
  I[2]=sum[2]*4./3.;

  f[0]=flux[0];
  f[1]=flux[Grid_Length-1];
  I[0]=(f[0]+f[1])/3.;

  integral=I[0]+I[1]+I[2];

//if(Grid_Length==0)
//  printf("Integral=%e\n",integral);

  return integral;
}



//Find interpolation polynomial
//given xa and ya, returns polynomial y and error dy
void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
  int i,m,ns=1;
  double den, dif, dift, ho, hp;
  double *c,*d;
  double ww;

  dif=fabs(x-xa[1]);
  c=vector(1,n);
  d=vector(1,n);
  for (i=1;i<=n;i++) {
      if ( (dift=fabs(x-xa[i])) < dif) {
	  ns=i;
	  dif=dift;
      }
      c[i]=ya[i];
      d[i]=ya[i];
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) {
      for (i=1;i<=n-m;i++)
	{
	  ho=xa[i]-x;
	  hp=xa[i+m]-x;
	  ww=c[i+1]-d[i];
	  if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
	  den=ww/den;
	  d[i]=hp*den;
	  c[i]=ho*den;
	}
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free_vector(d,1,n);
  free_vector(c,1,n);
}


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

#define NREND 1
#define FREE_ARG char*

double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NREND)*sizeof(double)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NREND;
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
  free((FREE_ARG) (v+nl-NREND));
}

void print_galaxy(char string[], int p, int halonr)
{
  int j;
/*	printf("%s Hnr=%d firstinFOF=%d prog=%d nestprog=%d Descendant=%d gal=%d Type=%d\n",
			string, Gal[p].HaloNr, Halo[halonr].FirstHaloInFOFgroup, Halo[halonr].FirstProgenitor,
			Halo[halonr].NextProgenitor, Halo[halonr].Descendant, p, Gal[p].Type);
	printf("     Mvir=%0.3e Vvir=%0.3e Hot=%0.3e Cold=%0.3e Eject=%0.3e disk=%0.3e bulge=%0.3e  GasDiskRadius=%0.3e\n",
			Gal[p].Mvir*1.e10, Gal[p].Vvir, Gal[p].HotGas*1.e10, Gal[p].ColdGas*1.e10, Gal[p].EjectedMass*1.e10,
			Gal[p].DiskMass*1.e10, Gal[p].BulgeMass*1.e10, Gal[p].GasDiskRadius);

	printf("     HotMetals=%0.3e ColdMetals=%0.3e diskMetals=%0.3e bulgeMetals=%0.3e\n",
			Gal[p].MetalsHotGas, Gal[p].MetalsColdGas,Gal[p].MetalsDiskMass, Gal[p].MetalsBulgeMass);


	printf("     x=%0.3f y=%0.3f z=%0.3f vx=%0.3f vy=%0.3f vz=%0.3f\n",
			Gal[p].Pos[0],Gal[p].Pos[1],Gal[p].Pos[2],Gal[p].Vel[0],Gal[p].Vel[1],Gal[p].Vel[2]);*/

  printf("%s Hnr=%d gal=%d snap=%d\n",string, Gal[p].HaloNr,p, Gal[p].SnapNum);
  printf(" Hot=%0.3e Cold=%0.3e Eject=%0.3e disk=%0.3e bulge=%0.3e  GasDiskRadius=%0.3e StellarDiskRadius=%0.3e BulgeSize=%0.3e\n",
	 Gal[p].HotGas*1.e10, Gal[p].ColdGas*1.e10, Gal[p].EjectedMass*1.e10,
	 Gal[p].DiskMass*1.e10, Gal[p].BulgeMass*1.e10, Gal[p].GasDiskRadius, Gal[p].StellarDiskRadius, Gal[p].BulgeSize);

  if(isnan(Gal[p].GasDiskRadius))
    exit(0);

}
