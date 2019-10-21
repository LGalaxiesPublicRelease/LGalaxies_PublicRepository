#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"
#include "mcmc_vars.h"
#ifdef ALL_SKY_LIGHTCONE
#include "lightcone.h"
#endif

/**@file init.c
 * @brief Sets up some unit conversion variables; converts SN and AGN feedback
 *        variables into internal units; and reads in input tables, including
 *        the desired output snapshots, the photometric and dust tables, the
 *        cooling functions and reionization tables.
 *
 *        <B>set_units()</B> - THIS IS FUNDAMENTAL TO UNDERSTAND THE UNITS IN
 *        THE CODE! sets ups some variables used to convert from physical to
 *        internal units (as UnitDensity_in_cgs). These are obtained from the
 *        three defined in input.par: UnitLength_in_cm, UnitMass_in_g and
 *        UnitVelocity_in_cm_per_s.
 *
 *        <B>read_output_snaps()</B> - reads in the list of output redshifts from
 *        file ./input/desired_output_redshifts.txt and converts them into snapsshots
 *        for the given cosmology.
 *
 *        <B>read_zlist()</B> - reads in 1/(z+1) from FileWithZList defined
 *        in ./input/input.par and creates a table with output
 *        redshift ZZ[] and ages Age[].
 *
 *        <B>read_file_nrs()</B> - Done if SPECIFYFILENR OFF - the dark matter files
 *        to read can be defined in a file, instead of being read sequentially.
 *        These are defined in FileNrDir, in input.par, and read into
 *        ListInputFilrNr[].
 *
 *
 *        <B>read_reionization()</B> - Reads in Reion_z[] and Reion_Mc[] from
 *        ./input/Mc.txt. These are used if UPDATEREIONIZATION ON to get Okamoto(2008)
 *        fitting parameters (Mc), instead of Kravtsov(2004)  for the Gnedin (2000)
 *        reionization formalism.
 *
 *
 *        <B>read_dust_tables()</B> - Reads in the dust extinction for the same bands
 *        read from the spectrophotometric tables both for the inter-galactic medium
 *        (**_Set_Ext_table.dat) and young birth clouds (**_Set_YSExt_table.dat).
 *        Detailed description at recipe_dust.c
 *
 *        <B>read_cooling_functions()</B> - calls the functions that read in the
 *        cooling functions defined in cool_func.c
 *
 *        In init.c, but called from other files, are the function to interpolate
 *        properties into different tables. They are all the same but interpolate
 *        into different properties (have different inputs):
 *
 *
 *        <B>find_interpolated_lum()</B> - Used by add_to_luminosities() to
 *        interpolates into the age and metallicity in the SSP tables.
 *
 *         <B>find_interpolate_reionization()</B> - Called from recipe_infall
 *        interpolates into the Reion_z[], Reion_Mc[] table, giving a value of Mc
 *        for a given redshift.
 *
 *        SuperNovae and AGN feedback parameters are converted into internal units:
 *
 *        \f$ AgnEfficiency = \frac{UnitMass_{\rm{g}}}{1.58e^{-26}UnitTime_{\rm{s}}}\f$
 *
 *        \f$ EnergySNcode = \frac{EnergySN}{UnitEnergy_{\rm{cgs}}} h; \f$
          \f$ EtaSNcode = EtaSN \frac{UnitMass_{\rm{g}}}{M_\odot h}. \f$

 *          */

//Needs to be moved to proto.h
void read_reionization(void);

/**@brief controlling recipe for init.c, calls functions to read in tables and
 *        defines some SN and AGN feedback parameters.
 *
 *        Calls set_units(), read_output_snaps(), read_zlist(), read_recgas(),
 *        read_file_nrs(), read_sfrz(), read_reionization(), read_dust_tables() and
 *        read_cooling_functions().
 *
 *        Converts EnergySN (->EnergySNcode), EtaSN (->EtaSNcode) and AgnEfficiency
 *        into internal units.
 *        */
void init(void)
{
  struct rlimit rlim;

  getrlimit(RLIMIT_CORE, &rlim);
  rlim.rlim_cur = RLIM_INFINITY;
  setrlimit(RLIMIT_CORE, &rlim);

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, 42);	/* start-up seed */

  set_units();
#ifdef GALAXYTREE
  ScaleFactor = pow(2, Hashbits) / BoxSize;
#endif

#ifdef DEBUG
  debug_galaxy_init();
#endif

/*#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN
  EnergySN*=1./RecycleFraction;
  EnergySNcode = EnergySN / UnitEnergy_in_cgs * Hubble_h;

  EnergySNII*=1./RecycleFraction;
  EnergySNIa*=1./RecycleFraction;
  EnergyAGB*=1./RecycleFraction;
  EnergySNIIcode = EnergySNII / UnitEnergy_in_cgs * Hubble_h;
  EnergySNIacode = EnergySNIa / UnitEnergy_in_cgs * Hubble_h;
  EnergyAGBcode = EnergyAGB / UnitEnergy_in_cgs * Hubble_h;
#else*/
  EnergySNcode = EnergySN / UnitEnergy_in_cgs * Hubble_h;
//#endif
  EtaSNcode = EtaSN * (UnitMass_in_g / SOLAR_MASS) / Hubble_h;


  //reads in the redshifts for the used Cosmology
  read_zlist();
  //reads in the redshifts for Original Cosmology
  read_zlist_original_cosm();

  //reads in the desired output snapshots
   read_output_snaps();


#ifdef SPECIFYFILENR
  /* read in the number of the files to be processed */
  read_file_nrs();
#endif

  if(ReionizationModel == 0)
    read_reionization();


  //Values of a for the beginning and end of reionization
  a0 = 1.0 / (1.0 + Reionization_z0);
  ar = 1.0 / (1.0 + Reionization_zr);
  
#ifdef H2_AND_RINGS
  //for(i=0; i<RNUM; i++) RingRadius[i]= pow(1.2,i-3)/1000;
  int i;
  for(i=0; i<RNUM; i++)  
    RingRadius[i]= 0.01*pow(2.0,i+1)/1000.; //The radii of the rings. (in unit of Mpc/h)
  //RingRadius[i]= 0.44*pow(1.5,i+1)/1000.; //The radii of the rings. (in unit of Mpc/h)
    // RingRadius[i]= 0.1*pow(1.7,i+1)/1000.; //The radii of the rings. (in unit of Mpc/h)


  for(i=0; i<RNUM; i++)
    {
      if(i==0)
	RingArea[i] = M_PI*RingRadius[i]*RingRadius[i];
      else
	RingArea[i] = M_PI*(RingRadius[i]*RingRadius[i]-RingRadius[i-1]*RingRadius[i-1]);
      /* with conversion factor from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2) -> *(1.e10/h)/(1.e6/h*1.e6/h)*/
      InverseRingArea[i] = 0.01 * Hubble_h / (RingArea[i]*WARM_PHASE_FACTOR);
    }

#endif

  read_cooling_functions();

//CREATE ARRAYS OF SFH TIME STRUCTURE:
#ifdef  STAR_FORMATION_HISTORY
  create_sfh_bins();
#endif

#ifndef MR_PLUS_MRII //If this option (to run with MCMC) is on, the tables are read later
#ifdef COMPUTE_SPECPHOT_PROPERTIES
//read in photometric tables
#ifdef PHOTTABLES_PRECOMPUTED
#ifdef MRII
  setup_LumTables_precomputed("MRII");
#else
#ifdef CATERPILLAR
  setup_LumTables_precomputed("CATERPILLAR");
#else
  setup_LumTables_precomputed("MR");
#endif
#endif
#endif //PHOTTABLES_PRECOMPUTED

#ifdef SPEC_PHOTABLES_ON_THE_FLY
  setup_Spec_LumTables_onthefly();
#endif
#endif
#endif

//READ IN THE YIELD TABLES, AND FORM NORMALISED YIELD ARRAYS:
#ifdef DETAILED_METALS_AND_MASS_RETURN
  read_yield_tables();
  init_integrated_yields();
  integrate_yields();
#endif

#ifdef H2_AND_RINGS
  if(H2FractionRecipe==0)
    init_H2fraction_KMT09();
  else
    if(H2FractionRecipe==1)
      init_H2fraction_KMT08();

  init_jump_index_H2Fraction();
#endif


#ifdef ALL_SKY_LIGHTCONE
  int nstep, nr, i;
  float previoustime, newtime, deltaT, time;

  for (nr = 0 ; nr < NCONES; nr++)
    lightcone_geometry(nr);
  read_distance_table();

  for(i = 1; i < Zlistlen; i++)
    for (nstep = 0; nstep < STEPS; nstep++)
    {
      previoustime = NumToTime(i-1);
      newtime = NumToTime(i);
      deltaT = previoustime - newtime;
      time = previoustime - (nstep + 0.5) * (deltaT / STEPS);
      LightCone_Time[i * STEPS + nstep] = time;
      LightCone_RR[i * STEPS + nstep] = Time_to_Distance(time);
//      printf("%d %d %g %g, z=%f\n",i, nstep, LightCone_RR[i * STEPS + nstep], LightCone_Time[i * STEPS + nstep], Distance_to_Redshift(Time_to_Distance(time)));
    }
#endif

}



/*@brief Reads in 1/(z+1) from FileWithZList defined
 *       in ./input/input.par for the list of output snapshots.
 *       Creates a table with redshift ZZ[] and ages Age[].*/
void read_zlist(void)
{
	int i;
  FILE *fd;
  char fname[1000];
  char sbuf[1000];
  //double dumb2, dumb;

  sprintf(fname, "%s", FileWithZList);

  if(!(fd = fopen(fname, "r")))
    {
      sprintf(sbuf,"can't read output list in file '%s'\n", fname);
      terminate(sbuf);
    }

  Zlistlen = 0;
  do
    {
      //if(fscanf(fd, "%d %lg %lf\n", &dumb, &ZZ[Zlistlen], &dumb2) == 3)
      if(fscanf(fd, "%lf\n", &AA[Zlistlen]) == 1)
	Zlistlen++;
      else
      	break;

    }
  while(Zlistlen < LastDarkMatterSnapShot+1);

  fclose(fd);

  for(i = 0; i < Zlistlen; i++)
     {
  	//convert AA[] into redshift - ZZ[]
       ZZ[i] = 1 / AA[i] - 1;
       //table with time in internal units (Mpc/Km/s/h)
       Age[i] = time_to_present(ZZ[i]);

     }

#ifndef MCMC
#ifdef PARALLEL
  if(ThisTask == 0)
    printf("found %d defined times in zlist.\n", Zlistlen);
#else
  printf("found %d defined times in zlist.\n", Zlistlen);
#endif
#endif

}


void read_zlist_new(void)
{
  int i, dummy_snap;
  double dummy_a;
  FILE *fd;
  char fname[1000];

  sprintf(fname, "%s", FileWithZList);
  if(!(fd = fopen(fname, "r")))
    {
  	char sbuf[1000];
  	sprintf(sbuf, "can't open file `%s'\n", fname);
  	terminate(sbuf);
    }

  Zlistlen = 0;
  do
    {
      if(fscanf(fd, "%d %lg %lg", &dummy_snap, &ZZ[Zlistlen], &dummy_a) == 3)
    	  Zlistlen++;
      else
	break;

    }
  while(Zlistlen < LastDarkMatterSnapShot+1);
  fclose(fd);

  for(i = 0; i < Zlistlen; i++)
  {
  	//convert redshift - ZZ[] into AA[]
  	AA[i] = 1/(ZZ[i] + 1);
  	//printf("z[%d]=%f\n",i,ZZ[i]);
  	//table with time in internal units (Mpc/Km/s/h)
  	if(ZZ[i]>=0.0)
  		Age[i] = time_to_present(ZZ[i]);
  	else
  		Age[i] = 0.0;
  		//break;
  }


#ifndef MCMC
#ifdef PARALLEL
  if(ThisTask == 0)
    printf("found %d defined times in zlist.\n", Zlistlen);
#else
  printf("found %d defined times in zlist.\n", Zlistlen);
#endif
#endif
}

void read_zlist_original_cosm(void)
{
  FILE *fd;
  char fname[1000];

  sprintf(fname, "%s", FileWithZList_OriginalCosm);
  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      terminate("in read_zlist_original_cosm");
    }

  Zlistlen = 0;
  do
    {
      if(fscanf(fd, " %lg ", &AA_OriginalCosm[Zlistlen]) == 1)
        {Zlistlen++;
	}
      else
        break;
    }
  while(Zlistlen < LastDarkMatterSnapShot+1);

  fclose(fd);

#ifndef MR_PLUS_MRII //option for MCMC
#ifdef PARALLEL
  if(ThisTask == 0)
    printf("found %d defined times in zlist.\n", Zlistlen);
#else
  printf("found %d defined times in zlist.\n", Zlistlen);
#endif
#endif
}


/**@brief Reads in the list of output snapshots from
 *        file /input/desired_output_snaps.txt*/
void read_output_snaps(void)
{
  int i, j;

#if !defined(GALAXYTREE) && !defined(ALL_SKY_LIGHTCONE) 
  char buf[1000];
  FILE *fd;

  LastSnapShotNr=0;

  sprintf(buf, "%s", FileWithOutputRedshifts);

  if(!(fd = fopen(buf, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "file `%s' not found.\n", buf);
      terminate(sbuf);
    }

  for(i = 0; i < NOUT; i++)
    {
      if(fscanf(fd, " %f ", &ListOutputRedshifts[i]) != 1)
        {
    	  char sbuf[1000];
    	  sprintf(sbuf, "I/O error in file '%s'\n", buf);
    	  terminate(sbuf);
        }

      //find the snapshot corresponding to the desired output redshift in ListOutputRedshifts[]
      for(j = 0; j < LastDarkMatterSnapShot+1; j++)
      	if(ListOutputRedshifts[i]>=ZZ[j])
      	{
      		if((ZZ[j-1]-ListOutputRedshifts[i])<(ListOutputRedshifts[i]-ZZ[j]) || ZZ[j]< 0.0)
      			ListOutputSnaps[i]=j-1;
    		  else
    		    ListOutputSnaps[i]=j;
    		  break;
      	}


#ifdef MCMC
      		if (ThisTask == 0 && CurrentMCMCStep==1)
#else
      		if (ThisTask == 0)
#endif
      			printf("requested z=%0.2f, available snap[%d] z=%f & snap[%d] z=%f, use snap[%d]\n",
      					ListOutputRedshifts[i], j-1, ZZ[j-1], j, ZZ[j], ListOutputSnaps[i]);


      //define the LastSnapShotNr as the highest snapshot need to be computed
      if(LastSnapShotNr<ListOutputSnaps[i])
    	LastSnapShotNr=ListOutputSnaps[i];
    }
  fclose(fd);

#else
  for(i = 0; i < NOUT; i++)
    ListOutputSnaps[i] = i;
  LastSnapShotNr=LastDarkMatterSnapShot;
#endif
}




/**@brief Sets up some variables used to convert from physical to internal
 *        units (as UnitDensity_in_cgs); These are obtained from the three
 *        defined in input.par: UnitLength_in_cm (cm to Mpc), UnitMass_in_g
 *        (g to 1e10Msun) and UnitVelocity_in_cm_per_s (cm/s to km/s).
 *
 *       As defined in input.par, \f$\rm{UnitLength}_{\rm{cm}}=
 *       3.08568\times 10^{24}\rm{cm}\f$, converts from cm into Mpc and
 *       \f$\rm{UnitVelocity}_{\rm{cm/s}}=10000\rm{cm/s}\f$, converts from
 *       cm/s to Km/s (cm to km). In set_units() \f$\rm{UnitTime}_{\rm{s}}\f$
 *       is derived from these two quantities:
 *
 *       \f$\frac{\rm{UnitLength}_{\rm{cm}}}{\rm{UnitVelocity}_{\rm{cm/s}}}
 *       =3.08568\times 10^{19}\rm{Mpc}~\rm{Km}^{-1}\rm{s}^{-1}\f$,
 *
 *       through the code \f$t_{\rm{dyn}}\f$ has internal units and its never
 *       converted (note that \f$t_{\rm{dyn}}\f$ has an h factor, as the code internal
 *       units which despite not being included is \f$\rm{UnitTime}_{\rm{s}}\f$ is
 *       included in the output of time_to_present() - so it is consistent).
 *
 *       \f$ \rm{UnitDensity}_{\rm{cgs}} =
 *       \frac{\rm{UnitMass}_{\rm{g}}}{\rm{UnitLength}_{\rm{cm}}^3}=6.769898\times 10^{-31}\f$,
 *       converts density in \f$\rm{g}~\rm{cm}^{-3}\f$ into internal units
 *       \f$(10^{10}M_{\odot}\rm{Mpc}^{-3})\f$
 *
 *        \f$ \rm{UnitPressure}_{\rm{cgs}} =
 *        \frac{\rm{UnitMass}_{\rm{g}}}{\rm{UnitLength}_{\rm{cm}} \times \rm{UnitTime}_{\rm{s}}^2}
 *        =6.769898\times 10^{-21}\f$, converts pressure in
 *        \f$\rm{g}~\rm{cm}^{-1}\rm{s}^{-2}\f$ into internal units
 *        \f$(10^{10}M_{\odot}~\rm{Mpc}^{-1}(Mpc/Mk/s) \f$
 *
 *       \f$ \rm{UnitCoolingRate}_{\rm{cgs}} =
 *       \frac{\rm{UnitPressure}_{\rm{cgs}}}{\rm{UnitTime}_{\rm{s}}}=2.193973\times 10^{-40}\f$,
 *       converts the cooling rate in \f$\rm{g}~\rm{cm}^{-1}\rm{s}^{-3}\f$ into
 *       internal units \f$(10^{10}M_{\odot}~\rm{Mpc}^{-1}(Mpc/Mk/s)^{-3}) \f$
 *
 *      \f$ \rm{UnitEnergy}_{\rm{cgs}} =
 *       \frac{\rm{UnitMass}_{\rm{g}} \times \rm{UnitLength}_{\rm{cm}}^2}{\rm{UnitTime}_{\rm{s}}^2}
 *       =1.989000\times 10^{53}\f$, converts energy in
 *       \f$\rm{g}~\rm{cm}^2\rm{s}^{-2}\f$ into internal units
 *       \f$(10^{10}M_{\odot}~\rm{Mpc}^{2}(Mpc/Mk/s)^{-2})\f$
 *
 *       \f$ \rm{Hubble} = \rm{HUBBLE} \times \rm{UnitTime}_{\rm{s}}=100.0001\f$, where
 *       \f$\rm{HUBBLE}=3.2407789\times 10^{-18} h~\rm{s}^{-1}\f$, is the hubble
 *       constante in \f$(h~\rm{Km}~\rm{s}^{-1}\rm{Mpc}^{-1})\f$.
 *
 *       Then define a constant:
 *       \f$ \rm{RhoCrit} = \frac{3\times \rm{Hubble}^2}{8\times \rm{M}\_\rm{PI} \times G}=
 *       27.75505~h^2~10^{10}M_{\odot}\rm{Mpc}^{-3}\f$,
 *
 *       */

void set_units(void)
{


	// SEC_PER_MEGAYEAR   3.155e13
	// SEC_PER_YEAR       3.155e7

  //UnitLength_in_cm & UnitVelocity_in_cm_per_s; defined at input.par
  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitTime_in_Megayears = UnitTime_in_s / SEC_PER_MEGAYEAR;
  UnitTime_in_years = 1e6*UnitTime_in_Megayears;


  //gravity in internal units
  G = GRAVITY / pow3(UnitLength_in_cm) * UnitMass_in_g * pow2(UnitTime_in_s);//43.00708

  //converts g.cm^-3 into internal units (1e10Msun Mpc^-3)
  UnitDensity_in_cgs = UnitMass_in_g / pow3(UnitLength_in_cm);//6.769898e-31

  //converts g.cm^-1s^-2 into internal units (10^10Msun.Mpc^-1(Mpc/Km/s)^-2) \f$
  UnitPressure_in_cgs = UnitMass_in_g / UnitLength_in_cm / pow2(UnitTime_in_s);//6.769898e-21

  //converts g.cm^-1.s^-3 into internal units (10^10Msun.Mpc^-1(Mpc/Km/s)^-3)
  UnitCoolingRate_in_cgs = UnitPressure_in_cgs / UnitTime_in_s;//2.193973e-40

  //converts g.cm^2.s^-2 into internal units (10^10Msun.Mpc^2(Mpc/Km/s)^-2)
  UnitEnergy_in_cgs = UnitMass_in_g * pow2(UnitLength_in_cm) / pow2(UnitTime_in_s);//1.989000e+53


  /* converts the Hubble constant from h.s^-1 into h.Km.s^-1.Mpc-1 */
  // Would make much more sense to define this including the h100 factor.
  Hubble = HUBBLE * UnitTime_in_s;//100.000

  /* compute a quantitity */
  //TODO only used by halomodel
  RhoCrit = 3 * Hubble * Hubble / (8 * M_PI * G);//27.75505 (h^2.10^10Msun.Mpc^-3)

}


#ifdef SPECIFYFILENR

void read_file_nrs(void)
{
  int i;
  char buf[1000];
  FILE *fd;
  sprintf(buf, "%s", FileNrDir);
  if(!(fd = fopen(buf, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "file `%s' not found.\n", buf);
      terminate(sbuf);
    }

  for(i = 0; i < 111; i++) //only 111files in ../input/filenrdir.txt are read in
    {
      if(fscanf(fd, " %d ", &ListInputFilrNr[i]) != 1)
	{
	  char sbuf[1000];
	  sprintf(sbuf, "I/O error in file '%s'\n", buf);
	  terminate(sbuf);
	}
    }
  fclose(fd);
  
}


#endif




void read_reionization(void)
{
  FILE *fd;
  int p;
  float dumb;

  if(!(fd = fopen(McFile, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "file `%s' not found.\n", McFile);
      terminate(sbuf);
    }

  for(p = 0; p < 45; p++)
    {
      fscanf(fd, "%f", &dumb);
      fscanf(fd, "%f", &Reion_z[p]);
      fscanf(fd, "%f", &Reion_Mc[p]);

    }
  Reion_z[45] = Reion_z[44];
  Reion_Mc[45] = Reion_Mc[44];


  fclose(fd);


}

void find_interpolate_reionization(double zcurr, int *tabindex, double *f1, double *f2)
{
  double frac;
  int idx;

  if(zcurr > Reion_z[0])	/* redshift is higher than the largest z in the table, take the 1st entry */
    {
      *tabindex = 0;
      *f1 = 1.;
      *f2 = 0.;
    }
  else if(zcurr <= Reion_z[44])	/* redshift smaller than the smallest rho in the table, take the last entry */
    {
      *tabindex = 44;
      *f1 = 1.;
      *f2 = 0.;
    }
  else
    {
      idx = 0;
      while(Reion_z[idx + 1] > zcurr)
	idx++;

      frac = (-zcurr + Reion_z[idx]) / (Reion_z[idx] - Reion_z[idx + 1]);
      *f1 = 1 - frac;
      *f2 = frac;
      *tabindex = idx;
    }

}
