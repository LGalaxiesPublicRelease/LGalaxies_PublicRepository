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
  double SpinParameter, dgas, Vmax;

  if(DiskRadiusModel == 0 || DiskRadiusModel == 1)
    {
      if (Gal[p].Type == 0)
	Vmax=Gal[p].Vmax;
      else
	Vmax=Gal[p].InfallVmax;

    if(Halo[halonr].Spin[0]==0 && Halo[halonr].Spin[1]==0 && Halo[halonr].Spin[2]==0)
    	dgas = Gal[p].Rvir / 10.0;
    else
    	dgas = 3.0 * sqrt(Halo[halonr].Spin[0] * Halo[halonr].Spin[0] +
			  Halo[halonr].Spin[1] * Halo[halonr].Spin[1] +
			  Halo[halonr].Spin[2] * Halo[halonr].Spin[2] ) / 2.0 / Vmax;

  	return dgas;

  }
  else
    /*  simpler prescription */
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

double get_gas_disk_radius(int p)
{
  double Vmax, Radius;

#ifndef H2_AND_RINGS
  if (Gal[p].Type == 0)
    Vmax=Gal[p].Vmax;
  else
    Vmax=Gal[p].InfallVmax;

  //if the spin is not available for the ColdGas
  if((Gal[p].ColdGasSpin[0]==0 && Gal[p].ColdGasSpin[1]==0 && Gal[p].ColdGasSpin[2]==0) || Gal[p].ColdGas==0.)
    Radius = Gal[p].Rvir / 10.0;
  else
    Radius = 3.0 * sqrt(Gal[p].ColdGasSpin[0] * Gal[p].ColdGasSpin[0] +
				      Gal[p].ColdGasSpin[1] * Gal[p].ColdGasSpin[1] +
				      Gal[p].ColdGasSpin[2] * Gal[p].ColdGasSpin[2] ) / 2.0 / Vmax;
#else
  int jj;

  if(Gal[p].ColdGas==0)
    Radius=0.;
  else if(Gal[p].ColdGas<1.0e-6)
    Radius=RingRadius[0]/2.;
  else
    {
      Radius=0.5*RingRadius[0]*Gal[p].ColdGasRings[0];
      for(jj=1;jj<RNUM;jj++)
	Radius+=(0.5*(RingRadius[jj-1]+RingRadius[jj])*Gal[p].ColdGasRings[jj]);
      Radius=3.0*Radius/Gal[p].ColdGas/2.0;      //2.0=mean radius/scale length for exponential disk
    }
  //else
  //  Radius = Gal[p].Rvir / 10.0;
#endif

  return Radius;
}

/** @brief Updates the stellar disk radius.
 *
 *  The stellar disk is assumed to be thin, in centrifugal equilibrium and to have
 *  an exponential density profile:
 *
 *  \f$ \Sigma(R_{\star})=
 *      \Sigma_{\star\rm{0}}e^{-\frac{R_{\star}}{R_{\rm{\star,d}}}}, \f$
 *
 *  where \f$R_{\rm{\star,d}}\f$ is the scale length of the gas disk and rd=get_initial_disk_radius(Gal[t].HaloNr, t)/3.;
 *  \f$\Sigma_{\star\rm{0}}\f$ is the corresponding central surface density.
 *
 *  Assuming a flat circular velocity curve (galaxy with a negligible self-gravity)
 *  in an isothermal dark matter halo, the gas disk scale length is given by:
 *
 *  \f$ R_{\rm{\star,d}}=\frac{J_{\star}/M_{\star}}{2V_{\rm{max}}}\f$,
 *
 *  assuming that the maximum circular velocity of satellite galaxies does not
 *  change after infall (inner dark matter regions are compact and don't change). */

double get_stellar_disk_radius(int p)
{
  double Vmax, Radius;

#ifndef H2_AND_RINGS
  if (Gal[p].Type == 0)
    Vmax=Gal[p].Vmax;
  else
    Vmax=Gal[p].InfallVmax;

  //if the spin is not available for the ColdGas
  if((Gal[p].DiskSpin[0]==0 && Gal[p].DiskSpin[1]==0 && Gal[p].DiskSpin[2]==0) || Gal[p].DiskMass==0.)
    Radius = Gal[p].Rvir / 10.0;
  else
    Radius = 3.0 * sqrt(Gal[p].DiskSpin[0] * Gal[p].DiskSpin[0] +
			Gal[p].DiskSpin[1] * Gal[p].DiskSpin[1] +
			Gal[p].DiskSpin[2] * Gal[p].DiskSpin[2] ) / 2.0 / Vmax;

#else
    int jj;
    if(Gal[p].DiskMass==0.)
      Radius=0.;
    else if (Gal[p].DiskMass<1.0e-6)
      Radius=RingRadius[0]/2.;
    else
      {
	Radius=0.5*RingRadius[0]*Gal[p].DiskMassRings[0];
	for(jj=1;jj<RNUM;jj++)
	  Radius+=(0.5*(RingRadius[jj-1]+RingRadius[jj])*Gal[p].DiskMassRings[jj]);
	Radius=3.0*Radius/Gal[p].DiskMass/2.0;      //2.0=mean radius/scale length for exponential disk
      }
#endif

    return Radius;
}




/** @brief Calculates the half mass radius of galaxies */
double half_mass_radius(int p, int do_ColdGas, int do_DiskMass, int do_BulgeMass)
{
  double r=0., rb, M;
  int ii;
#ifndef H2_AND_RINGS
  double Mdisk, rmax, Mgas, Mbulge, rd, rgd, dr, totmass;
  int N = 1000;
  #define RADIUS_RMIN 5e-7
  #define RADIUS_N 100
#else
  double BulgeMassRings[RNUM], massRings;
#endif



#ifndef H2_AND_RINGS
  rgd = Gal[p].ColdGasRadius/3.;
  rd = Gal[p].DiskRadius/3.;
  rb = Gal[p].BulgeSize;

  Mgas = Gal[p].ColdGas;
  Mdisk = Gal[p].DiskMass;
  Mbulge = Gal[p].BulgeMass;

  totmass = 0.;
  rmax=0.;

  if(do_ColdGas==1)
    {
      totmass+=Mgas;
      rmax=max(rmax,1.68*rgd);
    }
  if(do_DiskMass==1)
    {
      totmass+=Mdisk;
      rmax=max(rmax,1.68*rd);
    }
  if(do_BulgeMass==1)
    {
      totmass+=Mbulge;
      rmax=max(rmax,rb);
    }


  //rmax=max(rb,1.68*max(rd,rgd));
  if (rmax < 2.*RADIUS_RMIN)
  	return(rmax);
  dr=(rmax-RADIUS_RMIN)/(float)RADIUS_N;

    // increases the search radius until it encompasses half the total mass taking
    // into account the stellar disk, stellar bulge and cold gas disk.
  ii = 0;
  do {
      r = (RADIUS_RMIN) + (ii+0.5)* dr;

      M=0.;

      if(do_ColdGas==1)
	if(Mgas>0.)
	  M += Mgas*diskmass(r/rgd);
      if(do_DiskMass==1)
	if(Mdisk>0.)
	  M += Mdisk*diskmass(r/rd);

      // due to a bug in the code these functions were called when Mass=0 leading to nans being returned and no disruption for pure disks
      // or pure bulges. That bug is reproduced here to allow reconstruction of these erroneous results
#if !defined(GUO10) && !defined(GUO13) && !defined(HENRIQUES13)
  	if(Mbulge>0.)
#endif
  	  if(do_BulgeMass==1)
  	    M += Mbulge*bulgemass(r/rb);

  	ii++;
  	if(ii > N) terminate ("couldn't find half mass radius");
  }
  while(M < 0.5*totmass);

#else //if H2_AND_RINGS

  M=0.;
  if(do_ColdGas==1)
    M+=0.5*Gal[p].ColdGas;
  if(do_DiskMass==1)
    M+=0.5*Gal[p].DiskMass;
  if(do_BulgeMass==1)
    M+=0.5*Gal[p].BulgeMass;

#ifndef RINGS_IN_BULGES
  rb=Gal[p].BulgeSize;

  //bulge mass radius distribution according to Jaffe's radius. In each ring, mass=M_bulge*[1/(1+r_in/r_b)-1/(1+r_out/r_b)]
  if(Gal[p].BulgeSize<1.0e-6)
    for(ii=0;ii<RNUM;ii++)
      BulgeMassRings[ii]=0.0;
  else
    {
      //here mass(r)=mass*(1-1/(1+x)) is used which is equivalente to  mass(r)=mass*x/(1+x)
      BulgeMassRings[0]=Gal[p].BulgeMass*(1-1/(1+RingRadius[0]/rb));
      for(ii=1;ii<RNUM;ii++)
	BulgeMassRings[ii]=Gal[p].BulgeMass*(1/(1+RingRadius[ii-1]/rb)-1/(1+RingRadius[ii]/rb));
    }
#else
  for(ii=0;ii<RNUM;ii++)
    BulgeMassRings[ii]=Gal[p].BulgeMassRings[ii];
#endif

  if(M<1.0e-6)
    r=0.5*RingRadius[0];
  else
    {
      for(ii=0;ii<RNUM;ii++)
	{
	  massRings=0.;
	  if(do_ColdGas==1)
	    massRings+=Gal[p].ColdGasRings[ii];
	  if(do_DiskMass==1)
	    massRings+=Gal[p].DiskMassRings[ii];
	  if(do_BulgeMass==1)
	    massRings+=BulgeMassRings[ii];
	  //massRings=Gal[p].ColdGasRings[ii]+Gal[p].DiskMassRings[ii]+BulgeMassRings[ii];
	  if(M>massRings)
	    M-=massRings;
	  else
	    break;
	}
      if(ii==RNUM)
	r=RingRadius[RNUM-1];
      else
	{
	  //weight r by the remaining mass over the mass in current ring
	  if(ii==0)
	    r=M/massRings*RingRadius[ii];
	  else
	    r=M/massRings*RingRadius[ii]+(1-M/massRings)*RingRadius[ii-1];
	}
    }


#endif //H2_AND_RINGS


  return (r);
}


double isothermal_mass(double Mvir, double Rvir, double dr)
{
	return Mvir/Rvir * dr;
}


/** @brief Returns the mass of a disk within a given radius in units of the scale length
 *         Disk profile -> exponential */
double diskmass(double x)
{
  return 1.-(1.+x)*exp(-x);
}

/** @brief Returns the mass of a bulge at a certain radius.
 *         Bulge profile -> de Vaucouleurs type r^{1/4} law */

// The previous complicated expression seemed to be a long-winded way of saying that
// the density varies as 1/x^2(1+x)^2, leading to
double bulgemass(double x)
{
  return x/(1.+x);
}





#ifdef COMPUTE_SPECPHOT_PROPERTIES
//o->ObsMagDust[nlum]=ObsLumDiskDust+ObsLumBulgeDust;
/** @brief Calculates the half light radius of satellite galaxies on the Vband, NMAG=2 for 40 bands*/
double stellar_half_light_radius(struct GALAXY_OUTPUT *o)
{
  double r, rd, rb, L;
  int ii;
  double Ldisk, Lbulge, totL;
#ifndef OUTPUT_RINGS
  double rmax, dr;
  int N = 1000;
  #define SAT_RADIUS_RMIN 5e-7
  #define SAT_RADIUS_N 100
#else
  double BulgeMassRings[RNUM], LightRings;
#endif

  r=0.;
  rd=o->DiskRadius/3.;
  rb=o->BulgeSize;
  Ldisk=mag_to_lum(o->MagDust[2])-mag_to_lum(o->MagBulge[2]);
  Lbulge = mag_to_lum(o->MagBulge[2]);
  totL = mag_to_lum(o->MagDust[2]);

#ifndef OUTPUT_RINGS
  rmax=max(rb,1.68*rd);
  if (rmax < 2.*SAT_RADIUS_RMIN)
  	return(rmax);
  dr=(rmax-SAT_RADIUS_RMIN)/(float)SAT_RADIUS_N;

    // increases the search radius until it encompasses half the total mass taking
    // into account the stellar disk, stellar bulge and cold gas disk.
  ii = 0;
  do {
      L=0.;
      // Not sure that we need the 0.5 here - it's all a matter of definition
  	r = (SAT_RADIUS_RMIN) + (ii+0.5)* dr;
  	if(Ldisk>0.)
  	  L += Ldisk*diskmass(r/rd);

  	if(Lbulge>0.)
  	  L += Lbulge*bulgemass(r/rb);

  	ii++;
  	if(ii > N) terminate ("couldn't find half mass radius");
  }
  while(L < 0.5*totL);

#else //if OUTPUT_RINGS

  L=0.5*(totL);   //to find half mass radius

#ifndef RINGS_IN_BULGES
  rb=o->BulgeSize;
  //bulge mass radius distribution according to Jaffe's radius. In each ring, mass=M_bulge*[1/(1+r_in/r_b)-1/(1+r_out/r_b)]
  if(o->BulgeSize<1.0e-6)
    for(ii=0;ii<RNUM;ii++)
      BulgeMassRings[ii]=0.0;
  else
    {
      BulgeMassRings[0]=o->BulgeMass*(1-1/(1+RingRadius[0]/rb));
      for(ii=1;ii<RNUM;ii++)
	BulgeMassRings[ii]=o->BulgeMass*(1/(1+RingRadius[ii-1]/rb)-1/(1+RingRadius[ii]/rb));
    }
#else
  for(ii=1;ii<RNUM;ii++)
  	BulgeMassRings[ii]=o->BulgeMassRings[ii];
#endif //RINGS_IN_BULGES

  if(L<1.0e-6)
    r=0.5*RingRadius[0];
  else
    {
      for(ii=0;ii<RNUM;ii++)
	{
	  LightRings=(o->DiskMassRings[ii]/o->DiskMass)*Ldisk+(BulgeMassRings[ii]/o->BulgeMass)*Lbulge;  //total mass a each ring
	  if(L>LightRings)
	    L-=LightRings;
	  else
	    break;
	}
      if(ii==RNUM)
	r=RingRadius[RNUM-1];
      else
	{
	  if(ii==0)
	    r=L/LightRings*RingRadius[ii];
	  else
	    r=L/LightRings*RingRadius[ii]+(1-L/LightRings)*RingRadius[ii-1];
	}
    }
#endif //OUTPUT_RINGS


  return (r);
}
#endif



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

      Gal[p].ColdGasSpin[j] = Halo[halonr].Spin[j];
      Gal[p].DiskSpin[j] = Halo[halonr].Spin[j];
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
  //Gal[p].ReheatedGas = 0.0;
  Gal[p].EjectedMass = 0.0;
#ifdef EXCESS_MASS
  Gal[p].ExcessMass = 0.0;
#endif
  Gal[p].ICM = 0.0;

#ifdef TRACK_MASSGROWTH_CHANNELS
  Gal[p].MassFromInSitu = 0.0;
  Gal[p].MassFromMergers = 0.0;
  Gal[p].MassFromBursts = 0.0;
#endif
#ifdef TRACK_BURST
  Gal[p].BurstMass = 0.0;
#endif
  if (BlackHoleGrowth ==0)
    Gal[p].BlackHoleMass = 0.0;
  else if (BlackHoleGrowth ==1)
    Gal[p].BlackHoleMass = BlackHoleSeedMass;
  Gal[p].BlackHoleGas = 0.0;
  /*ram pressure*/
  Gal[p].HotRadius=Gal[p].Rvir;
#ifdef GALAXYTREE
  Gal[p].DisruptOn = 0;
#endif

  Gal[p].MetalsColdGas = metals_init();
  Gal[p].MetalsDiskMass = metals_init();
  Gal[p].MetalsBulgeMass = metals_init();
  Gal[p].MetalsHotGas = metals_init();
  //Gal[p].MetalsReheatedGas = metals_init();
  Gal[p].MetalsEjectedMass = metals_init();
#ifdef EXCESS_MASS
  Gal[p].MetalsExcessMass = metals_init();
#endif
  Gal[p].MetalsICM = metals_init();
#ifdef METALS_SELF
  Gal[p].MetalsHotGasSelf = metals_init();
#endif


#ifdef H2_AND_RINGS
  for(j=0;j<RNUM;j++)
    {
      Gal[p].ColdGasRings[j]=0.0;
      Gal[p].MetalsColdGasRings[j]=metals_init();
      Gal[p].DiskMassRings[j]=0.0;
      Gal[p].MetalsDiskMassRings[j]=metals_init();
#ifdef RINGS_IN_BULGES
      Gal[p].BulgeMassRings[j]=0.0;
      Gal[p].MetalsBulgeMassRings[j]=metals_init();
#endif
    }
#endif
  //inclination defined as the angle between galaxy spin and the z-axis
#ifdef COMPUTE_SPECPHOT_PROPERTIES
  Gal[p].CosInclination = 0.0;
#endif

  Gal[p].PrimordialAccretionRate = 0.0;
  Gal[p].CoolingRate = 0.0;
  Gal[p].CoolingRate_beforeAGN = 0.0;
  Gal[p].CoolingRadius = 0.0;
  Gal[p].CoolingGas = 0.0;
  Gal[p].QuasarAccretionRate=0.0;
  Gal[p].RadioAccretionRate=0.0;
  Gal[p].AGNheatingFromCentral = 0.0;
#ifdef TRACK_NMERGERS
  Gal[p].NMajorMergers = 0.0;
  Gal[p].NMinorMergers = 0.0;
#endif
  Gal[p].Sfr = 0.0;
  Gal[p].SfrBulge = 0.0;
#ifdef H2_AND_RINGS
 for(j=0;j<RNUM;j++) Gal[p].SfrRings[j]=0;
#endif


  Gal[p].StarMerge=0.0;

  Gal[p].XrayLum = 0.0;

  Gal[p].ColdGasRadius = get_initial_disk_radius(halonr, p);
  Gal[p].DiskRadius = Gal[p].ColdGasRadius;
  Gal[p].BulgeSize = 0.0;

#ifndef HT09_DISRUPTION
  Gal[p].OriMergTime = 0.0;
  Gal[p].MergTime = 0.0;
#else
  Gal[p].OriMergRadius = 0.0;
  Gal[p].MergRadius = 0.0;
#endif

#ifdef TRACK_SPLASHBACKS
  Gal[p].flagSplashBack=0;
  Gal[p].TimeSinceSplashBack=0.;
#endif

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
  int ll;
  for(ll=0;ll<NUM_ELEMENTS;ll++)
    {
      Gal[p].ColdGas_elements[ll]     = 0.;
      Gal[p].DiskMass_elements[ll]    = 0.;
      Gal[p].BulgeMass_elements[ll]   = 0.;
      Gal[p].HotGas_elements[ll]      = 0.;
      //Gal[p].ReheatedGas_elements[ll] = 0.;
      Gal[p].EjectedMass_elements[ll] = 0.;
#ifdef EXCESS_MASS
      Gal[p].ExcessMass_elements[ll]  = 0.;
#endif
      Gal[p].ICM_elements[ll]         = 0.;
#ifdef H2_AND_RINGS
      for(j=0;j<RNUM;j++)
	{
	  Gal[p].ColdGasRings_elements[j][ll]  = 0.;
	  Gal[p].DiskMassRings_elements[j][ll] = 0.;
#ifdef RINGS_IN_BULGES
	  Gal[p].BulgeMassRings_elements[j][ll] = 0.;
#endif
	}
#endif
    }
#endif
}

/*TODO take away magnitudes and work with luminositites*/

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
  double f1, f2, fmet1, fmet2, LuminosityToAdd; 
#ifdef OUTPUT_MOMAF_INPUTS
  double dLuminosityToAdd;
#endif
  double X1, age, tbc;
 	int N_AgeBins=1, ii;
  double upper_time;

  //TODO define elsewhere and maybe make the 10. an input parameter?
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

    //SF in random time in step
  	//double rand;
  	//rand=ran1(&mu_seed);
  	//printf("ran=%f\n",rand/100.*3.);
  	//time=upper_time-dt*rand/100.*3.;  	
  	//put the the time at the end of the bin (divide bin by 5 and put in 5th bin) 	
  	//time=upper_time-3*dt/(5.)-dt/(5.)/2.;


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
						                    f2 * LumTables[j][metindex + 1][zindex + 1][tabindex +	1]));
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



#ifdef HT09_DISRUPTION
void sub_to_luminosities(int p, float RemainFract)
{
  int outputbin, j;

#ifdef OUTPUT_REST_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++)
    {
	  for(j = 0; j < NMAG; j++)
		{
		  Gal[p].Lum[j][outputbin] *= RemainFract;
		  Gal[p].YLum[j][outputbin] *= RemainFract;
		  Gal[p].LumBulge[j][outputbin] *= RemainFract;
		  Gal[p].YLumBulge[j][outputbin] *= RemainFract;

		}
    }

#endif //OUTPUT_REST_MAGS

#ifdef COMPUTE_OBS_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++)
    {
	  for(j = 0; j < NMAG; j++)
		{
		  Gal[p].ObsLum[j][outputbin] *= RemainFract;
		  Gal[p].ObsYLum[j][outputbin] *= RemainFract;
		  Gal[p].ObsLumBulge[j][outputbin] *= RemainFract;
		  Gal[p].ObsYLumBulge[j][outputbin] *= RemainFract;
#ifdef OUTPUT_MOMAF_INPUTS
		  Gal[p].dObsLum[j][outputbin] *= RemainFract;
		  Gal[p].dObsYLum[j][outputbin] *= RemainFract;
		  Gal[p].dObsLumBulge[j][outputbin] *= RemainFract;
		  Gal[p].dObsYLumBulge[j][outputbin] *= RemainFract;
#endif
		}
    }
#endif // COMPUTE_OBS_MAGS

}
#endif //HT09_DISRUPTION



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

double mag_to_lum(double mag)
{
  if(mag < 99.0)
    return pow(10,-1.*(mag/2.5));
  else
    return 0.0;
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

#ifdef TRACK_SPLASHBACKS
  if(Gal[ngal].flagSplashBack==1)
    Gal[ngal].TimeSinceSplashBack+= NumToTime(Gal[ngal].SnapNum-1)-NumToTime(Gal[ngal].SnapNum);
#endif

}


/**@brief Updates properties of type 1 galaxies.
 *
 * If MERGE01 = 1, then a dynamical friction decay time scale is calculated
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

#ifdef TRACK_SPLASHBACKS
  Gal[ngal].flagSplashBack=0;
  Gal[ngal].TimeSinceSplashBack=0.;
#endif

#ifdef MERGE01

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
#ifndef HT09_DISRUPTION
      	if (Halo[Halo[halonr].FirstHaloInFOFgroup].FirstProgenitor == -1 )
      		Gal[ngal].MergTime = estimate_merging_time(prog,Halo[halonr].FirstHaloInFOFgroup,ngal);
      	else
      		Gal[ngal].MergTime = estimate_merging_time(prog,Halo[Halo[halonr].FirstHaloInFOFgroup].FirstProgenitor,ngal);
      	Gal[ngal].MergTime -= NumToTime(Halo[halonr].SnapNum) - NumToTime(Halo[prog].SnapNum);
	  	 //to calculate the position of type 2
      	Gal[ngal].OriMergTime=Gal[ngal].MergTime;
      	Gal[ngal].OriMvir = get_virial_mass(prog);
      	Gal[ngal].OriRvir = get_virial_radius(prog);
#else
      	int central_halonr;
      	if (Halo[Halo[halonr].FirstHaloInFOFgroup].FirstProgenitor == -1 )
      		central_halonr=Halo[halonr].FirstHaloInFOFgroup;
      	else
      		central_halonr=Halo[Halo[halonr].FirstHaloInFOFgroup].FirstProgenitor;

      	Gal[ngal].MergRadius = get_merging_radius (prog, central_halonr, ngal);
      	Gal[ngal].OriMergRadius = Gal[ngal].MergRadius;
      	Gal[ngal].OriMergmass = get_virial_mass(prog);
#endif
      }
    }
  }
#endif
  
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

  mass_checks(ngal,"model_misc.c",__LINE__);

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

#ifdef TRACK_SPLASHBACKS
  Gal[ngal].flagSplashBack=0;
  Gal[ngal].TimeSinceSplashBack=0.;
#endif

  if(HotGasOnType2Galaxies==0)
    {
      Gal[ngal].HotGas = 0.0;
      Gal[ngal].MetalsHotGas=metals_init();
#ifdef DETAILED_METALS_AND_MASS_RETURN
#ifdef INDIVIDUAL_ELEMENTS
      int kk;
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[ngal].HotGas_elements[kk]=0.;
#endif
#endif
      Gal[ngal].HotRadius = 0.0;
    }

  /* Estimate remaining merging timescale. */
  if (Gal[ngal].MergeOn == 0)
  {
  	//if central galaxy has no progenitor
  	if (mostmassive == -1)
    		mostmassive = halonr;
#ifndef HT09_DISRUPTION
  	Gal[ngal].MergTime = estimate_merging_time(prog,mostmassive,ngal);
  	Gal[ngal].MergTime -= NumToTime(Halo[halonr].SnapNum) - NumToTime(Halo[prog].SnapNum);
  	//to calculate the position of type 2
  	Gal[ngal].OriMergTime=Gal[ngal].MergTime;
  	Gal[ngal].OriMvir = get_virial_mass(prog);
  	Gal[ngal].OriRvir = get_virial_radius(prog);
#else
  	Gal[ngal].MergRadius = get_merging_radius (prog, mostmassive, ngal);
  	Gal[ngal].OriMergRadius = Gal[ngal].MergRadius;
  	Gal[ngal].OriMergmass=get_virial_mass(prog);
#endif
  }
  mass_checks(ngal,"model_misc.c",__LINE__);
}

void transfer_material(int p, char cp[], int q, char cq[], double fraction, char call_function[], int call_line)
  {
  /* Transfers a fraction of component cq of galaxy q onto component cp of galaxy p.
   * cp and cq must each be one of
   *
   * If -DTRACK_BURST is set then can also specify BurstMass as an option.  This is
   * a little different in that it is not a separate component, so it should only
   * be transferred if both cq and cp are BurstMass
   *
   */
  float Mass;
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals Metals;
#ifdef INDIVIDUAL_ELEMENTS
  int kk;
  float Yield[NUM_ELEMENTS];
#endif
#else
  float Metals;
#endif //DETAILED_METALS_AND_MASS_RETURN

#ifdef STAR_FORMATION_HISTORY
  int ii;
  float sfh_Mass[SFH_NBIN];
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals sfh_Metals[SFH_NBIN];
 #ifdef INDIVIDUAL_ELEMENTS
  float sfh_Elements[SFH_NBIN][NUM_ELEMENTS];
 #endif
#else
   float sfh_Metals[SFH_NBIN];
#endif
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
   int ll, outputbin;
#ifdef OUTPUT_REST_MAGS
  float Lum[NMAG][NOUT];
  float YLum[NMAG][NOUT];
#endif
#ifdef COMPUTE_OBS_MAGS
  float ObsLum[NMAG][NOUT];
  float ObsYLum[NMAG][NOUT];
#ifdef OUTPUT_MOMAF_INPUTS
  float dObsLum[NMAG][NOUT];
  float dObsYLum[NMAG][NOUT];
#endif
#endif
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

  /* Sanity checks */
  if (fraction > 1.000001) {//1.000001
      char sbuf[1000];
      sprintf(sbuf, "\nparent call from: %s, line %d \ntransfer_material: fraction>1\nfraction = %.11f\nFrom '%s' to '%s\n",
	      call_function, call_line, fraction,cq, cp);
      terminate(sbuf);
  }


#ifdef STAR_FORMATION_HISTORY
  if (Gal[p].sfh_ibin != Gal[q].sfh_ibin) {
      printf("\nparent call from: %s, line %d \n*** transfer_material: inconsistent itimes ***\n",call_function, call_line);
      for(ii=0;ii<SFH_NBIN;ii++)
	printf("Bin[%d] time_1=%e dt_1=%e Nbins_1=%d time_2=%e dt_2=%e Nbins_2=%d\n",ii,
	       Gal[p].sfh_t[ii],Gal[p].sfh_dt[ii],Gal[p].sfh_Nbins[ii],
	       Gal[q].sfh_t[ii],Gal[q].sfh_dt[ii],Gal[q].sfh_Nbins[ii]);
      exit(1);
  }
#endif
#ifdef TRACK_BURST
  if ((strcmp(cp,"BurstMass")==0 && !strcmp(cq,"BurstMass")==0) ||
      (strcmp(cq,"BurstMass")==0 && !strcmp(cp,"BurstMass")==0)) {
    terminate("\n*** transfer_material: used incorrectly with BurstMass ***\n");
  }
#endif

#ifdef H2_AND_RINGS
  if(strcmp(cq,"ColdGas")==0 || strcmp(cp,"ColdGas")==0)
    {
      char sbuf[1000];
      sprintf(sbuf, "\nparent call from: %s, line %d \n*** ColdGas transfered with transfer material. Must use transfer rings with H2_AND_RINGS",call_function, call_line);
      terminate(sbuf);
    }
  if(strcmp(cq,"DiskMass")==0 || strcmp(cp,"DiskMass")==0)
    {
      char sbuf[1000];
      sprintf(sbuf, "\nparent call from: %s, line %d \n*** DiskMass transfered with transfer material. Must use transfer rings with H2_AND_RINGS",call_function, call_line);
      terminate(sbuf);
    }
#endif



  //Initialize arrays to contain mass to transfer
  Mass = 0.;
  Metals=metals_add(metals_init(),metals_init(),1.);
#ifdef INDIVIDUAL_ELEMENTS
  for(kk=0;kk<NUM_ELEMENTS;kk++)
    Yield[kk] = 0.;
#endif

#ifdef STAR_FORMATION_HISTORY
  for (ii=0; ii<=Gal[q].sfh_ibin; ii++)
    {
      sfh_Mass[ii]=0.;
      sfh_Metals[ii]=metals_add(metals_init(),metals_init(),1.);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	sfh_Elements[ii][kk]=0.;
#endif
    }
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++)
    {
      for(ll = 0; ll < NMAG; ll++)
	{
#ifdef OUTPUT_REST_MAGS
	  Lum[ll][outputbin]=0.;
	  YLum[ll][outputbin]=0.;
#endif
#ifdef COMPUTE_OBS_MAGS
	  ObsLum[ll][outputbin]=0.;
	  ObsYLum[ll][outputbin]=0.;
#ifdef OUTPUT_MOMAF_INPUTS
	  dObsLum[ll][outputbin]=0.;
	  dObsYLum[ll][outputbin]=0.;
#endif
#endif
	}
    }
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES





  //MASS AND METALS TO BE TRANSFERED
  if (strcmp(cq,"ColdGas")==0)
    {
      Mass = fraction*Gal[q].ColdGas;
      Metals = metals_add(metals_init(),Gal[q].MetalsColdGas,fraction);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Yield[kk] = Gal[q].ColdGas_elements[kk]*fraction;
#endif
      //if there is SF, gas goes to stars into the last sfh bin
#ifdef STAR_FORMATION_HISTORY
      for (ii=0; ii<Gal[q].sfh_ibin; ii++)
	{
	  sfh_Mass[ii]=0.;
	  sfh_Metals[ii]=metals_add(metals_init(),Gal[q].sfh_MetalsDiskMass[ii],0);
#ifdef INDIVIDUAL_ELEMENTS
	  for(kk=0;kk<NUM_ELEMENTS;kk++)
	    sfh_Elements[ii][kk]=Gal[q].sfh_DiskMass_elements[ii][kk]*0.;
#endif
	}
      sfh_Mass[Gal[q].sfh_ibin]=fraction*Gal[q].ColdGas;
      sfh_Metals[Gal[q].sfh_ibin]=metals_add(metals_init(),Gal[q].MetalsColdGas,fraction);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	sfh_Elements[Gal[q].sfh_ibin][kk]=Gal[q].ColdGas_elements[kk]*fraction;
#endif
#endif
    }

  else if (strcmp(cq,"HotGas")==0)
    {
      Mass=fraction*Gal[q].HotGas;
      Metals = metals_add(metals_init(),Gal[q].MetalsHotGas,fraction);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Yield[kk] = Gal[q].HotGas_elements[kk]*fraction;
#endif
    }

  /*else if (strcmp(cq,"ReheatedGas")==0)
      {
        Mass=fraction*Gal[q].ReheatedGas;
        Metals = metals_add(metals_init(),Gal[q].MetalsReheatedGas,fraction);
  #ifdef INDIVIDUAL_ELEMENTS
        for(kk=0;kk<NUM_ELEMENTS;kk++)
  	Yield[kk] = Gal[q].ReheatedGas_elements[kk]*fraction;
  #endif
      }*/

  else if (strcmp(cq,"EjectedMass")==0)
    {
      Mass=fraction*Gal[q].EjectedMass;
      Metals = metals_add(metals_init(),Gal[q].MetalsEjectedMass,fraction);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Yield[kk] = Gal[q].EjectedMass_elements[kk]*fraction;
#endif
    }

#ifdef EXCESS_MASS
  else if (strcmp(cq,"ExcessMass")==0)
    {
      Mass=fraction*Gal[q].ExcessMass;
      Metals = metals_add(metals_init(),Gal[q].MetalsExcessMass,fraction);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Yield[kk] = Gal[q].ExcessMass_elements[kk]*fraction;
#endif
    }
#endif

  else if (strcmp(cq,"DiskMass")==0)
    {
      Mass = fraction*Gal[q].DiskMass;
      Metals=metals_add(metals_init(),Gal[q].MetalsDiskMass,fraction);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Yield[kk] = Gal[q].DiskMass_elements[kk]*fraction;
#endif

#ifdef STAR_FORMATION_HISTORY
      for (ii=0; ii<=Gal[q].sfh_ibin; ii++)
	{
	  sfh_Mass[ii]=fraction*Gal[q].sfh_DiskMass[ii];
	  sfh_Metals[ii]=metals_add(metals_init(),Gal[q].sfh_MetalsDiskMass[ii],fraction);
#ifdef INDIVIDUAL_ELEMENTS
	  for(kk=0;kk<NUM_ELEMENTS;kk++)
	    sfh_Elements[ii][kk]=Gal[q].sfh_DiskMass_elements[ii][kk]*fraction;
#endif
	}
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
      for(outputbin = 0; outputbin < NOUT; outputbin++)
	{
	  for(ll = 0; ll < NMAG; ll++)
	    {
#ifdef OUTPUT_REST_MAGS
	      Lum[ll][outputbin]=fraction*(Gal[q].Lum[ll][outputbin]-Gal[q].LumBulge[ll][outputbin]);
	      YLum[ll][outputbin]=fraction*(Gal[q].YLum[ll][outputbin]-Gal[q].YLumBulge[ll][outputbin]);
#endif
#ifdef COMPUTE_OBS_MAGS
	      ObsLum[ll][outputbin]=fraction*(Gal[q].ObsLum[ll][outputbin]-Gal[q].ObsLumBulge[ll][outputbin]);
	      ObsYLum[ll][outputbin]=fraction*(Gal[q].ObsYLum[ll][outputbin]-Gal[q].ObsYLumBulge[ll][outputbin]);
#ifdef OUTPUT_MOMAF_INPUTS
	      dObsLum[ll][outputbin]=fraction*(Gal[q].dObsLum[ll][outputbin]-Gal[q].dObsLumBulge[ll][outputbin]);
	      dObsYLum[ll][outputbin]=fraction*(Gal[q].dObsYLum[ll][outputbin]-Gal[q].dObsYLumBulge[ll][outputbin]);
#endif
#endif
	    }
	}
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES
    }

  else if (strcmp(cq,"BulgeMass")==0)
    {
      Mass=fraction*Gal[q].BulgeMass;
      Metals=metals_add(metals_init(),Gal[q].MetalsBulgeMass,fraction);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Yield[kk] = Gal[q].BulgeMass_elements[kk]*fraction;
#endif

#ifdef STAR_FORMATION_HISTORY
      for (ii=0; ii<=Gal[q].sfh_ibin; ii++)
	{
	  sfh_Mass[ii]=fraction*Gal[q].sfh_BulgeMass[ii];
	  sfh_Metals[ii]=metals_add(metals_init(),Gal[q].sfh_MetalsBulgeMass[ii],fraction);
#ifdef INDIVIDUAL_ELEMENTS
	  for(kk=0;kk<NUM_ELEMENTS;kk++)
	    sfh_Elements[ii][kk]=Gal[q].sfh_BulgeMass_elements[ii][kk]*fraction;
#endif
	}
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
      for(outputbin = 0; outputbin < NOUT; outputbin++)
	{
	  for(ll = 0; ll < NMAG; ll++)
	    {
#ifdef OUTPUT_REST_MAGS
	      Lum[ll][outputbin]=fraction*Gal[q].LumBulge[ll][outputbin];
	      YLum[ll][outputbin]=fraction*Gal[q].YLumBulge[ll][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
	      ObsLum[ll][outputbin]=fraction*Gal[q].ObsLumBulge[ll][outputbin];
	      ObsYLum[ll][outputbin]=fraction*Gal[q].ObsYLumBulge[ll][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	      dObsLum[ll][outputbin]=fraction*Gal[q].dObsLumBulge[ll][outputbin];
	      dObsYLum[ll][outputbin]=fraction*Gal[q].dObsYLumBulge[ll][outputbin];
#endif
#endif
	    }
	}
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES
    }

  else if (strcmp(cq,"ICM")==0)
    {
      Mass=fraction*Gal[q].ICM;
      Metals=metals_add(metals_init(),Gal[q].MetalsICM,fraction);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Yield[kk] = Gal[q].ICM_elements[kk]*fraction;
#endif

#ifdef STAR_FORMATION_HISTORY
      for (ii=0; ii<=Gal[q].sfh_ibin; ii++)
	{
	  sfh_Mass[ii]=fraction*Gal[q].sfh_ICM[ii];
	  sfh_Metals[ii]=metals_add(metals_init(),Gal[q].sfh_MetalsICM[ii],fraction);
#ifdef INDIVIDUAL_ELEMENTS
	  for(kk=0;kk<NUM_ELEMENTS;kk++)
	    sfh_Elements[ii][kk]=Gal[q].sfh_ICM_elements[ii][kk]*fraction;
#endif
	}
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
#ifdef ICL
      for(outputbin = 0; outputbin < NOUT; outputbin++)
	{
	  for(ll = 0; ll < NMAG; ll++)
	    {
#ifdef OUTPUT_REST_MAGS
	      Lum[ll][outputbin]=fraction*Gal[q].ICLLum[ll][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
	      ObsLum[ll][outputbin]=fraction*Gal[q].ObsICLLum[ll][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	      dObsLum[ll][outputbin]=fraction*Gal[q].dObsICLLum[ll][outputbin];
#endif
#endif
	    }
	}
#endif
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES
    }

#ifdef TRACK_BURST
  else if (strcmp(cq,"BurstMass")==0)
    {
      Mass = fraction*Gal[q].BurstMass;
#ifdef STAR_FORMATION_HISTORY
      for (ii=0; ii<=Gal[q].sfh_ibin; ii++) sfh_Mass[ii]=fraction*Gal[q].sfh_BurstMass[ii];
#endif
    }
#endif
  else {
    printf("\nparent call from: %s, line %d\nUnknown component type %s in call to transfer_material\n",call_function, call_line,cq);
    exit(1);
  }






  //Add to galaxy p
  if (strcmp(cp,"ColdGas")==0)
    {
      Gal[p].ColdGas += Mass;
      Gal[p].MetalsColdGas = metals_add(Gal[p].MetalsColdGas,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[p].ColdGas_elements[kk] += Yield[kk];
#endif
    }

  else if (strcmp(cp,"HotGas")==0)
    {
      Gal[p].HotGas += Mass;
      Gal[p].MetalsHotGas = metals_add(Gal[p].MetalsHotGas,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[p].HotGas_elements[kk] += Yield[kk];
#endif
#ifdef METALS_SELF
      if (p==q) Gal[p].MetalsHotGasSelf = metals_add(Gal[p].MetalsHotGasSelf,Metals,1.);
#endif
    }

  /*else if (strcmp(cp,"ReheatedGas")==0)
      {
        Gal[p].ReheatedGas += Mass;
        Gal[p].MetalsReheatedGas = metals_add(Gal[p].MetalsReheatedGas,Metals,1.);
  #ifdef INDIVIDUAL_ELEMENTS
        for(kk=0;kk<NUM_ELEMENTS;kk++)
  	Gal[p].ReheatedGas_elements[kk] += Yield[kk];
  #endif
      }*/


  else if (strcmp(cp,"EjectedMass")==0)
    {
      Gal[p].EjectedMass += Mass;
      Gal[p].MetalsEjectedMass = metals_add(Gal[p].MetalsEjectedMass,Metals,1.);
  #ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[p].EjectedMass_elements[kk] += Yield[kk];
  #endif
    }

#ifdef EXCESS_MASS
  else if (strcmp(cp,"ExcessMass")==0)
    {
      Gal[p].ExcessMass += Mass;
      Gal[p].MetalsExcessMass = metals_add(Gal[p].MetalsExcessMass,Metals,1.);
  #ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[p].ExcessMass_elements[kk] += Yield[kk];
  #endif
    }
#endif

  else if (strcmp(cp,"BlackHoleMass")==0)
    {
      Gal[p].BlackHoleMass += Mass;
  /*    Gal[p].MetalsBlackHoleMass = metals_add(Gal[p].MetalsBlackHoleMass,Metals,1.);
  #ifdef INDIVIDUAL_ELEMENTS
  for(kk=0;kk<NUM_ELEMENTS;kk++)
      Gal[p].BlackHoleMass_elements[kk] += Yield[kk];
  #endif*/
    }

  else if (strcmp(cp,"BlackHoleGas")==0)
    {
      Gal[p].BlackHoleGas += Mass;
   /*   Gal[p].MetalsBlackHoleGas = metals_add(Gal[p].MetalsBlackHoleGas,Metals,1.);
  #ifdef INDIVIDUAL_ELEMENTS
  for(kk=0;kk<NUM_ELEMENTS;kk++)
        Gal[p].BlackHoleMass_elements[kk] += Yield[kk];
  #endif*/
    }

  else if (strcmp(cp,"DiskMass")==0)
    {
      Gal[p].DiskMass += Mass;
      Gal[p].MetalsDiskMass=metals_add(Gal[p].MetalsDiskMass,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[p].DiskMass_elements[kk] += Yield[kk];
#endif

#ifdef STAR_FORMATION_HISTORY
      for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
	{
	  Gal[p].sfh_DiskMass[ii] += sfh_Mass[ii];
	  Gal[p].sfh_MetalsDiskMass[ii]=metals_add(Gal[p].sfh_MetalsDiskMass[ii],sfh_Metals[ii],1.);
#ifdef INDIVIDUAL_ELEMENTS
	  for(kk=0;kk<NUM_ELEMENTS;kk++)
	    Gal[p].sfh_DiskMass_elements[ii][kk] += sfh_Elements[ii][kk];
#endif
	}
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
      for(outputbin = 0; outputbin < NOUT; outputbin++)
	{
	  for(ll = 0; ll < NMAG; ll++)
	    {
#ifdef OUTPUT_REST_MAGS
	      Gal[p].Lum[ll][outputbin]+=Lum[ll][outputbin];
	      Gal[p].YLum[ll][outputbin]+=YLum[ll][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
	      Gal[p].ObsLum[ll][outputbin]+=ObsLum[ll][outputbin];
	      Gal[p].ObsYLum[ll][outputbin]+=ObsYLum[ll][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	      Gal[p].dObsLum[ll][outputbin]+=dObsLum[ll][outputbin];
	      Gal[p].dObsYLum[ll][outputbin]+=dObsYLum[ll][outputbin];
#endif
#endif
	    }
	}
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES
    }

  else if (strcmp(cp,"BulgeMass")==0)
    {
      Gal[p].BulgeMass += Mass;
      Gal[p].MetalsBulgeMass=metals_add(Gal[p].MetalsBulgeMass,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[p].BulgeMass_elements[kk] += Yield[kk];
#endif

#ifdef STAR_FORMATION_HISTORY
      for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
	{
	  Gal[p].sfh_BulgeMass[ii] += sfh_Mass[ii];
	  Gal[p].sfh_MetalsBulgeMass[ii]=metals_add(Gal[p].sfh_MetalsBulgeMass[ii],sfh_Metals[ii],1.);
#ifdef INDIVIDUAL_ELEMENTS
	  for(kk=0;kk<NUM_ELEMENTS;kk++)
	    Gal[p].sfh_BulgeMass_elements[ii][kk] += sfh_Elements[ii][kk];
#endif
	}
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
      for(outputbin = 0; outputbin < NOUT; outputbin++)
	{
	  for(ll = 0; ll < NMAG; ll++)
	    {
#ifdef OUTPUT_REST_MAGS
	      Gal[p].Lum[ll][outputbin]+=Lum[ll][outputbin];
	      Gal[p].YLum[ll][outputbin]+=YLum[ll][outputbin];
	      Gal[p].LumBulge[ll][outputbin]+=Lum[ll][outputbin];
	      Gal[p].YLumBulge[ll][outputbin]+=YLum[ll][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
	      Gal[p].ObsLum[ll][outputbin]+=ObsLum[ll][outputbin];
	      Gal[p].ObsYLum[ll][outputbin]+=ObsYLum[ll][outputbin];
	      Gal[p].ObsLumBulge[ll][outputbin]+=ObsLum[ll][outputbin];
	      Gal[p].ObsYLumBulge[ll][outputbin]+=ObsYLum[ll][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	      Gal[p].dObsLum[ll][outputbin]+=dObsLum[ll][outputbin];
	      Gal[p].dObsYLum[ll][outputbin]+=dObsYLum[ll][outputbin];
	      Gal[p].dObsLumBulge[ll][outputbin]+=dObsLum[ll][outputbin];
	      Gal[p].dObsYLumBulge[ll][outputbin]+=dObsYLum[ll][outputbin];
#endif
#endif
	    }
	}
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES
    }

  else if (strcmp(cp,"ICM")==0)
    {
      Gal[p].ICM += Mass;
      Gal[p].MetalsICM=metals_add(Gal[p].MetalsICM,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[p].ICM_elements[kk] += Yield[kk];
#endif

#ifdef STAR_FORMATION_HISTORY
      for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
	{
	  Gal[p].sfh_ICM[ii] += sfh_Mass[ii];
	  Gal[p].sfh_MetalsICM[ii]=metals_add(Gal[p].sfh_MetalsICM[ii],sfh_Metals[ii],1.);
#ifdef INDIVIDUAL_ELEMENTS
	  for(kk=0;kk<NUM_ELEMENTS;kk++)
	    Gal[p].sfh_ICM_elements[ii][kk] += sfh_Elements[ii][kk];
#endif
	}
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
#ifdef ICL
      for(outputbin = 0; outputbin < NOUT; outputbin++)
	{
	  for(ll = 0; ll < NMAG; ll++)
	    {
#ifdef OUTPUT_REST_MAGS
	      Gal[p].ICLLum[ll][outputbin]+=Lum[ll][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
	      Gal[p].ObsICLLum[ll][outputbin]+=ObsLum[ll][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	      Gal[p].dObsICLLum[ll][outputbin]+=dObsLum[ll][outputbin];
#endif
#endif
	    }
	}
#endif
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES
  }

#ifdef TRACK_BURST
  else if (strcmp(cp,"BurstMass")==0)
    {
      Gal[p].BurstMass += Mass;
#ifdef STAR_FORMATION_HISTORY
      for (ii=0; ii<=Gal[p].sfh_ibin; ii++) Gal[p].sfh_BurstMass[ii] += sfh_Mass[ii];
#endif
    }
#endif

  else {
    printf("\nparent call from: %s, line %d\nUnknown component type %s in call to transfer_material\n",call_function, call_line,cp);
    exit(1);
  }





  //Subtract from galaxy q;
  if (strcmp(cq,"ColdGas")==0)
    {
      Gal[q].ColdGas -= Mass;
      Gal[q].MetalsColdGas = metals_add(Gal[q].MetalsColdGas,Metals,-1.);
 #ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[q].ColdGas_elements[kk] -= Yield[kk];
 #endif
    }

  else if (strcmp(cq,"HotGas")==0)
    {
      Gal[q].HotGas -= Mass;
      Gal[q].MetalsHotGas = metals_add(Gal[q].MetalsHotGas,Metals,-1.);
 #ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[q].HotGas_elements[kk] -= Yield[kk];
 #endif
 #ifdef METALS_SELF
      Gal[q].MetalsHotGasSelf = metals_add(Gal[q].MetalsHotGasSelf,Metals,-1.);
 #endif
    }

  /*else if (strcmp(cq,"ReheatedGas")==0)
      {
        Gal[q].ReheatedGas -= Mass;
        Gal[q].MetalsReheatedGas = metals_add(Gal[q].MetalsReheatedGas,Metals,-1.);
   #ifdef INDIVIDUAL_ELEMENTS
        for(kk=0;kk<NUM_ELEMENTS;kk++)
  	Gal[q].ReheatedGas_elements[kk] -= Yield[kk];
   #endif
      }*/

  else if (strcmp(cq,"EjectedMass")==0)
    {
      Gal[q].EjectedMass -= Mass;
      Gal[q].MetalsEjectedMass = metals_add(Gal[q].MetalsEjectedMass,Metals,-1.);
 #ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[q].EjectedMass_elements[kk] -= Yield[kk];
 #endif
    }

#ifdef EXCESS_MASS
  else if (strcmp(cq,"ExcessMass")==0)
    {
      Gal[q].ExcessMass -= Mass;
      Gal[q].MetalsExcessMass = metals_add(Gal[q].MetalsExcessMass,Metals,-1.);
 #ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[q].ExcessMass_elements[kk] -= Yield[kk];
 #endif
    }
#endif

  else if (strcmp(cq,"DiskMass")==0)
    {
      Gal[q].DiskMass -= Mass;
      Gal[q].MetalsDiskMass=metals_add(Gal[q].MetalsDiskMass,Metals,-1.);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[q].DiskMass_elements[kk] -= Yield[kk];
#endif

#ifdef STAR_FORMATION_HISTORY
      for (ii=0; ii<=Gal[q].sfh_ibin; ii++)
	{
	  Gal[q].sfh_DiskMass[ii] -= sfh_Mass[ii];
	  Gal[q].sfh_MetalsDiskMass[ii]=metals_add(Gal[q].sfh_MetalsDiskMass[ii],sfh_Metals[ii],-1.);
#ifdef INDIVIDUAL_ELEMENTS
	  for(kk=0;kk<NUM_ELEMENTS;kk++)
	    Gal[q].sfh_DiskMass_elements[ii][kk] -= sfh_Elements[ii][kk];
#endif
	}
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
      for(outputbin = 0; outputbin < NOUT; outputbin++)
	{
	  for(ll = 0; ll < NMAG; ll++)
	    {
#ifdef OUTPUT_REST_MAGS
	      Gal[q].Lum[ll][outputbin]-=Lum[ll][outputbin];
	      Gal[q].YLum[ll][outputbin]-=YLum[ll][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
	      Gal[q].ObsLum[ll][outputbin]-=ObsLum[ll][outputbin];
	      Gal[q].ObsYLum[ll][outputbin]-=ObsYLum[ll][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	      Gal[q].dObsLum[ll][outputbin]-=dObsLum[ll][outputbin];
	      Gal[q].dObsYLum[ll][outputbin]-=dObsYLum[ll][outputbin];
#endif
#endif
	    }
	}
#endif //POST_PROCESS_MAGS
#endif// COMPUTE_SPECPHOT_PROPERTIES
    }

  else if (strcmp(cq,"BulgeMass")==0)
    {
      Gal[q].BulgeMass -= Mass;
      Gal[q].MetalsBulgeMass=metals_add(Gal[q].MetalsBulgeMass,Metals,-1.);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[q].BulgeMass_elements[kk] -= Yield[kk];
#endif

#ifdef STAR_FORMATION_HISTORY
      for (ii=0; ii<=Gal[q].sfh_ibin; ii++)
	{
	  Gal[q].sfh_BulgeMass[ii] -= sfh_Mass[ii];
	  Gal[q].sfh_MetalsBulgeMass[ii]=metals_add(Gal[q].sfh_MetalsBulgeMass[ii],sfh_Metals[ii],-1.);
#ifdef INDIVIDUAL_ELEMENTS
	  for(kk=0;kk<NUM_ELEMENTS;kk++)
	    Gal[q].sfh_BulgeMass_elements[ii][kk] -= sfh_Elements[ii][kk];
#endif
	}
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
      for(outputbin = 0; outputbin < NOUT; outputbin++)
	{
	  for(ll = 0; ll < NMAG; ll++)
	    {
#ifdef OUTPUT_REST_MAGS
	      Gal[q].Lum[ll][outputbin]-=Lum[ll][outputbin];
	      Gal[q].YLum[ll][outputbin]-=YLum[ll][outputbin];
	      Gal[q].LumBulge[ll][outputbin]-=Lum[ll][outputbin];
	      Gal[q].YLumBulge[ll][outputbin]-=YLum[ll][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
	      Gal[q].ObsLum[ll][outputbin]-=ObsLum[ll][outputbin];
	      Gal[q].ObsYLum[ll][outputbin]-=ObsYLum[ll][outputbin];
	      Gal[q].ObsLumBulge[ll][outputbin]-=ObsLum[ll][outputbin];
	      Gal[q].ObsYLumBulge[ll][outputbin]-=ObsYLum[ll][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	      Gal[q].dObsLum[ll][outputbin]-=dObsLum[ll][outputbin];
	      Gal[q].dObsYLum[ll][outputbin]-=dObsYLum[ll][outputbin];
	      Gal[q].dObsLumBulge[ll][outputbin]-=dObsLum[ll][outputbin];
	      Gal[q].dObsYLumBulge[ll][outputbin]-=dObsYLum[ll][outputbin];
#endif
#endif
	    }
	}
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES
    }

  else if (strcmp(cq,"ICM")==0)
    {
      Gal[q].ICM -= Mass;
      Gal[q].MetalsICM=metals_add(Gal[q].MetalsICM,Metals,-1.);
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[q].ICM_elements[kk] -= Yield[kk];
#endif

#ifdef STAR_FORMATION_HISTORY
      for (ii=0; ii<=Gal[q].sfh_ibin; ii++)
	{
	  Gal[q].sfh_ICM[ii] -= sfh_Mass[ii];
	  Gal[q].sfh_MetalsICM[ii]=metals_add(Gal[q].sfh_MetalsICM[ii],sfh_Metals[ii],-1.);
#ifdef INDIVIDUAL_ELEMENTS
	  for(kk=0;kk<NUM_ELEMENTS;kk++)
	    Gal[q].sfh_ICM_elements[ii][kk] -= sfh_Elements[ii][kk];
#endif
	}
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
#ifdef ICL
      for(outputbin = 0; outputbin < NOUT; outputbin++)
	{
	  for(ll = 0; ll < NMAG; ll++)
	    {
#ifdef OUTPUT_REST_MAGS
	      Gal[q].ICLLum[ll][outputbin]-=Lum[ll][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
	      Gal[q].ObsICLLum[ll][outputbin]-=ObsLum[ll][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	      Gal[q].dObsICLLum[ll][outputbin]-=dObsLum[ll][outputbin];
#endif
#endif
	    }
	}
#endif
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES
    }

#ifdef TRACK_BURST
  else if (strcmp(cq,"BurstMass")==0)
    {
      Gal[q].BurstMass -=0.;
#ifdef STAR_FORMATION_HISTORY
      for (i=0; i<=Gal[q].sfh_ibin; i++)
	Gal[q].sfh_BurstMass[i] -= sfh_Mass[i];
#endif
    }
#endif
  else {
    printf("\nparent call from: %s, line %d\nUnknown component type %s in call to transfer_material\n",call_function, call_line,cq);
    exit(1);
  }





/*************************************************
 *
 * set SFH and to 0 if total is 0
 *
 */
  //set metals, elements and SFH to 0, if total=0

  //galaxy p cold gas
 /* if(Gal[p].ColdGas == 0.)
    {
      Gal[p].MetalsColdGas = metals_init();
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[p].ColdGas_elements[kk] = 0.;
#endif
    }
  //galaxy q cold gas
  if(Gal[q].ColdGas == 0.)
    {
      Gal[q].MetalsColdGas = metals_init();
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[q].ColdGas_elements[kk] = 0.;
#endif
    }

  //galaxy p disk mass
  if(Gal[p].DiskMass == 0.)
    {
      Gal[p].MetalsDiskMass = metals_init();
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[p].DiskMass_elements[kk] = 0.;
#endif
    }
  //galaxy q disk mass
  if(Gal[q].DiskMass == 0.)
    {
      Gal[q].MetalsDiskMass = metals_init();
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[q].DiskMass_elements[kk] = 0.;
#endif
    }

  //galaxy p bulge mass
  if(Gal[p].BulgeMass == 0.)
    {
      Gal[p].MetalsBulgeMass = metals_init();
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[p].BulgeMass_elements[kk] = 0.;
#endif
    }
  //galaxy p bulge mass
  if(Gal[q].BulgeMass == 0.)
    {
      Gal[q].MetalsBulgeMass = metals_init();
#ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[q].BulgeMass_elements[kk] = 0.;
#endif
    }*/



#ifdef CHECK_SFH_AND_RINGS
  //set SFH to 0
#ifdef STAR_FORMATION_HISTORY
  if(Gal[p].DiskMass == 0. || Gal[q].DiskMass == 0. || Gal[p].BulgeMass == 0. || Gal[q].BulgeMass == 0.)
    for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
      {
	//galaxy q disk mass
	if(Gal[p].DiskMass == 0.)
	  {
	    Gal[p].sfh_DiskMass[ii] = 0.;
	    Gal[p].sfh_MetalsDiskMass[ii]=metals_init();
#ifdef INDIVIDUAL_ELEMENTS
	    for(kk=0;kk<NUM_ELEMENTS;kk++)
	      Gal[p].sfh_DiskMass_elements[ii][kk]=0.;
#endif
	  }
	//galaxy q disk mass
	if(Gal[q].DiskMass == 0.)
	  {
	    Gal[q].sfh_DiskMass[ii] = 0.;
	    Gal[q].sfh_MetalsDiskMass[ii]=metals_init();
#ifdef INDIVIDUAL_ELEMENTS
	    for(kk=0;kk<NUM_ELEMENTS;kk++)
	      Gal[q].sfh_DiskMass_elements[ii][kk]=0.;
#endif
	  }

	//galaxy p bulge mass
	if(Gal[p].BulgeMass == 0.)
	  {
	    Gal[p].sfh_BulgeMass[ii] = 0.;
	    Gal[p].sfh_MetalsBulgeMass[ii]=metals_init();
#ifdef INDIVIDUAL_ELEMENTS
	    for(kk=0;kk<NUM_ELEMENTS;kk++)
	      Gal[p].sfh_BulgeMass_elements[ii][kk]=0.;
#endif
	  }
	//galaxy q bulge mass
	if(Gal[q].BulgeMass == 0.)
	  {
	    Gal[q].sfh_BulgeMass[ii] = 0.;
	    Gal[q].sfh_MetalsBulgeMass[ii]=metals_init();
#ifdef INDIVIDUAL_ELEMENTS
	    for(kk=0;kk<NUM_ELEMENTS;kk++)
	      Gal[q].sfh_BulgeMass_elements[ii][kk]=0.;
#endif
	  }

      }
#endif
#endif //CHECK_SFH_AND_RINGS




  return;
}





#ifdef H2_AND_RINGS
void transfer_material_with_rings(int p, char cp[], int q, char cq[], double fractionRings[], char call_function[], int call_line)
{
  //if H2_AND_RINGS the material in cold gas rings also needs to be transferred
  // there are only rings in the cold gas and Diskmass

  //Variables for TOTAL quantities
  double Mass, fraction;
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals Metals;
#ifdef INDIVIDUAL_ELEMENTS
  int kk;
  double Yield[NUM_ELEMENTS];
#endif
#else
  double Metals;
#endif //DETAILED_METALS_AND_MASS_RETURN


  //Variables for SFH
#ifdef STAR_FORMATION_HISTORY
  int ii;
  double sfh_Mass[SFH_NBIN];
  double sfh_MassRings[RNUM][SFH_NBIN];
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals sfh_Metals[SFH_NBIN];
#ifdef INDIVIDUAL_ELEMENTS
  double sfh_Elements[SFH_NBIN][NUM_ELEMENTS];
#endif
#else
  double sfh_Metals[SFH_NBIN];
#endif
#endif

//Variables for LUMINOSITIES
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  int ll, outputbin;
#ifdef OUTPUT_REST_MAGS
  float Lum[NMAG][NOUT];
  float YLum[NMAG][NOUT];
#endif
#ifdef COMPUTE_OBS_MAGS
  float ObsLum[NMAG][NOUT];
  float ObsYLum[NMAG][NOUT];
#ifdef OUTPUT_MOMAF_INPUTS
  float dObsLum[NMAG][NOUT];
  float dObsYLum[NMAG][NOUT];
#endif
#endif
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES


  //Variables for RINGS
  int jj;
  double MassRings[RNUM];
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals MetalsRings[RNUM];
#ifdef INDIVIDUAL_ELEMENTS
  double YieldRings[RNUM][NUM_ELEMENTS];
#endif
#else
  double MetalsRings[RNUM];
#endif //DETAILED_METALS_AND_MASS_RETURN




  //Initialize arrays to contain mass to transfer
   Mass = 0.;
   Metals=metals_add(metals_init(),metals_init(),1.);
#ifdef INDIVIDUAL_ELEMENTS
   for(kk=0;kk<NUM_ELEMENTS;kk++)
     Yield[kk] = 0.;
#endif

#ifdef STAR_FORMATION_HISTORY
   for (ii=0; ii<=Gal[q].sfh_ibin; ii++)
     {
       sfh_Mass[ii]=0.;
       for (jj=0;jj<RNUM;jj++)
	 sfh_MassRings[jj][ii]=0.;
       sfh_Metals[ii]=metals_add(metals_init(),metals_init(),1.);
#ifdef INDIVIDUAL_ELEMENTS
       for(kk=0;kk<NUM_ELEMENTS;kk++)
	 sfh_Elements[ii][kk]=0.;
#endif
     }
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
   for(outputbin = 0; outputbin < NOUT; outputbin++)
     {
       for(ll = 0; ll < NMAG; ll++)
 	{
#ifdef OUTPUT_REST_MAGS
 	  Lum[ll][outputbin]=0.;
 	  YLum[ll][outputbin]=0.;
#endif
#ifdef COMPUTE_OBS_MAGS
 	  ObsLum[ll][outputbin]=0.;
 	  ObsYLum[ll][outputbin]=0.;
#ifdef OUTPUT_MOMAF_INPUTS
 	  dObsLum[ll][outputbin]=0.;
 	  dObsYLum[ll][outputbin]=0.;
#endif
#endif
 	}
     }
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES




   //MASS AND METALS TO BE TRANSFERED
   if (strcmp(cq,"ColdGas")==0)
     {
       for (jj=0;jj<RNUM;jj++)
	 {
	   Mass += fractionRings[jj]*Gal[q].ColdGasRings[jj];
	   MassRings[jj] = fractionRings[jj]*Gal[q].ColdGasRings[jj];
	   Metals = metals_add(Metals,Gal[q].MetalsColdGasRings[jj],fractionRings[jj]);
	   MetalsRings[jj] = metals_add(metals_init(),Gal[q].MetalsColdGasRings[jj],fractionRings[jj]);
#ifdef INDIVIDUAL_ELEMENTS
	   for(kk=0;kk<NUM_ELEMENTS;kk++)
	     {
	       Yield[kk] += Gal[q].ColdGasRings_elements[jj][kk]*fractionRings[jj];
	       YieldRings[jj][kk] = Gal[q].ColdGasRings_elements[jj][kk]*fractionRings[jj];
	     }
#endif
	 }
	   //if there is SF, gas goes to stars into the last sfh bin
#ifdef STAR_FORMATION_HISTORY
       for (ii=0; ii<Gal[q].sfh_ibin; ii++)
	 {
	   sfh_Mass[ii]=0.;
	   for (jj=0;jj<RNUM;jj++)
	     sfh_MassRings[jj][ii]=0.;
	   sfh_Metals[ii]=metals_add(metals_init(),Gal[q].sfh_MetalsDiskMass[ii],0.);
#ifdef INDIVIDUAL_ELEMENTS
	   for(kk=0;kk<NUM_ELEMENTS;kk++)
	     sfh_Elements[ii][kk]=0.;
#endif
	 }
       for (jj=0;jj<RNUM;jj++)
	 {
	   sfh_Mass[Gal[q].sfh_ibin]+=fractionRings[jj]*Gal[q].ColdGasRings[jj];
	   sfh_MassRings[jj][Gal[q].sfh_ibin]=fractionRings[jj]*Gal[q].ColdGasRings[jj];
	   sfh_Metals[Gal[q].sfh_ibin]=metals_add(sfh_Metals[Gal[q].sfh_ibin],Gal[q].MetalsColdGasRings[jj],fractionRings[jj]);
#ifdef INDIVIDUAL_ELEMENTS
	   for(kk=0;kk<NUM_ELEMENTS;kk++)
	     sfh_Elements[Gal[q].sfh_ibin][kk]+=Gal[q].ColdGasRings_elements[jj][kk]*fractionRings[jj];
#endif
	 }
#endif //STAR_FORMATION_HISTORY
     }

   else if (strcmp(cq,"HotGas")==0)
     {
       for (jj=0;jj<RNUM;jj++)
	 {   	 	   
	   Mass += fractionRings[jj]*Gal[q].HotGas;
   	   MassRings[jj] = fractionRings[jj]*Gal[q].HotGas;
   	   Metals = metals_add(Metals,Gal[q].MetalsHotGas,fractionRings[jj]);
   	   MetalsRings[jj] = metals_add(metals_init(),Gal[q].MetalsHotGas,fractionRings[jj]);
   #ifdef INDIVIDUAL_ELEMENTS
   	   for(kk=0;kk<NUM_ELEMENTS;kk++)
   	     {
   	       Yield[kk] += Gal[q].HotGas_elements[kk]*(fractionRings[jj]);
   	       YieldRings[jj][kk] = Gal[q].HotGas_elements[kk]*(fractionRings[jj]);
   	     }
   #endif
	 }
     }

   /*else if (strcmp(cq,"ReheatedGas")==0)
        {
          for (jj=0;jj<RNUM;jj++)
   	 {   	 
   	   Mass += fractionRings[jj]*Gal[q].ReheatedGas;
      	   MassRings[jj] = fractionRings[jj]*Gal[q].ReheatedGas;
      	   Metals = metals_add(Metals,Gal[q].MetalsReheatedGas,fractionRings[jj]);
      	   MetalsRings[jj] = metals_add(metals_init(),Gal[q].MetalsReheatedGas,fractionRings[jj]);
#ifdef INDIVIDUAL_ELEMENTS
      	   for(kk=0;kk<NUM_ELEMENTS;kk++)
      	     {
      	       Yield[kk] += Gal[q].ReheatedGas_elements[kk]*(fractionRings[jj]);
      	       YieldRings[jj][kk] = Gal[q].ReheatedGas_elements[kk]*(fractionRings[jj]);
      	     }
#endif
   	 }
        }*/

   else if (strcmp(cq,"DiskMass")==0)
     {
       for (jj=0;jj<RNUM;jj++)
     	 {
	   Mass += fractionRings[jj]*Gal[q].DiskMassRings[jj];
	   MassRings[jj]=fractionRings[jj]*Gal[q].DiskMassRings[jj];
	   Metals=metals_add(Metals,Gal[q].MetalsDiskMassRings[jj],fractionRings[jj]);
	   MetalsRings[jj]=metals_add(metals_init(),Gal[q].MetalsDiskMassRings[jj],fractionRings[jj]);
#ifdef INDIVIDUAL_ELEMENTS
	   for(kk=0;kk<NUM_ELEMENTS;kk++)
	     {
	       Yield[kk] += Gal[q].DiskMassRings_elements[jj][kk]*fractionRings[jj];
	       YieldRings[jj][kk] = Gal[q].DiskMassRings_elements[jj][kk]*fractionRings[jj];
	     }
	   for (ii=0; ii<=Gal[q].sfh_ibin; ii++)
	     {
	       sfh_Mass[ii]+=fractionRings[jj]*Gal[q].sfh_DiskMassRings[jj][ii];
	       sfh_MassRings[jj][ii]=fractionRings[jj]*Gal[q].sfh_DiskMassRings[jj][ii];
	     }
#endif
     	 }

#ifdef STAR_FORMATION_HISTORY
       for (ii=0; ii<=Gal[q].sfh_ibin; ii++)
	 {
	   if(Gal[q].sfh_DiskMass[ii]>0.)
	     {
	       fraction=sfh_Mass[ii]/Gal[q].sfh_DiskMass[ii];
	       sfh_Metals[ii]=metals_add(metals_init(),Gal[q].sfh_MetalsDiskMass[ii],fraction);
#ifdef INDIVIDUAL_ELEMENTS
	       for(kk=0;kk<NUM_ELEMENTS;kk++)
		 sfh_Elements[ii][kk]=Gal[q].sfh_DiskMass_elements[ii][kk]*fraction;
#endif
	     }
	 }
#endif


#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
       for(outputbin = 0; outputbin < NOUT; outputbin++)
	 {
	   for(ll = 0; ll < NMAG; ll++)
	     {
#ifdef OUTPUT_REST_MAGS
	       Lum[ll][outputbin]=fraction*(Gal[q].Lum[ll][outputbin]-Gal[q].LumBulge[ll][outputbin]);
	       YLum[ll][outputbin]=fraction*(Gal[q].YLum[ll][outputbin]-Gal[q].YLumBulge[ll][outputbin]);
#endif
#ifdef COMPUTE_OBS_MAGS
	       ObsLum[ll][outputbin]=fraction*(Gal[q].ObsLum[ll][outputbin]-Gal[q].ObsLumBulge[ll][outputbin]);
	       ObsYLum[ll][outputbin]=fraction*(Gal[q].ObsYLum[ll][outputbin]-Gal[q].ObsYLumBulge[ll][outputbin]);
#ifdef OUTPUT_MOMAF_INPUTS
	       dObsLum[ll][outputbin]=fraction*(Gal[q].dObsLum[ll][outputbin]-Gal[q].dObsLumBulge[ll][outputbin]);
	       dObsYLum[ll][outputbin]=fraction*(Gal[q].dObsYLum[ll][outputbin]-Gal[q].dObsYLumBulge[ll][outputbin]);
#endif
#endif
	     }
	 }
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

     }


#ifdef RINGS_IN_BULGES
   else if (strcmp(cq,"BulgeMass")==0)
     {
       for (jj=0;jj<RNUM;jj++)
	 {
	   Mass += fractionRings[jj]*Gal[q].BulgeMassRings[jj];
	   MassRings[jj]=fractionRings[jj]*Gal[q].BulgeMassRings[jj];
	   Metals=metals_add(Metals,Gal[q].MetalsBulgeMassRings[jj],fractionRings[jj]);
	   MetalsRings[jj]=metals_add(metals_init(),Gal[q].MetalsBulgeMassRings[jj],fractionRings[jj]);
#ifdef INDIVIDUAL_ELEMENTS
	   for(kk=0;kk<NUM_ELEMENTS;kk++)
	     {
	       Yield[kk] += Gal[q].BulgeMassRings_elements[jj][kk]*fractionRings[jj];
	       YieldRings[jj][kk] = Gal[q].BulgeMassRings_elements[jj][kk]*fractionRings[jj];
	     }
	   for (ii=0; ii<=Gal[q].sfh_ibin; ii++)
	     {
	       sfh_Mass[ii]+=fractionRings[jj]*Gal[q].sfh_BulgeMassRings[jj][ii];
	       sfh_MassRings[jj][ii]=fractionRings[jj]*Gal[q].sfh_BulgeMassRings[jj][ii];
	     }
#endif
	 }

#ifdef STAR_FORMATION_HISTORY
       for (ii=0; ii<=Gal[q].sfh_ibin; ii++)
	 {
	   if(Gal[q].sfh_BulgeMass[ii]>0.)
	     {
	       fraction=sfh_Mass[ii]/Gal[q].sfh_BulgeMass[ii];
	       sfh_Metals[ii]=metals_add(metals_init(),Gal[q].sfh_MetalsBulgeMass[ii],fraction);
#ifdef INDIVIDUAL_ELEMENTS
	       for(kk=0;kk<NUM_ELEMENTS;kk++)
		 sfh_Elements[ii][kk]=Gal[q].sfh_BulgeMass_elements[ii][kk]*fraction;
#endif
	     }
	 }
#endif


#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
       for(outputbin = 0; outputbin < NOUT; outputbin++)
	 {
	   for(ll = 0; ll < NMAG; ll++)
	     {
#ifdef OUTPUT_REST_MAGS
	       Lum[ll][outputbin]=fraction*Gal[q].LumBulge[ll][outputbin];
	       YLum[ll][outputbin]=fraction*Gal[q].YLumBulge[ll][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
	       ObsLum[ll][outputbin]=fraction*Gal[q].ObsLumBulge[ll][outputbin];
	       ObsYLum[ll][outputbin]=fraction*Gal[q].ObsYLumBulge[ll][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	       dObsLum[ll][outputbin]=fraction*Gal[q].dObsLumBulge[ll][outputbin];
	       dObsYLum[ll][outputbin]=fraction*Gal[q].dObsYLumBulge[ll][outputbin];
#endif
#endif
	     }
	 }
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

     }
#endif

   else
     {
       printf("\nparent call from: %s, line %d\nUnknown component type %s in call to transfer_material\n",call_function, call_line,cq);
       exit(1);
     }




  //Add to galaxy p
   if (strcmp(cp,"ColdGas")==0)
     {
       Gal[p].ColdGas += Mass;
       Gal[p].MetalsColdGas = metals_add(Gal[p].MetalsColdGas,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
       for(kk=0;kk<NUM_ELEMENTS;kk++)
	 Gal[p].ColdGas_elements[kk] += Yield[kk];
#endif
       //RINGS
       for (jj=0;jj<RNUM;jj++)
	 {
	   Gal[p].ColdGasRings[jj] += MassRings[jj];
	   Gal[p].MetalsColdGasRings[jj] = metals_add(Gal[p].MetalsColdGasRings[jj],MetalsRings[jj],1.);
#ifdef INDIVIDUAL_ELEMENTS
	   for(kk=0;kk<NUM_ELEMENTS;kk++)
	     Gal[p].ColdGasRings_elements[jj][kk] += YieldRings[jj][kk];
#endif
	 }
     }

   else if (strcmp(cp,"HotGas")==0)
     {
       Gal[p].HotGas += Mass;
       Gal[p].MetalsHotGas = metals_add(Gal[p].MetalsHotGas,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
       for(kk=0;kk<NUM_ELEMENTS;kk++)
	 Gal[p].HotGas_elements[kk] += Yield[kk];
#endif
#ifdef METALS_SELF
       if (p==q) Gal[p].MetalsHotGasSelf = metals_add(Gal[p].MetalsHotGasSelf,Metals,1.);
#endif
     }

   /*else if (strcmp(cp,"ReheatedGas")==0)
       {
         Gal[p].ReheatedGas += Mass;
         Gal[p].MetalsReheatedGas = metals_add(Gal[p].MetalsReheatedGas,Metals,1.);
  #ifdef INDIVIDUAL_ELEMENTS
         for(kk=0;kk<NUM_ELEMENTS;kk++)
           Gal[p].ReheatedGas_elements[kk] += Yield[kk];
  #endif
       }*/

   else if (strcmp(cp,"BlackHoleMass")==0)
     {
       Gal[p].BlackHoleMass += Mass;
       /*    Gal[p].MetalsBlackHoleMass = metals_add(Gal[p].MetalsBlackHoleMass,Metals,1.);
    #ifdef INDIVIDUAL_ELEMENTS
      for(kk=0;kk<NUM_ELEMENTS;kk++)
        Gal[p].BlackHoleMass_elements[kk] += Yield[kk];
    #endif*/
     }

   else if (strcmp(cp,"BlackHoleGas")==0)
     {
       Gal[p].BlackHoleGas += Mass;
       /*   Gal[p].MetalsBlackHoleGas = metals_add(Gal[p].MetalsBlackHoleGas,Metals,1.);
    #ifdef INDIVIDUAL_ELEMENTS
    for(kk=0;kk<NUM_ELEMENTS;kk++)
          Gal[p].BlackHoleMass_elements[kk] += Yield[kk];
    #endif*/
     }

   else if (strcmp(cp,"DiskMass")==0)
     {
       Gal[p].DiskMass += Mass;
       Gal[p].MetalsDiskMass=metals_add(Gal[p].MetalsDiskMass,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
       for(kk=0;kk<NUM_ELEMENTS;kk++)
	 Gal[p].DiskMass_elements[kk] += Yield[kk];
#endif
       //RINGS
       for(jj=0;jj<RNUM;jj++)
	 {
	   Gal[p].DiskMassRings[jj] += MassRings[jj];
	   Gal[p].MetalsDiskMassRings[jj]=metals_add(Gal[p].MetalsDiskMassRings[jj],MetalsRings[jj],1.);
#ifdef INDIVIDUAL_ELEMENTS
	   for(kk=0;kk<NUM_ELEMENTS;kk++)
	     Gal[p].DiskMassRings_elements[jj][kk] += YieldRings[jj][kk];
#endif
	 }

#ifdef STAR_FORMATION_HISTORY
       for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
	 {
	   Gal[p].sfh_DiskMass[ii] += sfh_Mass[ii];
	   for(jj=0;jj<RNUM;jj++)
	     Gal[p].sfh_DiskMassRings[jj][ii] += sfh_MassRings[jj][ii];
	   Gal[p].sfh_MetalsDiskMass[ii]=metals_add(Gal[p].sfh_MetalsDiskMass[ii],sfh_Metals[ii],1.);
#ifdef INDIVIDUAL_ELEMENTS
	   for(kk=0;kk<NUM_ELEMENTS;kk++)
	     Gal[p].sfh_DiskMass_elements[ii][kk] += sfh_Elements[ii][kk];
#endif
	 }
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
       for(outputbin = 0; outputbin < NOUT; outputbin++)
	 {
	   for(ll = 0; ll < NMAG; ll++)
	     {
#ifdef OUTPUT_REST_MAGS
	       Gal[p].Lum[ll][outputbin]+=Lum[ll][outputbin];
	       Gal[p].YLum[ll][outputbin]+=YLum[ll][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
	       Gal[p].ObsLum[ll][outputbin]+=ObsLum[ll][outputbin]; // Mass += fractionRings[jj]/RNUM*Gal[q].HotGas;
	       Gal[p].ObsYLum[ll][outputbin]+=ObsYLum[ll][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	       Gal[p].dObsLum[ll][outputbin]+=dObsLum[ll][outputbin];
	       Gal[p].dObsYLum[ll][outputbin]+=dObsYLum[ll][outputbin];
#endif
#endif
	     }
	 }
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES
     }

   else if (strcmp(cp,"BulgeMass")==0)
     {
       Gal[p].BulgeMass += Mass;
       Gal[p].MetalsBulgeMass=metals_add(Gal[p].MetalsBulgeMass,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
       for(kk=0;kk<NUM_ELEMENTS;kk++)
	 Gal[p].BulgeMass_elements[kk] += Yield[kk];
#endif

#ifdef RINGS_IN_BULGES
       //RINGS
       for(jj=0;jj<RNUM;jj++)
      	 {
      	   Gal[p].BulgeMassRings[jj] += MassRings[jj];
      	   Gal[p].MetalsBulgeMassRings[jj]=metals_add(Gal[p].MetalsBulgeMassRings[jj],MetalsRings[jj],1.);
#ifdef INDIVIDUAL_ELEMENTS
      	   for(kk=0;kk<NUM_ELEMENTS;kk++)
      	     Gal[p].BulgeMassRings_elements[jj][kk] += YieldRings[jj][kk];
#endif
      	 }
#endif

#ifdef STAR_FORMATION_HISTORY
       for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
	 {
	   Gal[p].sfh_BulgeMass[ii] += sfh_Mass[ii];
#ifdef RINGS_IN_BULGES
	   for(jj=0;jj<RNUM;jj++)
	     Gal[p].sfh_BulgeMassRings[jj][ii] += sfh_MassRings[jj][ii];
#endif
	   Gal[p].sfh_MetalsBulgeMass[ii]=metals_add(Gal[p].sfh_MetalsBulgeMass[ii],sfh_Metals[ii],1.);
#ifdef INDIVIDUAL_ELEMENTS
	   for(kk=0;kk<NUM_ELEMENTS;kk++)
	     Gal[p].sfh_BulgeMass_elements[ii][kk] += sfh_Elements[ii][kk];
#endif
	 }
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
       for(outputbin = 0; outputbin < NOUT; outputbin++)
	 {
	   for(ll = 0; ll < NMAG; ll++)
	     {
#ifdef OUTPUT_REST_MAGS
	       Gal[p].Lum[ll][outputbin]+=Lum[ll][outputbin];
	       Gal[p].YLum[ll][outputbin]+=YLum[ll][outputbin];
	       Gal[p].LumBulge[ll][outputbin]+=Lum[ll][outputbin];
	       Gal[p].YLumBulge[ll][outputbin]+=YLum[ll][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
	       Gal[p].ObsLum[ll][outputbin]+=ObsLum[ll][outputbin];
	       Gal[p].ObsYLum[ll][outputbin]+=ObsYLum[ll][outputbin];
	       Gal[p].ObsLumBulge[ll][outputbin]+=ObsLum[ll][outputbin];
	       Gal[p].ObsYLumBulge[ll][outputbin]+=ObsYLum[ll][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	       Gal[p].dObsLum[ll][outputbin]+=dObsLum[ll][outputbin];
	       Gal[p].dObsYLum[ll][outputbin]+=dObsYLum[ll][outputbin];
	       Gal[p].dObsLumBulge[ll][outputbin]+=dObsLum[ll][outputbin];
	       Gal[p].dObsYLumBulge[ll][outputbin]+=dObsYLum[ll][outputbin];
#endif
#endif
	     }
	 }
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES
     }

   else if (strcmp(cp,"ICM")==0)
     {
       Gal[p].ICM += Mass;
       Gal[p].MetalsICM=metals_add(Gal[p].MetalsICM,Metals,1.);
#ifdef INDIVIDUAL_ELEMENTS
       for(kk=0;kk<NUM_ELEMENTS;kk++)
	 Gal[p].ICM_elements[kk] += Yield[kk];
#endif

#ifdef STAR_FORMATION_HISTORY
       for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
	 {
	   Gal[p].sfh_ICM[ii] += sfh_Mass[ii];
	   Gal[p].sfh_MetalsICM[ii]=metals_add(Gal[p].sfh_MetalsICM[ii],sfh_Metals[ii],1.);
#ifdef INDIVIDUAL_ELEMENTS
	   for(kk=0;kk<NUM_ELEMENTS;kk++)
	     Gal[p].sfh_ICM_elements[ii][kk] += sfh_Elements[ii][kk];
#endif
	 }
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
#ifdef ICL
       for(outputbin = 0; outputbin < NOUT; outputbin++)
	 {
	   for(ll = 0; ll < NMAG; ll++)
	     {
#ifdef OUTPUT_REST_MAGS
	       Gal[p].ICLLum[ll][outputbin]+=Lum[ll][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
	       Gal[p].ObsICLLum[ll][outputbin]+=ObsLum[ll][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	       Gal[p].dObsICLLum[ll][outputbin]+=dObsLum[ll][outputbin];
#endif
#endif
	     }
	 }
#endif
#endif //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES
     }




  //Subtract from galaxy q;
   if (strcmp(cq,"ColdGas")==0)
     {
       Gal[q].ColdGas -= Mass;
       Gal[q].MetalsColdGas = metals_add(Gal[q].MetalsColdGas,Metals,-1.);
#ifdef INDIVIDUAL_ELEMENTS
       for(kk=0;kk<NUM_ELEMENTS;kk++)
	 Gal[q].ColdGas_elements[kk] -= Yield[kk];
#endif
       for (jj=0;jj<RNUM;jj++)
	 {
	   Gal[q].ColdGasRings[jj] -= MassRings[jj];
	   Gal[q].MetalsColdGasRings[jj] = metals_add(Gal[q].MetalsColdGasRings[jj],MetalsRings[jj],-1.);
#ifdef INDIVIDUAL_ELEMENTS
	   for(kk=0;kk<NUM_ELEMENTS;kk++)
	     Gal[q].ColdGasRings_elements[jj][kk] -= YieldRings[jj][kk];
#endif
	 }
     }

   else if (strcmp(cq,"HotGas")==0)
     {
        Gal[q].HotGas -= Mass;
        Gal[q].MetalsHotGas = metals_add(Gal[q].MetalsHotGas,Metals,-1.);
   #ifdef INDIVIDUAL_ELEMENTS
        for(kk=0;kk<NUM_ELEMENTS;kk++)
          Gal[q].HotGas_elements[kk] -= Yield[kk];
   #endif
   #ifdef METALS_SELF
        Gal[q].MetalsHotGasSelf = metals_add(Gal[q].MetalsHotGasSelf,Metals,-1.);
   #endif
     }

   /*else if (strcmp(cq,"ReheatedGas")==0)
       {
          Gal[q].ReheatedGas -= Mass;
          Gal[q].MetalsReheatedGas = metals_add(Gal[q].MetalsReheatedGas,Metals,-1.);
     #ifdef INDIVIDUAL_ELEMENTS
          for(kk=0;kk<NUM_ELEMENTS;kk++)
            Gal[q].ReheatedGas_elements[kk] -= Yield[kk];
     #endif
       }*/

   else if (strcmp(cq,"DiskMass")==0)
     {
       Gal[q].DiskMass -= Mass;
       Gal[q].MetalsDiskMass=metals_add(Gal[q].MetalsDiskMass,Metals,-1.);
#ifdef INDIVIDUAL_ELEMENTS
       for(kk=0;kk<NUM_ELEMENTS;kk++)
	 Gal[q].DiskMass_elements[kk] -= Yield[kk];
#endif
       //RINGS
       for(jj=0;jj<RNUM;jj++)
	 {
	   Gal[q].DiskMassRings[jj] -= MassRings[jj];
	   Gal[q].MetalsDiskMassRings[jj]=metals_add(Gal[q].MetalsDiskMassRings[jj],MetalsRings[jj],-1.);
#ifdef INDIVIDUAL_ELEMENTS
	   for(kk=0;kk<NUM_ELEMENTS;kk++)
	     Gal[q].DiskMassRings_elements[jj][kk] -= YieldRings[jj][kk];
#endif
	 }

#ifdef STAR_FORMATION_HISTORY
       for (ii=0; ii<=Gal[q].sfh_ibin; ii++)
	 {
	   Gal[q].sfh_DiskMass[ii] -= sfh_Mass[ii];
	   for(jj=0;jj<RNUM;jj++)
	     Gal[q].sfh_DiskMassRings[jj][ii] -= sfh_MassRings[jj][ii];
	   Gal[q].sfh_MetalsDiskMass[ii]=metals_add(Gal[q].sfh_MetalsDiskMass[ii],sfh_Metals[ii],-1.);
#ifdef INDIVIDUAL_ELEMENTS
	   for(kk=0;kk<NUM_ELEMENTS;kk++)
	     Gal[q].sfh_DiskMass_elements[ii][kk] -= sfh_Elements[ii][kk];
#endif
	 }
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
       for(outputbin = 0; outputbin < NOUT; outputbin++)
	 {
	   for(ll = 0; ll < NMAG; ll++)
	     {
#ifdef OUTPUT_REST_MAGS
	       Gal[q].Lum[ll][outputbin]-=Lum[ll][outputbin];
	       Gal[q].YLum[ll][outputbin]-=YLum[ll][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
	       Gal[q].ObsLum[ll][outputbin]-=ObsLum[ll][outputbin];
	       Gal[q].ObsYLum[ll][outputbin]-=ObsYLum[ll][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	       Gal[q].dObsLum[ll][outputbin]-=dObsLum[ll][outputbin];
	       Gal[q].dObsYLum[ll][outputbin]-=dObsYLum[ll][outputbin];
#endif
#endif
	     }
	 }
#endif //POST_PROCESS_MAGS
#endif// COMPUTE_SPECPHOT_PROPERTIES
     }

   else if (strcmp(cq,"BulgeMass")==0)
     {
       Gal[q].BulgeMass -= Mass;
       Gal[q].MetalsBulgeMass=metals_add(Gal[q].MetalsBulgeMass,Metals,-1.);
#ifdef INDIVIDUAL_ELEMENTS
       for(kk=0;kk<NUM_ELEMENTS;kk++) // Mass += fractionRings[jj]/RNUM*Gal[q].HotGas;
	 Gal[q].BulgeMass_elements[kk] -= Yield[kk];
#endif

#ifdef RINGS_IN_BULGES
       //RINGS
       for(jj=0;jj<RNUM;jj++)
	 {
   	   Gal[q].BulgeMassRings[jj] -= MassRings[jj];
   	   Gal[q].MetalsBulgeMassRings[jj]=metals_add(Gal[q].MetalsBulgeMassRings[jj],MetalsRings[jj],-1.);
#ifdef INDIVIDUAL_ELEMENTS
   	   for(kk=0;kk<NUM_ELEMENTS;kk++)
   	     Gal[q].BulgeMassRings_elements[jj][kk] -= YieldRings[jj][kk];
#endif
   	 }
#endif

#ifdef STAR_FORMATION_HISTORY
       for (ii=0; ii<=Gal[q].sfh_ibin; ii++)
   	 {
   	   Gal[q].sfh_BulgeMass[ii] -= sfh_Mass[ii];
#ifdef RINGS_IN_BULGES
   	   for(jj=0;jj<RNUM;jj++)
   	     Gal[q].sfh_BulgeMassRings[jj][ii] -= sfh_MassRings[jj][ii];
#endif
   	   Gal[q].sfh_MetalsBulgeMass[ii]=metals_add(Gal[q].sfh_MetalsBulgeMass[ii],sfh_Metals[ii],-1.);
#ifdef INDIVIDUAL_ELEMENTS
   	   for(kk=0;kk<NUM_ELEMENTS;kk++)
   	     Gal[q].sfh_BulgeMass_elements[ii][kk] -= sfh_Elements[ii][kk];
#endif
   	 }
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
       for(outputbin = 0; outputbin < NOUT; outputbin++)
	 {
	   for(ll = 0; ll < NMAG; ll++)
	     {
#ifdef OUTPUT_REST_MAGS
   	       Gal[q].Lum[ll][outputbin]-=Lum[ll][outputbin];
   	       Gal[q].YLum[ll][outputbin]-=YLum[ll][outputbin];
   	       Gal[q].LumBulge[ll][outputbin]-=Lum[ll][outputbin];
   	       Gal[q].YLumBulge[ll][outputbin]-=YLum[ll][outputbin];
#endif
#ifdef COMPUTE_OBS_MAGS
   	       Gal[q].ObsLum[ll][outputbin]-=ObsLum[ll][outputbin];
   	       Gal[q].ObsYLum[ll][outputbin]-=ObsYLum[ll][outputbin];
   	       Gal[q].ObsLumBulge[ll][outputbin]-=ObsLum[ll][outputbin];
   	       Gal[q].ObsYLumBulge[ll][outputbin]-=ObsYLum[ll][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
   	       Gal[q].dObsLum[ll][outputbin]-=dObsLum[ll][outputbin];
   	       Gal[q].dObsYLum[ll][outputbin]-=dObsYLum[ll][outputbin];
   	       Gal[q].dObsLumBulge[ll][outputbin]-=dObsLum[ll][outputbin];
   	       Gal[q].dObsYLumBulge[ll][outputbin]-=dObsYLum[ll][outputbin];
#endif
#endif
	     }
	 }
#endif //POST_PROCESS_MAGS
#endif// COMPUTE_SPECPHOT_PROPERTIES
     }










#ifdef CHECK_SFH_AND_RINGS

/*********************************************
 *
 * Ensure that no values became negative due
 * to precision
 *
 *********************************************/

#ifdef CHECK_NO_NEGATIVE_VALUES
  //cold gas galaxy p
  if(Gal[p].ColdGas < 0.0)
    if(Gal[p].ColdGas > -1e-6)
      {
	Gal[p].ColdGas = 0.;
	Gal[p].MetalsColdGas = metals_init();
#ifdef INDIVIDUAL_ELEMENTS
	for(kk=0;kk<NUM_ELEMENTS;kk++)
	  Gal[q].ColdGas_elements[kk] = 0.;
#endif
#ifdef H2_AND_RINGS
	for (jj=0; jj<RNUM; jj++)
   	  {
   	    Gal[p].ColdGasRings[jj]=0.;
   	    Gal[p].MetalsColdGasRings[jj] = metals_init();;
#ifdef INDIVIDUAL_ELEMENTS
   	 for(kk=0;kk<NUM_ELEMENTS;kk++)
   	    Gal[p].ColdGasRings_elements[jj][kk] = 0.;
#endif
   	  }
#endif
      }

  //cold gas galaxy q
  if(Gal[q].ColdGas < 0.0)
    if(Gal[q].ColdGas > -1e-6)
      {
	Gal[q].ColdGas = 0.;
	Gal[q].MetalsColdGas = metals_init();
 #ifdef INDIVIDUAL_ELEMENTS
	for(kk=0;kk<NUM_ELEMENTS;kk++)
	  Gal[q].ColdGas_elements[kk] = 0.;
 #endif
 #ifdef H2_AND_RINGS
	for (jj=0; jj<RNUM; jj++)
 	  {
 	    Gal[q].ColdGasRings[jj]=0.;
 	    Gal[q].MetalsColdGasRings[jj] = metals_init();;
#ifdef INDIVIDUAL_ELEMENTS
 	   for(kk=0;kk<NUM_ELEMENTS;kk++)
 	     Gal[q].ColdGasRings_elements[jj][kk] = 0.;
#endif
 	  }
#endif
      }

  //disk mass galaxy p
  if(Gal[p].DiskMass < 0.0)
    if(Gal[p].DiskMass > -1e-6)
      {
  	Gal[p].DiskMass = 0.;
  	Gal[p].MetalsDiskMass = metals_init();
#ifdef INDIVIDUAL_ELEMENTS
  	for(kk=0;kk<NUM_ELEMENTS;kk++)
  	Gal[q].DiskMass_elements[kk] = 0.;
#endif
#ifdef H2_AND_RINGS
  	for (jj=0; jj<RNUM; jj++)
  	  {
  	    Gal[p].DiskMassRings[jj]=0.;
  	    Gal[p].MetalsDiskMassRings[jj] = metals_init();;
#ifdef INDIVIDUAL_ELEMENTS
  	  for(kk=0;kk<NUM_ELEMENTS;kk++)
  	    Gal[p].DiskMassRings_elements[jj][kk] = 0.;
#endif
  	  }
#endif
        }

  //disk mass galaxy q
  if(Gal[q].DiskMass < 0.0)
    if(Gal[q].DiskMass > -1e-6)
      {
	Gal[q].DiskMass = 0.;
	Gal[q].MetalsDiskMass = metals_init();
#ifdef INDIVIDUAL_ELEMENTS
	for(kk=0;kk<NUM_ELEMENTS;kk++)
	  Gal[q].DiskMass_elements[kk] = 0.;
#endif
#ifdef H2_AND_RINGS
	for (jj=0; jj<RNUM; jj++)
	  {
	    Gal[q].DiskMassRings[jj]=0.; // Mass += fractionRings[jj]/RNUM*Gal[q].HotGas;
	    Gal[q].MetalsDiskMassRings[jj] = metals_init();;
#ifdef INDIVIDUAL_ELEMENTS
	    for(kk=0;kk<NUM_ELEMENTS;kk++)
	    Gal[q].DiskMassRings_elements[jj][kk] = 0.;
#endif
	  }
#endif
      }



  //bulge mass galaxy p
  if(Gal[p].BulgeMass < 0.0)
    if(Gal[p].BulgeMass > -1e-6)
      {
	Gal[p].BulgeMass = 0.;
	Gal[p].MetalsBulgeMass = metals_init();
#ifdef INDIVIDUAL_ELEMENTS
    	for(kk=0;kk<NUM_ELEMENTS;kk++)
    	  Gal[q].BulgeMass_elements[kk] = 0.;
#endif
#ifdef H2_AND_RINGS
#ifdef RINGS_IN_BULGES
    	for (jj=0; jj<RNUM; jj++)
    	  {
    	    Gal[p].BulgeMassRings[jj]=0.;
    	    Gal[p].MetalsBulgeMassRings[jj] = metals_init();;
#ifdef INDIVIDUAL_ELEMENTS
    	    for(kk=0;kk<NUM_ELEMENTS;kk++)
    	      Gal[p].BulgeMassRings_elements[jj][kk] = 0.;
#endif
    	  }
#endif
#endif
      }

  //Bulge mass galaxy q
  if(Gal[q].BulgeMass < 0.0)
    if(Gal[q].BulgeMass > -1e-6)
      {
  	Gal[q].BulgeMass = 0.;
  	Gal[q].MetalsBulgeMass = metals_init();
#ifdef INDIVIDUAL_ELEMENTS
  	for(kk=0;kk<NUM_ELEMENTS;kk++)
  	  Gal[q].BulgeMass_elements[kk] = 0.;
#endif
#ifdef H2_AND_RINGS
#ifdef RINGS_IN_BULGES
  	for (jj=0; jj<RNUM; jj++)
  	  {
  	    Gal[q].BulgeMassRings[jj]=0.;
  	    Gal[q].MetalsBulgeMassRings[jj] = metals_init();;
#ifdef INDIVIDUAL_ELEMENTS
  	    for(kk=0;kk<NUM_ELEMENTS;kk++)
  	      Gal[q].BulgeMassRings_elements[jj][kk] = 0.;
#endif
  	  }
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
   if(Gal[p].ColdGas == 0. || Gal[q].ColdGas == 0. || Gal[p].DiskMass == 0. || Gal[q].DiskMass == 0. || Gal[p].BulgeMass == 0. || Gal[q].BulgeMass == 0.)
     for (jj=0; jj<RNUM; jj++)
       {
 	//galaxy p cold gas
 	if(Gal[p].ColdGas == 0.)
 	  {
 	    Gal[p].ColdGasRings[jj]=0.;
 	    Gal[p].MetalsColdGasRings[jj] = metals_init();;
 #ifdef INDIVIDUAL_ELEMENTS
 	   for(kk=0;kk<NUM_ELEMENTS;kk++)
 	    Gal[p].ColdGasRings_elements[jj][kk] = 0.;
 #endif
 	  }
 	//galaxy q cold gas
 	if(Gal[q].ColdGas == 0.)
 	  {
 	    Gal[q].ColdGasRings[jj]=0.;
 	    Gal[q].MetalsColdGasRings[jj] = metals_init();;
 #ifdef INDIVIDUAL_ELEMENTS
 	   for(kk=0;kk<NUM_ELEMENTS;kk++)
 	    Gal[q].ColdGasRings_elements[jj][kk] = 0.;
 #endif
 	  }
 	//galaxy p disk mass
 	if(Gal[p].DiskMass == 0.)
 	  {
 	    Gal[p].DiskMassRings[jj]=0.;
 	    Gal[p].MetalsDiskMassRings[jj] = metals_init();;
 #ifdef INDIVIDUAL_ELEMENTS
 	   for(kk=0;kk<NUM_ELEMENTS;kk++)
 	    Gal[p].DiskMassRings_elements[jj][kk] = 0.;
 #endif
 	  }
 	//galaxy q disk mass
 	if(Gal[q].DiskMass == 0.)
 	  {
 	    Gal[q].DiskMassRings[jj]=0.;
 	    Gal[q].MetalsDiskMassRings[jj] = metals_init();;
 #ifdef INDIVIDUAL_ELEMENTS
 	   for(kk=0;kk<NUM_ELEMENTS;kk++)
 	    Gal[q].DiskMassRings_elements[jj][kk] = 0.;
 #endif
 	  }

#ifdef RINGS_IN_BULGES
 	//galaxy p Bulge mass
 	if(Gal[p].BulgeMass == 0.)
 	  {
 	    Gal[p].BulgeMassRings[jj]=0.;
 	    Gal[p].MetalsBulgeMassRings[jj] = metals_init();;
#ifdef INDIVIDUAL_ELEMENTS
 	    for(kk=0;kk<NUM_ELEMENTS;kk++)
 	      Gal[p].BulgeMassRings_elements[jj][kk] = 0.;
#endif
 	  }
 	//galaxy q Bulge mass
 	if(Gal[q].BulgeMass == 0.)
 	  {
 	    Gal[q].BulgeMassRings[jj]=0.;
 	    Gal[q].MetalsBulgeMassRings[jj] = metals_init();;
#ifdef INDIVIDUAL_ELEMENTS
 	    for(kk=0;kk<NUM_ELEMENTS;kk++)
 	      Gal[q].BulgeMassRings_elements[jj][kk] = 0.;
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

  float sum_rings, sum_rings_metals, ring_sum_minus_tot;
  float fract, fract_metals;

  //galaxy p cold gas
  if(Gal[p].ColdGas>0.)
    {
      sum_rings=0.;
      sum_rings_metals=0;
      for (jj=0; jj<RNUM; jj++)
	{
	  sum_rings+=Gal[p].ColdGasRings[jj];
	  sum_rings_metals+= metals_total(Gal[p].MetalsColdGasRings[jj]);
	} // Mass += fractionRings[jj]/RNUM*Gal[q].HotGas;

      if(sum_rings>0.)
	fract=Gal[p].ColdGas/sum_rings;
      else
	fract=0.;

      /*if(fract>0. && Gal[p].ColdGas>1e-6 && sum_rings>1e-6 && (fract>1.3 || fract<0.7))
      	{
      	  printf("fract=%0.5f %s %d\n",fract,call_function, call_line);
      	  terminate("Total and Ring Sum for Metals too different ColdGas");
      	}*/

      if(sum_rings_metals>0.)
      	fract_metals=metals_total(Gal[p].MetalsColdGas)/sum_rings_metals;
      else
	fract_metals=0.;

     /* if(fract_metals>0. && metals_total(Gal[p].MetalsColdGas)>1e-5 && sum_rings_metals>1e-5 && (fract_metals>1.3 || fract_metals<0.7))
      	{ // Mass += fractionRings[jj]/RNUM*Gal[q].HotGas;
      	  printf("ColdGasTotal=%e ColdGasRings=%e\n MetalsTotal=%e MetalsSum=%e fract=%0.5f %s %d\n",
		 Gal[p].ColdGas, sum_rings, metals_total(Gal[p].MetalsColdGas), sum_rings_metals, fract_metals,call_function, call_line);
      	  terminate("Total and Ring Sum for Metals too different ColdGas");
      	}*/

      if(sum_rings/Gal[p].ColdGas !=1.)
	{
	  for (jj=0;jj<RNUM;jj++)
	    {
	      Gal[p].ColdGasRings[jj] *= fract;
	      Gal[p].MetalsColdGasRings[jj] = metals_add(metals_init(),Gal[p].MetalsColdGasRings[jj],fract_metals);
/*#ifdef INDIVIDUAL_ELEMENTS
	      for(kk=0;kk<NUM_ELEMENTS;kk++)
		Gal[p].ColdGasRings_elements[jj][kk] = Gal[p].ColdGasRings_elements[jj][kk]*fract_metals;
#endif*/
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

  //galaxy q cold gas
  if(Gal[q].ColdGas >0.)
    {
      sum_rings=0.;
      sum_rings_metals=0;
      for (jj=0; jj<RNUM; jj++)
	{
	  sum_rings+=Gal[q].ColdGasRings[jj];
	  sum_rings_metals+= metals_total(Gal[q].MetalsColdGasRings[jj]);
     	}

      if(sum_rings>0.)
    	fract=Gal[q].ColdGas/sum_rings;
      else
    	fract=0.;
      /*if(fract>0. && Gal[q].ColdGas>1e-6 && sum_rings>1e-6  && (fract>1.3 || fract<0.7))
      	{
      	  printf("fract=%0.5f %s %d\n",fract,call_function, call_line);
      	  terminate("Total and Ring Sum for Metals too different ColdGas");
      	}*/

      if(sum_rings_metals>0.)
	fract_metals=metals_total(Gal[q].MetalsColdGas)/sum_rings_metals;
      else
	fract_metals=0.;
     /* if(fract_metals>0. && metals_total(Gal[q].MetalsColdGas)>1e-5 && sum_rings_metals>1e-5 && (fract_metals>1.3 || fract_metals<0.7))
	{
	  printf("ColdGasTotal=%e ColdGasRings=%e\n MetalsTotal=%e MetalsSum=%e fract=%0.5f %s %d\n",
	 		 Gal[q].ColdGas, sum_rings, metals_total(Gal[q].MetalsColdGas), sum_rings_metals, fract_metals,call_function, call_line);
	  terminate("Total and Ring Sum for Metals too different ColdGas");
	}*/
 // Mass += fractionRings[jj]/RNUM*Gal[q].HotGas;
      if(sum_rings/Gal[q].ColdGas !=1.)
	{
	  for (jj=0;jj<RNUM;jj++)
	    {
	      Gal[q].ColdGasRings[jj] *= fract;
	      Gal[q].MetalsColdGasRings[jj] = metals_add(metals_init(),Gal[q].MetalsColdGasRings[jj],fract_metals);
 /*#ifdef INDIVIDUAL_ELEMENTS
	      for(kk=0;kk<NUM_ELEMENTS;kk++)
		Gal[q].ColdGas_elementsRings[jj][kk] = Gal[q].ColdGas_elementsRings[jj][kk]*fract_metals;
 #endif*/
	    }

	  ring_sum_minus_tot=-Gal[q].ColdGas;
	  for (jj=0; jj<RNUM; jj++)
	    ring_sum_minus_tot+=Gal[q].ColdGasRings[jj];
	  if((ring_sum_minus_tot < -1e-4 && ring_sum_minus_tot < -1e-4*Gal[q].ColdGas) ||
	      (ring_sum_minus_tot >  1e-4 && ring_sum_minus_tot >  1e-4*Gal[q].ColdGas))
	    {
	      printf("ring_sum_minus_tot=%0.10e 1e-5*ColdGas=%0.10e\n", ring_sum_minus_tot, 1e-5*Gal[q].ColdGas);
	      terminate("");
	    }
	}
    }

   //galaxy p disk mass
  if(Gal[p].DiskMass>0.)
    {
      sum_rings=0.;
      sum_rings_metals=0;
      for (jj=0; jj<RNUM; jj++)
	{
	  sum_rings+=Gal[p].DiskMassRings[jj];
	  sum_rings_metals+= metals_total(Gal[p].MetalsDiskMassRings[jj]);
     	}

      if(sum_rings>0.)
	fract=Gal[p].DiskMass/sum_rings;
      else
	fract=0.;
      /*if(fract>0. && Gal[p].DiskMass>1e-6 && sum_rings>1e-6  && (fract>1.3 || fract<0.7))
      	{
      	  printf("fract=%0.5f %s %d\n",fract,call_function, call_line);
      	  terminate("Total and Ring Sum for Metals too different Disk Mass");
      	}*/

      if(sum_rings_metals>0.)
      	fract_metals=metals_total(Gal[p].MetalsDiskMass)/sum_rings_metals;
            else
      	fract_metals=0.;
     /* if(fract_metals>0. && metals_total(Gal[p].MetalsDiskMass)>1e-5 && sum_rings_metals>1e-5 && (fract_metals>1.3 || fract_metals<0.7))
      	{
	  for(jj=0;jj<RNUM;jj++)
	  	 printf("[%d] DiskMass=%e MetalsDisk=%e\n",jj,Gal[p].DiskMassRings[jj],Gal[p].MetalsDiskMassRings[jj]);

	  printf("ColdGasToal=%e ColdGasRings=%e\n MetalsTotal=%e MetalsSum=%e fract=%0.5f %s %d\n",
	 		 Gal[p].DiskMass, sum_rings, metals_total(Gal[p].MetalsDiskMass), sum_rings_metals, fract_metals,call_function, call_line);
      	  terminate("Total and Ring Sum for Metals too different DiskMass");
      	}*/

      if(sum_rings/Gal[p].DiskMass !=1.)
	{
	  for (jj=0;jj<RNUM;jj++)
	    {
	      Gal[p].DiskMassRings[jj] *= fract;
	      Gal[p].MetalsDiskMassRings[jj] = metals_add(metals_init(),Gal[p].MetalsDiskMassRings[jj],fract_metals);
/*#ifdef INDIVIDUAL_ELEMENTS
	      for(kk=0;kk<NUM_ELEMENTS;kk++)
		Gal[p].DiskMassRings_elements[jj][kk] = Gal[p].DiskMassRings_elements[jj][kk]*fract_metals;
#endif*/
	    }

	  ring_sum_minus_tot=-Gal[p].DiskMass;
	  for (jj=0; jj<RNUM; jj++)
	    ring_sum_minus_tot+=Gal[p].DiskMassRings[jj];
	  if((ring_sum_minus_tot < -1e-5 && ring_sum_minus_tot < -1e-5*Gal[p].DiskMass) ||
	      (ring_sum_minus_tot >  1e-5 && ring_sum_minus_tot >  1e-5*Gal[p].DiskMass))
	    {
	      printf("ring_sum_minus_tot=%0.10e 1e-5*ColdGas=%0.10e\n", ring_sum_minus_tot, 1e-5*Gal[p].DiskMass);
	      terminate("");
	    }
	}
    }



   //galaxy q disk mass
  if(Gal[q].DiskMass>0.)
    {
      sum_rings=0.;
      sum_rings_metals=0;
      for (jj=0; jj<RNUM; jj++)
	{
	  sum_rings+=Gal[q].DiskMassRings[jj];
	  sum_rings_metals+= metals_total(Gal[q].MetalsDiskMassRings[jj]);
     	}

      if(sum_rings>0.)
     	fract=Gal[q].DiskMass/sum_rings;
      else
     	fract=0.;
      /*if(fract>0. && Gal[q].DiskMass>1e-6 && sum_rings>1e-6   && (fract>1.3 || fract<0.7))
      	{
      	  printf("fract=%0.5f %s %d\n",fract,call_function, call_line);
      	  terminate("Total and Ring Sum for Metals too different DiskMass");
      	}*/

      if(sum_rings_metals>0.)
	fract_metals=metals_total(Gal[q].MetalsDiskMass)/sum_rings_metals;
      else
	fract_metals=0.;
    /*  if(fract_metals>0. && metals_total(Gal[q].MetalsDiskMass)>1e-5 && sum_rings_metals>1e-5 && (fract_metals>1.3 || fract_metals<0.7))
	{
	  printf("ColdGasToal=%e ColdGasRings=%e\n MetalsTotal=%e MetalsSum=%e fract=%0.5f %s %d\n",
	  	 		 Gal[q].DiskMass, sum_rings, metals_total(Gal[q].MetalsDiskMass), sum_rings_metals, fract_metals,call_function, call_line);
	  terminate("Total and Ring Sum for Metals too different DiskMass");
	}*/


      if(sum_rings/Gal[q].DiskMass !=1.)
	{
	  for (jj=0;jj<RNUM;jj++)
	    {
	      Gal[q].DiskMassRings[jj] *= fract;
	      Gal[q].MetalsDiskMassRings[jj] = metals_add(metals_init(),Gal[q].MetalsDiskMassRings[jj],fract_metals);
/*#ifdef INDIVIDUAL_ELEMENTS
	      for(kk=0;kk<NUM_ELEMENTS;kk++)
		Gal[q].DiskMassRings_elements[jj][kk] = Gal[q].DiskMassRings_elements[jj][kk]*fract_metals;
#endif*/
	    }

	  ring_sum_minus_tot=-Gal[q].DiskMass;
	  for (jj=0; jj<RNUM; jj++)
	    ring_sum_minus_tot+=Gal[q].DiskMassRings[jj];
	  if((ring_sum_minus_tot < -1e-5 && ring_sum_minus_tot < -1e-5*Gal[q].DiskMass) ||
	      (ring_sum_minus_tot >  1e-5 && ring_sum_minus_tot >  1e-5*Gal[q].DiskMass))
	    {
	      printf("ring_sum_minus_tot=%0.10e 1e-5*ColdGas=%0.10e\n", ring_sum_minus_tot, 1e-5*Gal[q].DiskMass);
	      terminate("");
	    }
	}
    }





#ifdef RINGS_IN_BULGES
  //galaxy p Bulge mass
   if(Gal[p].BulgeMass>0.)
     {
       sum_rings=0.;
       sum_rings_metals=0;
       for (jj=0; jj<RNUM; jj++)
 	{
 	  sum_rings+=Gal[p].BulgeMassRings[jj];
 	  sum_rings_metals+= metals_total(Gal[p].MetalsBulgeMassRings[jj]);
      	}

       if(sum_rings>0.)
 	fract=Gal[p].BulgeMass/sum_rings;
       else
 	fract=0.;
       /*if(fract>0. && Gal[p].BulgeMass>1e-6 && sum_rings>1e-6  && (fract>1.3 || fract<0.7))
       	{
       	  printf("fract=%0.5f %s %d\n",fract,call_function, call_line);
       	  terminate("Total and Ring Sum for Metals too different Bulge Mass");
       	}*/

       if(sum_rings_metals>0.)
       	fract_metals=metals_total(Gal[p].MetalsBulgeMass)/sum_rings_metals;
             else
       	fract_metals=0.;
      /* if(fract_metals>0. && metals_total(Gal[p].MetalsBulgeMass)>1e-5 && sum_rings_metals>1e-5 && (fract_metals>1.3 || fract_metals<0.7))
       	{
 	  for(jj=0;jj<RNUM;jj++)
 	  	 printf("[%d] BulgeMass=%e MetalsBulge=%e\n",jj,Gal[p].BulgeMassRings[jj],Gal[p].MetalsBulgeMassRings[jj]);

 	  printf("ColdGasToal=%e ColdGasRings=%e\n MetalsTotal=%e MetalsSum=%e fract=%0.5f %s %d\n",
 	 		 Gal[p].BulgeMass, sum_rings, metals_total(Gal[p].MetalsBulgeMass), sum_rings_metals, fract_metals,call_function, call_line);
       	  terminate("Total and Ring Sum for Metals too different BulgeMass");
       	}*/

       if(sum_rings/Gal[p].BulgeMass !=1.)
 	{
 	  for (jj=0;jj<RNUM;jj++)
 	    {
 	      Gal[p].BulgeMassRings[jj] *= fract;
 	      Gal[p].MetalsBulgeMassRings[jj] = metals_add(metals_init(),Gal[p].MetalsBulgeMassRings[jj],fract_metals);
 /*#ifdef INDIVIDUAL_ELEMENTS
 	      for(kk=0;kk<NUM_ELEMENTS;kk++)
 		Gal[p].BulgeMassRings_elements[jj][kk] = Gal[p].BulgeMassRings_elements[jj][kk]*fract_metals;
 #endif*/
 	    }

 	  ring_sum_minus_tot=-Gal[p].BulgeMass;
 	  for (jj=0; jj<RNUM; jj++)
 	    ring_sum_minus_tot+=Gal[p].BulgeMassRings[jj];
 	  if((ring_sum_minus_tot < -1e-5 && ring_sum_minus_tot < -1e-5*Gal[p].BulgeMass) ||
 	      (ring_sum_minus_tot >  1e-5 && ring_sum_minus_tot >  1e-5*Gal[p].BulgeMass))
 	    {
 	      printf("ring_sum_minus_tot=%0.10e 1e-5*ColdGas=%0.10e\n", ring_sum_minus_tot, 1e-5*Gal[p].BulgeMass);
 	      terminate("");
 	    }
 	}
     }



    //galaxy q Bulge mass
   if(Gal[q].BulgeMass>0.)
     {
       sum_rings=0.;
       sum_rings_metals=0;
       for (jj=0; jj<RNUM; jj++)
 	{
 	  sum_rings+=Gal[q].BulgeMassRings[jj];
 	  sum_rings_metals+= metals_total(Gal[q].MetalsBulgeMassRings[jj]);
      	}

       if(sum_rings>0.)
      	fract=Gal[q].BulgeMass/sum_rings;
       else
      	fract=0.;
       /*if(fract>0. && Gal[q].BulgeMass>1e-6 && sum_rings>1e-6   && (fract>1.3 || fract<0.7))
       	{
       	  printf("fract=%0.5f %s %d\n",fract,call_function, call_line);
       	  terminate("Total and Ring Sum for Metals too different BulgeMass");
       	}*/

       if(sum_rings_metals>0.)
 	fract_metals=metals_total(Gal[q].MetalsBulgeMass)/sum_rings_metals;
       else
 	fract_metals=0.;
     /*  if(fract_metals>0. && metals_total(Gal[q].MetalsBulgeMass)>1e-5 && sum_rings_metals>1e-5 && (fract_metals>1.3 || fract_metals<0.7))
 	{
 	  printf("ColdGasToal=%e ColdGasRings=%e\n MetalsTotal=%e MetalsSum=%e fract=%0.5f %s %d\n",
 	  	 		 Gal[q].BulgeMass, sum_rings, metals_total(Gal[q].MetalsBulgeMass), sum_rings_metals, fract_metals,call_function, call_line);
 	  terminate("Total and Ring Sum for Metals too different BulgeMass");
 	}*/


       if(sum_rings/Gal[q].BulgeMass !=1.)
 	{
 	  for (jj=0;jj<RNUM;jj++)
 	    {
 	      Gal[q].BulgeMassRings[jj] *= fract;
 	      Gal[q].MetalsBulgeMassRings[jj] = metals_add(metals_init(),Gal[q].MetalsBulgeMassRings[jj],fract_metals);
 /*#ifdef INDIVIDUAL_ELEMENTS
 	      for(kk=0;kk<NUM_ELEMENTS;kk++)
 		Gal[q].BulgeMassRings_elements[jj][kk] = Gal[q].BulgeMassRings_elements[jj][kk]*fract_metals;
 #endif*/
 	    }

 	  ring_sum_minus_tot=-Gal[q].BulgeMass;
 	  for (jj=0; jj<RNUM; jj++)
 	    ring_sum_minus_tot+=Gal[q].BulgeMassRings[jj];
 	  if((ring_sum_minus_tot < -1e-5 && ring_sum_minus_tot < -1e-5*Gal[q].BulgeMass) ||
 	      (ring_sum_minus_tot >  1e-5 && ring_sum_minus_tot >  1e-5*Gal[q].BulgeMass))
 	    {
 	      printf("ring_sum_minus_tot=%0.10e 1e-5*ColdGas=%0.10e\n", ring_sum_minus_tot, 1e-5*Gal[q].BulgeMass);
 	      terminate("");
 	    }
 	}
     }
#endif  //RINGS_IN_BULGES

#endif //CHECK_SFH_AND_RINGS






} //end transfer_rings
#endif //H2_AND_RINGS



void mass_checks(int igal, char call_function[], int call_line) {

  /* Some sanity checks on the masses of different components. 
   * If due to rounding error, then apply a correction;
   * otherwise print error message and exit
   */
	//ROB: Should probably make some of these for the elements

  int i;

#ifndef MASS_CHECKS
  return;
#endif

#ifdef STAR_FORMATION_HISTORY
#ifndef DETAILED_METALS_AND_MASS_RETURN
  double sfh_sum=0.;
#endif
#endif

#ifdef H2_AND_RINGS
  double ring_sum_minus_tot=0.;
#endif

#ifdef INDIVIDUAL_ELEMENTS
  float elements_total;
  int kk;
#endif

  //check if the gas mass is less than 0
  if(Gal[igal].ColdGas < 0.0) {
    if (Gal[igal].ColdGas > -1e-6)
      Gal[igal].ColdGas = 0.;
    else {
      printf("\n*** Mass check error, called from: %s, line: %d, ColdGas < 0. ***\n",call_function, call_line);
      printf("                Gal[%d].ColdGas = %g\n",igal,Gal[igal].ColdGas);
      terminate("");
    }
  }

  //check if the mass in metals is less than 0
  if(metals_total(Gal[igal].MetalsColdGas) < 0.0) {
    if (metals_total(Gal[igal].MetalsColdGas) > -1e-6)
      Gal[igal].MetalsColdGas = metals_init();
    else {
      printf("\n*** Mass check error, called from: %s, line: %d, MetalsColdGas < 0. ***\n",call_function, call_line);
      printf("                Gal[%d].MetalsColdGas = %g\n",igal,metals_total(Gal[igal].MetalsColdGas));
#ifdef INDIVIDUAL_ELEMENTS
      elements_total=0.;
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	elements_total+=Gal[igal].ColdGas_elements[kk];
      printf("ColdGas = %f, Total_metal_ele = %f, Snapnum = %i\n", Gal[igal].ColdGas, (elements_total/1.0e10)*Hubble_h,  Gal[igal].SnapNum);
#endif
      terminate("");
    }
  }

  //check if the mass in metals is greater than the gas mass
  if(metals_total(Gal[igal].MetalsColdGas) > Gal[igal].ColdGas) {
    if (metals_total(Gal[igal].MetalsColdGas) < 1e-6)
      Gal[igal].MetalsColdGas = metals_add(metals_init(),Gal[igal].MetalsColdGas,
					   Gal[igal].ColdGas/metals_total(Gal[igal].MetalsColdGas));
    else {
      printf("\n*** Mass check error, called from: %s, line: %d, MetalsColdGas > ColdGas ***\n",call_function, call_line);
      printf("          Gal[%d].MetalsColdGas = %g\n",igal,metals_total(Gal[igal].MetalsColdGas));
      printf("                Gal[%d].ColdGas = %g\n",igal,Gal[igal].ColdGas);
      terminate("");
    }
  }

  if(Gal[igal].HotGas < 0.0) {
    if (Gal[igal].HotGas > -1e-6)
      Gal[igal].HotGas = 0.;
    else {
      printf("\n*** Mass check error, called from: %s, line: %d, HotGas < 0. ***\n",call_function, call_line);
      printf("                Gal[%d].HotGas = %g\n",igal,Gal[igal].HotGas);
      terminate("");
    }
  }

  if(Gal[igal].HotGas > 0.0 && Gal[igal].HotRadius==0)
    {
      printf("\n*** Mass check error, called from: %s, line: %d, HotGas > 0. HotGasRadius=0.***\n",call_function, call_line);
      printf("      Gal[%d].HotGas = %g\n",igal,Gal[igal].HotGas);
      printf("      Gal[%d].HotGasRadius = %g\n",igal,Gal[igal].HotRadius);
      terminate("");
    }

  if(metals_total(Gal[igal].MetalsHotGas) < 0.0) {
    if (metals_total(Gal[igal].MetalsHotGas) > -1e-6)
      Gal[igal].MetalsHotGas = metals_init();
    else {
      printf("\n*** Mass check error, called from: %s, line: %d, MetalsHotGas < 0. ***\n",call_function, call_line);
      printf("                Gal[%d].MetalsHotGas = %g\n",igal,metals_total(Gal[igal].MetalsHotGas));
#ifdef INDIVIDUAL_ELEMENTS
      elements_total=0.;
      for(kk=0;kk<NUM_ELEMENTS;kk++)
      	elements_total+=Gal[igal].HotGas_elements[kk];
      printf("HotGas = %f, Total_metal_ele = %f, Snapnum = %i\n", Gal[igal].HotGas, (elements_total/1.0e10)*Hubble_h,  Gal[igal].SnapNum);
#endif
      terminate("");
    }
  }

  if(metals_total(Gal[igal].MetalsHotGas) > Gal[igal].HotGas) {
    if (metals_total(Gal[igal].MetalsHotGas) < 1e-5)
      Gal[igal].MetalsHotGas = metals_add(metals_init(),Gal[igal].MetalsHotGas,
					   Gal[igal].HotGas/metals_total(Gal[igal].MetalsHotGas));
   else {
      printf("\n***  Mass check error, called from: %s, line: %d, MetalsHotGas > HotGas ***\n",call_function, call_line);
      printf("          Gal[%d].MetalsHotGas = %g\n",igal,metals_total(Gal[igal].MetalsHotGas));
      printf("                Gal[%d].HotGas = %g\n",igal,Gal[igal].HotGas);
      printf("          Gal[%d].MetalsHotGas = %.11f\n",igal,metals_total(Gal[igal].MetalsHotGas));
      printf("                Gal[%d].HotGas = %.11f\n",igal,Gal[igal].HotGas);
      printf("             Gal[%d].BulgeMass = %g\n",igal,Gal[igal].BulgeMass);
      printf("           Gal[%d].EjectedMass = %g\n",igal,Gal[igal].EjectedMass);
#ifdef EXCESS_MASS
      printf("            Gal[%d].ExcessMass = %g\n",igal,Gal[igal].EjectedMass);
#endif
      printf("                  Snapnum = %i\n",Gal[igal].SnapNum);
      terminate("");
    }
  }

  if(Gal[igal].EjectedMass < 0.0) {
    if (Gal[igal].EjectedMass > -1e-6)
      Gal[igal].EjectedMass = 0.;
    else {
      printf("\n*** Mass check error, called from: %s, line: %d, EjectedMass < 0. ***\n",call_function, call_line);
      printf("                Gal[%d].EjectedMass = %g\n",igal,Gal[igal].EjectedMass);
      terminate("");
    }
  }

  if(metals_total(Gal[igal].MetalsEjectedMass) < 0.0) {
    if (metals_total(Gal[igal].MetalsEjectedMass) > -1e-6)
      Gal[igal].MetalsEjectedMass = metals_init();
    else {
      printf("\n*** Mass check error, called from: %s, line: %d, MetalsEjectedMass < 0. ***\n",call_function, call_line);
      printf("                Gal[%d].MetalsEjectedMass = %g\n",igal,metals_total(Gal[igal].MetalsEjectedMass));
      terminate("");
    }
  }

  if(metals_total(Gal[igal].MetalsEjectedMass) > Gal[igal].EjectedMass) {
    if (metals_total(Gal[igal].MetalsEjectedMass) < 1e-6)
      Gal[igal].MetalsEjectedMass = metals_add(metals_init(),Gal[igal].MetalsEjectedMass,
					   Gal[igal].EjectedMass/metals_total(Gal[igal].MetalsEjectedMass));
    else {
      printf("\n*** Mass check error, called from: %s, line: %d, MetalsEjectedMass > EjectedMass ***\n",call_function, call_line);
      printf("          Gal[%d].MetalsEjectedMass = %g\n",igal,metals_total(Gal[igal].MetalsEjectedMass));
      printf("                Gal[%d].EjectedMass = %g\n",igal,Gal[igal].EjectedMass);
      terminate("");
    }
  }

#ifdef EXCESS_MASS
  if(Gal[igal].ExcessMass < 0.0) {
    if (Gal[igal].ExcessMass > -1e-6)
      Gal[igal].ExcessMass = 0.;
    else {
      printf("\n*** Mass check error, called from: %s, line: %d, ExcessMass < 0. ***\n",call_function, call_line);
      printf("                Gal[%d].ExcessMass = %g\n",igal,Gal[igal].ExcessMass);
      terminate("");
    }
  }

  if(metals_total(Gal[igal].MetalsExcessMass) < 0.0) {
    if (metals_total(Gal[igal].MetalsExcessMass) > -1e-6)
      Gal[igal].MetalsExcessMass = metals_init();
    else {
      printf("\n*** Mass check error, called from: %s, line: %d, MetalsExcessMass < 0. ***\n",call_function, call_line);
      printf("                Gal[%d].MetalsExcessMass = %g\n",igal,metals_total(Gal[igal].MetalsExcessMass));
      terminate("");
    }
  }

  if(metals_total(Gal[igal].MetalsExcessMass) > Gal[igal].ExcessMass) {
    if (metals_total(Gal[igal].MetalsExcessMass) < 1e-6)
      Gal[igal].MetalsExcessMass = metals_add(metals_init(),Gal[igal].MetalsExcessMass,
					   Gal[igal].ExcessMass/metals_total(Gal[igal].MetalsExcessMass));
    else {
      printf("\n*** Mass check error, called from: %s, line: %d, MetalsExcessMass > ExcessMass ***\n",call_function, call_line);
      printf("          Gal[%d].MetalsExcessMass = %g\n",igal,metals_total(Gal[igal].MetalsExcessMass));
      printf("                Gal[%d].ExcessMass = %g\n",igal,Gal[igal].ExcessMass);
      terminate("");
    }
  }
#endif

  if(Gal[igal].DiskMass < 0.0) {
    if (Gal[igal].DiskMass > -1e-6)
      Gal[igal].DiskMass = 0.;
    else {
      printf("\n*** Mass check error, called from: %s, line: %d, DiskMass < 0. ***\n",call_function, call_line);
      printf("                Gal[%d].DiskMass = %g\n",igal,Gal[igal].DiskMass);
      terminate("");
    }
  }

  if(Gal[igal].BulgeMass < 0.0) {
    if (Gal[igal].BulgeMass > -1e-6)
      Gal[igal].BulgeMass = 0.;
    else {
      printf("\n*** Mass check error, called from: %s, line: %d, BulgeMass < 0. ***\n",call_function, call_line);
      printf("                Gal[%d].BulgeMass = %g\n",igal,Gal[igal].BulgeMass);
      terminate("");
    }
  }

  if(Gal[igal].ICM < 0.0) {
    if (Gal[igal].ICM > -1e-6)
      Gal[igal].ICM = 0.;
    else {
      printf("\n*** Mass check error, called from: %s, line: %d, ICM < 0. ***\n",call_function, call_line);
      printf("                Gal[%d].ICM = %g\n",igal,Gal[igal].ICM);
      terminate("");
    }
  }



  //CHECK for NAN
  if(isnan(Gal[igal].HotRadius))
      {
        printf("\n*** Mass check error, called from: %s, line: %d, HotRadius = nan. ***\n",call_function, call_line);
        printf("                Gal[%d].HotRadius = %g\n",igal,Gal[igal].HotRadius);
        terminate("");
      }
   if(isnan(Gal[igal].HotGas))
       {
         printf("\n*** Mass check error, called from: %s, line: %d, HotGas = nan. ***\n",call_function, call_line);
         printf("                Gal[%d].HotGas = %g\n",igal,Gal[igal].HotGas);
         terminate("");
       }

   if(isnan(Gal[igal].CoolingGas))
             {
               printf("\n*** Mass check error, called from: %s, line: %d, CoolingGas = nan. ***\n",call_function, call_line);
               printf("                Gal[%d].CoolingGas = %g\n",igal,Gal[igal].CoolingGas);
               terminate("");
             }

   if(isnan(Gal[igal].ColdGas))
     {
       printf("\n*** Mass check error, called from: %s, line: %d, ColdGas = nan. ***\n",call_function, call_line);
       printf("                Gal[%d].ColdGas = %g\n",igal,Gal[igal].ColdGas);
       terminate("");
     }
   if(isnan(Gal[igal].ColdGasRadius))
      {
        printf("\n*** Mass check error, called from: %s, line: %d, ColdGasRadius = nan. ***\n",call_function, call_line);
        printf("                Gal[%d].ColdGasRadius = %g\n",igal,Gal[igal].ColdGasRadius);
        terminate("");
      }


   /*if(isnan(Gal[igal].))
          {
            printf("\n*** Mass check error, called from: %s, line: %d,  = nan. ***\n",call_function, call_line);
            printf("                Gal[%d]. = %g\n",igal,Gal[igal].);
            terminate("");
          }*/



#ifdef TRACK_BURST
  if(Gal[igal].BurstMass < 0.0) {
    if (Gal[igal].BurstMass > -1e-6)
      Gal[igal].BurstMass = 0.;
    else {
      printf("\n*** Mass check error, called from: %s, line: %d, BurstMass < 0. ***\n",call_function, call_line);
      printf("                Gal[%d].BurstMass = %g\n",igal,Gal[igal].BurstMass);
      terminate("");
    }
  }
#endif

/* If DETAILED_METALS_AND_MASS_RETURN, sfh stores accumulation of 'stars', not 'stars-recycFrac'.
 * Therefore, it's sum doesn't equal DiskMass any more.*/
#ifndef DETAILED_METALS_AND_MASS_RETURN
#ifdef STAR_FORMATION_HISTORY
  sfh_sum=-Gal[igal].DiskMass;
  //sfh_sum=0.;
  for (i=0; i<=Gal[igal].sfh_ibin; i++)
    {
      sfh_sum+=Gal[igal].sfh_DiskMass[i];
      if(Gal[igal].sfh_DiskMass[i]<0.0)
	{
	  if(Gal[igal].sfh_DiskMass[i]> -1e-6)
	    Gal[igal].sfh_DiskMass[i]=0.0;
	  else
	    {
	      char sbuf[1000];
	      sprintf(sbuf, "\n*** Mass check error, called from: %s, line: %d, sfh_DiskMass<0. (sfh_Disk=%e)*** \n",call_function, call_line, Gal[igal].sfh_DiskMass[i]);
	      terminate(sbuf);
	    }
	  }
    }

  if((sfh_sum < -1e-4 && sfh_sum < -1e-4*Gal[igal].DiskMass) ||
     (sfh_sum >  1e-4 && sfh_sum >  1e-4*Gal[igal].DiskMass))
 // if((sfh_sum/Gal[igal].DiskMass > 1.0+1.e-3 || sfh_sum/Gal[igal].DiskMass < 1.0-1.e-3) && Gal[igal].DiskMass>1.e-3 )
    {
      printf("                     sfh_sum = %g\n",sfh_sum);
      printf("                Gal[%d].DiskMass = %g\n",igal,Gal[igal].DiskMass);
      printf("            Gal[%d].sfh_DiskMass = %g\n",igal,sfh_sum+Gal[igal].DiskMass);
      //for (i=0; i<=Gal[igal].sfh_ibin; i++)
	//printf("sfh[%d]=%g\n",i,Gal[igal].sfh_DiskMass[i]);
      char sbuf[1000];
      sprintf(sbuf, "\n*** Mass check error, called from: %s, line: %d, Inconsistent sfh for DiskMass.*** \n",call_function, call_line);
      terminate(sbuf);
    }

  sfh_sum=-Gal[igal].BulgeMass;
  for (i=0; i<=Gal[igal].sfh_ibin; i++) sfh_sum+=Gal[igal].sfh_BulgeMass[i];
  if((sfh_sum < -1e-4 && sfh_sum < -1e-4*Gal[igal].BulgeMass) ||
     (sfh_sum >  1e-4 && sfh_sum >  1e-4*Gal[igal].BulgeMass)) {
    printf("                     sfh_sum = %g\n",sfh_sum);
    printf("                Gal[%d].BulgeMass = %g\n",igal,Gal[igal].BulgeMass);
    printf("            Gal[%d].sfh_BulgeMass = %g\n",igal,sfh_sum+Gal[igal].BulgeMass);
    char sbuf[1000];
    sprintf(sbuf, "\n*** Mass check error, called from: %s, line: %d, Inconsistent sfh for BulgeMass. ***\n",call_function, call_line);
    terminate(sbuf);
  }

  sfh_sum=-Gal[igal].ICM;
  for (i=0; i<=Gal[igal].sfh_ibin; i++) sfh_sum+=Gal[igal].sfh_ICM[i];
  if(sfh_sum < -1e-4 || sfh_sum > 1e-4) {
    printf("                     sfh_sum = %g\n",sfh_sum);
    printf("                Gal[%d].ICM = %g\n",igal,Gal[igal].ICM);
    printf("            Gal[%d].sfh_ICM = %g\n",igal,sfh_sum+Gal[igal].ICM);
    for (i=0; i<=Gal[igal].sfh_ibin; i++) 
      printf("%d %f\n",i,Gal[igal].sfh_ICM[i]);
    char sbuf[1000];
    sprintf(sbuf, "\n*** Mass check error, called from: %s, line: %d, Inconsistent sfh for ICM. ***\n",call_function, call_line);
    terminate(sbuf);
  }

#endif
#endif //DETAILED_ENRICHEMENT


#ifdef H2_AND_RINGS

  ring_sum_minus_tot=-Gal[igal].DiskMass;
  for (i=0; i<RNUM; i++)
    {
      ring_sum_minus_tot+=Gal[igal].DiskMassRings[i];
      if(Gal[igal].DiskMassRings[i]<0.0)
	{
	  if(Gal[igal].DiskMassRings[i]> -1e-6)
	    Gal[igal].DiskMassRings[i]=0.0;
	  else
	    {
	      char sbuf[1000];
	      sprintf(sbuf, "\n*** Mass check error, called from: %s, line: %d, DiskMassRings<0. (DiskMassRings=%e)*** \n",call_function, call_line,Gal[igal].DiskMassRings[i]);
	      terminate(sbuf);
	    }
	}
    }
  if((ring_sum_minus_tot < -1e-6 && ring_sum_minus_tot < -1e-6*Gal[igal].DiskMass) ||
      (ring_sum_minus_tot >  1e-6 && ring_sum_minus_tot >  1e-6*Gal[igal].DiskMass))
  // if((sfh_sum/Gal[igal].DiskMass > 1.0+1.e-3 || sfh_sum/Gal[igal].DiskMass < 1.0-1.e-3) && Gal[igal].DiskMass>1.e-3 )
     {
       printf("            ring_sum_minus_tot = %g\n",ring_sum_minus_tot);
       printf("            Gal[%d].DiskMass = %g\n",igal,Gal[igal].DiskMass);
       printf("            Gal[%d].DiskMassRings = %g\n",igal,ring_sum_minus_tot+Gal[igal].DiskMass);
       //for (i=0; i<=Gal[igal].sfh_ibin; i++)
 	//printf("sfh[%d]=%g\n",i,Gal[igal].sfh_DiskMass[i]);
       char sbuf[1000];
       sprintf(sbuf, "\n*** Mass check error, called from: %s, line: %d, Inconsistent ring_sum for DiskMass.*** \n",call_function, call_line);
       terminate(sbuf);
     }

  ring_sum_minus_tot=-Gal[igal].ColdGas;
  for (i=0; i<RNUM; i++)
    {
      ring_sum_minus_tot+=Gal[igal].ColdGasRings[i];
      if(Gal[igal].ColdGasRings[i]<0.0)
	{
	  if(Gal[igal].ColdGasRings[i]> -1e-6)
	    Gal[igal].ColdGasRings[i]=0.0;
	  else
	    {
	      char sbuf[1000];
	      sprintf(sbuf, "\n*** Mass check error, called from: %s, line: %d, ColdGasRings<0. (ColdGasRings=%e)*** \n",call_function, call_line,Gal[igal].ColdGasRings[i]);
	      terminate(sbuf);
	    }
	}
    }
  if((ring_sum_minus_tot < -1e-6 && ring_sum_minus_tot < -1e-6*Gal[igal].ColdGas) ||
      (ring_sum_minus_tot >  1e-6 && ring_sum_minus_tot >  1e-6*Gal[igal].ColdGas))
    {
      printf("            ring_sum_minus_tot = %g\n",ring_sum_minus_tot);
      printf("            Gal[%d].ColdGas = %g\n",igal,Gal[igal].ColdGas);
      printf("            Gal[%d].ColdGasRings = %g\n",igal,ring_sum_minus_tot+Gal[igal].ColdGas);
      char sbuf[1000];
      sprintf(sbuf, "\n*** Mass check error, called from: %s, line: %d, Inconsistent ring_sum for ColdGas.*** \n",call_function, call_line);
      terminate(sbuf);
    }

  for (i=0; i<RNUM; i++)
    if(isnan(Gal[igal].ColdGasRings[i]))
      {
	char sbuf[1000];
	sprintf(sbuf, "\n*** Mass check error, called from: %s, line: %d, ColdGasRings<0. (ColdGasRings=%e)*** \n",call_function, call_line,Gal[igal].ColdGasRings[i]);
	terminate(sbuf);
      }

#endif

#ifdef BULGESIZE_DEBUG
  if ((Gal[igal].BulgeMass > TINY_MASS && Gal[igal].BulgeSize < TINY_LENGTH) ||
       (Gal[igal].BulgeMass < TINY_MASS && Gal[igal].BulgeSize > TINY_LENGTH)) {
      printf("\n*** Bulge Mass and Bulge Size inconsistent ***\n");
      printf("Gal=%d BulgeMass=%g, BulgeSize=%g\n",igal, Gal[igal].BulgeMass,Gal[igal].BulgeSize);
      terminate("");
  }
#endif

  return;
}


double separation_gal(int p, int q) {

  /* Calculates the separation of galaxies p and q, allowing for wrapping */

  int i;
  double sep1,sep2;

  sep2=0.;
  for (i=0; i<3; i++)
    {
      sep1 =  wrap(Gal[p].Pos[i] - Gal[q].Pos[i], BoxSize);
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
  int nfiles, filenr;
#ifndef OVERWRITE_OUTPUT
  int file;
#endif
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


#ifndef OVERWRITE_OUTPUT
#ifdef SPECIFYFILENR
	  file = ListInputFilrNr[filenr];
#else
	  file=filenr;
#endif

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
	      if(j==NTask) // Mass += fractionRings[jj]/RNUM*Gal[q].HotGas;
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



void re_set_parameters(int snapnum)
{
	 if(snapnum<26)
   	    {
		  SfrEfficiency = 0.0594;
		  AgnEfficiency = 4.04e-5;
		  BlackHoleGrowthRate = 0.00924;
		  FeedbackReheatingEpsilon = 6.72;
		  ReheatPreVelocity = 807.; // Mass += fractionRings[jj]/RNUM*Gal[q].HotGas;
		  ReheatSlope = 0.372;
		  FeedbackEjectionEfficiency = 1.16;
		  EjectPreVelocity = 131.;
		  EjectSlope = 2.42;
		  ReIncorporationFactor = 0.000125;
		  //ReIncorporationFactor = 0.0001;
		  Yield = 0.103;
   	    }
   	  else if(snapnum<31)
   	    {
   	  	SfrEfficiency = 0.0126;
   		  AgnEfficiency = 1.92e-5;
   		  BlackHoleGrowthRate = 0.0165;
   		  FeedbackReheatingEpsilon = 1.84;
   		  ReheatPreVelocity = 443.;
   		  ReheatSlope = 0.465;
   		  FeedbackEjectionEfficiency = 2.33;
   		  EjectPreVelocity = 355.;
   		  EjectSlope = 1.08;
   		  ReIncorporationFactor = 0.000891;
   		  //ReIncorporationFactor = 0.0000000001;
   		  Yield = 0.0138;
   	    }
   	  else if(snapnum<38)
   	    {
   		  SfrEfficiency = 0.0479;
   		  AgnEfficiency = 1.62e-5;
   		  BlackHoleGrowthRate = 0.00451;
   		  FeedbackReheatingEpsilon = 8.34;
   		  ReheatPreVelocity = 392.0;
   		  ReheatSlope = 0.876;
   		  FeedbackEjectionEfficiency = 0.891;
   		  EjectPreVelocity = 418.;
   		  EjectSlope = 0.682;
   		  ReIncorporationFactor = 0.822;
   		  Yield = 0.0192;
   	    }
   	  else
   	    {
   		  SfrEfficiency = 0.054;
   		  AgnEfficiency = 3.25e-5;
   		  BlackHoleGrowthRate = 0.0133;
   		  FeedbackReheatingEpsilon = 4.81;
   		  ReheatPreVelocity = 573.;
   		  ReheatSlope = 0.125;
   		  FeedbackEjectionEfficiency = 0.153;
   		  EjectPreVelocity = 171.;
   		  EjectSlope = 1.71;
   		  ReIncorporationFactor = 0.624;
   		  Yield = 0.0834;
   	    }
}



/*
void re_set_parameters(int snapnum)
{
	 if(snapnum<25)
   	    {
		  SfrEfficiency = 0.058;
		  //AgnEfficiency = 1.5e-5;
		  //BlackHoleGrowthRate = 0.083;
		  //FeedbackReheatingEpsilon = 9.3;
		  //ReheatPreVelocity = 60.;
		  //ReheatSlope = 0.43;
		  //FeedbackEjectionEfficiency = 1.9;
		  //EjectPreVelocity = 20.;
		  //EjectSlope = 1.9;
		  //ReIncorporationFactor = 0.77;
		  //Yield = 0.072;
   	    }
   	  else if(snapnum<30)
   	    {
   		  SfrEfficiency = 0.016;
   		  //AgnEfficiency = 1.7e-5;
   		  //BlackHoleGrowthRate = 0.05;
   		  //FeedbackReheatingEpsilon = 8.1;
   		  //ReheatPreVelocity = 45.;
   		  //ReheatSlope = 1.0;
   		  //FeedbackEjectionEfficiency = 3.4;
   		  //EjectPreVelocity = 19.;
   		  //EjectSlope = 1.5;
   		  //ReIncorporationFactor = 0.77;
   		  //Yield = 0.046;
   	    }
   	  else if(snapnum<38)
   	    {
   		  SfrEfficiency = 0.031;
   		  //AgnEfficiency = 1.0e-5;
   		  //BlackHoleGrowthRate = 0.018;
   		  //FeedbackReheatingEpsilon = 8.4;
   		  //ReheatPreVelocity = 26.0;
   		  //ReheatSlope = 0.65;
   		  //FeedbackEjectionEfficiency = 4.3;
   		  //EjectPreVelocity = 18.;
   		  //EjectSlope = 1.9;
   		  //ReIncorporationFactor = 0.77;
   		  //Yield = 0.073;
   	    }
   	  else
   	    {
   		  SfrEfficiency = 0.019;
   		  //AgnEfficiency = 5.0e-6;
   		  //BlackHoleGrowthRate = 0.074;
   		  //FeedbackReheatingEpsilon = 8.5;
   		  //ReheatPreVelocity = 110.;
   		  //ReheatSlope = 0.39;
   		  //FeedbackEjectionEfficiency = 0.78;
   		  //EjectPreVelocity = 30.;
   		  //EjectSlope = 2.0;
   		  //ReIncorporationFactor = 0.77;
   		  //Yield = 0.072;
   	    }
}*/




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
		for (i=1;i<=n-m;i++) {
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


void my_qsort_r(double *x, double *y, int first, int last)
{
  int pivot,j,i;
  double temp_x,temp_y;

  if(first<last)
    {
      pivot=first;
      i=first;
      j=last;

      while(i<j)
	{
	  while(x[i]<=x[pivot] && i<last)
	    i++;
	  while(x[j]>x[pivot])
	    j--;
	  if(i<j)
	    {
	      temp_x=x[i];
	      temp_y=y[i];
	      x[i]=x[j];
	      y[i]=y[j];
	      x[j]=temp_x;
	      y[j]=temp_y;
	    }
	}

      temp_x=x[pivot];
      temp_y=y[pivot];
      x[pivot]=x[j];
      y[pivot]=y[j];
      x[j]=temp_x;
      y[j]=temp_y;
      my_qsort_r(x,y,first,j-1);
      my_qsort_r(x,y,j+1,last);
    }
}


#ifdef DEBUG
void debug_galaxy_init()
{
  char buf[1000];

  sprintf(buf, "%s/debug.txt", OutputDir);
  printf("opening debug file: %s\n",buf);

#ifdef DEBUG_PRINT
  if(!(FdGalDebug = fopen(buf, "w")))
    {
      char sbuf[1000];
      sprintf(sbuf, "can't open file `%s'\n", buf);
      terminate(sbuf);
    }
  for(ii=0;ii<NDebugProps;ii++)
   fprintf(FdGalDebug,"%s\n", DebugProperties[ii]);
#endif

#ifdef DEBUG_READ_AND_CHECK
  if(!(FdGalDebug = fopen(buf, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "can't open file `%s'\n", buf);
      terminate(sbuf);
    }

  char check_prop_name[100], sbuf[1000];
  for(ii=0;ii<NDebugProps;ii++)
    {
      //read what is in the output file and compare with the names
      //in DebugProperties, adding a '\n' included in the output file
      fgets (check_prop_name, 100, FdGalDebug);
      sprintf(sbuf, "%s\n", DebugProperties[ii]);
      if (strcmp(sbuf,check_prop_name)!=0)
	{
	  printf("DebugProperties[%d]='%s', prop_name in debug file='%s'",ii,sbuf,check_prop_name);
	  terminate("");
	}
    }
#endif


}

/*char* deblank(char* input)
{
    int i,j;
    char *output=input;

    for (i = 0, j = 0; i<strlen(input); i++,j++)
      {
	if (input[i]!=' ')
	    output[j]=input[i];
	else
	  j--;

      }

    output[j]=0;
    return output;
}*/

void debug_galaxy(int p, char call_function[], int call_line)
{
#ifdef DEBUG_PRINT
  fprintf(FdGalDebug, "%d %d %d\n", Gal[p].HaloNr, p, Gal[p].SnapNum);
  fprintf(FdGalDebug, "%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f ",
	  Gal[p].HotGas, Gal[p].ColdGas, Gal[p].EjectedMass, Gal[p].DiskMass, Gal[p].BulgeMass,
	  Gal[p].DiskRadius, Gal[p].BulgeSize,
	  Gal[p].MetalsHotGas, Gal[p].MetalsColdGas,Gal[p].MetalsDiskMass, Gal[p].MetalsBulgeMass);

#ifdef H2_AND_RINGS
  int jj;
  for(jj=0;jj<RNUM;jj++)
    fprintf(FdGalDebug, "%0.20f ", Gal[p].ColdGasRings[jj]);
  fprintf(FdGalDebug, "\n");
#endif

  fprintf(FdGalDebug, "\n");
#endif

#ifdef DEBUG_PRINT_TO_CONSOLE
  printf("%s Hnr=%d gal=%d snap=%d\n",call_function, Gal[p].HaloNr, p, Gal[p].SnapNum);
  printf(" Hot=%0.3e Cold=%0.3e Eject=%0.3e disk=%0.3e bulge=%0.3e  ColdGasRadius=%0.3e DiskRadius=%0.3e BulgeSize=%0.3e\n",
	 Gal[p].HotGas*1.e10, Gal[p].ColdGas*1.e10, Gal[p].EjectedMass*1.e10,
	 Gal[p].DiskMass*1.e10, Gal[p].BulgeMass*1.e10, Gal[p].ColdGasRadius, Gal[p].DiskRadius, Gal[p].BulgeSize);

  printf(" HotMetals=%0.3e ColdMetals=%0.3e diskMetals=%0.3e bulgeMetals=%0.3e StellarMetallicity=%0.2f\n",
	 Gal[p].MetalsHotGas*1.e10, Gal[p].MetalsColdGas*1.e10,Gal[p].MetalsDiskMass*1.e10, Gal[p].MetalsBulgeMass*1.e10,
	 log10(Gal[p].MetalsDiskMass+Gal[p].MetalsBulgeMass)/((Gal[p].DiskMass+Gal[p].BulgeMass)*0.02));
#endif

#ifdef DEBUG_READ_AND_CHECK
  double CheckProps[NDebugProps];
  int halonr, galp, snapnum, ii;

  fscanf(FdGalDebug, "%d %d %d\n", &halonr, &galp, &snapnum);
  //printf("%d %d %d\n", halonr, galp, snapnum);
  if (halonr!=Gal[p].HaloNr)
     printf("error debug halonr=%d Gal[p].HaloNr=%d\n", halonr, Gal[p].HaloNr);
  if (galp!=p)
      printf("error debug galp=%d p=%d\n", galp, p);
  if (snapnum!=Gal[p].SnapNum)
      printf("error debug snapnum=%d Gal[p].SnapNum=%d\n", snapnum, Gal[p].SnapNum);
  if ((halonr!=Gal[p].HaloNr) || (galp!=p) || (snapnum!=Gal[p].SnapNum))
    {
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  for (ii=0;ii<NDebugProps;ii++)
    fscanf(FdGalDebug, "%lf", &CheckProps[ii]);

  char buf1[100], buf2[100];

  sprintf(buf1, "%0.20f", CheckProps[0]);
  sprintf(buf2, "%0.20f", Gal[p].HotGas);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("HotGas='%s' HotGas in file='%s'\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[1]);
  sprintf(buf2, "%0.20f", Gal[p].ColdGas);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("ColdGas=%s ColdGas in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[2]);
  sprintf(buf2, "%0.20f", Gal[p].EjectedMass);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("EjectedMass=%s EjectedMass in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[3]);
  sprintf(buf2, "%0.20f", Gal[p].DiskMass);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("DiskMass=%s DiskMass in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[4]);
  sprintf(buf2, "%0.20f", Gal[p].BulgeMass);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("BulgeMass=%s BulgeMass in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[5]);
  sprintf(buf2, "%0.20f", Gal[p].DiskRadius);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("DiskRadius=%s DiskRadius in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[6]);
  sprintf(buf2, "%0.20f", Gal[p].BulgeSize);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("BulgeSize=%s BulgeSize in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[7]);
  sprintf(buf2, "%0.20f", Gal[p].MetalsHotGas);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("MetalsHotGas=%s MetalsHotGas in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[8]);
  sprintf(buf2, "%0.20f", Gal[p].MetalsColdGas);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("MetalsColdGas=%s MetalsColdGas in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[9]);
  sprintf(buf2, "%0.20f", Gal[p].MetalsDiskMass);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("MetalsDiskMass=%s MetalsDiskMass in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[10]);
  sprintf(buf2, "%0.20f", Gal[p].MetalsBulgeMass);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("MetalsBulgeMass=%s MetalsBulgeMass in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

#ifdef H2_AND_RINGS
  sprintf(buf1, "%0.20f", CheckProps[11]);
  sprintf(buf2, "%0.20f", Gal[p].ColdGasRings[0]);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("ColdGasRings0=%s ColdGasRings0 in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[12]);
  sprintf(buf2, "%0.20f", Gal[p].ColdGasRings[1]);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("ColdGasRings1=%s ColdGasRings1 in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[13]);
  sprintf(buf2, "%0.20f", Gal[p].ColdGasRings[2]);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("ColdGasRings2=%s ColdGasRings2 in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[14]);
  sprintf(buf2, "%0.20f", Gal[p].ColdGasRings[3]);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("ColdGasRings3=%s ColdGasRings3 in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[15]);
  sprintf(buf2, "%0.20f", Gal[p].ColdGasRings[4]);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("ColdGasRings4=%s ColdGasRings4 in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[16]);
  sprintf(buf2, "%0.20f", Gal[p].ColdGasRings[5]);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("ColdGasRings5=%s ColdGasRings5 in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[17]);
  sprintf(buf2, "%0.20f", Gal[p].ColdGasRings[6]);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("ColdGasRings6=%s ColdGasRings6 in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[18]);
  sprintf(buf2, "%0.20f", Gal[p].ColdGasRings[7]);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("ColdGasRings7=%s ColdGasRings7 in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[19]);
  sprintf(buf2, "%0.20f", Gal[p].ColdGasRings[8]);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("ColdGasRings8=%s ColdGasRings8 in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[20]);
  sprintf(buf2, "%0.20f", Gal[p].ColdGasRings[9]);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("ColdGasRings9=%s ColdGasRings9 in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[21]);
  sprintf(buf2, "%0.20f", Gal[p].ColdGasRings[10]);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("ColdGasRings10=%s ColdGasRings10 in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }

  sprintf(buf1, "%0.20f", CheckProps[22]);
  sprintf(buf2, "%0.20f", Gal[p].ColdGasRings[11]);
  if(strcmp(buf1,buf2)!=0)
    {
      printf("ColdGasRings11=%s ColdGasRings11 in file=%s\n",buf2,buf1);
      printf("%d %d %d\n", halonr, galp, snapnum);
      printf("called from %s, line %d\n", call_function, call_line);
      terminate("");
    }
#endif
#endif
  if(isnan(Gal[p].ColdGasRadius))
    terminate("ColdGasRadius=nan");

}
#endif

#ifdef H2_AND_RINGS
void print_rings(char string[], int p)
{
  int jj;

  //if(Gal[p].Vvir>200. && Gal[p].Vvir<235. && Gal[p].Type==0 && Gal[p].SnapNum>40)
  if(HaloIDs[Gal[p].HaloNr].FirstHaloInFOFgroup > 40000000012495 && HaloIDs[Gal[p].HaloNr].FirstHaloInFOFgroup < 40000000012757 && Gal[p].Type==0)
    {
      //printf("haloid=%lld mainleaf=%lld\n",HaloIDs[Gal[p].HaloNr].HaloID, HaloIDs[Gal[p].HaloNr].LastProgenitor);
      printf("%s SNAP=%d H2  starsRings DiskMassRings\n", string, Gal[p].SnapNum);
      for(jj=0;jj<RNUM;jj++)
	printf("[%d] H=%0.2e H2=%0.2e DiskMassRings=%0.2e\n", jj, Gal[p].ColdGasRings[jj]*1.e10,
	       Gal[p].ColdGasRings[jj]*Gal[p].H2fractionRings[jj]*1.e10, Gal[p].DiskMassRings[jj]*1.e10);
    }
}

void print_check_rings_and_total(char string[], int p)
{
  double sum_rings=0.;
  int jj;

  for(jj=0;jj<RNUM;jj++)
    sum_rings+=Gal[p].ColdGasRings[jj];

  printf("%s total=%0.7e ring_sum=%0.7e\n",string, Gal[p].ColdGas*1.e10,sum_rings*1.e10);
}
#endif //H2_AND_RINGS


int string_length(char *s)
{
  int c = 0;

  while(*(s+c))
    c++;

  return c;
}


