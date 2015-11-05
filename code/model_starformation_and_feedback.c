#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @file recipe_starformation_and_feedback.c
 *  @brief recipe_starformation_and_feedback.c computes the amount of stars
 *         formed from the cold gas, the amount of gas reheated from cold to hot
 *         and the amount of gas ejected from hot to external.
 *
 * The routine is divided in two parts, star formation and SN feedback, with a
 * number of different implementations controlled by input parameters.
 *
 *
 *  0 -\f$M_{\rm{crit}}=3.8\times 10^9
 *     \left(\frac{V_{\rm{max}}}{200\,\rm{km s}^{-1}}\right)
 *     \left(\frac{r_{\rm{disk}}}{10\,\rm{kpc}}\right)M_{\odot}\f$
 *     (Eq. 16 Guo2010) (StarFormationModel = 0), \n
 *        - same as 1 but using \f$V_{\rm{max}}\f$ or \f$V_{\rm{max,infall}}\f$
 *          instead of \f$V_{\rm{vir}}\f$ and allowing SF in satellites. *

 *
 * There are 2 options for the <B>SN Feedback Recipe</B>:
 *
 * 0 - \f$\epsilon_{\rm{disk}}=\epsilon
 *      \biggl[0.5+\left(\frac{V_{\rm{max}}}{70km/s}\right)^{-\beta_1}\biggr]\f$,
 *     \f$\epsilon_{\rm{halo}}=\eta
 *      \biggl[0.5+\left(\frac{V_{\rm{max}}}{70km/s}\right)^{-\beta_2}\biggr]\f$
 *     (Eqs. 19 & 21 Guo2010)(FeedbackEjectionModel = 2)
 *     same as FeedbackEjectionModel = 1 * Vmax dependence.
 *
 * Also, Guo2010 alowed for type 1 satellite to have gas cycles and receive
 * gas from their own satellites when these are outside Rvir of the type 0.
 * */


/** @brief Main recipe, calculates the fraction of cold gas turned into stars due
  *        to star formation; the fraction of mass instantaneously recycled and
  *        returned to the cold gas; the fraction of gas reheated from cold to hot,
  *        ejected from hot to external and returned from ejected to hot due to
  *        SN feedback.   */
void starformation(int p, int centralgal, double time, double dt, int nstep)
{
	/*! Variables: reff-Rdisk, tdyn=Rdisk/Vmax, strdot=Mstar_dot, stars=strdot*dt*/
  double tdyn, strdot=0., stars, cold_crit, metallicitySF;
  
  if(Gal[p].Type == 0)
    {
      tdyn = Gal[p].GasDiskRadius / Gal[p].Vmax;
      cold_crit = SfrColdCrit * Gal[p].Vmax/200. * Gal[p].GasDiskRadius*100.;
    }
  else
    {
      tdyn = Gal[p].GasDiskRadius / Gal[p].InfallVmax;
      cold_crit = SfrColdCrit * Gal[p].InfallVmax/200. * Gal[p].GasDiskRadius*100.;
    }

  //standard star formation law (Croton2006, Delucia2007, Guo2010)
  if(StarFormationModel == 0)
    {
      if(Gal[p].ColdGas > cold_crit)
	strdot = SfrEfficiency * (Gal[p].ColdGas - cold_crit) / tdyn;
      else
	strdot = 0.0;
    }

  /*if(StarFormationModel == 1)
  {
  	strdot = ALTERNATIVE STAR FORMATION LAW

  }*/

  /* Note that Units of dynamical time are Mpc/Km/s - no conversion on dt needed */
  stars = strdot * dt;
  if(stars < 0.0)
    terminate("***error stars<0.0***\n");

  /*  update the star formation rate */
  /*Sfr=stars/(dt*steps)=strdot*dt/(dt*steps)=strdot/steps -> average over the STEPS*/
  Gal[p].Sfr += stars / (dt * STEPS);

  mass_checks("recipe_starform #1",p);
  mass_checks("recipe_starform #1.1",centralgal);

  /* update for star formation
   * updates Mcold, StellarMass, MetalsMcold and MetalsStellarMass
   * in Guo2010 case updates the stellar spin -> hardwired, not an option */

  /* Store the value of the metallicity of the cold phase when SF occurs */
  if (Gal[p].ColdGas > 0.)
  	metallicitySF= metals_total(Gal[p].MetalsColdGas)/Gal[p].ColdGas;
  else
    metallicitySF=0.;
 

  if (stars > 0.)
  	update_stars_due_to_reheat(p, centralgal, &stars);

  mass_checks("recipe_starform #2",p);
  mass_checks("recipe_starform #2.1",centralgal);


  // update_from_star_formation can only be called
  // after SD_feeedback recipe since stars need to be re_set once the reheated mass is known
  // (star formation and feedback share the same fraction of cold gas)
  if (stars > 0.)
    update_from_star_formation(p, stars, false, nstep); // false indicates not a burst

  update_massweightage(p, stars, time);

  if (stars > 0.)
    SN_feedback(p, centralgal, stars, "ColdGas");


#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  /*  Update the luminosities due to the stars formed */
  if (stars > 0.0)
    add_to_luminosities(p, stars, time, dt, metallicitySF);
#endif //NDEF POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

  if(Gal[p].DiskMass > 0.0)
    check_disk_instability(p);

  if (DiskRadiusModel== 0)
    get_stellar_disk_radius(p);

}






void update_stars_due_to_reheat(int p, int centralgal, double *stars)
{
  double MergeCentralVvir=0.;
  double fac;
  double CentralVvir=0.;
  double reheated_mass=0., ejected_mass=0.;
  /* SN FEEDBACK RECIPES */

  /* In Guo2010 type 1s can eject, reincorporate gas and get gas from their
   * own satellites (is not sent to the type 0 galaxy as in Delucia2007),
   * for gas flow computations:
   * If satellite is inside Rvir of main halo, Vvir of main halo used
   * If it is outside, the Vvir of its central subhalo is used. */

  //REHEAT
  CentralVvir = Gal[centralgal].Vvir; // main halo Vvir
  MergeCentralVvir = Gal[Gal[p].CentralGal].Vvir; //central subhalo Vvir

  // Feedback depends on the circular velocity of the host halo
  // Guo2010 - eq 18 & 19
  if(FeedbackReheatingModel == 0)
    {
      if (Gal[Gal[p].CentralGal].Type == 0)
	reheated_mass = FeedbackReheatingEpsilon * *stars *
	(.5+1./pow(Gal[Gal[p].CentralGal].Vmax/ReheatPreVelocity,ReheatSlope));
      else
	reheated_mass = FeedbackReheatingEpsilon * *stars *
	(.5+1./pow(Gal[Gal[p].CentralGal].InfallVmax/ReheatPreVelocity,ReheatSlope));

      if (reheated_mass * Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir > *stars * (EtaSNcode * EnergySNcode))
	reheated_mass = *stars * (EtaSNcode * EnergySNcode) / (Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir);
    }
  /*else if(FeedbackReheatingModel == 1)
    {
      reheated_mass =  ALTERNATIVE Reheating LAW ;
    }*/

  if((*stars + reheated_mass) > Gal[p].ColdGas)
    {
      fac = Gal[p].ColdGas / (*stars + reheated_mass);
      *stars *= fac;
      reheated_mass *= fac;
    }

}



/** @brief Updates the different components due to star formation: mass
  *        and metals in stars and cold gas and stellar spin. */
//void update_from_star_formation(int p, double time, double stars, double metallicity)
void update_from_star_formation(int p, double stars, bool flag_burst, int nstep)
{
  int i;
  double fraction;
  double stars_nett=0.;

  if(Gal[p].ColdGas <= 0. || stars <= 0.) {
    printf("update_from_star_formation: Gal[p].ColdGas <= 0. || stars <= 0.\n");
    exit(0);
  }

  stars_nett=(1 - RecycleFraction) * stars;
  /* Update the Stellar Spin when forming stars */
  if (Gal[p].DiskMass+stars_nett > 1.e-8)
    for (i = 0; i < 3; i++)
      Gal[p].StellarSpin[i]=((Gal[p].StellarSpin[i])*(Gal[p].DiskMass)+stars_nett*Gal[p].GasSpin[i])/(Gal[p].DiskMass+stars_nett);

    /*  Update Gas and Metals from star formation */
  mass_checks("update_from_star_formation #0",p);

  fraction=stars_nett/Gal[p].ColdGas;

#ifdef STAR_FORMATION_HISTORY
  Gal[p].sfh_DiskMass[Gal[p].sfh_ibin]+=stars_nett; //ROB: Add amount of stars formed to SFH history of the disk. (NOTE: ALL SF OCCURS IN THE DISK. sfh_BulgeMass only increases when stars are transferred to the bulge before they explode)
  Gal[p].sfh_MetalsDiskMass[Gal[p].sfh_ibin] = metals_add(Gal[p].sfh_MetalsDiskMass[Gal[p].sfh_ibin],Gal[p].MetalsColdGas,fraction);
#ifdef TRACK_BURST
  if (flag_burst) Gal[p].sfh_BurstMass[Gal[p].sfh_ibin]+=stars_nett;
#endif
#endif //STAR_FORMATION_HISTORY

  Gal[p].MetalsDiskMass=metals_add(Gal[p].MetalsDiskMass,Gal[p].MetalsColdGas,fraction);
  Gal[p].MetalsColdGas=metals_add(Gal[p].MetalsColdGas,Gal[p].MetalsColdGas,-fraction);

  Gal[p].DiskMass += stars_nett;
  Gal[p].ColdGas -= stars_nett;
#ifdef TRACK_BURST
  if (flag_burst) Gal[p].BurstMass+=stars_nett;
#endif

  mass_checks("update_from_star_formation #1",p);

  /* Formation of new metals - instantaneous recycling approximation - only SNII
   * Also recompute the metallicity of the cold phase.*/
  Gal[p].MetalsColdGas += Yield * stars;

  if (DiskRadiusModel == 0)
    get_stellar_disk_radius(p);

}

/* there are two modes for supernova feedback corresponding to when the mass returning
 * by dying stars is returned to the cold gas - reheat and ejection; and when the mass
 * is returned to the hot gas - onle ejection.*/
void SN_feedback(int p, int centralgal, double stars, char feedback_location[])
{
  double CentralVvir, MergeCentralVvir=0., EjectVmax, EjectVvir, SN_Energy, Reheat_Energy, fac;
  double reheated_mass=0., ejected_mass=0.;
  /* SN FEEDBACK MODEL */

  /* In Guo2010 type 1s can eject, reincorporate gas and get gas from their
   * own satellites (is not sent to the type 0 galaxy as in Delucia2007),
   * for gas flow computations:
   * If satellite is inside Rvir of main halo, Vvir of main halo used
   * If it is outside, the Vvir of its central subhalo is used. */

  if (strcmp(feedback_location,"HotGas")==0)
    reheated_mass = 0.;
  else
    if (strcmp(feedback_location,"ColdGas")==0)
	{
	CentralVvir = Gal[centralgal].Vvir; // main halo Vvir
	MergeCentralVvir = Gal[Gal[p].CentralGal].Vvir; //central subhalo Vvir

	mass_checks("recipe_starform #0",p);
	mass_checks("recipe_starform #0.1",centralgal);

	// Feedback depends on the circular velocity of the host halo
	// Guo2010 - eq 18 & 19
	if(FeedbackReheatingModel == 0)
	  {
	    if (Gal[Gal[p].CentralGal].Type == 0)
	      reheated_mass = FeedbackReheatingEpsilon * stars *
	      (.5+1./pow(Gal[Gal[p].CentralGal].Vmax/ReheatPreVelocity,ReheatSlope));
	    else
	      reheated_mass = FeedbackReheatingEpsilon * stars *
	      (.5+1./pow(Gal[Gal[p].CentralGal].InfallVmax/ReheatPreVelocity,ReheatSlope));

	    if (reheated_mass * Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir > stars * (EtaSNcode * EnergySNcode))
	      reheated_mass = stars * (EtaSNcode * EnergySNcode) / (Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir);
	  }
	/*else if(FeedbackReheatingModel == 1)
	  {
	    reheated_mass =  ALTERNATIVE Reheating LAW ;
	  }*/

	if(reheated_mass > Gal[p].ColdGas)
	  reheated_mass = Gal[p].ColdGas;

	}// end if feedback_location


  /* Determine ejection (for FeedbackEjectionModel 2 we have the dependence on Vmax)
   * Guo2010 - eq 22
   * Note that satellites can now retain gas and have their own gas cycle*/

  if (Gal[Gal[p].CentralGal].Type == 0)
    {
      EjectVmax=Gal[centralgal].Vmax;
      EjectVvir=Gal[centralgal].Vvir;// main halo Vvir
    }
  else
    {
      EjectVmax=Gal[Gal[p].CentralGal].InfallVmax;
      EjectVvir=Gal[Gal[p].CentralGal].Vvir; //central subhalo Vvir
    }

  if(FeedbackEjectionModel == 0)
    {
      ejected_mass = (FeedbackEjectionEfficiency* (EtaSNcode * EnergySNcode) * stars *
	      min(1./FeedbackEjectionEfficiency, .5+1/pow(EjectVmax/EjectPreVelocity,EjectSlope)) -
	      reheated_mass*EjectVvir*EjectVvir) /(EjectVvir*EjectVvir);
    }
  else if(FeedbackEjectionModel == 1)//the ejected material is assumed to have V_SN
    {
      SN_Energy = .5 * stars * (EtaSNcode * EnergySNcode);
      Reheat_Energy = .5 * reheated_mass * EjectVvir * EjectVvir;

      ejected_mass = (SN_Energy - Reheat_Energy)/(0.5 * FeedbackEjectionEfficiency*(EtaSNcode * EnergySNcode));

      //if VSN^2<Vvir^2 nothing is ejected
      if(FeedbackEjectionEfficiency*(EtaSNcode * EnergySNcode)<EjectVvir*EjectVvir)
	ejected_mass =0.0;
    }

  // Finished calculating mass exchanges, so just check that none are negative
  if (reheated_mass < 0.0) reheated_mass = 0.0;
  if (ejected_mass < 0.0) ejected_mass = 0.0;

  /* Update For Feedback */
  /* update cold, hot, ejected gas fractions and respective metallicities
   * there are a number of changes introduced by Guo2010 concerning where
   * the gas ends up */

  //ejected_mass = 0.01*Gal[centralgal].HotGas;
  if (reheated_mass + ejected_mass > 0.)
    update_from_feedback(p, centralgal, reheated_mass, ejected_mass);
}


/** @brief Updates cold, hot and external gas components due to SN
 *         reheating and ejection. */
void update_from_feedback(int p, int centralgal, double reheated_mass, double ejected_mass)
{
  double dis=0.;
  double massremain;
  double fraction;
  int merger_centre;

  //mass_checks("update_from_feedback #1",p);

  if(Gal[p].ColdGas > 0.)
  {
      //REHEAT
      // if galaxy is a type 1 or a type 2 orbiting a type 1 with hot gas being striped,
      //some of the reheated and ejected masses goes to the type 0 and some stays in the type 1

      if(Gal[p].Type ==0)
	{
	  transfer_gas(p,"Hot",p,"Cold",((float)reheated_mass)/Gal[p].ColdGas,"update_from_feedback", __LINE__);
	}
      else
	{
	  if(Gal[p].Type ==1)
	    merger_centre=centralgal;
	  else if(Gal[p].Type ==2)
	    merger_centre=Gal[p].CentralGal;

	  dis=separation_gal(centralgal,Gal[p].CentralGal)/(1+ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);

	  //compute share of reheated mass
	  if (dis<Gal[centralgal].Rvir && Gal[Gal[p].CentralGal].Type == 1)
	    {
	      //mass that remains on type1 (the rest goes to type 0) for reheat - massremain, for eject - ejected mass
	      massremain=reheated_mass*Gal[p].HotRadius/Gal[p].Rvir;
	      ejected_mass = ejected_mass*Gal[p].HotRadius/Gal[p].Rvir;

	      if (massremain > reheated_mass)
		massremain = reheated_mass;
	    }
	  else
	    massremain=reheated_mass;

	  //needed due to precision issues, since we first remove massremain and then (reheated_mass-massremain)
	  //from the satellite into the type 0 and type 1 the fraction might not add up on the second call
	  //since Gal[p].ColdGas is a float and reheated_mass & massremain are doubles
	  if((massremain + reheated_mass)>Gal[p].ColdGas)
	    massremain=Gal[p].ColdGas-reheated_mass;

	  //transfer massremain
	  transfer_gas(Gal[p].CentralGal,"Hot",p,"Cold",massremain/Gal[p].ColdGas,"update_from_feedback", __LINE__);

	  //transfer reheated_mass-massremain from galaxy to the type 0
	  if (reheated_mass > massremain)
	    if(Gal[p].ColdGas > 0.) //if the reheat to itself, left cold gas below limit do not reheat to central
	      transfer_gas(centralgal,"Hot",p,"Cold",(reheated_mass-massremain)/Gal[p].ColdGas,"update_from_feedback", __LINE__);
	}//types

  }//if(Gal[p].ColdGas > 0.)

  mass_checks("update_from_feedback #2",p);

  //DO EJECTION OF GAS
  if (Gal[Gal[p].CentralGal].HotGas > 0.)
    {
      if (ejected_mass > Gal[Gal[p].CentralGal].HotGas)
	ejected_mass = Gal[Gal[p].CentralGal].HotGas;  //either eject own gas or merger_centre gas for ttype 2's

      fraction=((float)ejected_mass)/Gal[Gal[p].CentralGal].HotGas;

      if (Gal[Gal[p].CentralGal].Type == 1)
	{
	  /* If type 1, or type 2 orbiting type 1 near type 0 */
	  if (FateOfSatellitesGas == 0)
	    transfer_gas(Gal[p].CentralGal,"Ejected",Gal[p].CentralGal,"Hot",fraction,"update_from_feedback", __LINE__);
	  else if (FateOfSatellitesGas == 1)
	    {
	      if (dis < Gal[centralgal].Rvir)
		transfer_gas(centralgal,"Hot",Gal[p].CentralGal,"Hot",fraction,"update_from_feedback", __LINE__);
	      else
		transfer_gas(Gal[p].CentralGal,"Ejected",Gal[p].CentralGal,"Hot",fraction,"update_from_feedback", __LINE__);
	    }
	}
      else // If galaxy type 0 or type 2 merging into type 0
	transfer_gas(centralgal,"Ejected",Gal[p].CentralGal,"Hot",fraction,"update_from_feedback", __LINE__);

    }//(Gal[Gal[p].CentralGal].HotGas > 0.)

}

//Age in Mpc/Km/s/h - code units
void update_massweightage(int p, double stars, double time)
{
  int outputbin;
  double age;

  for(outputbin = 0; outputbin < NOUT; outputbin++)
    {
      age = time - NumToTime(ListOutputSnaps[outputbin]);
      Gal[p].MassWeightAge[outputbin] += age * stars * (1. - RecycleFraction);
    }
}

/** @brief Checks for disk stability using the
 *         Mo, Mao & White (1998) criteria */
void check_disk_instability(int p)
{

  double Mcrit, fraction, stars, diskmass;

/** @brief Calculates the stability of the stellar disk as discussed
 *         in Mo, Mao & White (1998). For unstable stars, the required
 *         amount is transfered to the bulge to make the disk stable again.
 *         Mass, metals and luminosities updated. After Guo2010 the bulge
 *         size is followed and needs to be updated.
 *         Eq 34 & 35 in Guo2010 are used. */


  if(DiskInstabilityModel == 0)
    {
      /* check stellar disk -> eq 34 Guo2010*/
      if (Gal[p].Type != 0)
	Mcrit = Gal[p].InfallVmax * Gal[p].InfallVmax * Gal[p].StellarDiskRadius / G;
      else
	Mcrit = Gal[p].Vmax * Gal[p].Vmax * Gal[p].StellarDiskRadius / G;
      diskmass = Gal[p].DiskMass;
      stars = diskmass - Mcrit;
      fraction = stars / diskmass;
    }
  else if (DiskInstabilityModel == 1)
    stars = 0.;

  /* add excess stars to the bulge */
  if(stars > 0.0)
    {
      /* to calculate the bulge size */
      update_bulge_from_disk(p,stars);
      transfer_stars(p,"Bulge",p,"Disk",fraction);

      if(BHGrowthInDiskInstabilityModel == 1)
	if(Gal[p].ColdGas > 0.)
	  {
	    Gal[p].BlackHoleMass += Gal[p].ColdGas*fraction;
	    Gal[p].ColdGas -= Gal[p].ColdGas*fraction;
	  }

#ifndef POST_PROCESS_MAGS
      double Lumdisk;
      int outputbin, j;
#ifdef OUTPUT_REST_MAGS
      for(outputbin = 0; outputbin < NOUT; outputbin++)
	{
	  for(j = 0; j < NMAG; j++)
	    {
	      Lumdisk = Gal[p].Lum[j][outputbin]-Gal[p].LumBulge[j][outputbin];
	      Gal[p].LumBulge[j][outputbin] += fraction * Lumdisk;
	      Lumdisk = Gal[p].YLum[j][outputbin]-Gal[p].YLumBulge[j][outputbin];
	      Gal[p].YLumBulge[j][outputbin] += fraction * Lumdisk;
	    }
	}
#endif
#ifdef COMPUTE_OBS_MAGS
      for(outputbin = 0; outputbin < NOUT; outputbin++)
	{
	  for(j = 0; j < NMAG; j++)
	    {
	      Lumdisk = Gal[p].ObsLum[j][outputbin]-Gal[p].ObsLumBulge[j][outputbin];
	      Gal[p].ObsLumBulge[j][outputbin] += fraction * Lumdisk;
	      Lumdisk = Gal[p].ObsYLum[j][outputbin]-Gal[p].ObsYLumBulge[j][outputbin];
	      Gal[p].ObsYLumBulge[j][outputbin] += fraction * Lumdisk;
#ifdef OUTPUT_MOMAF_INPUTS
	      Lumdisk = Gal[p].dObsLum[j][outputbin]-Gal[p].dObsLumBulge[j][outputbin];
	      Gal[p].dObsLumBulge[j][outputbin] += fraction * Lumdisk;
	      Lumdisk = Gal[p].dObsYLum[j][outputbin]-Gal[p].dObsYLumBulge[j][outputbin];
	      Gal[p].dObsYLumBulge[j][outputbin] += fraction * Lumdisk;
#endif
	    }
	}
#endif
#endif

      if ((Gal[p].BulgeMass > 1e-9 && Gal[p].BulgeSize == 0.0)||
	  (Gal[p].BulgeMass == 0.0 && Gal[p].BulgeSize >1e-9))
	{
	  char sbuf[1000];
	  sprintf(sbuf, "bulgesize wrong in diskinstablility.c \n");
	  terminate(sbuf);
	}
    }// if(stars > 0.0)
  
}


/** @brief Introduced in Guo2010 to track the change in size of bulges
 *         after their growth due to disk instabilities. */

void update_bulge_from_disk(int p, double stars)
{      
  double bulgesize, diskmass, fint, massfrac, orisize;
  int  j;

/** @brief Updates bulge from disk instability -> stars represents the mass
  *        transfered to the bulge, which occupies a size in the bulge equal
  *        to the occupied in the disk. */


  orisize=Gal[p].BulgeSize; //remove, not used
  diskmass=(Gal[p].DiskMass);

  /* alpha_inter=2.0/C=0.5 (alpha larger than in mergers since
   * the existing and newly formed bulges are concentric)*/
  fint=4.0;

  /* update the stellardisk spin due to the angular momentum transfer
   * from disk to bulge changing the specific angular momentum for disk stars.
   * This should be done on the main routine, as this is update bulge.*/
  massfrac=stars/diskmass;
  for (j = 0; j <3 ; j++)
  {
  	if(massfrac==1) //everything transferred to the bulge
  		Gal[p].StellarSpin[j]=0;
  	else
  		Gal[p].StellarSpin[j]=Gal[p].StellarSpin[j]/(1-massfrac);

  }

  /* update disksize done, disk mass is automatically given by total-bulge*/

//GET BULGE SIZE - Eq. 35 in Guo2010
  /* if previous Bulge Mass = 0
     -> bulge size is given directly from newly formed bulge */
  if(Gal[p].BulgeMass <1.e-9) {
    /* size of newly formed bulge, which consists of the stellar mass
     * transfered from the disk. This is calculated using bulge_from_disk
     * which receives Delta_M/DiskMass and returns Rb/Rd. From eq 35 and
     * since DiskMass=2PISigma(Rd)^2 we see that Delta_M/DiskMass=1-(1+Rb/Rd)*exp(-Rb/Rd),
     * so function bulge_from_disk avoids calculating the slow "ln" function */
    bulgesize=bulge_from_disk(stars/diskmass)*Gal[p].StellarDiskRadius/3.;
    Gal[p].BulgeSize=bulgesize;

  }      
  else {
    bulgesize=bulge_from_disk(stars/diskmass)*Gal[p].StellarDiskRadius/3.;
    /* combine the old with newly formed bulge and calculate the
     * bulge size assuming energy conservation as for mergers but
     * using alpha=2. - eq 33 */
    Gal[p].BulgeSize=(Gal[p].BulgeMass+stars)*(Gal[p].BulgeMass+stars)/
      (Gal[p].BulgeMass*Gal[p].BulgeMass/Gal[p].BulgeSize+stars*stars/bulgesize+fint*Gal[p].BulgeMass*stars/(Gal[p].BulgeSize+bulgesize));
  }

  if((Gal[p].BulgeMass + stars > 1.e-9 && Gal[p].BulgeSize == 0.0)
     || (Gal[p].BulgeMass + stars == 0 && Gal[p].BulgeSize > 1.e-9))
  {
  	char sbuf[1000];
  	printf("GasDiskMass=%e GasDiskSize=%e \nStellarDiskMass=%e StellarDiskSize=%e \nBulgeMass=%e Bulgesize=%e\n",
  			Gal[p].ColdGas, Gal[p].GasDiskRadius, Gal[p].DiskMass, Gal[p].StellarDiskRadius,
  			Gal[p].BulgeMass, Gal[p].BulgeSize);

  	printf("TransferSize=%e, OriBulgeSize=%e\n", bulgesize, orisize);

  	sprintf(sbuf, "bulgesize or mass wrong in disk instablility");
      terminate(sbuf);
    }


}


/** @brief Calculates the size of the disk that contains the
 *         mass transfered to the bulge. */
double bulge_from_disk(double frac)
{
  double x1,x2,x0,value;
/** @brief Calculates the size of the disk that contains the
 *         mass transfered to the bulge. The bulge is assumed
 *         to form with the same size. avoid doing "ln" from eq 35*/
  x1=0.0;
  x2=1.;
  while ((func_size(x2,frac) * func_size(x1,frac))>0) {
    x1=x2;
    x2=x2*2;
  }
  x0=x1+(x2-x1)/2.;
  value=func_size(x0,frac);
  if (value < 0) 
    value = -value;

  while(value>0.00001) {
    if(func_size(x0,frac)*func_size(x2,frac)>0)
      x2=x0;
    else
      x1=x0;
    x0=x1+(x2-x1)/2.;
    value=func_size(x0,frac);
    if (value < 0) 
      value = -value;
  }
    
  return x0;
}


double func_size(double x, double a)
{
  return  exp(-x)*(1+x)-(1-a);
}  





