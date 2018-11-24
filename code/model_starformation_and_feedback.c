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
 *        - same as 0 but using \f$V_{\rm{max}}\f$ or \f$V_{\rm{max,infall}}\f$
 *          instead of \f$V_{\rm{vir}}\f$ and allowing SF in satellites. *

 *
 * There are 2 options for the <B>SN Feedback Recipe</B>:
 *
 * 0 - \f$\epsilon_{\rm{disk}}=\epsilon
 *      \biggl[0.5+\left(\frac{V_{\rm{max}}}{70km/s}\right)^{-\beta_1}\biggr]\f$,
 *     \f$\epsilon_{\rm{halo}}=\eta
 *      \biggl[0.5+\left(\frac{V_{\rm{max}}}{70km/s}\right)^{-\beta_2}\biggr]\f$
 *     (Eqs. 19 & 21 Guo2010)(FeedbackReheatingModel = 0)
 *
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
  /* Variables: reff-Rdisk, tdyn=Rdisk/Vmax, strdot=Mstar_dot, stars=strdot*dt*/
  double tdyn, strdot=0., stars, cold_crit;
  int ii;
#ifdef H2_AND_RINGS
  double strdotr[RNUM], starsRings[RNUM];
  double sfe, cold_crit_rate, SigmaGas, SigmaGasratio;
  //double sfe, vmax, cold_crit_rate, SigmaGas, SigmaGasratio, sigmah50=0.0;
  int j;
#endif
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  double metallicitySF;
#endif
#endif
  double Vmax, gas_radius;

  if(Gal[p].Type == 0)
    Vmax=Gal[p].Vmax;
  else
    Vmax=Gal[p].InfallVmax;

  gas_radius = Gal[p].ColdGasRadius;
  //the h factor in Gal[p].ColdGasRadius goes into cold_crit
  tdyn = gas_radius / Vmax;
  cold_crit = SfrColdCrit * Vmax/200. * gas_radius*100.;

  //standard star formation law (Croton2006, Delucia2007, Guo2010, Henriques2015)
  if(StarFormationModel == 0)
    {
      if(Gal[p].ColdGas > cold_crit)
	strdot = SfrEfficiency * (Gal[p].ColdGas - cold_crit) / tdyn;
      else
	strdot = 0.0;
    }

#ifdef H2_AND_RINGS
  update_h2fraction(p);

  sfe=SfrEfficiency*UnitTime_in_years/Hubble_h; //convert from yr-1 into code units of time // the unit of sfe here is (km/s)/(Mpc/h)=1.022e-3 h^{+1} /Gyr
  //sfe=SfrEfficiency*UnitTime_in_years/Hubble_h/pow((1+ZZ[Gal[p].SnapNum]),1.0); //convert from yr-1 into code units of time

  if(SFRtdyn==1)
    sfe= (sfe/tdyn)/UnitTime_in_years*Hubble_h; // for star formation rate proportional to 1/t_dyn
    //sfe= sfe/1.8/tdyn; // for star formation rate proportional to 1/t_dyn

  for(j=0;j<RNUM;j++)
    {
      if(StarFormationModel == 0)
	  strdotr[j]=strdot*Gal[p].ColdGasRings[j]/Gal[p].ColdGas;
      else if(StarFormationModel == 2)
	{
	  if(Gal[p].Type == 0)
	    cold_crit_rate = SfrColdCrit * Gal[p].Vmax/200. * Gal[p].ColdGasRadius/Gal[p].ColdGas*100.;
	  else
	    cold_crit_rate = SfrColdCrit * Gal[p].InfallVmax/200. * Gal[p].ColdGasRadius/Gal[p].ColdGas*100.;

	  if(cold_crit_rate < 1 && cold_crit_rate>=0)
	    strdotr[j] = sfe * Gal[p].ColdGasRings[j] * (1 - cold_crit_rate);
	  else strdotr[j] = 0.0;
       	}

      else if(StarFormationModel == 3) /*The star formation law in Krumholz et al. 2009*/
	{
	  double SigmaGas0=85.0, SF_Law_pow=0.33;

	  if(j==0)
	    SigmaGas = Gal[p].ColdGasRings[j] / (M_PI* RingRadius[j]*RingRadius[j])/WARM_PHASE_FACTOR*Clumpingfactor;
	  else
	    SigmaGas = Gal[p].ColdGasRings[j] / (M_PI*(RingRadius[j]*RingRadius[j]-RingRadius[j-1]*RingRadius[j-1]))/WARM_PHASE_FACTOR*Clumpingfactor;

	  /* convert from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2) */
	  SigmaGas=SigmaGas*0.01*Hubble_h;
	  SigmaGasratio=pow(SigmaGas/SigmaGas0,SF_Law_pow);

	  if(SigmaGasratio<1.0 && SigmaGasratio>0.0)
	    strdotr[j] = sfe/SigmaGasratio * Gal[p].ColdGasRings[j]*Gal[p].H2fractionRings[j]/WARM_PHASE_FACTOR ;
	  //Only cold H2 component is proportional to star formation rate.
	  else
	    if(SigmaGasratio>=1.0)
	      strdotr[j] = sfe*SigmaGasratio * Gal[p].ColdGasRings[j]*Gal[p].H2fractionRings[j]/WARM_PHASE_FACTOR ;
	    else strdotr[j]=0.0;
	}

      else if(StarFormationModel == 4)	/*The star formation law in Fu et al. 2010*/
	{
	  //double sigma_H2, N_sf=1.0, sigma2_crit=70, area;

	  if(Gal[p].H2fractionRings[j]>=0.0)
	      strdotr[j] = sfe * Gal[p].ColdGasRings[j]*Gal[p].H2fractionRings[j] / WARM_PHASE_FACTOR;
	  else strdotr[j]=0.0;
	}

      else  strdotr[j] = 0.0;
    }

  for (j=0,strdot=0;j<RNUM;j++)
    strdot+=strdotr[j];

#endif // H2_AND_RINGS



  /* Note that Units of dynamical time are Mpc/Km/s - no conversion on dt needed
   * be mentioned 3.06e19 to 3.15e19 */

  if(strdot < 0.0)
    strdot =0.;
    //terminate("***error stars<0.0***\n");

  stars = strdot * dt;

#ifdef H2_AND_RINGS
  for(j=0;j<RNUM;j++)
    {
      if(strdotr[j] < 0.0)
	strdotr[j] =0.;
      starsRings[j] = strdotr[j] * dt;
      if(starsRings[j] <0.0) starsRings[j] =0.0;
    }
#endif

  //otherwise this check is done inside update_stars_due_to_reheat for stars+reheatedmass!>coldgas
#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN
  if(stars > Gal[p].ColdGas)
  	stars = Gal[p].ColdGas;
#ifdef H2_AND_RINGS
  for(j=0;j<RNUM;j++)
    if(starsRings[j]>Gal[p].ColdGasRings[j])
      starsRings[j]=Gal[p].ColdGasRings[j];
#endif
#endif

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
  mass_checks(centralgal,"model_starformation_and_feedback.c",__LINE__);

  /* Store the value of the metallicity of the cold phase when SF occurs */

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  metallicitySF=0.;
  if (Gal[p].ColdGas > 0.)
    {
      for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	metallicitySF+= Gal[p].MetalsColdGas[ii];
      metallicitySF/=Gal[p].ColdGas;

    }
#endif
#endif

//if FEEDBACK_COUPLED_WITH_MASS_RETURN feedback happens only when stars die,
//there is no need to balance it with SF
#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
  if (stars > 0.)
#ifndef H2_AND_RINGS
    update_stars_due_to_reheat(p, centralgal, &stars);
#else
    update_stars_due_to_reheat(p, centralgal, &stars, starsRings);
#endif
#endif //FEEDBACK_COUPLED_WITH_MASS_RETURN


/* if update_stars_due_to_reheat is commented out, uncomment this lines:
   if(stars > Gal[p].ColdGas)
    stars = Gal[p].ColdGas;
 #ifdef H2_AND_RINGS
   for(j=0;j<RNUM;j++)
     if(starsRings[j]>Gal[p].ColdGasRings[j])
       starsRings[j]=Gal[p].ColdGasRings[j];
 #endif*/


  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
  mass_checks(centralgal,"model_starformation_and_feedback.c",__LINE__);

  /*  update the star formation rate */
   /*Sfr=stars/(dt*steps)=strdot*dt/(dt*steps)=strdot/steps -> average over the STEPS*/
   Gal[p].Sfr += stars / (dt * STEPS);
 #ifdef H2_AND_RINGS
   for(j=0;j<RNUM;j++)
     Gal[p].SfrRings[j] += starsRings[j] / (dt * STEPS);
 #endif

  // update_from_star_formation can only be called
  // after SD_feeedback recipe since stars need to be re_set once the reheated mass is known
  // (star formation and feedback share the same fraction of cold gas)
  if (stars > 0.)
 #ifndef H2_AND_RINGS
    update_from_star_formation(p, stars, "insitu", nstep); // false indicates not a burst
 #else
    update_from_star_formation(p, stars, starsRings, "insitu", nstep); // false indicates not a burst
 #endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  //  Update the luminosities due to the stars formed
  if (stars > 0.0)
    add_to_luminosities(p, stars, time, dt, metallicitySF);
#endif
#endif

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
  mass_checks(centralgal,"model_starformation_and_feedback.c",__LINE__);

  update_massweightage(p, stars, time);

#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
  /* ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN feedback is only called when stars die,
   * inside DETAILED_METALS_AND_MASS_RETURN */
  if (stars > 0.)
#ifndef H2_AND_RINGS
    SN_feedback(p, centralgal, stars, "ColdGas");
#else
    SN_feedback(p, centralgal, stars, starsRings, "ColdGas");
#endif
#endif

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
  mass_checks(centralgal,"model_starformation_and_feedback.c",__LINE__);


  if(DiskInstabilityModel==0)
    if(Gal[p].DiskMass > 0.0)
      check_disk_instability(p,dt);

  if (DiskRadiusModel == 0)
    Gal[p].DiskRadius=get_stellar_disk_radius(p);

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
  mass_checks(centralgal,"model_starformation_and_feedback.c",__LINE__);

}




#ifndef H2_AND_RINGS
void update_stars_due_to_reheat(int p, int centralgal, double *stars)
#else
void update_stars_due_to_reheat(int p, int centralgal, double *stars, double starsRings[])
#endif
{
  double reheated_mass, frac, Radius_low=0.,totmetals;
  int ii;

#ifndef H2_AND_RINGS
  totmetals=0.;
  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
    totmetals+= Gal[p].MetalsColdGas[ii];
  reheated_mass=compute_SN_reheat(p, centralgal, *stars, Gal[p].ColdGas, totmetals, Radius_low, Gal[p].ColdGasRadius);
  if((*stars + reheated_mass) > Gal[p].ColdGas)
    {
      frac = Gal[p].ColdGas / (*stars + reheated_mass);
      *stars *= frac;
    }
#else
  int jj;
  for(jj=0;jj<RNUM;jj++)
    {
      if(jj>0)
    	Radius_low=RingRadius[jj-1];

      totmetals=0.;
      for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	totmetals+= Gal[p].MetalsColdGasRings[jj][ii];

      reheated_mass=compute_SN_reheat(p, centralgal, starsRings[jj], Gal[p].ColdGasRings[jj], totmetals, Radius_low, RingRadius[jj]);
      if((starsRings[jj] + reheated_mass) > Gal[p].ColdGasRings[jj])
	{
	  frac = Gal[p].ColdGasRings[jj] / (starsRings[jj] + reheated_mass);
  	  starsRings[jj] *= frac;
  	}
       *stars+=starsRings[jj];
    }
#endif
}




/** @brief Updates the different components due to star formation: mass
  *        and metals in stars and cold gas and stellar spin. */
//void update_from_star_formation(int p, double time, double stars, double metallicity)
#ifndef H2_AND_RINGS
void update_from_star_formation(int p, double stars, char type_of_event[], int nstep)
#else
void update_from_star_formation(int p, double stars, double starsRings[], char type_of_event[], int nstep)
#endif
{
  int ii;
  double stars_to_add=0., NonRecycledFraction=0.;
#ifdef H2_AND_RINGS
  int jj;
  double stars_to_addr[RNUM], fractionRings[RNUM];
#else
  double fraction;
#endif
  if(Gal[p].ColdGas <= 0. || stars <= 0.) {
    printf("Gal[p].ColdGas <= 0. || stars <= 0., Coldgas=%0.5e stars=%0.5e, in function update_from_star_formation, model_starformation_and_feedback.c line:%d\n",Gal[p].ColdGas, stars,__LINE__);
    exit(0);
  }

  /* If DETAILED_METALS_AND_MASS_RETURN, no longer an assumed instantaneous
   * recycled fraction. Mass is returned over time via SNe and AGB winds.
   * Update the Stellar Spin when forming stars */
#ifndef DETAILED_METALS_AND_MASS_RETURN
  NonRecycledFraction=(1 - RecycleFraction);
#else
  NonRecycledFraction=1.;
#endif

#ifndef H2_AND_RINGS
  stars_to_add=NonRecycledFraction*stars;
#else
  for(jj=0;jj<RNUM;jj++)
    {
      stars_to_addr[jj]=NonRecycledFraction * starsRings[jj];
      stars_to_add+=stars_to_addr[jj];
    }
#endif //H2_AND_RINGS

  if (Gal[p].DiskMass+stars_to_add > 1.e-8)
    for (ii = 0; ii < 3; ii++)
      Gal[p].DiskSpin[ii]=((Gal[p].DiskSpin[ii])*(Gal[p].DiskMass)+stars_to_add*Gal[p].ColdGasSpin[ii])/(Gal[p].DiskMass+stars_to_add);

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);


#ifdef H2_AND_RINGS
  /*  Update Gas and Metals from star formation */
  for(jj=0;jj<RNUM;jj++)
   if(Gal[p].ColdGasRings[jj]>0.)
     fractionRings[jj]=stars_to_addr[jj]/Gal[p].ColdGasRings[jj];
   else
     fractionRings[jj]=0.;
  transfer_material_with_rings(p,"DiskMass",p,"ColdGas",fractionRings,"model_starformation_and_feedback.c", __LINE__);
#else
  fraction=stars_to_add/Gal[p].ColdGas;
  transfer_material(p,"DiskMass",p,"ColdGas",fraction,"model_starformation_and_feedback.c", __LINE__);
#endif




#ifdef TRACK_MASSGROWTH_CHANNELS
  //for this calculation we want just the long lived mass and
  //take the instantaneous recycling aproximation even
  //for the detailed chemical enrichment because it is not
  //possible to know which component to eject mass from afterwards
  double long_lived_mass;
  long_lived_mass=stars_to_add;
#ifdef DETAILED_METALS_AND_MASS_RETURN
  long_lived_mass*=(1 - RecycleFraction);
#endif

  if(strcmp(type_of_event,"insitu")==0)
      {
        Gal[p].MassFromInSitu+=long_lived_mass;
#ifdef STAR_FORMATION_HISTORY
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
        Gal[p].sfh_MassFromInSitu[Gal[p].sfh_ibin]+=long_lived_mass;
#endif
#endif
      }

  if(strcmp(type_of_event,"merger")==0)
    {
      Gal[p].MassFromBursts+=long_lived_mass;
#ifdef STAR_FORMATION_HISTORY
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
      Gal[p].sfh_MassFromBursts[Gal[p].sfh_ibin]+=long_lived_mass;
#endif
#endif
    }
#endif //TRACK_MASSGROWTH_CHANNELS

#ifdef TRACK_BURST
  if(strcmp(type_of_event,"merger")==0)
    {
      Gal[p].BurstMass+=stars_to_add;
#ifdef STAR_FORMATION_HISTORY
      Gal[p].sfh_BurstMass[Gal[p].sfh_ibin]+=stars_to_add;
#endif
    }
#endif

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

  if (FeedbackReheatingModel == 0 || FeedbackReheatingModel == 1)
    {
 /* stars (instead of star_to_add) used because the Yield is defined as a
  * fraction of all stars formed, not just long lived */
#ifdef DETAILED_METALS_AND_MASS_RETURN
#ifdef METALS_SELF
      Gal[p].MetalsHotGasSelf.str.type2 += Yield * FracZtoHot * stars;
#endif
#else //IF NOT DETAILED_METALS_AND_MASS_RETURN
      //This part is not used if OPT+=DELAYED_ENRICHMENT_AND MASS_RETURN as yield
      //and recycling fraction are not fixed:
#ifndef H2_AND_RINGS
      Gal[p].MetalsColdGas[0] += Yield * (1.-FracZtoHot) * stars;
#else
      for(jj=0;jj<RNUM;jj++)
	{
	  Gal[p].MetalsColdGasRings[jj][0] += Yield* (1.-FracZtoHot) * starsRings[jj];
	  Gal[p].MetalsColdGas[0] += Yield* (1.-FracZtoHot) * starsRings[jj];
	}
#endif
      Gal[Gal[p].CentralGal].MetalsHotGas[0] += Yield * FracZtoHot * stars;
      Gal[Gal[p].CentralGal].HotGas += Yield * FracZtoHot * stars;
#ifdef METALS_SELF
      Gal[p].MetalsHotGasSelf[0] += Yield * FracZtoHot * stars;
#endif
#endif //DETAILED_METALS_AND_MASS_RETURN
    }

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

  if (DiskRadiusModel == 0)
    Gal[p].DiskRadius = get_stellar_disk_radius(p);


}












/* there are two modes for supernova feedback corresponding to when the mass returning
 * by dying stars is returned to the cold gas - reheat and ejection; and when the mass
 * is returned to the hot gas - onle ejection.*/
#ifndef H2_AND_RINGS
void SN_feedback(int p, int centralgal, double stars, char feedback_location[])
#else
void SN_feedback(int p, int centralgal, double stars, double starsRings[], char feedback_location[])
#endif
{
  double EjectVmax, EjectVvir, SN_Energy, Reheat_Energy, ReScaled_EnergySNcode;
  double reheated_mass=0., ejected_mass=0., totmetals;
  double Radius_low;
  int ii;
  /* SN FEEDBACK RECIPES */
#ifdef H2_AND_RINGS
  int jj;
  double reheated_massr[RNUM];
#endif

  Radius_low=0.;
  //REHEAT
#ifndef H2_AND_RINGS
  //when FEEDBACK_COUPLED_WITH_MASS_RETURN some mass goes into HOTGAS and does not produce reheating
  if (strcmp(feedback_location,"HotGas")==0)
    reheated_mass=0;
  else
    {
      totmetals=0.;
      for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
     	totmetals+= Gal[p].MetalsColdGas[ii];
      reheated_mass=compute_SN_reheat(p, centralgal, stars, Gal[p].ColdGas, totmetals, Radius_low, Gal[p].ColdGasRadius);
    }
#else
  reheated_mass=0.0;
  //stars=0;
  for(jj=0;jj<RNUM;jj++)
    {
      if(jj>0)
	Radius_low=RingRadius[jj-1];

      //when FEEDBACK_COUPLED_WITH_MASS_RETURN some mass goes into HOTGAS and does not produce reheating
      if (strcmp(feedback_location,"HotGas")==0)
	  reheated_massr[jj]=0.;
      else
	{
	  totmetals=0.;
	  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	    totmetals+= Gal[p].MetalsColdGasRings[jj][ii];
	  reheated_massr[jj]=compute_SN_reheat(p, centralgal, starsRings[jj], Gal[p].ColdGasRings[jj], totmetals, Radius_low, RingRadius[jj]);
	}

      reheated_mass+=reheated_massr[jj];
      //stars+=starsRings[jj];
    }

  //reheated_mass > Gal[p].ColdGas might happen due to precision
  if(reheated_mass > Gal[p].ColdGas)
     reheated_mass = Gal[p].ColdGas;

#endif


  // Determine ejection (for FeedbackEjectionModel 0 we have the dependence on Vmax) Guo2010 - eq 22
  //EJECT
  if(FeedbackEagleScaling == 1)
     {
       //ReScaled_EnergySNcode=EnergySNcode*EAGLE2015_rescale_of_EnergySN (ColdGas, MetalsColdGas, Radius_low, Radius_high);
       //ReScaled_EnergySNcode=EnergySNcode/pow((1 + ZZ[Halo[Gal[p].HaloNr].SnapNum]),4.0);
       double SigmaGas;
       SigmaGas = Gal[p].ColdGas / (M_PI*(Gal[p].ColdGasRadius*Gal[p].ColdGasRadius-Radius_low*Gal[p].ColdGasRadius));
       // convert from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2)
       SigmaGas=SigmaGas*0.01*Hubble_h;

       ReScaled_EnergySNcode=EnergySNcode*SigmaGas;

       if(ReScaled_EnergySNcode>EnergySNcode)
          ReScaled_EnergySNcode=EnergySNcode;
     }
   else
     ReScaled_EnergySNcode=EnergySNcode;


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
      ejected_mass = (FeedbackEjectionEfficiency* (EtaSNcode * ReScaled_EnergySNcode) * stars *
	  min(1./FeedbackEjectionEfficiency, .5+1/pow(EjectVmax/EjectPreVelocity,EjectSlope)) -
	  reheated_mass*EjectVvir*EjectVvir) /(EjectVvir*EjectVvir);
      //ejected_mass = EtaSNcode * ReScaled_EnergySNcode;
    }
  else if(FeedbackEjectionModel == 1)//the ejected material is assumed to have V_SN
    {
      SN_Energy = .5 * stars * (EtaSNcode * ReScaled_EnergySNcode);
      Reheat_Energy = .5 * reheated_mass * EjectVvir * EjectVvir;

      ejected_mass = (SN_Energy - Reheat_Energy)/(0.5 * FeedbackEjectionEfficiency*(EtaSNcode * ReScaled_EnergySNcode));

      //if VSN^2<Vvir^2 nothing is ejected
      if(FeedbackEjectionEfficiency*(EtaSNcode * ReScaled_EnergySNcode)<EjectVvir*EjectVvir)
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
    {
#ifndef H2_AND_RINGS
      update_from_feedback(p, centralgal, reheated_mass, ejected_mass);
#else
      update_from_feedback(p, centralgal, reheated_mass, ejected_mass,  reheated_massr);
#endif
    }

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

}







/*#ifdef H2_AND_RINGS
double compute_SN_reheat(int p, int centralgal, double stars, double ColdGas, int jj)
#else
double compute_SN_reheat(int p, int centralgal, double stars, double ColdGas)
#endif*/
double compute_SN_reheat(int p, int centralgal, double stars, double ColdGas, double MetalsColdGas, double Radius_low, double Radius_high)
{
  double reheated_mass=0.;
  double ReScaled_EnergySNcode;
  double MergeCentralVvir=0.;

  /* In Guo2010 type 1s can eject, reincorporate gas and get gas from their
   * own satellites (is not sent to the type 0 galaxy as in Delucia2007),
   * for gas flow computations:
   * If satellite is inside Rvir of main halo, Vvir of main halo used
   * If it is outside, the Vvir of its central subhalo is used. */

  if(FeedbackEagleScaling == 1)
    {
      //ReScaled_EnergySNcode=EnergySNcode*EAGLE2015_rescale_of_EnergySN (ColdGas, MetalsColdGas, Radius_low, Radius_high);
      //ReScaled_EnergySNcode=EnergySNcode/pow((1 + ZZ[Halo[Gal[p].HaloNr].SnapNum]),4.0);

      double SigmaGas;
#ifdef H2_AND_RINGS
      SigmaGas = ColdGas/(M_PI *(Radius_high*Radius_high-Radius_low*Radius_low))/WARM_PHASE_FACTOR*Clumpingfactor;
#else
      SigmaGas = ColdGas / (M_PI*(Radius_high*Radius_high-Radius_low*Radius_low));
#endif
      // convert from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2)
      SigmaGas=SigmaGas*0.01*Hubble_h;

      ReScaled_EnergySNcode=EnergySNcode*SigmaGas;

      if(ReScaled_EnergySNcode>EnergySNcode)
         ReScaled_EnergySNcode=EnergySNcode;
    }
  else
    ReScaled_EnergySNcode=EnergySNcode;

  //REHEAT
  if(ColdGas>0.)
    {
  // Feedback depends on the circular velocity of the host halo
  // Guo2010 - eq 18 & 19
      if(FeedbackReheatingModel == 0)
	{
	  if (Gal[Gal[p].CentralGal].Type == 0)
	    reheated_mass = FeedbackReheatingEpsilon * stars * (.5+1./pow(Gal[Gal[p].CentralGal].Vmax/ReheatPreVelocity,ReheatSlope));
	  else
	    reheated_mass = FeedbackReheatingEpsilon * stars * (.5+1./pow(Gal[Gal[p].CentralGal].InfallVmax/ReheatPreVelocity,ReheatSlope));

	  if(FeedbackReheatingDeansityScaling==1)
	    {
	      double SigmaGas;
#ifdef H2_AND_RINGS
	      SigmaGas = ColdGas/(M_PI *(Radius_high*Radius_high-Radius_low*Radius_low))/WARM_PHASE_FACTOR*Clumpingfactor;
#else
	      SigmaGas = ColdGas / (M_PI*(Radius_high*Radius_high-Radius_low*Radius_low));
#endif
	      // convert from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2)
	      SigmaGas=SigmaGas*0.01*Hubble_h;

	      if(SigmaGas>0.)
		reheated_mass*=pow(SigmaGas*0.05,2);
	     // reheated_mass/=pow((1 + ZZ[Halo[Gal[p].HaloNr].SnapNum]),5.0); - much lower reheating, much higher metals
	     // reheated_mass*=pow((1 + ZZ[Halo[Gal[p].HaloNr].SnapNum]),2.0); // much higher reheating, much lower stellar metals, gas metals high-z unnafected
	    }

	  if (reheated_mass * Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir > stars * (EtaSNcode * ReScaled_EnergySNcode))
	    reheated_mass = stars * (EtaSNcode * ReScaled_EnergySNcode) / (Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir);
	}
    }
  else
    reheated_mass=0.;

  if(reheated_mass > ColdGas)
    reheated_mass = ColdGas;

 /*
      double rd, ringtot;
    int jj;
     rd=Gal[p].ColdGasRadius;
    ringtot=1.-(1+RingRadius[RNUM-1]/rd)/exp(RingRadius[RNUM-1]/rd);
    reheated_massr[0]=(1-(1+RingRadius[0]/rd)/exp(RingRadius[0]/rd))/ringtot*reheated_mass;
    //aux_mass=reheated_massr[0];
    for(jj=1; jj<RNUM; jj++)
      {
      reheated_massr[jj]= ((1+RingRadius[jj-1]/rd)/exp(RingRadius[jj-1]/rd)-(1+RingRadius[jj]/rd)/exp(RingRadius[jj]/rd))/ringtot*reheated_mass;
      //aux_mass+=reheated_massr[jj];
      }*/

  return reheated_mass;



}

double EAGLE2015_rescale_of_EnergySN (double ColdGas, double MetalsColdGas, double Radius_low, double Radius_high)
{
double fth_min=0.3, fth_max=3.0, fth;
double metallicity_Z, z_solar=0.0127;
double nH_birth, nH_grams=1.6737236e-24, nH_0=0.67;
double n=2*log(10);

metallicity_Z=MetalsColdGas/ColdGas;
nH_birth=ColdGas*UnitMass_in_g/nH_grams/(4./3.*M_PI*UnitLength_in_cm*(Radius_high*Radius_high-Radius_low*Radius_low));

fth=fth_min+(fth_max-fth_min)/( 1 + pow(metallicity_Z/(0.1*z_solar),n) * pow(nH_birth/nH_0,-1.*n) );

return  fth;
}

/** @brief Updates cold, hot and external gas components due to SN
 *         reheating and ejection. */
#ifndef H2_AND_RINGS
void update_from_feedback(int p, int centralgal, double reheated_mass, double ejected_mass)
#else
void update_from_feedback(int p, int centralgal, double reheated_mass, double ejected_mass, double reheated_massr[])
#endif
{
  double dis=0.;
  double MassRemain=0.;
  double fraction;
  int merger_centre=0;
#ifdef H2_AND_RINGS
  double fractionRings[RNUM], tmpfractionRings[RNUM], MassRemainRings[RNUM];
  int jj;
#endif

  if(Gal[p].ColdGas > 0.) {
    /* REHEAT if galaxy is a type 1 or a type 2 orbiting a type 1 with hot gas
     * being stripped, some of the reheated and ejected masses goes to the type
     * 0 and some stays in the type 1 */
#ifdef H2_AND_RINGS
    for(jj=0;jj<RNUM;jj++)
      if(Gal[p].ColdGasRings[jj]>0.)
	fractionRings[jj]=reheated_massr[jj]/(Gal[p].ColdGasRings[jj]);
      else
	fractionRings[jj]=0.;
#endif

    if(Gal[p].Type ==0)
      {
#ifdef H2_AND_RINGS
	transfer_material_with_rings(p,"HotGas",p,"ColdGas",fractionRings,"model_starformation_and_feedback.c", __LINE__);
	//transfer_material_with_rings(p,"ReheatedGas",p,"ColdGas",fractionRings,"model_starformation_and_feedback.c", __LINE__);
#else
	transfer_material(p,"HotGas",p,"ColdGas",((float)reheated_mass)/((float)Gal[p].ColdGas),"model_starformation_and_feedback.c", __LINE__);
	//transfer_material(p,"ReheatedGas",p,"ColdGas",((float)reheated_mass)/((float)Gal[p].ColdGas),"model_starformation_and_feedback.c", __LINE__);
#endif
      }

    //For satellite galaxies compute how much gas stays in the galaxy and how much goes to central companion
    else if(Gal[p].Type <3)
      {
	if(Gal[p].Type ==1)
	  merger_centre=centralgal;
	else if(Gal[p].Type ==2)
	  merger_centre=Gal[p].CentralGal;

	if(HotGasOnType2Galaxies==0)
	  //if no hot gas in type 2's, share gas between  0 and 1
	  dis=separation_gal(centralgal,Gal[p].CentralGal)/(1+ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);
	else if(HotGasOnType2Galaxies==1)
	  //if hot gas in type 2's, share gas between itself and mergercentre
	  dis=separation_gal(merger_centre,p)/(1+ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);

	//compute share of reheated mass
	if ((dis<Gal[centralgal].Rvir && Gal[Gal[p].CentralGal].Type == 1 &&  HotGasOnType2Galaxies==0) ||
	    (dis<Gal[merger_centre].Rvir && HotGasOnType2Galaxies==1))
	  {
	    //mass that remains on type1 (the rest goes to type 0)
	    //for reheat - MassRemain, for eject - ejected_mass
	    MassRemain=reheated_mass*Gal[p].HotRadius/Gal[p].Rvir;
	    ejected_mass = ejected_mass*Gal[p].HotRadius/Gal[p].Rvir;
	    if (MassRemain > reheated_mass)
	      MassRemain = reheated_mass;
#ifdef H2_AND_RINGS
	    for(jj=0;jj<RNUM;jj++)
	      {
		MassRemainRings[jj]=reheated_massr[jj]*Gal[p].HotRadius/Gal[p].Rvir;
		if (MassRemainRings[jj] > reheated_massr[jj])
		  MassRemainRings[jj] = reheated_massr[jj];
	      }
#endif
	  }
	else
	  {
	    MassRemain=reheated_mass;
#ifdef H2_AND_RINGS
	    for(jj=0;jj<RNUM;jj++)
	      MassRemainRings[jj] = reheated_massr[jj];
#endif
	  }

      //needed due to precision issues, since we first remove MassRemain and
      //then (reheated_mass-MassRemain) from the satellite into the type 0 and
      //type 1 the fraction might not add up on the second call since
      //Gal[p].ColdGas is a float and reheated_mass & MassRemain are doubles
      if((MassRemain + reheated_mass)>Gal[p].ColdGas)
	MassRemain=Gal[p].ColdGas-reheated_mass;

#ifdef H2_AND_RINGS
      for(jj=0;jj<RNUM;jj++)
	if((MassRemainRings[jj] + reheated_massr[jj])>Gal[p].ColdGasRings[jj])
	  MassRemainRings[jj] = Gal[p].ColdGasRings[jj]-reheated_massr[jj];
#endif

      mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
      //transfer MassRemain

      if(reheated_mass>0.)
	{
#ifdef H2_AND_RINGS
	  //for(jj=0;jj<RNUM;jj++)
	  //  tmpfractionRings[jj]=fractionRings[jj]*(MassRemain/reheated_mass);
	  for(jj=0;jj<RNUM;jj++)
	    if(Gal[p].ColdGasRings[jj]>0.)
	      tmpfractionRings[jj]=MassRemainRings[jj]/(Gal[p].ColdGasRings[jj]);
	    else
	      tmpfractionRings[jj]=0.;

	  if(HotGasOnType2Galaxies==0)
	    {
	    //tranfer to itself if type 1, merger centre if type 2
	    if(Gal[p].CentralGal==p)
	      transfer_material_with_rings(Gal[p].CentralGal,"HotGas",p,"ColdGas", tmpfractionRings,"model_starformation_and_feedback.c", __LINE__);
	      //transfer_material_with_rings(Gal[p].CentralGal,"ReheatedGas",p,"ColdGas", tmpfractionRings,"model_starformation_and_feedback.c", __LINE__);
	    else
	      transfer_material_with_rings(Gal[p].CentralGal,"HotGas",p,"ColdGas", tmpfractionRings,"model_starformation_and_feedback.c", __LINE__);
	    }
	  else if(HotGasOnType2Galaxies==1)
	    //tranfer to itself
	    transfer_material_with_rings(p,"HotGas",p,"ColdGas",tmpfractionRings,"model_starformation_and_feedback.c", __LINE__);
	    //transfer_material_with_rings(p,"ReheatedGas",p,"ColdGas",tmpfractionRings,"model_starformation_and_feedback.c", __LINE__);
#else
	  if(HotGasOnType2Galaxies==0)
	    {
	    //tranfer to itself if type 1, merger centre if type 2
	      if(Gal[p].CentralGal==p)
		transfer_material(Gal[p].CentralGal,"HotGas",p,"ColdGas", MassRemain/Gal[p].ColdGas,"model_starformation_and_feedback.c", __LINE__);
	     // transfer_material(Gal[p].CentralGal,"ReheatedGas",p,"ColdGas", MassRemain/Gal[p].ColdGas,"model_starformation_and_feedback.c", __LINE__);
	      else
		transfer_material(Gal[p].CentralGal,"HotGas",p,"ColdGas", MassRemain/Gal[p].ColdGas,"model_starformation_and_feedback.c", __LINE__);
	    }
	  else if(HotGasOnType2Galaxies==1)
	    //tranfer to itself
	    transfer_material(p,"HotGas",p,"ColdGas",MassRemain/Gal[p].ColdGas,"model_starformation_and_feedback.c", __LINE__);
	    //transfer_material(p,"ReheatedGas",p,"ColdGas",MassRemain/Gal[p].ColdGas,"model_starformation_and_feedback.c", __LINE__);
#endif
	}

      mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

      //transfer reheated_mass-MassRemain from galaxy to the type 0
      if (reheated_mass > MassRemain)
	if(Gal[p].ColdGas > 0.)
	  {
	    //if the reheat to itself left cold gas below limit do not reheat
	    //to central
#ifdef H2_AND_RINGS
	    //cannot use tmpfractionRings defined from fractionRings since
	    // Gal[p].ColdGasRings has changed from MassRemain above
	    for(jj=0;jj<RNUM;jj++)
	      if(Gal[p].ColdGasRings[jj]>0.)
		fractionRings[jj]=(reheated_massr[jj]-MassRemainRings[jj])/Gal[p].ColdGasRings[jj];
		//fractionRings[jj]=(reheated_massr[jj]/Gal[p].ColdGasRings[jj])*	((reheated_mass-MassRemain)/reheated_mass);
	      else
		fractionRings[jj]=0.;


	  if(HotGasOnType2Galaxies==0)	 //tranfer to type 0
	    transfer_material_with_rings(centralgal,"HotGas",p,"ColdGas",fractionRings,"model_starformation_and_feedback.c", __LINE__);
	  else if(HotGasOnType2Galaxies==1)//tranfer to merger centre
	      transfer_material_with_rings(merger_centre,"HotGas",p,"ColdGas",fractionRings,"model_starformation_and_feedback.c", __LINE__);

#else
	  ///with rings ColdGas is a double and using (float) might cause
	  //(float)(reheated_mass-MassRemain)/Gal[p].ColdGas to be >1
	  if(HotGasOnType2Galaxies==0) //tranfer to type 0
	    transfer_material(centralgal,"HotGas",p,"ColdGas", (float)(reheated_mass-MassRemain)/Gal[p].ColdGas,"model_starformation_and_feedback.c", __LINE__);
	  else if(HotGasOnType2Galaxies==1) //tranfer to merger centre
	    transfer_material(merger_centre,"HotGas",p,"ColdGas", (float)(reheated_mass-MassRemain)/Gal[p].ColdGas,"model_starformation_and_feedback.c", __LINE__);
#endif //H2_AND_RINGS
	  }

      //} //Gal[p].Type !=2
      mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
    }//types

  }//if(Gal[p].ColdGas > 0.)

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

   //DO EJECTION OF GAS
  if ( (Gal[Gal[p].CentralGal].HotGas > 0. && HotGasOnType2Galaxies==0) ||
       (Gal[p].HotGas > 0. && HotGasOnType2Galaxies==1) ) {
    
    if(HotGasOnType2Galaxies==0) {
      if (ejected_mass > Gal[Gal[p].CentralGal].HotGas)
	//either eject own gas or merger_centre gas for ttype 2's
	ejected_mass = Gal[Gal[p].CentralGal].HotGas;
      fraction=ejected_mass/Gal[Gal[p].CentralGal].HotGas;
      
    }
    else if(HotGasOnType2Galaxies==1) {
      if (ejected_mass > Gal[p].HotGas && HotGasOnType2Galaxies==1)
	ejected_mass = Gal[p].HotGas;  //always eject own gas
      fraction=ejected_mass/Gal[p].HotGas;
    }

    if (Gal[Gal[p].CentralGal].Type == 1) {      /* If type 1, or type 2 orbiting type 1 near type 0 */

      if (FateOfSatellitesGas == 0) {
	if(HotGasOnType2Galaxies==0)
	  transfer_material(Gal[p].CentralGal,"EjectedMass",Gal[p].CentralGal,"HotGas",fraction,"model_starformation_and_feedback.c", __LINE__);
	else if(HotGasOnType2Galaxies==1)
	  transfer_material(Gal[p].CentralGal,"EjectedMass",p,"HotGas", fraction,"model_starformation_and_feedback.c", __LINE__);
      }
      else if (FateOfSatellitesGas == 1) {
	if (dis < Gal[centralgal].Rvir)
	  transfer_material(centralgal,"HotGas",Gal[p].CentralGal,"HotGas",fraction,"model_starformation_and_feedback.c", __LINE__);
	else
	  transfer_material(Gal[p].CentralGal,"EjectedMass",Gal[p].CentralGal,"HotGas",fraction,"model_starformation_and_feedback.c", __LINE__);
      }
    }
    else {
      // If galaxy type 0 or type 2 merging into type 0
      if(HotGasOnType2Galaxies==0)
	transfer_material(centralgal,"EjectedMass",Gal[p].CentralGal,"HotGas",fraction,"model_starformation_and_feedback.c", __LINE__);
      else if(HotGasOnType2Galaxies==1)
	transfer_material(centralgal,"EjectedMass",p,"HotGas",fraction,"model_starformation_and_feedback.c", __LINE__);
    }

  }//(Gal[Gal[p].CentralGal].HotGas > 0.)

#ifdef H2_AND_RINGS
  update_h2fraction(p);
#endif

  mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

}


//Age in Mpc/Km/s/h - code units
void update_massweightage(int p, double stars, double time)
{
  int outputbin;
  double age;

  for(outputbin = 0; outputbin < NOUT; outputbin++)
    {
      age = time - NumToTime(ListOutputSnaps[outputbin]);
#ifdef DETAILED_METALS_AND_MASS_RETURN
      Gal[p].MassWeightAge[outputbin] += age * stars;
#else
      Gal[p].MassWeightAge[outputbin] += age * stars * (1. - RecycleFraction);
#endif
    }
}





/** @brief Checks for disk stability using the
 *         Mo, Mao & White (1998) criteria */

void check_disk_instability(int p, double dt)
{

  double Mcrit, fraction, stars, diskmass;
#ifdef H2_AND_RINGS
  double rstar, vmax;
  int j;
#endif
/** @brief Calculates the stability of the stellar disk as discussed
 *         in Mo, Mao & White (1998). For unstable stars, the required
 *         amount is transfered to the bulge to make the disk stable again.
 *         Mass, metals and luminosities updated. After Guo2010 the bulge
 *         size is followed and needs to be updated.
 *         Eq 34 & 35 in Guo2010 are used. */

  diskmass = Gal[p].DiskMass;

#ifndef H2_AND_RINGS
  /* check stellar disk -> eq 34 Guo2010*/
  if (Gal[p].Type != 0)
    Mcrit = Gal[p].InfallVmax * Gal[p].InfallVmax * Gal[p].DiskRadius / G;
  else
    Mcrit = Gal[p].Vmax * Gal[p].Vmax * Gal[p].DiskRadius / G;
#else
   if(diskmass<1.0e-6)
     rstar=0.5*RingRadius[0];
   else
     {
       rstar=0.5*RingRadius[0]*Gal[p].DiskMassRings[0];
       for(j=1;j<RNUM;j++)
	 rstar+=0.5*(RingRadius[j-1]+RingRadius[j])*Gal[p].DiskMassRings[j];
       rstar=rstar/diskmass/2.0;      //2.0=mean radius/scale length for exponential disk
     }

   if (Gal[p].Type != 0)
     vmax=Gal[p].InfallVmax;
   else
     vmax=Gal[p].Vmax;

   Mcrit = vmax * vmax * rstar / G;
#endif

   mass_checks(p,"model_starformation_and_feedback.c",__LINE__);

   stars = diskmass - Mcrit;
   fraction = stars / diskmass;

   /* add excess stars to the bulge */
   if(stars > 0.0)
     {
       /* to calculate the bulge size */
       update_bulgesize_from_disk_instability(p,stars);


       //The bulge will be formed in the same place as the disk was, so the disk rings
       //are transferred directly into bulge rings
#ifdef H2_AND_RINGS
       double dmass, fractionRings[RNUM];
       for(j=0;j<RNUM;j++)
	 fractionRings[j]=0.;

       dmass=stars;
       j=0; //avoid non-definded j if dmass<1e-6
      // if(dmass>1.0e-6)
	 for(j=0;j<RNUM;j++)
	   {
	     //mass is transfered first from the inner rings
	     //until the necessary mass is achieved
	     if(dmass>Gal[p].DiskMassRings[j])
	       {
		 dmass-=Gal[p].DiskMassRings[j];
		 fractionRings[j]=1.;
	       }
	     else break;
	   }

       //check needed in case there is a ring with 0 mass in the middle
       if(Gal[p].DiskMassRings[j]>0.)
	 fractionRings[j]=dmass/Gal[p].DiskMassRings[j];
       else
	 fractionRings[j]=0.;


       //for(j=0;j<RNUM;j++)
	// fractionRings[j]=dmass/Gal[p].DiskMass;
       transfer_material_with_rings(p,"BulgeMass",p,"DiskMass",fractionRings,"model_starformation_and_feedback.c", __LINE__);
       mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
#else
       transfer_material(p,"BulgeMass",p,"DiskMass",fraction, "model_starformation_and_feedback.c", __LINE__);
#endif

       if(BHGrowthInDiskInstabilityModel == 1)
	 if(Gal[p].ColdGas > 0.)
	   {
#ifdef H2_AND_RINGS

	     for(j=0;j<RNUM;j++)
	       {
		 fractionRings[j]*=0.1*BlackHoleGrowthRate  / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir)));
		 fractionRings[j] = min(1.0,fractionRings[j]);
	       }


	    // for(j=0;j<RNUM;j++)
	     //  fractionRings[j]*=BlackHoleGrowthRate*Gal[p].BlackHoleMass/(Gal[p].DiskMass+Gal[p].BulgeMass);

	     transfer_material_with_rings(p,"BlackHoleMass",p,"ColdGas",fractionRings,"model_starformation_and_feedback.c", __LINE__);
	     mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
#else
	     fraction*=0.1*BlackHoleGrowthRate  / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir)));
	     transfer_material(p,"BlackHoleMass",p,"ColdGas",fraction, "model_starformation_and_feedback.c", __LINE__);
#endif
	     Gal[p].QuasarAccretionRate += fraction*Gal[p].ColdGas / (dt*STEPS);

	   }


#ifdef BULGESIZE_DEBUG
       mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
       if ((Gal[p].BulgeMass > TINY_MASS && Gal[p].BulgeSize < TINY_LENGTH)||
	   (Gal[p].BulgeMass < TINY_MASS && Gal[p].BulgeSize > TINY_LENGTH)) {
	   printf("BulgeMass=%g BulgeSize=%g\n",Gal[p].BulgeMass,Gal[p].BulgeSize);
	   terminate("bulgesize wrong in disk instablility\n");
       }
#endif

	/*mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
       if ((Gal[p].BulgeMass > 1e-9 && Gal[p].BulgeSize == 0.0)||
       (Gal[p].BulgeMass == 0.0 && Gal[p].BulgeSize >1e-9))
	 {
	   char sbuf[1000];
	   sprintf(sbuf, "bulgesize wrong in disk instablility.c \n");
	   printf("BulgeMass=%g BulgeSize=%g\n",Gal[p].BulgeMass,Gal[p].BulgeSize);
	   terminate(sbuf);
	 }*/

     }// if(stars > 0.0)
   mass_checks(p,"model_starformation_and_feedback.c",__LINE__);
} //end check_disk_instability


/** @brief Introduced in Guo2010 to track the change in size of bulges
 *         after their growth due to disk instabilities. */

void update_bulgesize_from_disk_instability(int p, double stars)
{      
  double bulgesize, diskmass, fint, massfrac;
  int  j;


/** @brief Updates bulge from disk instability -> stars represents the mass
  *        transfered to the bulge, which occupies a size in the bulge equal
  *        to that occupied in the disk. */


   // alpha_inter=2.0/C=0.5 (alpha larger than in mergers since
   // the existing and newly formed bulges are concentric)
   fint=4.0;

   // disk instabilities are assumed to NOT remove angular momentum.
   // Since the mass of the disk changes, the spin is changed by the same
   // amount to keep angular momentum constant


  diskmass=Gal[p].DiskMass;
  massfrac=stars/diskmass;
  for (j = 0; j <3 ; j++)
    {
      if(massfrac==1) //everything transferred to the bulge
	Gal[p].DiskSpin[j]=0;
      else
	Gal[p].DiskSpin[j]=Gal[p].DiskSpin[j]/(1-massfrac);
    }
  if (DiskRadiusModel == 0)
    Gal[p].DiskRadius = get_stellar_disk_radius(p);



#ifndef H2_AND_RINGS

#ifdef BULGESIZE_DEBUG
  double orisize;
  orisize=Gal[p].BulgeSize;
#endif




  //GET BULGE SIZE - Eq. 35 in Guo2010
  /* size of newly formed bulge, which consists of the stellar mass transfered
   * from the disk. This is calculated using bulge_from_disk which receives
   * Delta_M/DiskMass and returns Rb/Rd. From eq 35 and since
   * DiskMass=2PISigma(Rd)^2 we see that
   * Delta_M/DiskMass=1-(1+Rb/Rd)*exp(-Rb/Rd), so function bulge_from_disk
   * avoids calculating the slow "ln" function */
  bulgesize=bulge_from_disk(stars/diskmass)*Gal[p].DiskRadius/3.;
  if(Gal[p].BulgeMass < TINY_MASS) {
      /* if previous Bulge Mass = 0
       * -> bulge size is given directly from newly formed bulge */
      Gal[p].BulgeSize=bulgesize;
  } else {
      /* combine the old with newly formed bulge and calculate the
       * bulge size assuming energy conservation as for mergers but
       * using alpha=2. - eq 33 */
      Gal[p].BulgeSize=(Gal[p].BulgeMass+stars)*(Gal[p].BulgeMass+stars) /
	  (Gal[p].BulgeMass*Gal[p].BulgeMass/Gal[p].BulgeSize + stars*stars/bulgesize + fint*Gal[p].BulgeMass*stars/(Gal[p].BulgeSize+bulgesize));
  }

  /* Added by PAT to see if the cause of low bulgesizes could be propagation of
   * tightly-bound bulges.  These originate in unfeasibly small disks. */
#define BULGESIZE_MIN 1e-4
  Gal[p].BulgeSize=max(Gal[p].BulgeSize,BULGESIZE_MIN);

#ifdef BULGESIZE_DEBUG
  if((Gal[p].BulgeMass + stars > TINY_MASS && Gal[p].BulgeSize < TINY_LENGTH)
     || (Gal[p].BulgeMass + stars < TINY_MASS && Gal[p].BulgeSize > TINY_LENGTH)) {
      printf("Original DiskMass=%e, DiskSize=%e\nOriginal BulgeMass=%e, BulgeSize=%e\nTransferred stars=%e, bulgesize=%e\nFinal BulgeMass=%e, BulgeSize=%e\n",
	     Gal[p].DiskMass, Gal[p].DiskRadius, Gal[p].BulgeMass, orisize, stars, bulgesize, Gal[p].BulgeMass+stars, Gal[p].BulgeSize);
      terminate("bulgesize or mass wrong in disk instablility");
    }
#endif
  
#else //H2_AND_RINGS

  /*size of new formed bulge, which consist of the stellar mass trasfered from the disk*/
  /*combine the old bulge with the new materials and caculate the bulge size assuming energy conservation */
  diskmass=stars;
  //j=0;
 // if(diskmass>1.0e-6)
  //  {
      for(j=0;j<RNUM;j++)
  	{
	  //mass is transfered first from the inner rings first
	  //until the necessary mass is achieved
	  if(diskmass>Gal[p].DiskMassRings[j])
	      diskmass-=Gal[p].DiskMassRings[j];
	  else break;
  	}
      if(j==RNUM)
	bulgesize=RingRadius[RNUM-1];
      else
	{
	  if(j==0)
	    bulgesize=diskmass/Gal[p].DiskMassRings[j]*RingRadius[j];
	  else
	    bulgesize=diskmass/Gal[p].DiskMassRings[j]*RingRadius[j]+(1-diskmass/Gal[p].DiskMassRings[j])*RingRadius[j-1];
	}
  //  }
  //else bulgesize=0.5*RingRadius[0];

  if(Gal[p].BulgeMass <1.e-9)
    Gal[p].BulgeSize=bulgesize;
  else
    Gal[p].BulgeSize=(Gal[p].BulgeMass+stars)*(Gal[p].BulgeMass+stars)/
    (Gal[p].BulgeMass*Gal[p].BulgeMass/Gal[p].BulgeSize+stars*stars/bulgesize+fint*Gal[p].BulgeMass*stars/(Gal[p].BulgeSize+bulgesize));

    /*if ((Gal[p].BulgeMass+stars > 1.e-8 && Gal[p].BulgeSize == 0.0)||(Gal[p].BulgeMass+stars == 0 && Gal[p].BulgeSize >1.e-8))
      {
        printf("bulgesize wrong in disk instablility. Diskmass %f, bulgemass %f, bulgesize %f, coldgas %f, masstransfer %f transsize %f\n",
	       Gal[p].DiskMass, Gal[p].BulgeMass, Gal[p].BulgeSize, Gal[p].ColdGas, stars, bulgesize);
        exit(0);
      }*/

#endif //H2_AND_RINGS

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

