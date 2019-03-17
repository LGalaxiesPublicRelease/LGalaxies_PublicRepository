#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @file recipe_cooling.c
 *  @brief recipe_cooling.c calculates the amount of mass that cools
 *  from the hot to the cold phase at each timestep (this fraction
 *  is then reduced by AGN heating) and updates the hot and cold gas
 *  fractions accordingly.
 *
 * This recipe calculates the amount of mass that cools from the hot
 * to the cold phase at each timestep. Two different infalling regimes
 * are assumed depending on the redshift and mass of the halos. This is
 * the recipe proposed by (White & Rees, 1978) and that has since been
 * used in most SA.
 *
 * Before only central galaxies could cool gas, so R_hot was just set
 * equal to Rvir. After Guo2010, satellites can have cooling, but since
 * they are loosing dark matter, Rvir is not a good approximation of R_Hot
 * so Gal.HotRadius was introduced.
 *
 * -> At early times and for low mass halos, the cooling radius can be
 * larger then the virial radius. In this case, despite the infalling
 * gas being shock heated to the virial temperature, it condenses within
 * a halo dynamical time and a quasi-static atmosphere cannot form.
 * Single line on the code, with no calculations needed, just that
 * \f$t_{\rm{dyn,h}}=R_{\rm{vir}}/V_{\rm{vir}}\f$
 * and therefore
 * \f$\dot{M}_{\rm{cool}}=0.5\frac{m_{\rm{hot}}}{t_{\rm{dyn}}}
 * =0.5m_{\rm{hot}}\frac{V_{\rm{vir}}}{R_{\rm{vir}}}\f$
 *
 *
 * -> For massive halos and at late times, the cooling radius lies within
 * the virial radius and the infalling gas does form a quasi-static hot
 * atmosphere that extends to the virial radius. This gas can cool at
 * later times and its accretion into central regions is modeled through
 * a cooling flow.
 *
 * \f$\dot{M}_{\rm{cool}}=
 *  m_{\rm{hot}}\frac{r_{\rm{cool}}}{R_{\rm{vir}}}\frac{1}{t_{\rm{cool}}}\f$
 *  eq5 (no 0.5 factor...)
 *
 *
 * In both cases the condensed gas settles into a central disk where
 * star formation occurs. The cooling rate is given by the cooling-AGN
 * heating.
 *
 * BH quiet accretion rate:
 *  \f$\dot{m}_{\rm BH,R}=k_{\rm AGN}\left(\frac{m_{\rm
    BH}}{10^8\,M_{\odot}}\right)\left(\frac{f_{\rm
    hot}}{0.1}\right)\left(\frac{V_{\rm
    vir}}{200\,\mathrm{km\,s}^{-1}}\right)^3\f$

   Luminosity from quiet accretion:
   \f$L_{\rm BH}=\eta\,\dot{m}_{\rm BH,R}\,c^2,\f$

   Corresponding reduction in cooling:
    \f$\dot{m}'_{\rm{cool}}=
    \dot{m}_{\rm{cool}}-\frac{L_{BH}}{\frac{1}{2}V^2_{\rm{vir}}}\f$
 *
 *
*/

/** @brief main cooling recipe, where the cooling rates are calculated (in the absence of any heating). */

void compute_cooling(int p, double dt, int ngal)
{
  double Vvir, Rvir, Mvir, x, lambda, tdyn_halo, rcool, temp, HotGas, MetalsHotGas, HotRadius;
  double coolingGas, logZ, rho_rcool, rhorsq0;
  /* Note that mu could be determined from the metallicity of the gas;
   * here we assume that it is fixed.  
   * The value is chosen to match the HWT15 version of the code */
  const double MU=0.59265;
  //MU_RATIO is 2*(n_e/n)*(n_i/n)/(3*mu)
  const double MU_RATIO=0.28086;
  double mumH=MU*PROTONMASS;
#ifdef BETAPROF
  /* Ratio of virial radius to core radius of beta profile.
   * For now fix by hand; later can free it up as a parameter. */
  const double X=10.;
  const double atanX=atan(X);
  //LAMBDA_RATIO is (n_e/n)*(n_i/n)
  const double LAMBDA_RATIO=0.25;
  double a, f, dt_ratio, tau_ratio, teq_ratio, tauCoolP, fg, fg0, Mvir;
#endif

  mass_checks(p,"cooling.c",__LINE__);

  HotGas = Gal[p].HotGas;
  if(HotGas > TINY_MASS) {
      MetalsHotGas = metals_total(Gal[p].MetalsHotGas);
      Mvir = Gal[p].Mvir;
      Vvir = Gal[p].Vvir;
      Rvir = Gal[p].Rvir;
      // As implemented in Hen15, dynamical time is defined at the edge of the halo
      tdyn_halo = Rvir / Vvir; // tdyn_halo = dynamical time at edge of halo

    /* temp -> Temperature of the Gas in Kelvin, obtained from
     * hydrostatic equilibrium k_B*T=0.5*mu*m_H*(Vc)^2 assuming Vvir~Vc.
     * Note that the constant below evaluatess to 36.345 (listed as 35.9 in old code)*/
      temp = 0.5*mumH*pow2(UnitVelocity_in_cm_per_s*Vvir)/BOLTZMANN;

      if(MetalsHotGas > 0)
	  logZ = log10(MetalsHotGas / HotGas);
      else
	  logZ = -10.0;

      // Cooling rate (/n^2)
      lambda = get_metaldependent_cooling_rate(log10(temp), logZ); // erg cm^3 / s

      // Uncertain if Gal[p].HotRadius is defined for centrals, hence:
      if (Gal[p].Type == 0)
	  HotRadius = Rvir;
      else 
	  HotRadius = Gal[p].HotRadius;
    
#ifdef BETAPROF
      a = Rvir/X; // Core radius of beta profile.
      f = X - atanX; // Fractional mass cf power law solution with same outer density.
      // eq 26 & 4 of BTH18
      //LAMBDA_RATIO is to convert from n_e n_i lambda to n^2 lambda;
      // The Hubble factor is because the units in L-Galaxies do not include them
      tauCoolP = 20.*GRAVITY*pow2(mumH*UnitLength_in_cm*Rvir)/(Hubble_h*LAMBDA_RATIO*lambda*UnitTime_in_s);
      tauCoolP *= pow2(f)/pow3(X);
      fg0 = HotGas/Mvir;
      dt_ratio=dt/tdyn_halo;
      tau_ratio=tdyn_halo*fg0/tauCoolP;
      if (tau_ratio <=1)
	  fg=fg0/(1+tau_ratio*dt_ratio);
      else {
	  teq_ratio=log(tau_ratio);
	  if (dt_ratio <= teq_ratio)
	      fg=fg0*exp(-dt_ratio);
	  else
	      fg=fg0/(tau_ratio*(1+(dt_ratio-teq_ratio)));
      }
      coolingGas = (fg0-fg)*Mvir;
      //printf("coolingGas=%g\n",coolingGas);
      // Use instantaneous X-ray luminosity, integrated out to infinity eq 23 of BTH18
      // Use Hot gas fraction before cooling as that corresponds more closely to emitted X-rays
      //Gal[p].XrayLum = log10(LAMBDA_RATIO*lambda*pow2(HotGas*UnitMass_in_g/(Hubble_h*4.*f*mumH))/pow3(a*UnitLength_in_cm/Hubble_h));
      // Use time-averaged X-ray luminosity, Mdot * H, eq 24 of BTH18; Hubble factors cancel here.
      Gal[p].XrayLum = log10(2.5*BOLTZMANN*temp*coolingGas*UnitMass_in_g/(mumH*dt*UnitTime_in_s));
#else
      /* An isothermal density profile for the hot gas is assumed here.
       * The density is rhorsq0 /r^2. */
      rhorsq0 = HotGas / (4 * M_PI * HotRadius);
      /* eq. 3 and 4 Guo2010 or  eq 3-5 of Hen15.  Note that these do not agree with each other
       * or with the expression implemented below! */
      x = lambda / (PROTONMASS * BOLTZMANN * temp); // units cm^3 / g s
      x *= (UnitDensity_in_cgs * UnitTime_in_s);  // now in internal units (apart from h factor).
      /* NOTE: because code units are actually UnitMass, UnitLength & UnitTime divided by h,
       * the above expression for x is actually missing a factor h. */
      rcool = sqrt(MU_RATIO * rhorsq0 * x * tdyn_halo);
    
      // This presumably keeps track of the largest value of the CoolingRadius that this halo has reached.
      // It seems an odd thing to track.  As far as I can tell, it is never used but is written out as a
      // diagnistic after dividing by sqrt(Hubble_h).
      if (Gal[p].CoolingRadius < rcool)
	  Gal[p].CoolingRadius = rcool;
      
      //if Hotradius is used, when galaxies become type 1's there will be a discontinuity in the cooling
      if(rcool > Rvir) // INFALL DOMINATED REGIME
	  //coolingGas = HotGas; - Delucia 2007
	  /*comes in to keep the continuity (Delucia2004) */
	  //put h
	  coolingGas = HotGas * Vvir * dt / HotRadius;
      else // HOT PHASE REGIME -> TODO - WHERE IS THE 0.5 factor of eq 5
	  /*coolingGas = (HotGas / Rvir) * (rcool / tdyn_halo) * dt * 0.5; */
	  coolingGas = (HotGas / HotRadius) * (rcool / tdyn_halo) * dt ;
#endif //BETAPROF
        
      // Photoionizing background
      if (log10(temp) < 4.0)
	  coolingGas = 0.;
      // Sanity checks
      if(coolingGas > HotGas)
	  coolingGas = HotGas;
      else if(coolingGas < 0.0)
	  coolingGas = 0.0;      
    
  }
  else {
      coolingGas = 0.0;
  }

  Gal[p].CoolingGas = coolingGas;

#ifndef BETAPROF
  /* Calculate X-ray luminosity here before adding in heating.  Basically undoes cooling calculation */
  /* Need to take the log or else get an overflow in cgs units on output. */
  Gal[p].XrayLum = log10(1.5*(coolingGas*UnitMass_in_g)/(dt*UnitTime_in_s)*temp*BOLTZMANN/(MU*PROTONMASS)) ;
#endif

  mass_checks(p,"cooling.c",__LINE__);

}   


#ifdef AGN_FEEDBACK

/** @brief calculates the energy released by black holes due to passive accretion,
 * that will be used to reduced the cooling.*/


/** @brief do_AGN_heating calculates the amount of energy released by
 * black holes due to passive accretion, Which is then used to reduce
 * the cooling.
 * This should be extended to include all forms of cooling.
 * Parameters:
 *    + AGN_disk_epsilon - fraction of mass flowing through the accretion disk that is converted to energy.
 *    + AGN_disk_feedback_efficiency - fraction of energy that goes into feedback, subject to:
 *    + AGN_disk_eddington_limit - upper limit on fraction of energy going into feedback.
 *    + AGN_accretion_mode - 0: none (i.e. Croton); 1: hot; 2: cooling; 3: cold. 
 *    + AGN_accretion_efficiency - normalisation of AGN accretion rate
 *    + AGN_accretion_power - BH power in AGN accretion rate formula (normalised to 1e8 Msun/h)
 *    + AGN_tau_accrete - accretion rate limiter (multiple of appropriate dynamical time)
 *    + AGN_tau_reheat - cold gas reheating limiter
 *    + AGN_tau_halo - hot gas ejection limiter
 */

#include <stdbool.h>

void do_AGN_heating(double dt, int ngal)
{
    bool inFoFHalo;
    double AGNaccreted, AGNcoeff, fraction, EddRate;
    double dist, HotGas, HotRadius;
    double CoolingGas;
    int p, p_CentralGal, FoFCentralGal;

    double Vmax, Vmax_halo, Vmax_FoFhalo, alpha, H_disk, R_disk, c_s, AccRateQ, AccRateR, AccRate;
    double MaxFactor;
    double Lum,LEdd,Lum_feedback,E_feedback,Reheat,Tdyn_halo,Tdyn_FoFhalo,t_disk;
    double AGN_disk_epsilon, epsilonBH, epsilon_feedback, AGN_disk_eddington_limit;
    double AGN_tau_accrete, tau_reheat, tau_halo, tau_FoFhalo;
    double c2=pow2(C/UnitVelocity_in_cm_per_s); // c^2 - Speed of light squared in code units

    // Effective feedback efficiency measured in terms of mass accreted onto BH.
    epsilonBH=AGN_disk_epsilon/(1.-AGN_disk_epsilon);

    for (p = 0; p < ngal; p++)
	if(Gal[p].Type == 0)
	    FoFCentralGal=p;

    for (p = 0; p < ngal; p++) {
	Gal[p].CoolingRate_beforeAGN += Gal[p].CoolingGas / (dt*STEPS); // averaged over all steps
	HotGas = Gal[p].HotGas;
	HotRadius = Gal[p].HotRadius;
	CoolingGas = Gal[p].CoolingGas;

	// Determine some halo properties that we will use later
	// Properties of this halo
	// I think that Vmax gives the best measure of the current halo escape velocity [PAT]?
	Vmax = Gal[p].Vmax;
	R_disk = Gal[p].GasDiskRadius/3.;  // exponential scale length of disk
	t_disk = R_disk/Vmax;
	// Dynamical time of halo containing this galaxy
	p_CentralGal = Gal[p].CentralGal;
	Tdyn_halo = Gal[p_CentralGal].Rvir/Gal[p_CentralGal].Vvir;
	Vmax_halo = Gal[p_CentralGal].Vmax;
	// Dynamical time of FOF halo
	Tdyn_FoFhalo = Gal[FoFCentralGal].Rvir/Gal[FoFCentralGal].Vvir;
	Vmax_FoFhalo = Gal[FoFCentralGal].Vmax;
	// Flag whether galaxy is in main FOF halo or not.
	// Distance of galaxy at centre of this halo to that at the centre of the FOF group:
	dist=separation_gal(FoFCentralGal,Gal[p].CentralGal)/(1+ZZ[Halo[Gal[FoFCentralGal].HaloNr].SnapNum]);
	if (Gal[p].Type==1 && dist < Gal[FoFCentralGal].Rvir)
	    inFoFHalo = true;
	else
	    inFoFHalo=false;

	// Quasar mode accretion rate
	AccRateQ = Gal[p].BlackHoleAccrete/dt;

	// Radio mode accretion rate
	if (AGN_accretion_mode == 0) {
	    terminate("Croton BH radio mode not yet implemented\n");
	}
	else if (AGN_accretion_mode == 1) {
	    AccRateR = AGN_accretion_efficiency * pow(Gal[p].BlackHoleMass/0.01,AGN_accretion_power) 
		* Gal[p].HotGas / Tdyn_halo;
	    AGNaccreted=min( AccRateR*dt, Gal[p].HotGas*(1-exp(-dt/(AGN_tau_accrete*Tdyn_halo))) );
	}
	else if (AGN_accretion_mode == 2) {
	    AccRateR = AGN_accretion_efficiency * pow(Gal[p].BlackHoleMass/0.01,AGN_accretion_power) 
		* Gal[p].CoolingGas / Tdyn_halo;
	    AGNaccreted=min( AccRateR*dt, Gal[p].CoolingGas*(1-exp(-dt/(AGN_tau_accrete*Tdyn_halo))) );
	}
	else if (AGN_accretion_mode == 3) {
	    H_disk = 250.e-6*Hubble_h; // Half-thickness of disk is 250 pc
	    c_s = 10.;  // Sound speed / vel disp. of disk in km/s
#ifdef BHACCVMAX
	    AccRateR = AGN_accretion_efficiency * pow(G*Gal[p].BlackHoleMass/(R_disk*Vmax),2)/(H_disk*Vmax) * Gal[p].ColdGas;
#else
	    AccRateR = AGN_accretion_efficiency * pow(G*Gal[p].BlackHoleMass/(R_disk*c_s),2)/(H_disk*c_s) * Gal[p].ColdGas;
#endif
	    AGNaccreted=min( AccRateR*dt, Gal[p].ColdGas*(1-exp(-dt/(AGN_tau_accrete*t_disk))) );
	}
	else 
	    terminate("Invalid value for AGN_accretion_mode\n");
#ifdef DEBUG
	if (Gal[p].Type==0 && Gal[p].BlackHoleMass > 1 && Gal[p].BlackHoleMass < 1.01) 
	    printf("BHMass, ColdGas, HotGas, MaxFactorAccrete = %12.3g, %12.3g, %12.3g, %12.3g\n",Gal[p].BlackHoleMass,Gal[p].ColdGas,Gal[p].HotGas,1-exp(-dt/(AGN_tau_accrete*t_disk)));
#endif

	// Add in accretion from quasar mode
	AGNaccreted = min(Gal[p].ColdGas, AGNaccreted + Gal[p].BlackHoleAccrete);
	AccRate=AGNaccreted/dt;
	// Limit to Eddington rate (but only once the BH exists) - Replace with seed!
	// Eddington luminosity in code units:
	LEdd = 1.3e48 * Gal[p].BlackHoleMass / (UnitEnergy_in_cgs / UnitTime_in_s);
	EddRate = LEdd/(epsilonBH*c2);
	if (Gal[p].BlackHoleMass > TINY_MASS)
	    AGNaccreted = min(AccRate,EddRate)*dt;
	else
	    AGNaccreted = AccRate*dt;
	AccRate=AGNaccreted/dt;

#ifdef DEBUG
	if (AGNaccreted > 1e10) {
	    //printf("BlackHoleMass = %g\n",Gal[p].BlackHoleMass);
	    //printf("EDDrate = %g\n",EDDrate);
	    //printf("AGN_accretion_efficiency = %g\n",AGN_accretion_efficiency);
	    //printf("ColdGas = %g\n",Gal[p].ColdGas);
	    //printf("AccRateR = %g\n",AccRateR);
	    printf("%12.3g %12.3g %12.3g %12.3g\n",AGNaccreted,Gal[p].BlackHoleMass,Gal[p].ColdGas,Gal[p].BlackHoleAccrete);
	    //printf("UnitTime_in_yr=%12.3g\n",UnitTime_in_s/3.16e7);
	    //exit(1);
	}
#endif

	Gal[p].LumAGN = 0.;
	if (AGNaccreted > 0.) { //------------------------------------------------------------------

	// Do Mass and metal transfers
   	Gal[p].BlackHoleMass += AGNaccreted;
	// Remove metals from ColdGas.  This does really work, honest.
	if (Gal[p].ColdGas > 0.)
	    fraction=AGNaccreted/Gal[p].ColdGas;
	else
	    fraction = 0.;
	Gal[p].MetalsColdGas = metals_add(Gal[p].MetalsColdGas, Gal[p].MetalsColdGas,-fraction);
	Gal[p].ColdGas -= AGNaccreted;
	
	// Set AGN power available for feedback and hence mass that can be reheated.
	Lum = epsilonBH*c2*AccRate;
	// Some fraction of Lum goes into radiated AGN power and some fraction AGN_disk_feedback_efficiency goes into feedback.
	// The fraction that goes into feedback is limited to AGN_disk_eddington_limit times the Eddington rate.
	Lum_feedback = min(AGN_disk_feedback_efficiency*Lum, AGN_disk_eddington_limit*LEdd);
	E_feedback = Lum_feedback*dt;
	Gal[p].LumAGN = Lum-Lum_feedback;

	// First reheat cold gas, up to the maximum allowed
	if (Gal[p].ColdGas > 0.) {
	    MaxFactor = 1.-exp(-dt/(AGN_tau_reheat*t_disk));
	    Reheat=min(E_feedback/(0.5*Vmax*Vmax),Gal[p].ColdGas*MaxFactor);
	    E_feedback-=Reheat*(0.5*Vmax*Vmax);
	    fraction=Reheat/Gal[p].ColdGas;
	    if((Gal[p].Type == 0) || (Gal[p].Type == 1)) 
		transfer_gas(p,"Hot",p,"Cold",fraction,"do_AGN_heating", __LINE__);
	    else
		transfer_gas(p_CentralGal,"Hot",p,"Cold",fraction,"do_AGN_heating", __LINE__);
	}

	// Next use any residual heating to suppress cooling, if any.
	if (E_feedback > 0. && Gal[p_CentralGal].CoolingGas > 0.) {
	    Reheat=min(E_feedback/(0.5*Vmax_halo*Vmax_halo),Gal[p_CentralGal].CoolingGas);
	    E_feedback-=Reheat*0.5*Vmax_halo*Vmax_halo;
	    Gal[p_CentralGal].CoolingGas-=Reheat;
	}
	
	// Eject gas from your own halo
	if (E_feedback >0. && Gal[p_CentralGal].HotGas > 0.) {
	    MaxFactor = 1.-exp(-dt/(AGN_tau_halo*Tdyn_halo));
	    Reheat=min(E_feedback/(0.5*Vmax_halo*Vmax_halo),Gal[p_CentralGal].HotGas*MaxFactor);
	    E_feedback-=Reheat*0.5*Vmax_halo*Vmax_halo;
	    fraction=Reheat/Gal[p_CentralGal].HotGas;
	    transfer_gas(p_CentralGal,"AGNfeedback",p_CentralGal,"Hot",fraction,"do_AGN_heating", __LINE__);
	}

	// If a Type 1 within a Type 0 halo, can also eject gas from the Type 0 halo
	if (inFoFHalo && E_feedback >0. && Gal[FoFCentralGal].HotGas > 0.) {
	    MaxFactor = 1.-exp(-dt/(AGN_tau_halo*Tdyn_FoFhalo));
	    Reheat=min(E_feedback/(0.5*Vmax_FoFhalo*Vmax_FoFhalo),Gal[FoFCentralGal].HotGas*MaxFactor);
	    E_feedback-=Reheat*0.5*Vmax_FoFhalo*Vmax_FoFhalo;
	    fraction=Reheat/Gal[FoFCentralGal].HotGas;
	    transfer_gas(FoFCentralGal,"AGNfeedback",FoFCentralGal,"Hot",fraction,"do_AGN_heating", __LINE__);
	}

	if(Gal[p].CoolingGas < 0.0)
	    Gal[p].CoolingGas = 0.0;
	if(Gal[p].HotGas < 0.0)
	    Gal[p].HotGas = 0.0;
	if(Gal[p].ColdGas < 0.0)
	    Gal[p].ColdGas = 0.0;
	if(Gal[p_CentralGal].CoolingGas < 0.0)
	    Gal[p_CentralGal].CoolingGas = 0.0;
	if(Gal[FoFCentralGal].CoolingGas < 0.0)
	    Gal[FoFCentralGal].CoolingGas = 0.0;
	if(Gal[p_CentralGal].HotGas < 0.0)
	    Gal[p_CentralGal].HotGas = 0.0;
	if(Gal[FoFCentralGal].HotGas < 0.0)
	    Gal[FoFCentralGal].HotGas = 0.0;

	} //---------------------------------------------------------------------------------------

	mass_checks("cooling_recipe #2.",p);
    }
}

#else // AGN_FEEDBACK

/** @brief calculates the energy released by black holes due to passive accretion,
  * that will be used to reduced the cooling.*/


/** @brief do_AGN_heating calculates the amount of energy released by
  * black holes due to passive accretion, Which is then used to reduce
  * the cooling.
  *
  * There is one parameter, AgnEfficiency, which is the efficiency of AGN
  * passive accretion and consequently of cooling flow reheating, Eqs. 10,
  * 11 & 12 in Croton 2006. There is a AGNRadioModeModel =2 (empirical)
  * with options 3 and 4 representing Bondi-Hoyle and cold cloud accretion.
  * The three should be identical and the use of empirical avoids people
  * shouting about duty-cycles being inconsistent with Bondi-Hoyle & etc.
  */

// AGN accretion and heating is assumed to go from the gas from the type 0
// galaxy at the centre of the FOF group to the most massive black hole
// inside Rvir. If the most massive black hole is in a type 1 this produces
// all the heating affecting the type 0 galaxy and will add to the accretion
// occuring on the type 1
// this is done to account for the fact that the centres of FOF groups can
// switch between galaxies. As a result a very small galaxy without a black hole
// might be assigned the centre of a cluster leading to huge cooling. It is therefore
// not necessary to do the same correction for satellites of subhalos.


void do_AGN_heating(double dt, int ngal, int FOF_centralgal)
{
  double AGNrate;  // The mass accretion rate onto the BH in code units
  double AGNheating, AGNaccreted, AGNcoeff, fraction, EDDrate, FreeFallRadius;
  double dist, HotGas, HotRadius, Rvir, Vvir, Mvir;
  double LeftOverEnergy, CoolingGas;
  int p;

  for (p = 0; p < ngal; p++) {
      Gal[p].CoolingRate_beforeAGN += Gal[p].CoolingGas / (dt*STEPS);
      AGNrate=0.;
      LeftOverEnergy = 0.;
      HotGas = Gal[p].HotGas;
      HotRadius = Gal[p].HotRadius;
      CoolingGas = Gal[p].CoolingGas;
      Mvir = Gal[p].Mvir;
      Rvir = Gal[p].Rvir;
      Vvir = Gal[p].Vvir;
      
      mass_checks(p,"cooling.c",__LINE__);
      
      if(HotGas > 0.0) {
	  
	  if(AGNRadioModeModel == 0) {
	      // Henriques 2015.  Proportional to MBH * MHotgas.  The factor of 10 here is just from
	      // the way things are normalised in the paper.
	      AGNrate = AgnEfficiency * (UnitTime_in_s*SOLAR_MASS)/(UnitMass_in_g*SEC_PER_YEAR)
		  * (Gal[p].BlackHoleMass/Hubble_h) * (HotGas/Hubble_h) * 10.;
	  }
	  
	  else if(AGNRadioModeModel == 1) {
	      //empirical (standard) accretion recipe - Eq. 10 in Croton 2006
	      AGNrate = AgnEfficiency / (UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS)
		  * (Gal[p].BlackHoleMass / 0.01) * pow3(Vvir / 200.0)
		  * ((HotGas / HotRadius * Rvir / Mvir) / 0.1);
	  }
	  
	  else if(AGNRadioModeModel == 2 || AGNRadioModeModel == 3) {
	      double x, lambda, temp, logZ, MetalsHotGas;
	      const double MU=0.59265;
	      //MU_RATIO is 2*(n_e/n)*(n_i/n)/(3*mu)
	      const double MU_RATIO=0.28086;
	      double mumH=MU*PROTONMASS;

	      // The following repeats a calculation in compute_cooling().  Surely we can be more efficient?
	      temp = 0.5*mumH*UnitVelocity_in_cm_per_s*UnitVelocity_in_cm_per_s/BOLTZMANN * Vvir*Vvir;
	      MetalsHotGas = metals_total(Gal[p].MetalsHotGas);
	      if(MetalsHotGas > 0)
		logZ = log10(MetalsHotGas / HotGas);
	      else
		  logZ = -10.0;
	      lambda = get_metaldependent_cooling_rate(log10(temp), logZ);
	      x = PROTONMASS * BOLTZMANN * temp / lambda; // now this has units sec g/cm^3
	      x /= (UnitDensity_in_cgs * UnitTime_in_s);  // now in internal units

	      /* Bondi-Hoyle accretion recipe -- efficiency = 0.15
	       * Eq. 29 in Croton 2006 */
	      // Actually this does NOT agree with that paper - no need for magic numbers.  PAT.
	      if(AGNRadioModeModel == 2) {
		  AGNrate = (2.5 * M_PI * G) * (0.75 * 0.6 * x) * Gal[p].BlackHoleMass * 0.15;
	      }
	      /* Cold cloud accretion recipe -- trigger: Rff = 50 Rdisk,
	       * and accretion rate = 0.01% cooling rate.  Eq. 25 in Croton 2006 */
	      else if(AGNRadioModeModel == 3) {
		  FreeFallRadius = HotGas / (6.0 * 0.6 * x * Rvir * Vvir) / HotRadius * Rvir;
		  if(Gal[p].BlackHoleMass > 0.0 && FreeFallRadius < Gal[p].ColdGasRadius * 50.0)
		    AGNrate = 0.0001 * CoolingGas / dt;
		  else
		    AGNrate = 0.0;
	      }
	  }

	  /* Eddington accretion rate */
	  /* Note that this assumes an efficiency of 50%
	   * - it ignores the e/(1-e) factor in L = e/(1-e) Mdot c^2 */
	  EDDrate = 1.3e48 * Gal[p].BlackHoleMass / (UnitEnergy_in_cgs / UnitTime_in_s) / 9e10;

	  /* accretion onto BH is always limited by the Eddington rate */
	  if(AGNrate > EDDrate) AGNrate = EDDrate;

	  /*  accreted mass onto black hole the value of dt puts an h factor into AGNaccreted as required for code units */
	  AGNaccreted = AGNrate * dt;

	  /* cannot accrete more mass than is available! */
	  if(AGNaccreted > HotGas) AGNaccreted = HotGas;

	  /*  coefficient to heat the cooling gas back to the virial temperature of the halo */
	  /*  1.34e5 = sqrt(2*eta*c^2), eta=0.1 (standard efficiency) and c in km/s
	   *  Eqs. 11 & 12 in Croton 2006 */
	  AGNcoeff = (1.34e5 / Vvir) * (1.34e5 / Vvir);

	  /*  cooling mass that can be suppressed from AGN heating */
	  AGNheating = AGNcoeff * AGNaccreted;

	  if(AGNRadioModeModel == 0 && Gal[p].Type==1) {
	      dist=separation_gal(p,FOF_centralgal);
	      if(dist < Gal[FOF_centralgal].Rvir) {
		  if(AGNheating > (Gal[p].CoolingGas + Gal[FOF_centralgal].CoolingGas)) {
		      AGNheating = (Gal[p].CoolingGas + Gal[FOF_centralgal].CoolingGas);
		      AGNaccreted = (Gal[p].CoolingGas + Gal[FOF_centralgal].CoolingGas) / AGNcoeff;
		  }
		  if(AGNheating > Gal[p].CoolingGas)
		      LeftOverEnergy = AGNheating - Gal[p].CoolingGas;
	      }
	  }
	  else
	      if(AGNheating > Gal[p].CoolingGas)
		  AGNaccreted = Gal[p].CoolingGas / AGNcoeff;

	  /* limit heating to cooling rate */
	  if(AGNheating > Gal[p].CoolingGas)
	      AGNheating = Gal[p].CoolingGas;

	  mass_checks(p,"cooling.c",__LINE__);

	  /*  accreted mass onto black hole */
	  Gal[p].BlackHoleMass += AGNaccreted; //ROB: transfer_mass functions should be used here
	  Gal[p].RadioAccretionRate += AGNaccreted / (dt*STEPS);
	  fraction=AGNaccreted/Gal[p].HotGas;
	  Gal[p].HotGas -= AGNaccreted;
	  Gal[p].MetalsHotGas = metals_add(Gal[p].MetalsHotGas,Gal[p].MetalsHotGas, -fraction);
	  mass_checks(p,"cooling.c",__LINE__);
#ifdef INDIVIDUAL_ELEMENTS
	  int kk;
	  for(kk=0;kk<NUM_ELEMENTS;kk++)
	    Gal[p].HotGas_elements[kk] *= (1-fraction);
#endif
#ifdef METALS_SELF
	  Gal[p].MetalsHotGasSelf = 	metals_add(Gal[p].MetalsHotGasSelf,Gal[p].MetalsHotGasSelf,-fraction);
#endif	

      }
      else
	  AGNheating = 0.0;

      mass_checks(p,"cooling.c",__LINE__);

      Gal[p].CoolingGas -= AGNheating;
      if(Gal[p].CoolingGas < 0.0) Gal[p].CoolingGas = 0.0;
      Gal[p].CoolingRate += Gal[p].CoolingGas / (dt*STEPS);

      if(AGNRadioModeModel == 0 && LeftOverEnergy>0.) {
	  Gal[FOF_centralgal].CoolingGas -= LeftOverEnergy;
	  if(Gal[FOF_centralgal].CoolingGas < 0.0)
	      Gal[FOF_centralgal].CoolingGas = 0.0;
	  else
	      Gal[FOF_centralgal].CoolingRate -= LeftOverEnergy / (dt*STEPS);
	}

      mass_checks(p,"cooling.c",__LINE__);

  }
}

#endif //AGN_FEEDBACK

/** @brief updates the fractions of hot and cold gas due to cooling. */
 /** @brief cool_gas_onto_galaxy updates the fractions of hot and cold gas
    * due to cooling. This is done for the mass, metals and, after Guo2010,
    * spin components */

void cool_gas_onto_galaxy(int p, double dt)
{
  double fraction, Mdisk, Mcool;
#ifdef H2_AND_RINGS
  double ringtot, rd, fractionRings[RNUM];
  int jj;
#endif

  /** @brief cool_gas_onto_galaxy updates the fractions of hot and cold gas
    * due to cooling. This is done for the mass, metals and, after Guo2010,
    * spin components */
  Mdisk=Gal[p].ColdGas;


  Mcool=Gal[p].CoolingGas;
  if (Mcool>Gal[p].HotGas)
    Mcool = Gal[p].HotGas;

  /*  add the fraction 1/STEPS of the total cooling gas to the cold disk */
  if(Mcool > 0.0)
    {
      fraction=((float)Mcool)/Gal[p].HotGas;

      double new_radius, halospinpar,coldspinpar;
      /*halospinpar=sqrt(Halo[Gal[p].HaloNr].Spin[0] * Halo[Gal[p].HaloNr].Spin[0] +
		       Halo[Gal[p].HaloNr].Spin[1] * Halo[Gal[p].HaloNr].Spin[1] +
		       Halo[Gal[p].HaloNr].Spin[2] * Halo[Gal[p].HaloNr].Spin[2] );*/


#ifdef H2_AND_RINGS
      rd=get_initial_disk_radius(Gal[p].HaloNr, p)/3.;
      ringtot=1.-(1+RingRadius[RNUM-1]/rd)/exp(RingRadius[RNUM-1]/rd);

      //normalized fraction of mass in each ring * fraction of mass to transfer = fraction in each ring
      fractionRings[0]=((1-(1+RingRadius[0]/rd)/exp(RingRadius[0]/rd))/ringtot)*fraction;
      for(jj=1; jj<RNUM; jj++)
     	fractionRings[jj]= (((1+RingRadius[jj-1]/rd)/exp(RingRadius[jj-1]/rd)-(1+RingRadius[jj]/rd)/exp(RingRadius[jj]/rd))/ringtot)*fraction;

      transfer_material_with_rings(p,"ColdGas",p,"HotGas", fractionRings,"cool_gas_onto_galaxy", __LINE__);

      /*if(Gal[p].ReheatedGas>0.)
	{
	  double Vmax;
	  if (Gal[p].Type == 0)
	    Vmax=Gal[p].Vmax;
	  else
	    Vmax=Gal[p].InfallVmax;

	  Mcool = Gal[p].ReheatedGas / (Gal[p].ColdGasRadius / Vmax) * dt;
	  //Mcool = Gal[p].ReheatedGas / (Gal[p].Rvir / Vmax) * dt;


	  if(Mcool>Gal[p].ReheatedGas)
	    Mcool=Gal[p].ReheatedGas;
	  fraction=Mcool/Gal[p].ReheatedGas;

	  //rd=get_initial_disk_radius(Gal[p].HaloNr, p)/3.;
	  rd=Gal[p].ColdGasRadius/3.;
	  ringtot=1.-(1+RingRadius[RNUM-1]/rd)/exp(RingRadius[RNUM-1]/rd);
	  //normalized fraction of mass in each ring * fraction of mass to transfer = fraction in each ring

	  fractionRings[0]=((1-(1+RingRadius[0]/rd)/exp(RingRadius[0]/rd))/ringtot)*fraction;
	  for(jj=1; jj<RNUM; jj++)
	    fractionRings[jj]= (((1+RingRadius[jj-1]/rd)/exp(RingRadius[jj-1]/rd)-(1+RingRadius[jj]/rd)/exp(RingRadius[jj]/rd))/ringtot)*fraction;

	  transfer_material_with_rings(p,"ColdGas",p,"ReheatedGas", fractionRings,"cool_gas_onto_galaxy", __LINE__);
	}*/

#else
      transfer_material(p,"ColdGas",p,"HotGas",fraction,"cool_gas_onto_galaxy", __LINE__);

      /*if(Gal[p].ReheatedGas>0.)
    	{
    	  double Vmax;
    	  if (Gal[p].Type == 0)
    	    Vmax=Gal[p].Vmax;
    	  else
    	    Vmax=Gal[p].InfallVmax;

    	  Mcool = Gal[p].ReheatedGas / (Gal[p].DiskRadius / Vmax) * dt;
    	  //Mcool = Gal[p].ReheatedGas / (Gal[p].Rvir / Vmax) * dt;

    	  if(Mcool>Gal[p].ReheatedGas);
    	    Mcool=Gal[p].ReheatedGas;
    	  fraction=Mcool/Gal[p].ReheatedGas;
    	  transfer_material(p,"ColdGas",p,"ReheatedGas",fraction,"cool_gas_onto_galaxy", __LINE__);
    	}*/
#endif

      if (DiskRadiusModel == 0)
    	{
	  int ii;
	  if (Gal[p].ColdGas > 0.0)
    	    {
    	      for (ii=0;ii<3;ii++)
    		Gal[p].ColdGasSpin[ii]=(Gal[p].ColdGasSpin[ii]*Mdisk+Halo[Gal[p].HaloNr].Spin[ii]*Mcool)/(Gal[p].ColdGas);
	      //Gal[p].ColdGasSpin[ii]=(Gal[p].ColdGasSpin[ii]*Mdisk+1./sqrt(3)*halospinpar*Mcool)/(Gal[p].ColdGas);
    	      Gal[p].ColdGasRadius=get_gas_disk_radius(p);
    	    }
    	   else
    	     {
    	       for (ii=0; ii<3; ii++)
    		 Gal[p].ColdGasSpin[ii]=0.0;
    	       Gal[p].ColdGasRadius = 0.;
    	     }
    	}

    }

}

