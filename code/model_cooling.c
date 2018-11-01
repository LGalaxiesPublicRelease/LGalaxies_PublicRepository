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

/** @brief main cooling recipe, where the cooling rates are calculated */

void compute_cooling(int p, double dt, int ngal)
{
  double Vvir, Rvir, x, lambda, tcool, rcool, temp, tot_hotMass, tot_metals=0., HotRadius;
  double coolingGas, logZ, rho_rcool, rho0;

  mass_checks(p,"cooling.c",__LINE__);

  tot_hotMass = Gal[p].HotGas;
  int ii;
  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
    tot_metals += Gal[p].MetalsHotGas[ii];

  if(tot_hotMass > 1.0e-6)
  {
  	/* TODO - Should Rvir be used at all in this recipe after Guo2010?
     * probably always HotRadius*/
    Vvir = Gal[p].Vvir;
    Rvir = Gal[p].Rvir;
    
      
    tcool = Rvir / Vvir; // tcool = t_dynamical = Rvir/Vvir


    /* temp -> Temperature of the Gas in Kelvin, obtained from
     * hidrostatic equilibrium KT=0.5*mu_p*(Vc)^2 assuming Vvir~Vc */
    temp = 35.9 * Vvir * Vvir;
      
    if (Gal[p].Type == 0)
      HotRadius = Gal[p].Rvir;
    else 
      HotRadius = Gal[p].HotRadius;
    
    if(tot_metals > 0)
      logZ = log10(tot_metals / tot_hotMass);
    else
      logZ = -10.0;

    //eq. 3 and 4 Guo2010
    lambda = get_metaldependent_cooling_rate(log10(temp), logZ);
    x = PROTONMASS * BOLTZMANN * temp / lambda; // now this has units sec g/cm^3
    x /= (UnitDensity_in_cgs * UnitTime_in_s);  // now in internal units
    rho_rcool = x / (0.28086 * tcool);
    /* an isothermal density profile for the hot gas is assumed here */
    rho0 = tot_hotMass / (4 * M_PI * HotRadius);
    rcool = sqrt(rho0 / rho_rcool);
    
    if (Gal[p].CoolingRadius < rcool)
      Gal[p].CoolingRadius = rcool;
      
    //if Hotradius is used, when galaxies become type 1's there will be a discontinuity in the cooling
    if(rcool > Rvir) // INFALL DOMINATED REGIME
      //coolingGas = tot_hotMass; - Delucia 2007
      /*comes in to keep the continuity (Delucia2004) */
      //put h
      coolingGas = tot_hotMass / (HotRadius / Vvir) * dt;
    else // HOT PHASE REGIME -> TODO - WHERE IS THE 0.5 factor of eq 5
      /*coolingGas = (tot_hotMass / Rvir) * (rcool / tcool) * dt * 0.5; */
      coolingGas = (tot_hotMass / HotRadius) * (rcool / tcool) * dt ;

    //Photoionizing background
    if (log10(temp) < 4.0)
      coolingGas = 0.;
    
    if(coolingGas > tot_hotMass)
      coolingGas = tot_hotMass;
    else if(coolingGas < 0.0)
      coolingGas = 0.0;      
    

  }
  else
  {
    coolingGas = 0.0;
  }

  Gal[p].CoolingGas = coolingGas;

  mass_checks(p,"cooling.c",__LINE__);

}   



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
//occuring on the type 1
// this is done to account for the fact that the centres of FOF groups can
//switch between galaxies. As a result a very small galaxy without a black hole
//might be assigned the centre of a cluster leading to huge cooling. It is therefore
//not necessary to do the same correction for satellites of subhalos.



void do_AGN_heating(double dt, int ngal, int FOF_centralgal)
{
  double AGNrate, AGNheating, AGNaccreted, AGNcoeff, fraction, EDDrate, FreeFallRadius;
  double dist, HotGas, HotRadius, Rvir, Vvir, Mvir;
  double LeftOverEnergy, CoolingGas;
  int p, ii;

  for (p = 0; p < ngal; p++)
    {
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

      if(HotGas > 0.0)
	{
	  if(AGNRadioModeModel == 0)
	    // Henriques 2015.  Proportional to MBH * MHotgas
	    {
	     // if(ZZ[Gal[p].SnapNum]<1.)
		  AGNrate =	AgnEfficiency * (UnitTime_in_s*SOLAR_MASS)/(UnitMass_in_g*SEC_PER_YEAR)
		    * Gal[p].BlackHoleMass/Hubble_h* (HotGas/Hubble_h) * 10.;
		//  AGNrate =	AgnEfficiency * (UnitTime_in_s*SOLAR_MASS)/(UnitMass_in_g*SEC_PER_YEAR)
		// 		    * Gal[p].BlackHoleMass/Hubble_h* 5000.;


	      ///((1+ZZ[Gal[p].SnapNum])*(1+ZZ[Gal[p].SnapNum])*(1+ZZ[Gal[p].SnapNum]));
	      /*AGNrate = AgnEfficiency * (UnitTime_in_s*SOLAR_MASS)/(UnitMass_in_g*SEC_PER_YEAR)
	       *          *pow(Gal[p].BlackHoleMass,0.5)*(Gal[p].HotGas * Gal[p].Rvir/Gal[p].HotRadius)*0.1;*/
	      /*else
		AGNrate =	0.;*/
	    }
	  else if(AGNRadioModeModel == 1)
	    {
	      //empirical (standard) accretion recipe - Eq. 10 in Croton 2006
	      AGNrate = AgnEfficiency / (UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS)
  		        * (Gal[p].BlackHoleMass / 0.01) * pow3(Vvir / 200.0)
			* ((HotGas / HotRadius * Rvir / Mvir) / 0.1);
	    }
	  else if(AGNRadioModeModel == 2 || AGNRadioModeModel == 3)
	    {
	      double x, lambda, temp, logZ, tot_metals=0.;

	      int ii;
	      for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	         tot_metals += Gal[p].MetalsHotGas[ii];

	      /* temp -> Temperature of the Gas in Kelvin, obtained from
	       * hidrostatic equilibrium KT=0.5*mu_p*(Vc)^2 assuming Vvir~Vc */
	      temp = 35.9 * Vvir * Vvir;
	      if(tot_metals > 0)
		logZ = log10(tot_metals / HotGas);
	      else
  				logZ = -10.0;
	      lambda = get_metaldependent_cooling_rate(log10(temp), logZ);
	      x = PROTONMASS * BOLTZMANN * temp / lambda; // now this has units sec g/cm^3
	      x /= (UnitDensity_in_cgs * UnitTime_in_s);  // now in internal units

	      /* Bondi-Hoyle accretion recipe -- efficiency = 0.15
	       * Eq. 29 in Croton 2006 */
	      if(AGNRadioModeModel == 2)
		{
		  AGNrate = (2.5 * M_PI * G) * (0.75 * 0.6 * x) * Gal[p].BlackHoleMass * 0.15;
		}
	      else if(AGNRadioModeModel == 3)
		{
		  /* Cold cloud accretion recipe -- trigger: Rff = 50 Rdisk,
		   * and accretion rate = 0.01% cooling rate
		   * Eq. 25 in Croton 2006 */
		  FreeFallRadius = HotGas / (6.0 * 0.6 * x * Rvir * Vvir) /	HotRadius * Rvir;
		  if(Gal[p].BlackHoleMass > 0.0 && FreeFallRadius < Gal[p].ColdGasRadius * 50.0)
		    AGNrate = 0.0001 * CoolingGas / dt;
		  else
		    AGNrate = 0.0;
		}
	    }

	  /* Eddington rate */
	  /* Note that this assumes an efficiency of 50%
	   * - it ignores the e/(1-e) factor in L = e/(1-e) Mdot c^2 */
	  EDDrate = 1.3e48 * Gal[p].BlackHoleMass / (UnitEnergy_in_cgs / UnitTime_in_s) / 9e10;

	  /* accretion onto BH is always limited by the Eddington rate */
	  if(AGNrate > EDDrate)
	    AGNrate = EDDrate;

	  /*  accreted mass onto black hole the value of dt puts an h factor into AGNaccreted as required for code units */
	  AGNaccreted = AGNrate * dt;

	  /* cannot accrete more mass than is available! */
	  if(AGNaccreted > HotGas)
	    AGNaccreted = HotGas;

	  /*  coefficient to heat the cooling gas back to the virial temperature of the halo */
	  /*  1.34e5 = sqrt(2*eta*c^2), eta=0.1 (standard efficiency) and c in km/s
	   *  Eqs. 11 & 12 in Croton 2006 */
	  AGNcoeff = (1.34e5 / Vvir) * (1.34e5 / Vvir);

	  /*  cooling mass that can be suppressed from AGN heating */
	  AGNheating = AGNcoeff * AGNaccreted;


	  if(AGNRadioModeModel == 0 && Gal[p].Type==1) 
	  {
	      dist=separation_gal(p,FOF_centralgal);
	      if(dist < Gal[FOF_centralgal].Rvir) 
	      {
		  if (AGNheating > (Gal[p].CoolingGas + Gal[FOF_centralgal].CoolingGas)) 
		  {
		      AGNheating = (Gal[p].CoolingGas + Gal[FOF_centralgal].CoolingGas);
		      AGNaccreted = (Gal[p].CoolingGas + Gal[FOF_centralgal].CoolingGas) / AGNcoeff;
		  }
		  if (AGNheating > Gal[p].CoolingGas) LeftOverEnergy = AGNheating - Gal[p].CoolingGas;
	      }
	      // Added the following line in as a fix because without it Type 1 BHs can grow ridiculously large.
	      // In actual fact, this whole section needs fully documenting as it is not clear what is supposed to be going on.
	      else 
		  if (AGNheating > Gal[p].CoolingGas) AGNaccreted = Gal[p].CoolingGas / AGNcoeff;
	  }
	  else
	      if (AGNheating > Gal[p].CoolingGas) AGNaccreted = Gal[p].CoolingGas / AGNcoeff;

	  /* limit heating to cooling rate */
	  if (AGNheating > Gal[p].CoolingGas) AGNheating = Gal[p].CoolingGas;

	  mass_checks(p,"cooling.c",__LINE__);

	  /*  accreted mass onto black hole */
	  Gal[p].BlackHoleMass += AGNaccreted; //ROB: transfer_mass functions should be used here
	  Gal[p].RadioAccretionRate += AGNaccreted / (dt*STEPS);
	  fraction=AGNaccreted/Gal[p].HotGas;
	  Gal[p].HotGas -= AGNaccreted;
	  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	    Gal[p].MetalsHotGas[ii] += (-fraction * Gal[p].MetalsHotGas[ii]);
	  mass_checks(p,"cooling.c",__LINE__);
#ifdef INDIVIDUAL_ELEMENTS
	  int kk;
	  for(kk=0;kk<NUM_ELEMENTS;kk++)
	    Gal[p].HotGas_elements[kk] *= (1-fraction);
#endif
#ifdef METALS_SELF
	  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	    Gal[p].MetalsHotGasSelf[ii] += (-fraction * Gal[p].MetalsHotGasSelf[ii]);
#endif	

	}
      else
	AGNheating = 0.0;

      mass_checks(p,"cooling.c",__LINE__);

      Gal[p].CoolingGas -= AGNheating;

      if(Gal[p].CoolingGas < 0.0)
	Gal[p].CoolingGas = 0.0;

      Gal[p].CoolingRate += Gal[p].CoolingGas / (dt*STEPS);


      if(AGNRadioModeModel == 0 && LeftOverEnergy>0.)
	{
	  Gal[FOF_centralgal].CoolingGas -= LeftOverEnergy;

	  if(Gal[FOF_centralgal].CoolingGas < 0.0)
	    Gal[FOF_centralgal].CoolingGas = 0.0;
	  else
	    Gal[FOF_centralgal].CoolingRate -= LeftOverEnergy / (dt*STEPS);
	}


      mass_checks(p,"cooling.c",__LINE__);

  }
}


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
      //determine the xray luminosity of any cooling gas in this snapshot (White & Frenk 1991 eq21)
      Gal[p].XrayLum = log10(2.5 * (Mcool / dt) * 6.31 * Gal[p].Vvir * Gal[p].Vvir) + 35.0;

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
  else //if(Mcool > 0.0)
    Gal[p].XrayLum = 0.0;
}

