/*  Copyright (C) <2016>  <L-Galaxies>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @file model_cooling.c
 *  @brief model_cooling.c calculates the amount of mass that cools
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
  double Vvir, Rvir, x, lambda, tcool, rcool, temp, tot_hotMass, tot_metals, HotRadius;
  double coolingGas, logZ, rho_rcool, rho0;

  mass_checks("cooling_recipe #1",p);

  tot_hotMass = Gal[p].HotGas;
  tot_metals = metals_total(Gal[p].MetalsHotGas);

  if(tot_hotMass > 1.0e-6)
  {

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
      coolingGas = tot_hotMass / (HotRadius / Vvir) * dt;
    else // HOT PHASE REGIME
      /*coolingGas = (tot_hotMass / Rvir) * (rcool / tcool) * dt */
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

  mass_checks("cooling_recipe #1.5",p);

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



void do_AGN_heating(double dt, int ngal)
{
  double AGNrate, AGNheating, AGNaccreted, AGNcoeff, fraction, EDDrate, FreeFallRadius;
  double dist, HotGas, HotRadius, Rvir, Vvir, Mvir;
  double LeftOverEnergy, CoolingGas, AGNAccretedFromCentral;
  int p, FoFCentralGal;

  if(AGNRadioModeModel == 0)
    {
      for (p = 0; p < ngal; p++)
	if(Gal[p].Type == 0)
	  FoFCentralGal=p;
    }

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

      if(HotGas > 0.0)
	{
	  if(AGNRadioModeModel == 0)
	    AGNrate = AgnEfficiency * (UnitTime_in_s*SOLAR_MASS)/(UNITMASS_IN_G*SEC_PER_YEAR)
	              * Gal[p].BlackHoleMass/Hubble_h * (HotGas/Hubble_h) * 10.;
	  else if(AGNRadioModeModel == 2)
	    {
	      //empirical (standard) accretion recipe - Eq. 10 in Croton 2006
	      AGNrate = AgnEfficiency / (UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS)
  	  		* (Gal[p].BlackHoleMass / 0.01) * pow3(Vvir / 200.0)
			* ((HotGas / HotRadius * Rvir / Mvir) / 0.1);
	    }
	  else if(AGNRadioModeModel == 3 || AGNRadioModeModel == 4)
	    {
	      double x, lambda, temp, logZ, tot_metals;

	      tot_metals = metals_total(Gal[p].MetalsHotGas);

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
	      if(AGNRadioModeModel == 3)
		AGNrate = (2.5 * M_PI * G) * (0.75 * 0.6 * x) * Gal[p].BlackHoleMass * 0.15;
	      else if(AGNRadioModeModel == 4)
		{
		  /* Cold cloud accretion recipe -- trigger: Rff = 50 Rdisk,
		   * and accretion rate = 0.01% cooling rate
		   * Eq. 25 in Croton 2006 */
		  FreeFallRadius = HotGas / (6.0 * 0.6 * x * Rvir * Vvir) /	HotRadius * Rvir;
		  if(Gal[p].BlackHoleMass > 0.0 && FreeFallRadius < Gal[p].GasDiskRadius * 50.0)
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
	      if(dist < Gal[FoFCentralGal].Rvir)
		{
		  if(AGNheating > (Gal[p].CoolingGas + Gal[FoFCentralGal].CoolingGas))
		    {
		      AGNheating = (Gal[p].CoolingGas + Gal[FoFCentralGal].CoolingGas);
		      AGNaccreted = (Gal[p].CoolingGas + Gal[FoFCentralGal].CoolingGas) / AGNcoeff;
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




	  /*  accreted mass onto black hole */
	  Gal[p].BlackHoleMass += AGNaccreted; //ROB: transfer_mass functions should be used here
	  Gal[p].RadioAccretionRate += AGNaccreted / (dt*STEPS);
	  fraction=AGNaccreted/Gal[p].HotGas;
	  Gal[p].HotGas -= AGNaccreted;
	  Gal[p].MetalsHotGas = metals_add(Gal[p].MetalsHotGas,Gal[p].MetalsHotGas, -fraction);

#ifdef INDIVIDUAL_ELEMENTS
	  Gal[p].HotGas_elements = elements_add(Gal[p].HotGas_elements,Gal[p].HotGas_elements,-fraction);
#endif
#ifdef METALS_SELF
	  Gal[p].MetalsHotGasSelf = 	metals_add(Gal[p].MetalsHotGasSelf,Gal[p].MetalsHotGasSelf,-fraction);
#endif	

	}
      else
	AGNheating = 0.0;


      Gal[p].CoolingGas -= AGNheating;

      if(Gal[p].CoolingGas < 0.0)
	Gal[p].CoolingGas = 0.0;

      Gal[p].CoolingRate += Gal[p].CoolingGas / (dt*STEPS);

      if(AGNRadioModeModel == 0 && LeftOverEnergy>0.)
  	{
	  Gal[FoFCentralGal].CoolingGas -= LeftOverEnergy;

	  if(Gal[FoFCentralGal].CoolingGas < 0.0)
	    Gal[FoFCentralGal].CoolingGas = 0.0;
	  else
	    Gal[FoFCentralGal].CoolingRate -= LeftOverEnergy / (dt*STEPS);
  	}


      mass_checks("cooling_recipe #2.",p);

  }
}


/** @brief updates the fractions of hot and cold gas due to cooling. */

void cool_gas_onto_galaxy(int p, double dt)
{
  double fraction,Mdisk,Mcool;
  int i;

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

    // We already know that 0<mcool<=Gal[p].HotGas
    fraction=((float)Mcool)/Gal[p].HotGas;
    transfer_gas(p,"Cold",p,"Hot",fraction,"cool_gas_onto_galaxy", __LINE__);

    if (DiskRadiusModel == 0)
    {
      if (Gal[p].ColdGas != 0.0)
      	for (i=0;i<3;i++)
      		Gal[p].GasSpin[i]=(Gal[p].GasSpin[i]*Mdisk+Gal[p].HaloSpin[i]*Mcool)/(Gal[p].ColdGas);
      get_gas_disk_radius(p);
    }
  }
  else
  	Gal[p].XrayLum = 0.0;
}



