#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

// TODO this could be in recipe_starformation_and_feedback

/** @file  recipe_reincorporation.c
 *  @brief recipe_reincorporation.c calculates the fraction of ejected gas that
 *         gets reincorporated into the hot fraction per timestep.
 *
 *         2 options are available to reincorporate the gas from the external
 *         reservoir:
 *           -\f$\dot{M}_{\rm{eject}}=
 *                -\gamma \left(\frac{M_{\rm{ejected}}}{t_{\rm{dyn,h}}}\right)\f$
 *         (Eq. 3 Delucia2004) (ReIncorporationRecipe == 1);
 *           -\f$\dot{M}_{\rm{eject}}=
 *                -\gamma \left(\frac{V_{\rm{vir}}}{\rm{220km/s}}\right)
 *                \left(\frac{M_{\rm{ejected}}}{t_{\rm{dyn,h}}}\right)\f$
 *         (Eq. 23 Guo2010) (ReIncorporationRecipe == 2)
 *
 *         Note that ReIncorporationRecipe == 1 doesn't correspond to Delucia2007
 *         anymore as Qi changed the place where the reincorporated gas ends.
 *         In delucia2007 only central galaxies reincorporated. Now satellites
 *         outside central Rvir can also do so.
 *
 **/
/** @brief reincorporates ejected gas back into the central galaxy hot halo */

void reincorporate_gas(int p, double dt)
{
  double reincorporated, fraction, reinc_time;

  reincorporated = 0.;

  mass_checks("reincorporate_gas #1",p);

  /* DeLucia2007 -> Mdot_eject=-gama_ej * M_ejected/tdyn */
  if(ReIncorporationRecipe == 0)
 	reincorporated = ReIncorporationFactor * Gal[p].EjectedMass / (Gal[p].Rvir / Gal[p].Vvir) * dt;

  /* Guo2010 -> Mdot_eject=-gama_ej * M_ejected/tdyn * Vvir/220 */
  else
 	if(ReIncorporationRecipe == 1)
      reincorporated = ReIncorporationFactor * Gal[p].EjectedMass / (Gal[p].Rvir / Gal[p].Vvir) * Gal[p].Vvir/220. *dt ;

   /* Henriques2012b Mdot_eject=-gama_ej*M_ejected*M_vir
    * Mvir should be in units of 1e12, but inside the
    * code Mvir is already in units of 1.e10*/
 	else
 	  if(ReIncorporationRecipe == 2)
 	  {
 	  	// Oppenheimer & Dave 2008
 	  	/* reinc_time= pow(Gal[p].Mvir/Hubble_h,-0.6)*ReIncorporationFactor/UnitTime_in_years;
 	  	   if(Gal[p].Mvir/Hubble_h > 500.)
 	  	   reinc_time= pow(500.,-0.6)*ReIncorporationFactor/UnitTime_in_years;
 	  	   reincorporated = Gal[p].EjectedMass / reinc_time * dt; */


 	  	reinc_time= (Hubble_h/Gal[p].Mvir)*(ReIncorporationFactor/UnitTime_in_years);
 	  	reincorporated = Gal[p].EjectedMass / reinc_time * dt;

 	    /*reinc_time= ReIncorporationFactor/UnitTime_in_years * pow((1+ZZ[Gal[p].SnapNum]),ReincZpower)
 	                * pow(Gal[p].Vvir/220.,-ReincVelocitypower) / (hubble_of_z(Gal[p].HaloNr));
 	            reincorporated = Gal[p].EjectedMass / reinc_time  *dt ;*/



 	  }

  if(FeedbackRecipe == 1)
  {
	  reincorporated = ReIncorporationFactor * Gal[p].EjectedMass /
	  	(Gal[p].Rvir * min(FeedbackEjectionEfficiency,1.)*sqrt(EtaSNcode * EnergySNcode)/(Gal[p].Vvir*Gal[p].Vvir))
	  	* Gal[p].Vvir/220. * 1.e-6* dt ;

	  //reincorporated = ReIncorporationFactor * Gal[p].EjectedMass / (Gal[p].Rvir / Gal[p].Vvir) *
	  	//                   min(FeedbackEjectionEfficiency,1.) * Gal[p].Vvir/220. *dt ;
	  //reincorporated = ReIncorporationFactor * Gal[p].EjectedMass / (Gal[p].Rvir / Gal[p].Vvir) * Gal[p].Vvir/220. *dt ;
	  //reincorporated = 0.;
  }

  if (reincorporated > Gal[p].EjectedMass)
    reincorporated = Gal[p].EjectedMass;
	
  mass_checks("reincorporate_gas #1.5",p);

  /*Update ejected and hot gas contents*/
  if (Gal[p].EjectedMass > 0.) {
    fraction=((float)reincorporated)/Gal[p].EjectedMass;
    transfer_gas(p,"Hot",p,"Ejected",fraction,"reincorporate_gas", __LINE__);
  }

  mass_checks("reincorporate_gas #2",p);

}

