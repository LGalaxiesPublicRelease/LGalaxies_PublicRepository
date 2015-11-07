#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @file  model_reincorporation.c
 *  @brief model_reincorporation.c calculates the fraction of ejected gas that
 *         gets reincorporated into the hot fraction per timestep.
 *
 *         23options are available to reincorporate the gas from the external
 *         reservoir:
 *           -\f$\dot{M}_{\rm{eject}}=
 *                -\gamma \left(\frac{M_{\rm{ejected}}}{t_{\rm{dyn,h}}}\right)\f$
 *         (Eq. 3 Delucia2004) (ReIncorporationModel == 2);
 *           -\f$\dot{M}_{\rm{eject}}=
 *                -\gamma \left(\frac{V_{\rm{vir}}}{\rm{220km/s}}\right)
 *                \left(\frac{M_{\rm{ejected}}}{t_{\rm{dyn,h}}}\right)\f$
 *         (Eq. 23 Guo2010) (ReIncorporationModel == 1)
 *
 **/
/** @brief reincorporates ejected gas back into the central galaxy hot halo */

void reincorporate_gas(int p, double dt)
{
  double reincorporated, fraction, reinc_time;

  reincorporated = 0.;

  mass_checks("reincorporate_gas #1",p);


  if(FeedbackEjectionModel == 0)
    {
      if(ReIncorporationModel == 0)
	{
	  reinc_time= (Hubble_h/Gal[p].Mvir)*(ReIncorporationFactor/UnitTime_in_years);
	  reincorporated = Gal[p].EjectedMass / reinc_time * dt;
	  /* Henriques2013 Mdot_eject=-gama_ej*M_ejected*M_vir Mvir should be in units of 1e12, but inside the
	   * code Mvir is already in units of 1.e10*/
	}
      else
	if(ReIncorporationModel == 1)
	  reincorporated = ReIncorporationFactor * Gal[p].EjectedMass / (Gal[p].Rvir / Gal[p].Vvir) * Gal[p].Vvir/220. *dt ;
      /* Guo2010 -> Mdot_eject=-gama_ej * M_ejected/tdyn * Vvir/220 */
	else
	  if(ReIncorporationModel == 2)
	    reincorporated = ReIncorporationFactor * Gal[p].EjectedMass / (Gal[p].Rvir / Gal[p].Vvir) * dt;
    }
  else if(FeedbackEjectionModel == 1)
    {
      reincorporated = ReIncorporationFactor * Gal[p].EjectedMass /
	  (Gal[p].Rvir * min(FeedbackEjectionEfficiency,1.)*sqrt(EtaSNcode * EnergySNcode)/(Gal[p].Vvir*Gal[p].Vvir))
	  * Gal[p].Vvir/220. * 1.e-6* dt ;
    }

  if (reincorporated > Gal[p].EjectedMass)
    reincorporated = Gal[p].EjectedMass;
	
  mass_checks("reincorporate_gas #1.5",p);

  /*Update ejected and hot gas contents*/
  if (Gal[p].EjectedMass > 0.)
    {
      fraction=((float)reincorporated)/Gal[p].EjectedMass;
      transfer_gas(p,"Hot",p,"Ejected",fraction,"reincorporate_gas", __LINE__);
    }

  mass_checks("reincorporate_gas #2",p);

}

