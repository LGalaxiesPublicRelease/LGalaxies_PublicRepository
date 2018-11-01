#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"


/**@file recipe_infall.c
 * @brief recipe_infall.c calculates the amount of gas that infalls
 *       into the galaxy hot gas component at each time step. This
 *       is derived from the baryonic fraction taking reionization
 *       into account.
 *
 *         */

double infall_recipe(int ngal)
{
    int FOF_centralgal;   // Galaxy at centre of FOF group
    int i;
    double tot_mass, reionization_modifier, infallingMass;
    double mvir;
#ifndef INFALL_UPDATE
    double dis, rvir;     // Mvir, Rvir of FOF_centralgal
#endif
    double zcurr;   // z+1
    
    /*  need to add up all the baryonic mass asociated with the full halo to check
     *  what baryonic fraction is missing/in excess. That will give the mass of gas
     *  that need to be added/subtracted to the hot phase, gas that infalled.*/
    tot_mass = 0.0;
    FOF_centralgal = Gal[0].FOFCentralGal;
    mvir = Gal[FOF_centralgal].Mvir;
    zcurr = ZZ[Halo[Gal[0].HaloNr].SnapNum];

    for(i = 0; i < ngal; i++) {    /* Loop over all galaxies in the FoF-halo */


	/* If galaxy is orbiting a galaxy inside Rvir of the type 0 it will contribute
	 * to the baryon sum */
	/* Note: Rob Yate's fix to excess baryon content is simply to comment
	 * out the following if statement so that we add up all the baryons in
	 * the FOF group. */
#ifndef INFALL_UPDATE
	/* dis is the separation of the galaxy which i orbits from the type 0 */
	dis=separation_gal(FOF_centralgal,Gal[i].CentralGal)/(zcurr+1);
	rvir = Gal[FOF_centralgal].Rvir;

	if ( dis < rvir ) {
#endif
	    tot_mass += Gal[i].DiskMass + Gal[i].BulgeMass + Gal[i].ICM + Gal[i].BlackHoleMass;
	    tot_mass += Gal[i].ColdGas + Gal[i].HotGas  + Gal[i].EjectedMass + Gal[i].BlackHoleGas;
	    //tot_mass += Gal[i].ColdGas + Gal[i].HotGas  + Gal[i].ReheatedGas + Gal[i].EjectedMass + Gal[i].BlackHoleGas;
#ifndef INFALL_UPDATE
	}
#endif
    }

    /* The logic here seems to be that the infalling mass is supposed to be
     * newly-accreted diffuse gas.  So we take the total mass and subtract all
     * the components that we already know about.  In principle, this could lead
     * to negative infall. */
    /* The reionization modifier is applied to the whole baryon mass, not just the 
     * diffuse component.  It is not obvious that this is the correct thing to do. */
    /* The baryonic fraction is conserved by adding/subtracting the infallingMass
     * calculated here to/from the hot gas of the central galaxy of the FOF
     * This is done in main.c where the infall recipe is called.
     * If ReionizationOn=1, the impact of reonization on the fraction of infalling
     * gas is computed, this is done using the Gnedin formalism with a choice
     * of fitting parameters to the formulas proposed by these authors.
     * If ReionizationOn=0, Okamoto2008 is used. In both cases reionization has
     * the effect of reducing the fraction of baryons that collapse into dark
     * matter halos, reducing the
     * amount of infalling gas. */
    if(ReionizationModel == 2)
	reionization_modifier = 1.0;
    else
	reionization_modifier = do_reionization(mvir, zcurr);

    infallingMass = reionization_modifier * BaryonFrac * mvir - tot_mass;

    return infallingMass;
  
}

double do_reionization(float Mvir, double Zcurr)
{
  double modifier, a, alpha;
  //Gnedin (2000)
  double f_of_a, a_on_a0, a_on_ar, Mfiltering, Mjeans, Mchar, mass_to_use;
  double Tvir, Vchar, omegaZ, xZ, deltacritZ, HubbleZ;
  //Okamoto (2008)
  double x0,x,delta_c,delta_0,tau,Mc;
  int tabindex; 
  double f1, f2;

  modifier=1.;
  if (ReionizationModel == 0) {
      /* reionization recipie described in Gnedin (2000), with the fitting
       * from Okamoto et al. 2008 -> Qi(2010)*/

      alpha = 2.0;
      a = 1. / (1 + Zcurr);
      a0 = 1.;

      x = - (1 - Omega) * a * a * a / (Omega + (1-Omega) * a * a * a);
      delta_c = (178 + 82 *x - 39 * x * x) / (1. + x);

      x0 = - (1 - Omega) * a0 * a0 * a0 / (Omega + (1-Omega) * a0 * a0 * a0);
      delta_0 = (178 + 82 *x0 - 39 * x0 * x0) / (1. + x0);

      tau = 0.73 * pow(Zcurr + 1,0.18)* exp(-pow(0.25 * Zcurr, 2.1));

      Mc = pow(tau  / (1 + Zcurr),3./2) * sqrt(delta_0 / delta_c);

      /* if use Okamoto et al. 2008*/

      find_interpolate_reionization(Zcurr, &tabindex, &f1, &f2);
      Mc = f1*log10(Reion_Mc[tabindex])+f2*log10(Reion_Mc[tabindex+1]);
      Mc = pow(10, Mc-10);

      modifier = pow(1 + (pow(2, alpha/3.) -1) * pow(Mc / Mvir, alpha), -3./alpha);
  }
  else if (ReionizationModel == 1)
    {
      /** reionization recipie described in Gnedin (2000), using the fitting */
      /*  formulas given by Kravtsov et al (2004) Appendix B, used after Delucia 2004*/

      /*  here are two parameters that Kravtsov et al keep fixed. */
      /*  alpha gives the best fit to the Gnedin data */
      alpha = 6.0;
      Tvir = 1e4;

      /*  calculate the filtering mass */

      a = 1.0 / (1.0 + Zcurr);
      a_on_a0 = a / a0;
      a_on_ar = a / ar;

      if(a <= a0)
	f_of_a = 3.0 * a / ((2.0 * alpha) * (5.0 + 2.0 * alpha)) * pow(a_on_a0, alpha);
      else if((a > a0) && (a < ar))
	f_of_a = (3.0 / a) * a0 * a0 *
	  (1.0 / (2.0 + alpha) - 2.0 * pow(a_on_a0, -0.5) / (5.0 + 2.0 * alpha)) +
	      a * a / 10.0 - (a0 * a0 / 10.0) * (5.0 - 4.0 * pow(a_on_a0, -0.5));
      else
	f_of_a = (3.0 / a) * (a0 * a0 * (1.0 / (2.0 + alpha) - 2.0 * pow(a_on_a0, -0.5) / (5.0 + 2.0 * alpha)) +
		  (ar * ar / 10.0) * (5.0 - 4.0 * pow(a_on_ar, -0.5)) - (a0 * a0 / 10.0) * (5.0 - 4.0 * pow(a_on_a0, -0.5)) +
		     a * ar / 3.0 - (ar * ar / 3.0) * (3.0 - 2.0 * pow(a_on_ar, -0.5)));
    
    /*  this is in units of 10^10Msun/h, note mu=0.59 and mu^-1.5 = 2.21 */
    Mjeans = 25.0 * pow(Omega, -0.5) * 2.21;
    Mfiltering = Mjeans * pow(f_of_a, 1.5);
    
    
    /*  calculate the characteristic mass coresponding to a halo temperature of 10^4K */
    Vchar = sqrt(Tvir / 36.0);
    omegaZ = Omega * (pow3(1.0 + Zcurr) / (Omega * pow3(1.0 + Zcurr) + OmegaLambda));
    xZ = omegaZ - 1.0;
    deltacritZ = 18.0 * M_PI * M_PI + 82.0 * xZ - 39.0 * xZ * xZ;
    HubbleZ = Hubble * sqrt(Omega * pow3(1.0 + Zcurr) + OmegaLambda);
    
    Mchar = Vchar * Vchar * Vchar / (G * HubbleZ * sqrt(0.5 * deltacritZ));

    /*  we use the maximum of Mfiltering and Mchar */
    mass_to_use = max(Mfiltering, Mchar);
    modifier = 1.0 / pow3(1.0 + 0.26 * (mass_to_use / Mvir));
  }
  else if (ReionizationModel ==2)
    {
      printf("Should not be called with this option\n");
      exit(0);
    }
  
  return modifier;

}


/* The gas that infalled is added to the hot gas of the central galaxy. 
 *
 * TODO: Note that the algorithm allows negative infall in order to balance 
 * the total baryon content.  In that case, should we decrement the metals
 * accordingly? - or change algorithm NOT to allow negative infall? */

void add_infall_to_hot(double infallingGas) {

    int FOF_centralgal;
    
    FOF_centralgal = Gal[0].FOFCentralGal;

#ifdef EXCESS_MASS
    if (InfallModel == 0) {
#endif

#if defined(GUO10) || defined(GUO13) || defined(HENRIQUES13)
//if infallingGas is negative set the limit to be removed at hotgas
  if(infallingGas<0. && -1.*infallingGas>Gal[FOF_centralgal].HotGas)
    infallingGas=-Gal[FOF_centralgal].HotGas;
#endif
  Gal[FOF_centralgal].HotGas += infallingGas;

#ifdef INDIVIDUAL_ELEMENTS
    Gal[FOF_centralgal].HotGas_elements[0] += 0.75*(infallingGas/Hubble_h)*1.0e10; //H
    Gal[FOF_centralgal].HotGas_elements[1] += 0.25*(infallingGas/Hubble_h)*1.0e10; //He
#endif

#ifdef EXCESS_MASS
    } else if (InfallModel == 1) {
	double fraction;
	if (infallingGas > 0.) {
	    double M_infalltoHot, newgas;
	    // By preference, we take the infalling gas from any excess component that already exists
	    newgas = infallingGas - Gal[FOF_centralgal].ExcessMass;
#ifdef INDIVIDUAL_ELEMENTS
	    // Extra hydrogen and helium added to ExcessMass (to be reincorporated below)
	    if (newgas > 0.) {
		Gal[FOF_centralgal].ExcessMass_elements[0] += 0.75*newgas*1.0e10/Hubble_h; //H
		Gal[FOF_centralgal].ExcessMass_elements[1] += 0.25*newgas*1.0e10/Hubble_h; //He
	    }
#endif
	    // Top up ExcessMass if required to reach infallingGas; no new elements are needed (except H, He done above)
	    Gal[FOF_centralgal].ExcessMass = max(Gal[FOF_centralgal].ExcessMass,infallingGas);
	    // We are going to transfer to the HotGas and Ejected phases in the same proportion as already exists
      	    M_infalltoHot = (Gal[FOF_centralgal].HotGas/(Gal[FOF_centralgal].HotGas + Gal[FOF_centralgal].EjectedMass)) * infallingGas;
	    // To minimise mixing, first transfer as much of M_infalltoHot as possible from Ejected to Hot
	    if (Gal[FOF_centralgal].EjectedMass > 0.) {
		fraction = min(M_infalltoHot/Gal[FOF_centralgal].EjectedMass,1.);
		M_infalltoHot -= fraction*Gal[FOF_centralgal].EjectedMass;
		transfer_material(FOF_centralgal,"HotGas",FOF_centralgal,"EjectedMass",fraction,"add_infall_to_hot", __LINE__);	    
	    }
	    // Next add the infallingGas to the Ejected component
	    fraction = infallingGas/Gal[FOF_centralgal].ExcessMass;
	    transfer_material(FOF_centralgal,"EjectedMass",FOF_centralgal,"ExcessMass",fraction,"add_infall_to_hot", __LINE__);
	    // Then transfer the rest of M_infalltoHot from the Ejected component to the HotGas component, if required
	    if (M_infalltoHot > 0.) {
		fraction = M_infalltoHot/Gal[FOF_centralgal].EjectedMass;
		// Fix for rounding error whilst still catching genuine bugs.
		if (fraction>1. && fraction<1.0001) fraction=1.;
		transfer_material(FOF_centralgal,"HotGas",FOF_centralgal,"EjectedMass",fraction,"add_infall_to_hot", __LINE__);
	    }
	} else {  // infallingGas <= 0.
	    double outflowingGas, M_outfromEjected, M_outfromHotGas;
	    // Can't expel more gas than we already have:
	    outflowingGas = min(-infallingGas,Gal[FOF_centralgal].HotGas+Gal[FOF_centralgal].EjectedMass);
	    if (outflowingGas > 0.) {
		// We are going to transfer from the HotGas and Ejected phases in the same proportion as already exists
		M_outfromEjected = (Gal[FOF_centralgal].EjectedMass/(Gal[FOF_centralgal].HotGas + Gal[FOF_centralgal].EjectedMass)) * outflowingGas;
		M_outfromHotGas = outflowingGas - M_outfromEjected;
		// To minimise mixing, first transfer as much of M_outfromEjected as possible from Ejected to Excess
		if (Gal[FOF_centralgal].EjectedMass > 0.) {
		    fraction = min(M_outfromEjected/Gal[FOF_centralgal].EjectedMass,1.);
		    M_outfromEjected -= fraction*Gal[FOF_centralgal].EjectedMass;
		    transfer_material(FOF_centralgal,"ExcessMass",FOF_centralgal,"EjectedMass",fraction,"add_infall_to_hot", __LINE__);
		}
		// Next move from HotGas to the Ejected component
		if (Gal[FOF_centralgal].HotGas > 0.) {
		    fraction = M_outfromHotGas/Gal[FOF_centralgal].HotGas;
		    transfer_material(FOF_centralgal,"EjectedMass",FOF_centralgal,"HotGas",fraction,"add_infall_to_hot", __LINE__);
		}
		// Then transfer the rest of M_outfromEjected from the Ejected component to the Excess component, if required
		if (M_outfromEjected > 0.) {
		    fraction = M_outfromEjected/Gal[FOF_centralgal].EjectedMass;
		    transfer_material(FOF_centralgal,"ExcessMass",FOF_centralgal,"EjectedMass",fraction,"add_infall_to_hot", __LINE__);
		}
	    }
	}
    }
#endif //EXCESS_MASS

}
