#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @file recipe_disrupt.c
 *  @brief recipe_disrupt.c checks if a type 2 satellite galaxy should
 *         or not be disrupted due to tidal forces
 *
 *  This routine takes into account the tidal effects that satellite
 *  galaxies experience while orbiting a central companion. Since the
 *  baryonic component is more compact and denser that the dark matter
 *  it assumes that only type 2 satellites are affected (those that have
 *  already lost their dark matter halo). It assumes that the disruption
 *  is complete and instantaneous with all the satellite material being
 *  transferred into the central galaxy.
 *
 *  The satellite is assumed to orbit a singular isothermal potential:
 *
 *  \f$\phi(R)=V^2_{\rm{vir}}\rm{ln} R\f$ (Eq. 28 Guo2010)
 *
 *  Assuming conservation of energy and angular momentum along the
 *  orbit, its pericentric distance (closest point to the centre along
 *  the orbit) can be estimated from:
 *
 *  \f$\left(\frac{R}{R_{\rm{peri}}}\right)^2=
 *  \frac{lnR/R_{\rm{peri}}+\frac{1}{2}(V/V_{\rm{vir}})^2}
 *  {\frac{1}{2}(V_{\rm{t}}/V_{\rm{vir}})^2}\f$
 *  (Eq. 29 Guo2010).
 *
 *  The main halo density at this point is compared with the baryonic
 *  mass (cold gas + stellar) density of the satellite within its half
 *  mass radius. If
 *
 *  \f$ \frac{M_{\rm{DM,halo}}(R_{\rm{peri}})}{R^3_{\rm{peri}}}\equiv
 *  \rho_{\rm{DM,halo}}>
 *  \rho_{\rm{sat}}\equiv\frac{M_{\rm{sat}}}{R^3_{\rm{sat,half}}}\f$
 *  (Eq. 30 Guo2010)
 *
 *  the galaxy is disrupted.
 *
 *  TODO: shouldn't galaxy p be voided after disruption?
 *
 *  */


void disrupt(int p, int centralgal)
{
  double rho_sat, rho_cen;
  double cen_mass, r_sat, radius;
#ifdef H2_AND_RINGS
  int jj;
  double fractionRings[RNUM];
#endif
  /* If the main halo density at the pericentre (closest point in the orbit
   * to the central galaxy)is larger than the satellite's density at the
   * half mass radius, the satellite is completely disrupted. Note that since
   * the satellite is a type 2 the only mass components remaining and
   * contributing to the density are the cold gas and stellar mass. */

  //TODO If we are passing in centralgal then we should not set it here 
  centralgal=Gal[p].CentralGal;
 
  mass_checks(centralgal,"model_disrupt.c",__LINE__);
  mass_checks(p,"model_disrupt.c",__LINE__);

  /* Radius calculated at the peri-point */
  radius=peri_radius(p, centralgal);
  if (radius < 0) {
   terminate("must be wrong \n");
  }

  /* Calculate the density of the main central halo at radius (the peri-centre).
   * The tidal forces are caused by the dark matter of the main halo, hence Mvir
   * is used. Assume isothermal. */
  cen_mass=Gal[centralgal].Mvir*radius/Gal[centralgal].Rvir;
  rho_cen=cen_mass/pow3(radius);

  if (Gal[p].DiskMass+Gal[p].BulgeMass>0)
  {
      // Calculate the rho according to the real geometry using all the baryonic material (ColdGas+DiskMass+BulgeMass)
      int do_ColdGas=1, do_DiskMass=1, do_BulgeMass=1;
      rho_sat=(Gal[p].DiskMass+Gal[p].BulgeMass+Gal[p].ColdGas)/pow3(half_mass_radius(p,do_ColdGas,do_DiskMass,do_BulgeMass));
  }
  else
    rho_sat=0.0;

  /* If density of the main halo is larger than that of the satellite baryonic
   * component, complete and instantaneous disruption is assumed. Galaxy becomes
   * a type 3 and all its material is transferred to the central galaxy. */
  if (rho_cen > rho_sat)
  {
  	Gal[p].Type = 3;
#ifdef GALAXYTREE
    int q;
    q = Gal[Gal[p].CentralGal].FirstProgGal;
    if (q >= 0)
    {
      // add progenitors of Gal[p] to the list of progentitors of Gal[p].CentralGal
      while (GalTree[q].NextProgGal >= 0)
      	q = GalTree[q].NextProgGal;
	
      GalTree[q].NextProgGal = Gal[p].FirstProgGal;
      
      if(GalTree[q].NextProgGal >= NGalTree)
	    {
	      printf("q=%d p=%d GalTree[q].NextProgGal=%d NGalTree=%d\n",
		     q, p, GalTree[q].NextProgGal, NGalTree);
	      terminate("problem");
	    }
    }

    if(q < 0)
    	terminate("this shouldn't happen");
	
     // TODO if !(q>=0) shouldn't we set
      // Gal[Gal[p].CentralGal].FirstProgGal to Gal[p].FirstProgGal  ??

    q = GalTree[q].NextProgGal;

    if(q < 0)
    	terminate("inconsistency");

    if(HaloGal[GalTree[q].HaloGalIndex].GalTreeIndex != q)
    	terminate("inconsistency");

    HaloGal[GalTree[q].HaloGalIndex].DisruptOn = 1;
#endif
    /* Put gas component to the central galaxy hot gas and stellar material into the ICM.
     * Note that the satellite should have no extended components. */
    // TODO Shouldn't the stars end up in the bulge? - this is a close merger

    //if BlackHoleDisruptGrowthRate>0 some disrupted mass goes into the black hole
    Gal[centralgal].BlackHoleMass+=BlackHoleDisruptGrowthRate*Gal[p].ColdGas;
    Gal[p].ColdGas-=BlackHoleDisruptGrowthRate*Gal[p].ColdGas;

#ifdef TRACK_MASSGROWTH_CHANNELS
  	Gal[p].MassFromInSitu = 0.;
  	Gal[p].MassFromMergers = 0.;
  	Gal[p].MassFromBursts = 0.;
#ifdef STAR_FORMATION_HISTORY
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
  	int ii;
  	for (ii=0; ii<=Gal[p].sfh_ibin; ii++) Gal[p].sfh_MassFromInSitu[ii] = 0.;
  	for (ii=0; ii<=Gal[p].sfh_ibin; ii++) Gal[p].sfh_MassFromMergers[ii] = 0.;
  	for (ii=0; ii<=Gal[p].sfh_ibin; ii++) Gal[p].sfh_MassFromBursts[ii] = 0.;
#endif
#endif
#endif

#ifdef TRACK_BURST
    /* Transfer burst component first */
    //transfer_material(centralgal,"BurstMass",p,"BurstMass",
	//	   (Gal[p].DiskMass+Gal[p].BulgeMass)/(Gal[p].DiskMass+Gal[p].BulgeMass+Gal[p].ICM),"disrupt.c", __LINE__);
#endif

#ifdef H2_AND_RINGS
    for (jj=0;jj<RNUM;jj++)
      fractionRings[jj]=1.;
    transfer_material_with_rings(centralgal,"ICM",p,"DiskMass",fractionRings,"disrupt.c", __LINE__);
    transfer_material_with_rings(centralgal,"ICM",p,"BulgeMass",fractionRings,"disrupt.c", __LINE__);
#else
    transfer_material(centralgal,"ICM",p,"BulgeMass",1.,"disrupt.c", __LINE__);
    transfer_material(centralgal,"ICM",p,"DiskMass",1.,"disrupt.c", __LINE__);
#endif

#ifdef H2_AND_RINGS
    transfer_material_with_rings(centralgal,"HotGas",p,"ColdGas",fractionRings,"disrupt.c", __LINE__);
#else
    transfer_material(centralgal,"HotGas",p,"ColdGas",1.,"disrupt.c", __LINE__);
#endif
    transfer_material(centralgal,"HotGas",p,"HotGas",1.,"disrupt.c", __LINE__);
    //transfer_material(centralgal,"ReheatedGas",p,"ReheatedGas",1.,"disrupt.c", __LINE__);

  } //if (rho_cen > rho_sat)
  mass_checks(centralgal,"model_disrupt.c",__LINE__);
  mass_checks(p,"model_disrupt.c",__LINE__);
  
}

/** @brief Calculates the distance of the satellite to the pericentre of the
  *        main dark matter halo. */

double peri_radius(int p, int centralgal)
{
  int i;
  double a, b, v[3], r[3], x, x0;
  for(i = 0; i < 3; i++)
    {
	  r[i] = wrap(Gal[p].Pos[i]-Gal[centralgal].Pos[i],BoxSize);
	  r[i] /= (1 + ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);
	  v[i] = Gal[p].Vel[i] - Gal[centralgal].Vel[i];
	  //r[i] = wrap(Gal[p].Pos_notupdated[i]-Gal[centralgal].Pos[i],BoxSize);
	  //r[i] /= (1 + ZZ[Halo[Gal[centralgal].HaloNr].SnapNum]);
	  //v[i] = Gal[p].Vel_notupdated[i] - Gal[centralgal].Vel[i];
    }

  b = 1 / 2. * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) / pow2(Gal[centralgal].Vvir);
  a =
    1 / 2. * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] -
	      pow2(r[0] * v[0] + r[1] * v[1] + r[2] * v[2])
		   / (r[0] * r[0] + r[1] * r[1] + r[2] * r[2])) / pow2(Gal[centralgal].Vvir);

  x = sqrt(b / a);
  x0 = 1000;
  while(fabs(x0 - x) >= 1.e-8)
    {
      x0 = x;
      x = sqrt((log(x0) + b) / a);
    }
  if(x == 0)
    {
      terminate("wrong in peri_radius \n");
    }

  return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])/x;
}



