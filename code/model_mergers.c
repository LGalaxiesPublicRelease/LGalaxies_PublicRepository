#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @file model_mergers.c
 *  @brief Calculates the merging time, the central galaxy (for type 1's),
 *         adds galaxies together, calculates SF from bursts and grows
 *         black holes.
 *
 *
 *       <B>set_merger_center</B> - calculates the central galaxy for type 1's,
 *       since type 1's can also merge. Therefore,
 *       they need a merger central galaxy and will also have a merger clock
 *       (needed for millennium two, since due to the high resolution, haloes
 *       are very difficult to disrupt and orbit around forever).
 *
 *
 *
 *       <B>estimate_merging_time</B> sets up a merger clock. Originally this
 *       was done only for type 2 galaxies. The positions of galaxies in the
 *       model are normally given by the position of their dark matter halo.
 *       However, when galaxies become satellites, their dark matter haloes
 *       are stripped by the central object to the point where there is none
 *       left. At this point, the dark matter of the satellite becomes part of
 *       the main halo, but the galaxy's position should continue to evolve
 *       due to the dynamical friction force caused by the dark matter around
 *       it. \n
 *       The code does keep track of the position of the most bounded particle
 *       when the satellite's halo was disrupted, but this is not used to track
 *       the galaxy position. Instead a clock is set, giving the time left until
 *       the satellite mergers with the central galaxy. Before, this was only
 *       done for type 2's (satellites that lost a halo). Guo2010 included the
 *       option to set the clock also for type 1's (satellites with
 *       a halo), since for the highest resolution millennium 2, their haloes
 *       can be very small and orbit forever around the central companion.\n
 *       This time is computed using the Chandrasekhar's formula for dynamical
 *       friction, as in Binney & Tremaine 1987:
 *
 *       \f$F_{\rm{df}}=
 *               -\frac{4\pi {\rm{G}}^2 m^2_{\rm{sat}} \ln(\Lambda) \rho B(x)}
 *                     {v^2_{\rm{rel}}}.\f$
 *
 *       Which gives (B&T Eq. 7.26):
 *
 *       \f$t_{df}\approx1.17\frac{V_{\rm{vir}}r^2_{\rm{sat}}}
 *                                {{\rm{G}}\,m_{\rm{sat}}\ln(\Lambda)},\f$
 *
 *       that is afterwards multiplied by 2 (after Delucia2007 to fit the
 *       data). When the merging time reaches zero, the satellite is assumed to
 *       merge with the central galaxy.
 *
 *
 *
 *       <B>deal_with_galaxy_merger</B> deals with the process, according to
 *       the mass fraction of the merger. Major if
 *       \f$M_{\rm{sat}}/M_{\rm{central}}>0.3\f$ and minor otherwise. It calls
 *       - add_galaxies_together - Add the cold and stellar phase of the merged
 *       galaxy to the central one. Also form a bulge at the central galaxy
 *       with the stars from the satellite in a minor merger if
 *       BulgeFormationInMinorMergersOn=1 (Major mergers are dealt later).
 *       - Then calls grow_black_hole - Grows black hole through accretion from
 *       cold gas during mergers (due to the instabilities triggered), as in
 *       Kauffmann & Haehnelt (2000). This is commonly referred as the quasar
 *       mode, main responsible for the black hole growth. After Croton2006 this
 *       mode is active even in minor mergers:
 *       \f$\Delta m_{\rm BH,Q}=M_{\rm{BH,min}}
 *          \frac{f_{\rm BH}(m_{\rm sat}/m_{\rm central})\,m_{\rm cold}}
 *           {1+(280\,\mathrm{km\,s}^{-1}/V_{\rm vir})^2}.\f$

 *       - Finally the burst of star formation due to the merger is treated.
 *           - If StarBurstModel = 0 (since Croton2006), the Somerville 2001
 *           model of bursts is used collisional_starburst_recipe(). The burst
 *           can happen for both major and minor mergers, with a fraction of
 *           the added cold gas from the satellite and central being consumed:
 *           \f$\dot{m}_{\star}^{\rm{burst}}
 *             = 0.56 \left(\frac{m_{\rm{sat}}}{m_{\rm{central}}}\right)^{0.7}
 *               m_{\rm{gas}}\f$.
 *           SN Feedback from starformation is computed and the sizes of bulge
 *           and disk followed.
 *
 *       - When a major merger occurs, the disk of both merging galaxies is
 *       completely destroyed to form a bulge. In either type of mergers, the
 *       bulge size is updated using Eq. 33 in Guo2010:
 *       \f$C\frac{GM^2_{\rm{new,bulge}}}{R_{\rm{new,bulge}}}=
 *          C\frac{GM^2_1}{R_1}+C\frac{GM^2_2}{R_2}+\alpha_{\rm{inter}}
 *          \frac{GM_1M_2}{R_1+R_2}\f$*/

/** @brief Calculates the central galaxies for type 1's. */

int set_merger_center(int fofhalo)
{
  /** @brief Get id of central galaxy, since type 1's
   *         can have their merger clock started before they become type 2's
   *         if M_star>M_vir. Introduced for millennium 2, where they can have
   *         very small masses, still be followed and never merge. At this
   *         moment the centre is still not known, reason why this function is
   *         needed. Also, if the type 1 merges, all the type 2's associated with
   *         it will need to know the "new" central galaxy they are merging into. */

  int prog, i, first_occupied, type, halonr, currentgal;
  double lenmax;
  i=0;

  halonr = fofhalo;

  //loop on all the haloes in current FOF group - to find a merger centre
  while(halonr >= 0)
  {
      lenmax = 0;
      first_occupied = Halo[halonr].FirstProgenitor;
      prog = Halo[halonr].FirstProgenitor;

      /* If the main progenitor of the current halo had no galaxies,
       * set first_occupied to the most massive progenitor. */
      if(prog >= 0)
	{
	  if(HaloAux[prog].NGalaxies == 0)
	    while(prog >= 0)
	      {
		for(i = 0, currentgal = HaloAux[prog].FirstGalaxy; i < HaloAux[prog].NGalaxies; i++)
		  {
		    type = HaloGal[currentgal].Type;

		    if(type == 0 || type == 1)
		      {
			if(Halo[prog].Len > lenmax)
			  {
			    lenmax = Halo[prog].Len;
			    first_occupied = prog;
			  }
		      }
		    currentgal = HaloGal[currentgal].NextGalaxy;
		  }
		prog = Halo[prog].NextProgenitor;
	      }
	}

      prog = Halo[halonr].FirstProgenitor;

      while(prog >= 0)//loop over all the progenitors
	{
    	  for(i = 0, currentgal = HaloAux[prog].FirstGalaxy; i < HaloAux[prog].NGalaxies; i++)//loop over all the galaxies in a given progenitor
    	    {
    	      type = HaloGal[currentgal].Type;

    	      if(type == 0 || type == 1) // the galaxy is a type 0 or 1?
    		if(prog == first_occupied) //is the main progenitor?
    		  if(halonr == Halo[halonr].FirstHaloInFOFgroup) //is the main halo?
    		    return currentgal;
    	      currentgal = HaloGal[currentgal].NextGalaxy;
    	    }
    	  prog = Halo[prog].NextProgenitor;
	}

      //if the halo has no galaxies, return 0
      if(i == 0)
	if(Halo[halonr].FirstHaloInFOFgroup == halonr)
	  return i;

      halonr = Halo[halonr].NextHaloInFOFgroup;
  }

  char sbuf[1000];
  sprintf(sbuf, "wrong in finding the central fof %d gal %d\n", fofhalo, currentgal);
  terminate(sbuf);
  return 0;
}

/** @brief Calculates the merging time whenever a galaxy becomes a satellite*/

double estimate_merging_time(int halonr, int mother_halonr, int p)
{
  int central_halonr;
  double coulomb, mergtime, SatelliteMass, SatelliteRadius, MotherHaloRvir;

  /** @brief Binney & Tremaine 1987 - 7.26 merging time for satellites due to
   *         dynamical friction. After Delucia2007 *2, shown to agree with
   *         Kolchin2008 simulations in Delucia2010. This is set when a galaxy
   *         becomes a type 2 or being a type 1 \f$M_{\rm{star}}>M_{\rm{vir}}\f$.
   *         In DeLucia2007 they could only merge into a type 0, now (after
   *         guo2010) they can merge into a type 1. */


  /*  recipe updated for more accurate merging time (see BT eq 7.26),
     now satellite radius at previous timestep is included */
  central_halonr = Halo[Halo[halonr].Descendant].FirstProgenitor;
  if(Gal[p].Type == 1)
    central_halonr=mother_halonr;
  if(central_halonr == halonr)
    {
      terminate("can't be...!\n");
    }


  coulomb = log(Halo[mother_halonr].Len / ((double) Halo[halonr].Len) + 1);

  /*  should include stellar+cold gas in SatelliteMass! */
  SatelliteMass = get_virial_mass(halonr)+(Gal[p].DiskMass+Gal[p].BulgeMass);

  SatelliteRadius = separation_halo(central_halonr,halonr)/(1 + ZZ[Halo[halonr].SnapNum]);

  int j;
  for (j = 0; j < 3; j++)
	Gal[p].DistanceToCentralGal[j] =  wrap(Halo[central_halonr].Pos[j] - Halo[halonr].Pos[j], BoxSize);


  MotherHaloRvir = get_virial_radius(mother_halonr);
  if(SatelliteRadius > MotherHaloRvir)
    SatelliteRadius = MotherHaloRvir;

  if(SatelliteMass > 0.0) {
    mergtime = 1.17 * SatelliteRadius * SatelliteRadius * get_virial_velocity(mother_halonr) /
               (coulomb * G * SatelliteMass); // Binney & Tremaine Eq.7.26

    /* change introduced by Delucia2007 to fit observations */
    mergtime = MergerTimeMultiplier*mergtime;
    //mergtime = 2.*mergtime;
  }
  else
    mergtime = -99999.9;

  return mergtime;

}

/** @brief Deals with all the physics triggered by mergers */

void deal_with_galaxy_merger(int p, int merger_centralgal, int centralgal, double time, double deltaT, int nstep)
{

/** @brief Deals with the physics triggered by mergers, according to the mass
 *         fraction of the merger \f$(M_{\rm{sat}}/M_{\rm{central}}><0.3)\f$.
 *         Add the cold and stellar phase of the satellite galaxy to the central
 *         one, form a bulge at the central galaxy with the stars from the
 *         satellite in a minor merger if BulgeFormationInMinorMergersOn=1.
 *         Grows black hole through accretion from cold gas "quasar mode".
 *         If StarBurstModel = 0, the Somerville 2001 model
 *         of bursts is used, SN Feedback from starformation is computed and
 *         the sizes of bulge and disk followed. When a major merger occurs,
 *         the disk of both merging galaxies is completely destroyed to form
 *         a bulge. New stars form of to the bulge*/

  double mi, ma, mass_ratio, Mcstar, Mcgas, Mcbulge, Mpstar, Mpgas,Mpbulge;
  double frac;
#ifdef GALAXYTREE
  int q;

  mass_checks("deal_with_galaxy_merger #0",p);
  mass_checks("deal_with_galaxy_merger #0",merger_centralgal);
  mass_checks("deal_with_galaxy_merger #0",centralgal);

  q = Gal[merger_centralgal].FirstProgGal;
  if(q >= 0)
    {
      while(GalTree[q].NextProgGal >= 0)
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
    terminate("may not happen");

  q = GalTree[q].NextProgGal;

  if(HaloGal[GalTree[q].HaloGalIndex].GalTreeIndex != q)
    terminate("inconsistency");

  HaloGal[GalTree[q].HaloGalIndex].MergeOn = 2;

  if(Gal[p].Type == 1)
    HaloGal[GalTree[q].HaloGalIndex].MergeOn = 3;
#endif


  /* flag galaxy as finished */
  Gal[p].Type = 3;

  /*  calculate mass ratio of merging galaxies */
  mi = Gal[p].DiskMass+Gal[p].BulgeMass+Gal[p].ColdGas;
  ma = Gal[merger_centralgal].DiskMass+Gal[merger_centralgal].BulgeMass+Gal[merger_centralgal].ColdGas;
  if(max(mi,ma) > 0.)
    mass_ratio = min(mi,ma) / max(mi,ma);
  else
    mass_ratio = 1.0;

  /* record the gas and stellar component  mass of merger central and satellite
   * galaxies the before merger */
  Mcstar=(Gal[merger_centralgal].DiskMass+Gal[merger_centralgal].BulgeMass);
  Mcbulge=Gal[merger_centralgal].BulgeMass;
  Mcgas=Gal[merger_centralgal].ColdGas;
  Mpstar=(Gal[p].DiskMass+Gal[p].BulgeMass);
  Mpbulge=Gal[p].BulgeMass;
  Mpgas=Gal[p].ColdGas;

  mass_checks("deal_with_galaxy_merger #1",p);
  mass_checks("deal_with_galaxy_merger #1",merger_centralgal);
  mass_checks("deal_with_galaxy_merger #1",centralgal);



  /* Add the cold and stellar phase of the merged galaxy to the central one.
     Also form a bulge if BulgeFormationInMinorMergersOn is set on (transfer stars
     from satellite disk to central bulge). In a major merger (dealt at the
     make_bulge_from_burst) the disk of the central (now made up of central and
     satellite will be moved to the bulge). Any new stars formed will go to the bulge */

  add_galaxies_together(merger_centralgal, p);

  mass_checks("deal_with_galaxy_merger #2",p);
  mass_checks("deal_with_galaxy_merger #2",merger_centralgal);
  mass_checks("deal_with_galaxy_merger #2",centralgal);

  /* grow black hole through accretion from cold disk during mergers, as in
   * Kauffmann & Haehnelt (2000) + minor mergers - Quasar Mode */
  if(AGNRadioModeModel != 5)
    grow_black_hole(merger_centralgal, mass_ratio, deltaT);

  mass_checks("deal_with_galaxy_merger #3",p);
  mass_checks("deal_with_galaxy_merger #3",merger_centralgal);
  mass_checks("deal_with_galaxy_merger #3",centralgal);

  if (StarBurstModel == 0)
    {
    /* Starburst as in Somerville 2001, with feedback computed inside. */
  	/* All star formation happens in the disk, but in a major merger this will then
  	 * be destroyed with everything moved to the bulge. */
      frac=collisional_starburst_recipe(mass_ratio, merger_centralgal, centralgal, time, deltaT);
      bulgesize_from_merger(mass_ratio,merger_centralgal,p,
			    Mcstar, Mcbulge, Mcgas, Mpstar, Mpbulge, Mpgas, frac);

      mass_checks("deal_with_galaxy_merger #3.5",p);
      mass_checks("deal_with_galaxy_merger #3.5",merger_centralgal);
      mass_checks("deal_with_galaxy_merger #3.5",centralgal);

      if(mass_ratio > ThreshMajorMerger)
	make_bulge_from_burst(merger_centralgal);
    }

  mass_checks("deal_with_galaxy_merger #4",p);
  mass_checks("deal_with_galaxy_merger #4",merger_centralgal);
  mass_checks("deal_with_galaxy_merger #4",centralgal);

  /* If we are in the presence of a minor merger, check disk stability (the disk
   * is completely destroyed in major mergers)*/
  if(mass_ratio < ThreshMajorMerger && (Gal[merger_centralgal].DiskMass+Gal[merger_centralgal].BulgeMass) > 0.0)
    check_disk_instability(merger_centralgal);

  if ((Gal[merger_centralgal].BulgeMass > 1.e-6 && Gal[merger_centralgal].BulgeSize == 0.0) ||
      (Gal[merger_centralgal].BulgeMass == 0.0 && Gal[merger_centralgal].BulgeSize >1.e-6)) {
  	char sbuf[1000];
  	sprintf(sbuf, "central: stellarmass %f, bulgemass %f, bulgesize %f, coldgas %f,gasdisk %f,stellardisk %f \n",
  			(Gal[merger_centralgal].DiskMass+Gal[merger_centralgal].BulgeMass),Gal[merger_centralgal].BulgeMass,
  			Gal[merger_centralgal].BulgeSize,Gal[merger_centralgal].ColdGas,Gal[merger_centralgal].GasDiskRadius,
  			Gal[merger_centralgal].StellarDiskRadius);
  	terminate(sbuf);
  }

  if (DiskRadiusModel == 0) {
    get_gas_disk_radius(merger_centralgal);
    get_stellar_disk_radius(merger_centralgal);
  }

  mass_checks("deal_with_galaxy_merger #5",p);
  mass_checks("deal_with_galaxy_merger #5",merger_centralgal);
  mass_checks("deal_with_galaxy_merger #5",centralgal);

}


/** @brief Grows black holes, through accretion from cold gas during mergers,
 *          as in Kauffmann & Haehnelt (2000) - Quasar Mode.  */

void grow_black_hole(int merger_centralgal, double mass_ratio, double deltaT)
{
  double BHaccrete, fraction;

  /** @brief Grows black hole through accretion from cold gas during mergers,
   *         as in Kauffmann & Haehnelt (2000).
   *                              in main.c */

  if(Gal[merger_centralgal].ColdGas > 0.0)
    {
    BHaccrete = BlackHoleGrowthRate * mass_ratio
      / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[merger_centralgal].Vvir))) * Gal[merger_centralgal].ColdGas;
    /* redshift dependent accretion, not published */
    /* BHaccrete = BlackHoleGrowthRate * (1.0 + ZZ[Halo[halonr].SnapNum]) * mass_ratio */

    /* cannot accrete more gas than is available! */
    if(BHaccrete > Gal[merger_centralgal].ColdGas)
      BHaccrete = Gal[merger_centralgal].ColdGas;

    fraction=BHaccrete/Gal[merger_centralgal].ColdGas;

    Gal[merger_centralgal].BlackHoleMass += BHaccrete;
    Gal[merger_centralgal].QuasarAccretionRate += BHaccrete / deltaT;

    Gal[merger_centralgal].ColdGas -= BHaccrete;
    Gal[merger_centralgal].MetalsColdGas=
      metals_add(Gal[merger_centralgal].MetalsColdGas, Gal[merger_centralgal].MetalsColdGas,-fraction);

  }
}

/** @brief Adds all the components of the satellite galaxy into its
 *         central companion. */

void add_galaxies_together(int t, int p)
{
  /** @brief All the components of the satellite galaxy are added to the
   *         correspondent component of the central galaxy. Cold gas spin
   *         is updated and a bulge is formed at the central galaxy, with
   *         the stars of the satellite if  BulgeFormationInMinorMergersOn=1.
   *         In case of a major merger, everything that was put in the disk of
   *         the central galaxy will be moved into the bulge
   */
  int outputbin, j;
  float tspin[3],tmass,pmass;

  /* t central, p satellite */

  mass_checks("add_galaxies_together #0",p);
  mass_checks("add_galaxies_together #0.1",t);

  /* angular momentum transfer between gas*/
  tmass= Gal[t].ColdGas;
  pmass= Gal[p].ColdGas;

  Gal[t].MergeSat +=(Gal[p].DiskMass+Gal[p].BulgeMass);
  Gal[p].MergeSat=0.;

  transfer_gas(t,"Cold",p,"Cold",1.,"add_galaxies_together", __LINE__);
  //transfer_gas(t,"Ejected",p,"Cold",1.,"add_galaxies_together", __LINE__);
  transfer_gas(t,"Hot",p,"Hot",1.,"add_galaxies_together", __LINE__);
  transfer_gas(t,"Ejected",p,"Ejected",1.,"add_galaxies_together", __LINE__);
#ifdef TRACK_BURST
    /* The whole burst component gets transferred */
  transfer_stars(t,"Burst",p,"Burst",1.);
#endif
  if(BulgeFormationInMinorMergersOn)
    transfer_stars(t,"Bulge",p,"Disk",1.);
  else
    transfer_stars(t,"Disk",p,"Disk",1.);
  transfer_stars(t,"Bulge",p,"Bulge",1.);
  transfer_stars(t,"ICM",p,"ICM",1.);

  Gal[t].BlackHoleMass += Gal[p].BlackHoleMass;
  Gal[p].BlackHoleMass=0.;
  Gal[t].StarMerge += Gal[p].StarMerge;
  Gal[p].StarMerge=0.;

  mass_checks("add_galaxies_together #1",p);
  mass_checks("add_galaxies_together #1.1",t);

  /*update the gas spin*/
  for(j=0;j<3;j++)
    tspin[j]=Gal[t].GasSpin[j]*tmass+Gal[t].HaloSpin[j]*pmass;
  if (Gal[t].ColdGas != 0)
    for (j=0;j<3;j++)
     Gal[t].GasSpin[j]=tspin[j]/(Gal[t].ColdGas);

  Gal[t].Sfr += Gal[p].Sfr;

  if(BulgeFormationInMinorMergersOn)
    Gal[t].SfrBulge += Gal[p].Sfr;

  Gal[t].QuasarAccretionRate += Gal[p].QuasarAccretionRate;
  Gal[t].RadioAccretionRate += Gal[p].RadioAccretionRate;

  for(outputbin = 0; outputbin < NOUT; outputbin++)
  	 Gal[t].MassWeightAge[outputbin] += Gal[p].MassWeightAge[outputbin];

#ifndef  POST_PROCESS_MAGS

/* Add the luminosities of the satellite and central galaxy */
#ifdef OUTPUT_REST_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++)
    {

    for(j = 0; j < NMAG; j++)
    {
      Gal[t].Lum[j][outputbin] += Gal[p].Lum[j][outputbin];
      Gal[t].YLum[j][outputbin] += Gal[p].YLum[j][outputbin];
#ifdef ICL
      Gal[t].ICLLum[j][outputbin] += Gal[p].ICLLum[j][outputbin];
#endif
    }
    if(BulgeFormationInMinorMergersOn)
    {
    	for(j = 0; j < NMAG; j++)
    	{
    		Gal[t].LumBulge[j][outputbin] += Gal[p].Lum[j][outputbin];
    		Gal[t].YLumBulge[j][outputbin] += Gal[p].YLum[j][outputbin];
      }
    }
    else
    {
      for(j = 0; j < NMAG; j++)
      {
      	Gal[t].LumBulge[j][outputbin]  += Gal[p].LumBulge[j][outputbin];
      	Gal[t].YLumBulge[j][outputbin] += Gal[p].YLumBulge[j][outputbin];
      }
    }
    }
#endif // OUTPUT_REST_MAGS

#ifdef COMPUTE_OBS_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++)
    {
      for(j = 0; j < NMAG; j++) {
      Gal[t].ObsLum[j][outputbin]   += Gal[p].ObsLum[j][outputbin];
      Gal[t].ObsYLum[j][outputbin]  += Gal[p].ObsYLum[j][outputbin];
#ifdef ICL
      Gal[t].ObsICL[j][outputbin]  += Gal[p].ObsICL[j][outputbin];
#endif

#ifdef OUTPUT_MOMAF_INPUTS
      Gal[t].dObsLum[j][outputbin] += Gal[p].dObsLum[j][outputbin];
      Gal[t].dObsYLum[j][outputbin] += Gal[p].dObsYLum[j][outputbin];
#ifdef ICL
      Gal[t].dObsICL[j][outputbin]  += Gal[p].dObsICL[j][outputbin];
#endif
#endif
    }
    if(BulgeFormationInMinorMergersOn) {
      for(j = 0; j < NMAG; j++) {
	Gal[t].ObsLumBulge[j][outputbin]   += Gal[p].ObsLum[j][outputbin];
	Gal[t].ObsYLumBulge[j][outputbin]  += Gal[p].ObsYLum[j][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	Gal[t].dObsLumBulge[j][outputbin]  += Gal[p].dObsLum[j][outputbin];
	Gal[t].dObsYLumBulge[j][outputbin] += Gal[p].dObsYLum[j][outputbin];
#endif
      }
    }
    else
    {
      for(j = 0; j < NMAG; j++) {
	Gal[t].ObsLumBulge[j][outputbin]   += Gal[p].ObsLumBulge[j][outputbin];
	Gal[t].ObsYLumBulge[j][outputbin]  += Gal[p].ObsYLumBulge[j][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
	Gal[t].dObsLumBulge[j][outputbin]  += Gal[p].dObsLumBulge[j][outputbin];
	Gal[t].dObsYLumBulge[j][outputbin] += Gal[p].dObsYLumBulge[j][outputbin];
#endif
      }
    }
  }
#endif //COMPUTE_OBS_MAGS
#endif  //POST_PROCESS_MAGS
}


/** @brief In a major merger, both disks are destroyed and all the mass transferred
 *         to the bulge. The galaxies have already been merged, so all we need to do here
 *         is transfer stars from disk to bulge. */

void make_bulge_from_burst(int p)
{
	int outputbin;
  /* generate bulge */
  transfer_stars(p,"Bulge",p,"Disk",1.);

  /*  update the star formation rate */
  Gal[p].SfrBulge  = Gal[p].Sfr;

#ifndef  POST_PROCESS_MAGS
  int j;
#ifdef OUTPUT_REST_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++) {
    for(j = 0; j < NMAG; j++) {
      Gal[p].LumBulge[j][outputbin]  = Gal[p].Lum[j][outputbin];
      Gal[p].YLumBulge[j][outputbin] = Gal[p].YLum[j][outputbin];
    }
  }
#endif
#ifdef COMPUTE_OBS_MAGS
  for(outputbin = 0; outputbin < NOUT; outputbin++) {
    for(j = 0; j < NMAG; j++) {
      Gal[p].ObsLumBulge[j][outputbin]   = Gal[p].ObsLum[j][outputbin];
      Gal[p].ObsYLumBulge[j][outputbin]  = Gal[p].ObsYLum[j][outputbin];
#ifdef OUTPUT_MOMAF_INPUTS
      Gal[p].dObsLumBulge[j][outputbin]  = Gal[p].dObsLum[j][outputbin];
      Gal[p].dObsYLumBulge[j][outputbin] = Gal[p].dObsYLum[j][outputbin];
#endif
    }
  }
#endif
#endif //POST_PROCESS_MAGS
}

/** @brief Merger burst recipe from Somerville 2001 (used after Croton2006) */

double collisional_starburst_recipe(double mass_ratio, int merger_centralgal, int centralgal,
				  double time, double deltaT)
{
  /** @brief If StarBurstModel = 1 (since Croton2006), the Somerville 2001
   *         model of bursts is used. The burst can happen for both major
   *         and minor mergers, with a fraction of the added cold gas from
   *         the satellite and central being consumed. SN Feedback from
   *         starformation is computed and the sizes of bulge and disk
   *         followed (not done for the other burst mode).*/

  double mstars, metallicitySF, eburst, Ggas;

  /* This is the major and minor merger starburst recipe of Somerville 2001.
   * The coefficients in eburst are taken from TJ Cox's PhD thesis and should
   * be more accurate then previous. */

  Ggas=Gal[merger_centralgal].ColdGas;

  /* the bursting fraction given the mass ratio */
  /* m_dot = 0.56*(m_sat/m_central)^0.7*m_gas */
  eburst = SfrBurstEfficiency * pow(mass_ratio, SfrBurstSlope);
  //eburst = 0.56 * pow(mass_ratio, 0.7);
  mstars = eburst * Gal[merger_centralgal].ColdGas;
  if(mstars < 0.0)
    mstars = 0.0;

  /*  update the star formation rate */
  Gal[merger_centralgal].Sfr += mstars / deltaT;

  /* Store the value of the metallicity of the cold phase when SF occurs.
   * Used to update luminosities below */
  metallicitySF = metals_total(Gal[merger_centralgal].MetalsColdGas)/Gal[merger_centralgal].ColdGas;

  if (mstars > 0.)
  	update_stars_due_to_reheat(merger_centralgal, centralgal, &mstars);

  // update_from_star_formation can only be called
  // after SD_feeedback recipe since stars need to be re_set once the reheated mass is known
  // (star formation and feedback share the same fraction of cold gas)
  int nstep=-1;
  if (mstars > 0.)
    update_from_star_formation(merger_centralgal, mstars, true, nstep); // true indicates starburst

	mass_checks("collisional_starburst_recipe #2",merger_centralgal);

  update_massweightage(merger_centralgal, mstars, time);

  if (mstars > 0.)
  	SN_feedback(merger_centralgal, centralgal, mstars, "ColdGas");


#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  /*  update the luminosities due to the stars formed */
  if (mstars > 0.0)
    add_to_luminosities(merger_centralgal, mstars, time, deltaT / STEPS, metallicitySF);
#endif //NDEF POST_PROCESS_MAGS
#endif

  if (Ggas > 0.)
    return mstars/Ggas;
  else
    return 0.0;

}

/** @brief Calculates the bulge size after a merger. */
void  bulgesize_from_merger(double mass_ratio,int merger_centralgal,int p,
			    double Mcstar,double Mcbulge,double Mcgas,
			    double Mpstar,double Mpbulge,double Mpgas,double frac)
{
  /** @brief For any type of merger calculates the new bulge size using
   *         Eq. 33 in Guo2010:
   *
   *         \f$C\frac{GM^2_{\rm{new,bulge}}}{R_{\rm{new,bulge}}}=
   *             C\frac{GM^2_1}{R_1}+C\frac{GM^2_2}{R_2}+		\
   *             alpha_{\rm{inter}}\frac{GM_1M_2}{R_1+R_2}\f$.
   *
   *         This implementation assumed that the new bulge occupies the
   *         same space as the components that formed it. */

  double Mc,Rc;
  double Mp,Rp;
  double fint,c;

  fint=0.5;
  c=0.5;

  /* calculate radius for the object that will form the new bulge - Rc and Rp */
  /* Minor Merger */
  if(mass_ratio < ThreshMajorMerger) {
    /* In a minor merger only the stars of the satellite galaxy are moved
     * to the bulge of the central galaxy, therefore only stellar
     * components are used to compute radius and masses. */
    frac=0.0;
    /* in a minor merger only consider the bulge mass of the central galaxy */
    Mc=Mcbulge;
    Rc=Gal[merger_centralgal].BulgeSize;
    /* and stellarmass of satellite*/
    Mp=Mpstar;
    if (Mp >0.0)
      Rp=(Gal[p].StellarDiskRadius/3.*1.68*(Mpstar-Mpbulge)+Gal[p].BulgeSize*Mpbulge)/Mpstar;
    else
      Rp=0.0;
  }
  /* Major Merger */
  else {
    /* on a major merger both stellar and gas (after a burst) components
     * from the final bulge and need to be considered */
    /* Mc = bulge mass + burst of central*/
    Mc=Mcstar+frac*Mcgas;
    if (Mc > 0.0)
      Rc=(Gal[merger_centralgal].StellarDiskRadius/3.*1.68*(Mcstar-Mcbulge)+Gal[merger_centralgal].BulgeSize*Mcbulge+Gal[merger_centralgal].GasDiskRadius*frac*Mcgas/3.*1.68)/(Mcgas*frac+Mcstar);
    else
      Rc=0.0;
    /* and satellite Mp */
    Mp=Mpstar+frac*Mpgas;
    if (Mp > 0.0)
      Rp=(Gal[p].StellarDiskRadius/3.*1.68*(Mpstar-Mpbulge)+Gal[p].BulgeSize*Mpbulge+Gal[p].GasDiskRadius*frac*Mpgas/3.*1.68)/(Mpgas*frac+Mpstar);
    else
      Rp=0.0;
  }

  if(Rp>0. && Rp<1.e-8)
  	Rp=1.e-8;
  if(Rc>0. && Rc<1.e-8)
    	Rc=1.e-8;
  /* If both original radius are bigger then 0 then this is Eq. 33 in Guo 2010
   * solved for R_new,bulge with all terms divided by G and C. */
  if(Rc >= 1.e-8 && Rp >= 1.e-8)
  	Gal[merger_centralgal].BulgeSize=(Mp+Mc)*(Mp+Mc)/(Mp*Mp/Rp+Mc*Mc/Rc+fint/c*Mp*Mc/(Rp+Rc));

  if(Rc >= 1.e-8 && Rp <= 1.e-8)
    Gal[merger_centralgal].BulgeSize=(Mp+Mc)*(Mp+Mc)/(Mc*Mc/Rc+fint/c*Mp*Mc/(Rp+Rc));

  if(Rc <= 1.e-8 && Rp <=1.e-8)
    Gal[merger_centralgal].BulgeSize=0.0;

  if(Rc <= 1.e-8 && Rp >= 1.e-8)
    Gal[merger_centralgal].BulgeSize=(Mp+Mc)*(Mp+Mc)/(Mp*Mp/Rp+fint/c*Mp*Mc/(Rp+Rc));


  if ((Mp+Mc > 0.0 && Gal[merger_centralgal].BulgeSize == 0.0 )||(Mp+Mc == 0.0 && Gal[merger_centralgal].BulgeSize> 0.0))
    {
      char sbuf[1000];
      printf("halonr=%d, merger_centralgal %d\n\n", Gal[merger_centralgal].HaloNr, merger_centralgal);

      printf("New Bulge Mass from Central (MC)=%e\n New Bulge Mass from Satellite (Mp)=%e\n NewBulge size=%e\n\n",
	     Mp,Mc,Gal[merger_centralgal].BulgeSize);

      if(mass_ratio < ThreshMajorMerger)
	{
	  printf("minor merger, new mass=original mass\n");
	  printf("New Bulge Mass From Central (Mc)   = Mcbulge = %f\n",Mcbulge);
	  printf("New Bulge Mass From Satellite (Mp) = Mpstar  = %f\n",Mpstar);
	}
      else
	{
	  printf("New Bulge From Central (Mc)   = Mcstar+frac*Mcgas = %f+%f*%f\n", Mcstar, frac, Mcgas);
	  printf("New Bulge From Satellite (Mp) = Mpstar+frac*Mpgas = %f+%f*%f\n", Mpstar, frac, Mpgas);
	}

      printf("BulgeSize from Central (Rc)=%e\nBulgeSize from Satellite (Rp)=%e\nmass ratio=%f\n\n",Rc,Rp, mass_ratio);

      printf("the following masses don't tell a lot because things have been merged already!!!\n");
      printf("    sat: BulgeMass=%0.7f, BulgeSize=%0.7f, GasMass=%0.7f, GasSize=%0.7f, DiskMass=%0.7f StellarSize=%0.7f \n",
	     Gal[p].BulgeMass, Gal[p].BulgeSize, Gal[p].ColdGas, Gal[p].GasDiskRadius,
	     Gal[p].DiskMass, Gal[p].StellarDiskRadius);
      printf(	"central: BulgeMass=%0.7f, BulgeSize=%0.7f, GasMass=%0.7f, GasSize=%0.7f, DiskMass=%0.7f StellarSize=%0.7f \n",
		Gal[merger_centralgal].BulgeMass,  Gal[merger_centralgal].BulgeSize, Gal[merger_centralgal].ColdGas,
		Gal[merger_centralgal].GasDiskRadius, Gal[merger_centralgal].DiskMass, Gal[merger_centralgal].StellarDiskRadius);

      sprintf(sbuf,"\n bulgesize wrong in merger");
      terminate(sbuf);
      exit(0);
  }

}
