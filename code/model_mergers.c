#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** @file recipe_mergers.c
 *  @brief Calculates the merging time, the central galaxy (for type 1's),
 *         adds galaxies together, calculates SF from bursts and grows
 *         black holes.
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
 *       MERGE01 option, that sets the clock also for type 1's (satellites with
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

void deal_with_galaxy_merger(int p, double time, double deltaT, int nstep)
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

    int FOF_centralgal, merger_centralgal;
    double mi, ma, mass_ratio, MgasCentral, MstarCentral, MbulgeCentral, MgasSat, MstarSat, MbulgeSat;
    double RgasCentral, RStellarDiskCentral, RgasSat, RStellarDiskSat;
    double frac;

    FOF_centralgal=Gal[p].FOFCentralGal;
    merger_centralgal=Gal[p].MergerCentralGal;

    mass_checks(p,"model_mergers.c",__LINE__);
    mass_checks(merger_centralgal,"model_mergers.c",__LINE__);
    mass_checks(FOF_centralgal,"model_mergers.c",__LINE__);

#ifdef GALAXYTREE
    int q;

    q = Gal[merger_centralgal].FirstProgGal;
    if(q >= 0) 
    {
		while(GalTree[q].NextProgGal >= 0) q = GalTree[q].NextProgGal;
		GalTree[q].NextProgGal = Gal[p].FirstProgGal;
		if(GalTree[q].NextProgGal >= NGalTree) 
		{
	    	printf("q=%d p=%d GalTree[q].NextProgGal=%d NGalTree=%d\n",
		  	 q, p, GalTree[q].NextProgGal, NGalTree);
	  	    terminate("problem");
		}
    }
    if(q < 0)
	terminate("q<0\n");
    q = GalTree[q].NextProgGal;

    if(HaloGal[GalTree[q].HaloGalIndex].GalTreeIndex != q)
	terminate("inconsistency");

    HaloGal[GalTree[q].HaloGalIndex].MergeOn = 2;

    if(Gal[p].Type == 1)
	HaloGal[GalTree[q].HaloGalIndex].MergeOn = 3;
#endif


  /*  calculate mass ratio of merging galaxies */
  mi = Gal[p].DiskMass+Gal[p].BulgeMass+Gal[p].ColdGas;
  ma = Gal[merger_centralgal].DiskMass+Gal[merger_centralgal].BulgeMass+Gal[merger_centralgal].ColdGas;
  if(max(mi,ma) > 0.)
    mass_ratio = min(mi,ma) / max(mi,ma);
  else
    mass_ratio = 1.0;

  /* record the gas and stellar component  mass of merger central and satellite
   * galaxies the before merger */
  MstarCentral=(Gal[merger_centralgal].DiskMass+Gal[merger_centralgal].BulgeMass);
  MbulgeCentral=Gal[merger_centralgal].BulgeMass;
  MgasCentral=Gal[merger_centralgal].ColdGas;
  MstarSat=(Gal[p].DiskMass+Gal[p].BulgeMass);
  MbulgeSat=Gal[p].BulgeMass;
  MgasSat=Gal[p].ColdGas;

  mass_checks(p,"model_mergers.c",__LINE__);
  mass_checks(merger_centralgal,"model_mergers.c",__LINE__);

  //compute the sizes before galaxies are merged to use in the calculation
  //of the resulting bulge size (the calculation is only needed here because
  //with H2_AND_RINGS the sizes of stellar and gaseous disks are computed
  //on the fly)
  RgasCentral=get_gas_disk_radius(merger_centralgal)/3.;
  RStellarDiskCentral=get_stellar_disk_radius(merger_centralgal)/3.;
  RgasSat=get_gas_disk_radius(p)/3.;
  RStellarDiskSat=get_stellar_disk_radius(p)/3.;


  /* Add the cold and stellar phase of the merged galaxy to the central one.
     Also form a bulge if BulgeFormationInMinorMergersOn is set on (transfer stars
     from satellite disk to central bulge). In a major merger (dealt at the
     make_bulge_from_burst) the disk of the central (now made up of central and
     satellite will be moved to the bulge). Any new stars formed will go to the bulge */

  add_galaxies_together(merger_centralgal, p, deltaT);

  mass_checks(p,"model_mergers.c",__LINE__);
  mass_checks(merger_centralgal,"model_mergers.c",__LINE__);

  /* grow black hole through accretion from cold disk during mergers, as in
   * Kauffmann & Haehnelt (2000) + minor mergers - Quasar Mode */
  if(AGNRadioModeModel <4)
    grow_black_hole(merger_centralgal, mass_ratio, deltaT);

  mass_checks(p,"model_mergers.c",__LINE__);
  mass_checks(merger_centralgal,"model_mergers.c",__LINE__);

  /* Starburst as in Somerville 2001, with feedback computed inside.
   * All star formation happens in the disk, but in a major merger this will then
   * be destroyed with everything moved to the bulge. */
  if (StarBurstModel == 0)
    {
      frac=collisional_starburst_recipe(mass_ratio, merger_centralgal, FOF_centralgal, time, deltaT);
      bulgesize_from_merger(mass_ratio,merger_centralgal,p, MstarCentral, MbulgeCentral, MgasCentral,
      			    MstarSat, MbulgeSat, MgasSat, frac, RgasCentral, RStellarDiskCentral, RgasSat, RStellarDiskSat);

      mass_checks(p,"model_mergers.c",__LINE__);
      mass_checks(merger_centralgal,"model_mergers.c",__LINE__);
      mass_checks(FOF_centralgal,"model_mergers.c",__LINE__);

      if(mass_ratio > ThreshMajorMerger)
	make_bulge_from_burst(merger_centralgal);
    }

#ifdef H2_AND_RINGS
  //Bulge Mass was added into the same place as the disk, it will now be redistributed
  //according to a Jaffe profile and after the new bulge size has been calculated
  distribute_bulge_material(merger_centralgal);
#endif

  mass_checks(p,"model_mergers.c",__LINE__);
  mass_checks(merger_centralgal,"model_mergers.c",__LINE__);
  mass_checks(FOF_centralgal,"model_mergers.c",__LINE__);

  /* If we are in the presence of a minor merger, check disk stability (the disk
   * is completely destroyed in major mergers)*/
  if(DiskInstabilityModel==0)
    {
    //if(mass_ratio < ThreshMajorMerger && (Gal[merger_centralgal].ColdGas) > 0.0)
    //  check_disk_instability_gas(merger_centralgal, deltaT/STEPS);
    if(mass_ratio < ThreshMajorMerger && (Gal[merger_centralgal].DiskMass+Gal[merger_centralgal].BulgeMass) > 0.0)
         check_disk_instability(merger_centralgal, deltaT/STEPS);
    }

  mass_checks(merger_centralgal,"model_mergers.c",__LINE__);

  /* Not supported option to shrink bulge sizes in gas rich mergers */
#ifdef SHRINK_IN_GAS_RICH_MERGER
  if (MgasCentral+MstarCentral+MgasSat+MstarSat > 1.e-8 ) {
    Gal[merger_centralgal].BulgeSize /= 1+pow((MgasCentral+MgasSat)/(MgasCentral+MstarCentral+MgasSat+MstarSat)/0.15,1.);
    // Don;t let bulge size shrink below 100 pc/h.
    if (Gal[merger_centralgal].BulgeSize < 1.e-4)
      Gal[merger_centralgal].BulgeSize = 1.e-4;
  }
#endif

 if ((Gal[merger_centralgal].BulgeMass > 1.e-6 && Gal[merger_centralgal].BulgeSize == 0.0) ||
      (Gal[merger_centralgal].BulgeMass == 0.0 && Gal[merger_centralgal].BulgeSize >1.e-6)) {
  	char sbuf[1000];
  	sprintf(sbuf, "2 central: stellarmass %0.10f, bulgemass %0.10f, bulgesize %0.10f, stellardisksize %0.10f \n",
  			(Gal[merger_centralgal].DiskMass+Gal[merger_centralgal].BulgeMass),Gal[merger_centralgal].BulgeMass,
  			Gal[merger_centralgal].BulgeSize, Gal[merger_centralgal].DiskRadius);
  	terminate(sbuf);
  }

  if (DiskRadiusModel == 0)
    {
      Gal[merger_centralgal].ColdGasRadius=get_gas_disk_radius(merger_centralgal);
      Gal[merger_centralgal].DiskRadius = get_stellar_disk_radius(merger_centralgal);
      // The following may not be necessary - depnds what happens in collisonal_starburst_recipe
      if (merger_centralgal != FOF_centralgal) 
	{
	  Gal[FOF_centralgal].ColdGasRadius=get_gas_disk_radius(FOF_centralgal);
	  Gal[FOF_centralgal].DiskRadius = get_stellar_disk_radius(FOF_centralgal);
	}
    }

  /* flag galaxy as finished */
  Gal[p].Type = 3;

  mass_checks(p,"model_mergers.c",__LINE__);
  mass_checks(merger_centralgal,"model_mergers.c",__LINE__);

}


/** @brief Grows black holes, through accretion from cold gas during mergers,
 *          as in Kauffmann & Haehnelt (2000) - Quasar Mode. */

void grow_black_hole(int merger_centralgal, double mass_ratio, double deltaT)
{
  double BHaccrete, fraction;
#ifdef H2_AND_RINGS
  int jj;
  double fractionRings[RNUM];
#endif
  /** @brief Grows black hole through accretion from cold gas during mergers,
   *         as in Kauffmann & Haehnelt (2000). In addition, black holes can grow
   *         during minor mergers.
   *         BlackHoleGrowth == 0 gives instantaneous accretion onto the black hole;
   *         BlackHoleGrowth == 1 instead feeds an accretion disk: accretion occurs
   *         in main.c */

  if(Gal[merger_centralgal].ColdGas > 0.0)
    {
      BHaccrete = BlackHoleGrowthRate * mass_ratio
	  / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[merger_centralgal].Vvir))) * Gal[merger_centralgal].ColdGas;
     // BHaccrete = BlackHoleGrowthRate * mass_ratio * Gal[merger_centralgal].ColdGas;

      // redshift dependent accretion, not published
      // BHaccrete = BlackHoleGrowthRate * (1.0 + ZZ[Halo[halonr].SnapNum]) * mass_ratio

      // cannot accrete more gas than is available
      if(BHaccrete > Gal[merger_centralgal].ColdGas)
	BHaccrete = Gal[merger_centralgal].ColdGas;

      fraction=BHaccrete/Gal[merger_centralgal].ColdGas;


     // print_rings("first call", merger_centralgal);
#ifdef H2_AND_RINGS
      for(jj=0;jj<RNUM;jj++)
       	  fractionRings[jj]=fraction;
     /* double mass_to_transfer;
      mass_to_transfer=BHaccrete;
      for(jj=0;jj<RNUM;jj++)
     	  fractionRings[jj]=0.;
      //transfer from inner circles first
      for(jj=0;jj<RNUM;jj++)
	{RStellarDiskSat
	  if(mass_to_transfer<Gal[merger_centralgal].ColdGasRings[jj])
	    {
	      fractionRings[jj]=mass_to_transfer/Gal[merger_centralgal].ColdGasRings[jj];
	      break;
	    }
	  else
	    {
	      fractionRings[jj]=1.;
	      mass_to_transfer-=Gal[merger_centralgal].ColdGasRings[jj];
	    }

	}*/
      //doesn't matter the destination, since there are no rings in BlackHoleMass or BlackHoleGas
      //but maybe there will be in the future
      if (BlackHoleGrowth == 0)
	{
	  transfer_material_with_rings(merger_centralgal,"BlackHoleMass",merger_centralgal,"ColdGas", fractionRings,"model_mergers.c", __LINE__);
	}
      else if (BlackHoleGrowth == 1)
	transfer_material_with_rings(merger_centralgal,"BlackHoleGas",merger_centralgal,"ColdGas", fractionRings,"model_mergers.c", __LINE__);
#else

      if (BlackHoleGrowth == 0)
      	{
	  transfer_material(merger_centralgal,"BlackHoleMass",merger_centralgal,"ColdGas",fraction,"model_mergers.c", __LINE__);
      	  //Gal[merger_centralgal].QuasarAccretionRate += BHaccrete / deltaT;
      	}
      else if (BlackHoleGrowth == 1)
	  transfer_material(merger_centralgal,"BlackHoleGas",merger_centralgal,"ColdGas",fraction,"model_mergers.c", __LINE__);
#endif //H2_AND_RINGS

      Gal[merger_centralgal].QuasarAccretionRate += BHaccrete / deltaT;
    }
}



/** @brief Adds all the components of the satellite galaxy into its
 *         central companion. */

void add_galaxies_together(int t, int p, double deltaT)
{
  /** @brief All the components of the satellite galaxy are added to the
   *         correspondent component of the central galaxy. Cold gas spin
   *         is updated and a bulge is formed at the central galaxy, with
   *         the stars of the satellite if  BulgeFormationInMinorMergersOn=1.
   *         In case of a major merger, everything that was put in the disk of
   *         the central galaxy will be moved into the bulge
   *
   * TODO Even though galaxy p has been set to type 3 (ie a non-galaxy), it would
   * make mass conservation more explicit to zero the properties of galaxy p after
   * the merger.
   * TODO Correct artificial diffusion of metals when BulgeFormationInMinorMergersOn=1. */
  int outputbin, j, ii;
  float tspin[3],tmass,pmass;
#ifdef TRACK_NMERGERS
  double mass_ratio, tmass_total, pmass_total;
#endif
  /* t central, p satellite */

  mass_checks(p,"model_mergers.c",__LINE__);
  mass_checks(t,"model_mergers.c",__LINE__);

  /* angular momentum transfer between gas*/
  tmass= Gal[t].ColdGas;
  pmass= Gal[p].ColdGas;

  //Gal[t].MergeSat +=(Gal[p].DiskMass+Gal[p].BulgeMass);
  //Gal[p].MergeSat=0.;

#ifdef TRACK_NMERGERS
  tmass_total = Gal[t].DiskMass+Gal[t].BulgeMass+Gal[t].ColdGas;
  pmass_total = Gal[p].DiskMass+Gal[p].BulgeMass+Gal[p].ColdGas;

  if(max(tmass_total,pmass_total) > 0.)
     mass_ratio = min(tmass_total,pmass_total) / max(tmass_total,pmass_total);
  else
    mass_ratio = 1.0;
#endif

#ifdef H2_AND_RINGS
  double ringtot, fractionRings[RNUM], rd;
  int jj;

  //distribute the satellite gas evenly through the rings and then assume the same
  //profile as for infall when it goes to the disk of the central galaxy
  rd=get_initial_disk_radius(Gal[t].HaloNr, t)/3.;
  ringtot=1-(1+RingRadius[RNUM-1]/rd)/exp(RingRadius[RNUM-1]/rd);
  //fractionRings * RNUM gives the fraction in a ring as a function of the mass in a ring
  //without * RNUM gives the fraction in a ring as a function of the total mass

  //normalized fraction of mass in each ring
  fractionRings[0]=(1-(1+RingRadius[0]/rd)/exp(RingRadius[0]/rd))/ringtot;
  for(j=1; j<RNUM; j++)
    fractionRings[j]=((1+RingRadius[j-1]/rd)/exp(RingRadius[j-1]/rd)-(1+RingRadius[j]/rd)/exp(RingRadius[j]/rd))/ringtot;

  for(j=0; j<RNUM; j++)
    if(fractionRings[j]<0.)
      fractionRings[j]=0.;


  for(j=0; j<RNUM; j++)
    {
      Gal[p].ColdGasRings[j]=Gal[p].ColdGas*fractionRings[j];
      for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
     	Gal[p].MetalsColdGasRings[j][ii] = (Gal[p].MetalsColdGas[ii] * fractionRings[j]);

#ifdef DETAILED_METALS_AND_MASS_RETURN
#ifdef INDIVIDUAL_ELEMENTS //All: [H][He][Cb][N][O][Ne][Mg][Si][S][Ca][Fe] or //Only [H][He][O][Mg][Fe]
      int kk;
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[p].ColdGasRings_elements[j][kk]=Gal[p].ColdGas_elements[kk]*fractionRings[j];
#endif
#endif
      fractionRings[j]=1.;
    }
  transfer_material_with_rings(t,"ColdGas",p,"ColdGas",fractionRings,"model_mergers.c", __LINE__);
#else
  transfer_material(t,"ColdGas",p,"ColdGas",1.,"model_mergers.c", __LINE__);
  //transfer_material(t,"EjectedMass",p,"ColdGas",1.,"model_mergers.c", __LINE__);
#endif
  transfer_material(t,"HotGas",p,"HotGas",1.,"model_mergers.c", __LINE__);
  //transfer_material(t,"ReheatedGas",p,"ReheatedGas",1.,"model_mergers.c", __LINE__);
  transfer_material(t,"EjectedMass",p,"EjectedMass",1.,"model_mergers.c", __LINE__); //TODO chose move to ejected or hot

  mass_checks(p,"model_mergers.c",__LINE__);
  mass_checks(t,"model_mergers.c",__LINE__);


#ifdef TRACK_MASSGROWTH_CHANNELS
  Gal[t].MassFromMergers += Gal[p].DiskMass + Gal[p].BulgeMass;

  Gal[p].MassFromInSitu = 0.;
  Gal[p].MassFromMergers = 0.;
  Gal[p].MassFromBursts = 0.;

#ifdef STAR_FORMATION_HISTORY
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
  for (ii=0; ii<=Gal[t].sfh_ibin; ii++) Gal[t].sfh_MassFromMergers[ii] += Gal[t].sfh_DiskMass[ii]+Gal[t].sfh_BulgeMass[ii];

  for (ii=0; ii<=Gal[p].sfh_ibin; ii++) Gal[p].sfh_MassFromInSitu[ii] = 0.;
  for (ii=0; ii<=Gal[p].sfh_ibin; ii++) Gal[p].sfh_MassFromMergers[ii] = 0.;
  for (ii=0; ii<=Gal[p].sfh_ibin; ii++) Gal[p].sfh_MassFromBursts[ii] = 0.;
#endif
#endif
#endif

#ifdef TRACK_BURST
    /* The whole burst component gets transferred */
  transfer_material(t,"BurstMass",p,"BurstMass",1.,"model_mergers.c", __LINE__);
#endif
#ifdef H2_AND_RINGS
  for(j=0; j<RNUM; j++)
      fractionRings[j]=1.;
#endif

  //Bulge occupies the same place as the disk it forms from, after the merger is finished the material will be
  // distrubted in a Jaffe profile
#ifdef H2_AND_RINGS
  if(BulgeFormationInMinorMergersOn)
    transfer_material_with_rings(t,"BulgeMass",p,"DiskMass",fractionRings,"model_mergers.c", __LINE__);
  else
    transfer_material_with_rings(t,"DiskMass",p,"DiskMass",fractionRings,"model_mergers.c", __LINE__);

  if(Gal[p].BulgeMass>0.)
    transfer_material_with_rings(t,"BulgeMass",p,"BulgeMass",fractionRings,"model_mergers.c", __LINE__);

#else
  if(BulgeFormationInMinorMergersOn)
    transfer_material(t,"BulgeMass",p,"DiskMass",1.,"model_mergers.c", __LINE__);
  else
    transfer_material(t,"DiskMass",p,"DiskMass",1.,"model_mergers.c", __LINE__);
  if(Gal[p].BulgeMass>0.)
    transfer_material(t,"BulgeMass",p,"BulgeMass",1.,"model_mergers.c", __LINE__);
#endif


  transfer_material(t,"ICM",p,"ICM",1.,"model_mergers.c", __LINE__);

  Gal[t].BlackHoleMass += Gal[p].BlackHoleMass;
  Gal[p].BlackHoleMass=0.;
  Gal[t].BlackHoleGas += Gal[p].BlackHoleGas;
  Gal[p].BlackHoleGas=0.;
  //Gal[t].StarMerge += Gal[p].StarMerge;
  //Gal[p].StarMerge=0.;

  mass_checks(p,"model_mergers.c",__LINE__);
  mass_checks(t,"model_mergers.c",__LINE__);

  /*update the gas spin - gasdiskradius updated in the end of deal_with_galaxy_mergers*/
  for(ii=0;ii<3;ii++)
    tspin[ii]=Gal[t].ColdGasSpin[ii]*tmass+Gal[t].HaloSpin[ii]*pmass;

  /*double halospinpar=sqrt(Halo[GaRStellarDiskSatl[t].HaloNr].Spin[0] * Halo[Gal[t].HaloNr].Spin[0] +
			   Halo[Gal[t].HaloNr].Spin[1] * Halo[Gal[t].HaloNr].Spin[1] +
			   Halo[Gal[t].HaloNr].Spin[2] * Halo[Gal[t].HaloNr].Spin[2] );
  for(ii=0;ii<3;ii++)
      tspin[ii]=Gal[t].ColdGasSpin[ii]*tmass+1./sqrt(3)*halospinpar*pmass;*/
  if (Gal[t].ColdGas != 0)
    for (ii=0;ii<3;ii++)
	  Gal[t].ColdGasSpin[ii]=tspin[ii]/(Gal[t].ColdGas);


  Gal[t].Sfr += Gal[p].Sfr;
  if(BulgeFormationInMinorMergersOn)
    Gal[t].SfrBulge += Gal[p].Sfr;

#ifdef TRACK_NMERGERS
  Gal[t].NMajorMergers += Gal[p].NMajorMergers;
  Gal[t].NMinorMergers += Gal[p].NMinorMergers;

  if(mass_ratio > ThreshMajorMerger)
    Gal[t].NMajorMergers += 1.;// / deltaT;
  else
    Gal[t].NMinorMergers += 1.;// / deltaT;
#endif

  /* Add the black hole accretion rates.  This makes little sense but is not
   * used if the superior BlackHoleGrowth==1 switch is on. */
  Gal[t].QuasarAccretionRate += Gal[p].QuasarAccretionRate;
  Gal[t].RadioAccretionRate += Gal[p].RadioAccretionRate;

  for(outputbin = 0; outputbin < NOUT; outputbin++)
  	 Gal[t].MassWeightAge[outputbin] += Gal[p].MassWeightAge[outputbin];

}


/** @brief In a major merger, both disks are destroyed and all the mass transferred
 *         to the bulge. The galaxies have already been merged, so all we need to do here
 *         is transfer stars from disk to bulge of the central galaxy. */
void make_bulge_from_burst(int p)
{
  /* generate bulge */
#ifdef H2_AND_RINGS
  int jj;
  double fractionRings[RNUM];
  for (jj=0;jj<RNUM;jj++)
    fractionRings[jj]=1.;
  transfer_material_with_rings(p,"BulgeMass",p,"DiskMass",fractionRings,"model_mergers.c", __LINE__);
#else
  transfer_material(p,"BulgeMass",p,"DiskMass",1.,"model_mergers.c", __LINE__);
#endif
  mass_checks(p,"model_mergers.c",__LINE__);

  /*  update the star formation rate */
  Gal[p].SfrBulge  = Gal[p].Sfr;
#ifdef H2_AND_RINGS
  for(jj=0;jj<RNUM;jj++) Gal[p].SfrRings[jj]=0;
#endif

}

/** @brief Merger burst recipe from Somerville 2001 (used after Croton2006) */

double collisional_starburst_recipe(double mass_ratio, int merger_centralgal, int centralgal,
				  double time, double deltaT)
{
  /** @brief If StarBurstModel = 0 (since Croton2006), the Somerville 2001
   *         model of bursts is used. The burst can happen for both major
   *         and minor mergers, with a fraction of the added cold gas from
   *         the satellite and central being consumed. SN Feedback from
   *         starformation is computed and the sizes of bulge and disk
   *         followed (not done for the other burst mode).*/

  double mstars, eburst, Ggas;
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  double metallicitySF=0.;
#endif
#endif
#ifdef H2_AND_RINGS
  double mstarsRings[RNUM];
  int j;
#endif
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
  /*double halospinpar=sqrt(Halo[GaRStellarDiskSatl[t].HaloNr].Spin[0] * Halo[Gal[t].HaloNr].Spin[0] +
			   Halo[Gal[t].HaloNr].Spin[1] * Halo[Gal[t].HaloNr].Spin[1] +
			   Halo[Gal[t].HaloNr].Spin[2] * Halo[Gal[t].HaloNr].Spin[2] );*/

#ifdef H2_AND_RINGS
  for(j=0;j<RNUM;j++)
    {
      mstarsRings[j] = eburst * Gal[merger_centralgal].ColdGasRings[j];
      if(mstarsRings[j] < 0.0) mstarsRings[j] = 0.0;
    }
#endif

#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN //otherwise there is another check inside SN_feedback
  if(mstars > Gal[merger_centralgal].ColdGas)
    mstars = Gal[merger_centralgal].ColdGas;
#ifdef H2_AND_RINGS
  for(j=0;j<RNUM;j++)
    if(mstarsRings[j] > Gal[merger_centralgal].ColdGasRings[j])
      mstarsRings[j] = Gal[merger_centralgal].ColdGasRings[j];
#endif
#endif

  /* Store the value of the metallicity of the cold phase when SF occurs.
   * Used to update luminosities below */
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  int ii;
  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
      metallicitySF += Gal[merger_centralgal].MetalsColdGas[ii];
  metallicitySF /= Gal[merger_centralgal].ColdGas;
#endif
#endif

//if FEEDBACK_COUPLED_WITH_MASS_RETURN feedback happens only when stars die,
//there is no need to balance it with SF
#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
  if (mstars > 0.)
#ifndef H2_AND_RINGS
    update_stars_due_to_reheat(merger_centralgal, centralgal, &mstars);
#else
    update_stars_due_to_reheat(merger_centralgal, centralgal, &mstars, mstarsRings);
#endif
#endif //FEEDBACK_COUPLED_WITH_MASS_RETURN
  /*double halospinpar=sqrt(Halo[GaRStellarDiskSatl[t].HaloNr].Spin[0] * Halo[Gal[t].HaloNr].Spin[0] +
			   Halo[Gal[t].HaloNr].Spin[1] * Halo[Gal[t].HaloNr].Spin[1] +
			   Halo[Gal[t].HaloNr].Spin[2] * Halo[Gal[t].HaloNr].Spin[2] );*/
 
  /*  update the star formation rate */
  Gal[merger_centralgal].Sfr += mstars / deltaT;
#ifdef H2_AND_RINGS
  for(j=0;j<RNUM;j++)
    Gal[merger_centralgal].SfrRings[j] += mstarsRings[j] / deltaT;
#endif
  mass_checks(merger_centralgal,"model_mergers.c",__LINE__);


  // update_from_star_formation can only be called
  // after SD_feeedback recipe since stars need to be re_set once the reheated mass is known
  // (star formation and feedback share the same fraction of cold gas)
  int nstep=-1;
  if (mstars > 0.)
#ifndef H2_AND_RINGS
    update_from_star_formation(merger_centralgal, mstars, "merger", nstep); // true indicates starburst
#else
    update_from_star_formation(merger_centralgal, mstars, mstarsRings, "merger", nstep); // true indicates starburst
#endif


  mass_checks(merger_centralgal,"model_mergers.c",__LINE__);

  update_massweightage(merger_centralgal, mstars, time);

#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
  if (mstars > 0.)
#ifndef H2_AND_RINGS
    SN_feedback(merger_centralgal, centralgal, mstars, "ColdGas");
#else
    SN_feedback(merger_centralgal, centralgal, mstars, mstarsRings, "ColdGas");
#endif
#endif




#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
  //update the luminosities due to the stars formed
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
void bulgesize_from_merger(double mass_ratio, int merger_centralgal, int p,
			   double MstarCentral,double MbulgeCentral,double MgasCentral,
			   double MstarSat, double MbulgeSat,double MgasSat, double frac,
			   double RgasCentral, double RStellarDiskCentral, double RgasSat, double RStellarDiskSat)
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

  double Mcen,Rcen;
  double Msat,Rsat;
  double fint=0.5;
  double c=0.5;

  /* calculate radius for the object that will form the new bulge - Rcen and Rsat */
  /* Minor Merger */
  if(mass_ratio < ThreshMajorMerger) {
      /* In a minor merger only the stars of the satellite galaxy are moved
       * to the bulge of the central galaxy, therefore only stellar
       * components are used to compute radius and masses. */
      /* In a minor merger only consider the bulge mass of the central galaxy */
      Mcen=MbulgeCentral;
      Rcen=Gal[merger_centralgal].BulgeSize;
      /* and stellarmass of satellite*/
      if(BulgeFormationInMinorMergersOn)
	  Msat=MstarSat;
      else
	  Msat=MbulgeSat;

      if (Msat >0.0)
	  Rsat=(RStellarDiskSat*1.68*(MstarSat-MbulgeSat)+Gal[p].BulgeSize*MbulgeSat)/MstarSat;
      else
	  Rsat=0.0;
  }
  /* Major Merger */
  else {
      /* on a major merger both stellar and gas (after a burst) components
       * form the final bulge and need to be considered */

      // Mcen = bulge mass + burst of central
      Mcen=MstarCentral+frac*MgasCentral;
      if (Mcen > 0.0)
	  Rcen=(RStellarDiskCentral*1.68*(MstarCentral-MbulgeCentral)+Gal[merger_centralgal].BulgeSize*MbulgeCentral+RgasCentral*frac*MgasCentral*1.68)/(MgasCentral*frac+MstarCentral);
      else
	  Rcen=0.0;
      // and satellite Msat
      Msat=MstarSat+frac*MgasSat;
      if (Msat > 0.0)
	  Rsat=(RStellarDiskSat*1.68*(MstarSat-MbulgeSat)+Gal[p].BulgeSize*MbulgeSat+RgasSat*frac*MgasSat*1.68)/(MgasSat*frac+MstarSat);
      else
	  Rsat=0.0;
  }

  float tiny;
  tiny=1.e-8;
  if(Rsat>0. && Rsat<tiny)
  	Rsat=tiny;
  if(Rcen>0. && Rcen<tiny)
    	Rcen=tiny;
  /* If both original radius are bigger then 0 then this is Eq. 33 in Guo 2010
   * solved for R_new,bulge with all terms divided by G and C. */

  // The problem with these conditions is that, once we have an erroneously small value of BulgeSize, we will keep it.
  if(Rcen >= tiny && Rsat >= tiny)
      Gal[merger_centralgal].BulgeSize=(Msat+Mcen)*(Msat+Mcen)/(Msat*Msat/Rsat+Mcen*Mcen/Rcen+fint/c*Msat*Mcen/(Rsat+Rcen));
  else if(Rcen >= tiny)
      Gal[merger_centralgal].BulgeSize=(Msat+Mcen)*(Msat+Mcen)/(Mcen*Mcen/Rcen+fint/c*Msat*Mcen/(Rsat+Rcen));
  else if(Rcen <= tiny && Rsat >= tiny)
      Gal[merger_centralgal].BulgeSize=(Msat+Mcen)*(Msat+Mcen)/(Msat*Msat/Rsat+fint/c*Msat*Mcen/(Rsat+Rcen));
  else
      Gal[merger_centralgal].BulgeSize=0.0;

  if ((Msat+Mcen > 0.0 && Gal[merger_centralgal].BulgeSize == 0.0 )||
      (Msat+Mcen == 0.0 && Gal[merger_centralgal].BulgeSize> 0.0)) {
      char sbuf[1000];
      printf("halonr=%d, merger_centralgal %d\n\n", Gal[merger_centralgal].HaloNr, merger_centralgal);
      printf("New Bulge Mass from Central (Mcen)=%e\n New Bulge Mass from Satellite (Msat)=%e\n NewBulge size=%e\n\n",
	     Msat,Mcen,Gal[merger_centralgal].BulgeSize);
      if(mass_ratio < ThreshMajorMerger) {
	  printf("minor merger, new mass=original mass\n");
	  printf("New Bulge Mass From Central (Mcen)   = MbulgeCentral = %f\n",MbulgeCentral);
	  printf("New Bulge Mass From Satellite (Msat) = MstarSat  = %f\n",MstarSat);
      }
      else {
	  printf("New Bulge From Central (Mcen)   = MstarCentral+frac*MgasCentral = %f+%f*%f\n", MstarCentral, frac, MgasCentral);
	  printf("New Bulge From Satellite (Msat) = MstarSat+frac*MgasSat = %f+%f*%f\n", MstarSat, frac, MgasSat);
      }
      printf("BulgeSize from Central (Rcen)=%e\nBulgeSize from Satellite (Rsat)=%e\nmass ratio=%f\n\n",Rcen,Rsat, mass_ratio);
      printf("the following masses don't tell a lot because things have been merged already!!!\n");
      printf("    sat: BulgeMass=%0.7f, BulgeSize=%0.7f, GasMass=%0.7f, GasSize=%0.7f, DiskMass=%0.7f StellarSize=%0.7f \n",
	     Gal[p].BulgeMass, Gal[p].BulgeSize, Gal[p].ColdGas, Gal[p].ColdGasRadius,
	     Gal[p].DiskMass, Gal[p].DiskRadius);
      printf(	"central: BulgeMass=%0.7f, BulgeSize=%0.7f, GasMass=%0.7f, GasSize=%0.7f, DiskMass=%0.7f StellarSize=%0.7f \n",
		Gal[merger_centralgal].BulgeMass,  Gal[merger_centralgal].BulgeSize, Gal[merger_centralgal].ColdGas,
		Gal[merger_centralgal].ColdGasRadius, Gal[merger_centralgal].DiskMass, Gal[merger_centralgal].DiskRadius);
      sprintf(sbuf,"\n bulgesize wrong in merger");
      terminate(sbuf);
      exit(0);
  }

}

#ifdef H2_AND_RINGS
  //Bulge Mass was added into the same place as the disk, it will now be redistributed
  //according to a Jaffe profile and after the new bulge size has been calculated
void distribute_bulge_material(int merger_centralgal)
{
  double rb=Gal[merger_centralgal].BulgeSize, TotMassInsideRings=0.,fractionRings[RNUM];
  int jj, ii;

  if(rb>0.)
    TotMassInsideRings=(RingRadius[RNUM-1]/rb)/(1+RingRadius[RNUM-1]/rb);

  if(TotMassInsideRings>0.)
    {
      fractionRings[0]=(RingRadius[0]/rb)/(1+RingRadius[0]/rb)/TotMassInsideRings;
      for(jj=1; jj<RNUM; jj++)
	fractionRings[jj]= ((RingRadius[jj]/rb)/(1+RingRadius[jj]/rb)-(RingRadius[jj-1]/rb)/(1+RingRadius[jj-1]/rb))/TotMassInsideRings;
    }
  else
    for(jj=0; jj<RNUM; jj++)
      fractionRings[jj]=0.;

  //RINGS
  for(jj=0;jj<RNUM;jj++)
    {
      Gal[merger_centralgal].BulgeMassRings[jj] = fractionRings[jj]*Gal[merger_centralgal].BulgeMass;
      for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	Gal[merger_centralgal].MetalsBulgeMassRings[jj][ii] = (Gal[merger_centralgal].MetalsBulgeMass[ii] * fractionRings[jj]);

#ifdef INDIVIDUAL_ELEMENTS
      int kk;
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[merger_centralgal].BulgeMassRings_elements[jj][kk] = fractionRings[jj]*Gal[merger_centralgal].BulgeMass_elements[kk];
#endif
#ifdef STAR_FORMATION_HISTORY
      for (ii=0; ii<=Gal[merger_centralgal].sfh_ibin; ii++)
	Gal[merger_centralgal].sfh_BulgeMassRings[jj][ii] = fractionRings[jj]*Gal[merger_centralgal].sfh_BulgeMass[ii];
#endif
    }
}
#endif //H2_AND_RINGS
