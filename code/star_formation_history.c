/*
 * star_formation_history.c
 *
 * Routine to track the star-formation history of galaxies.
 * Keeps bins of the size of the current step of the semi-analytic.
 * As a galaxy ages, bins are merged in factors of two.
 *
 *  Created on: Oct 26, 2010
 *      Author: Peter Thomas (petert)
 *
 *  Changed on: April 09, 2012
 *          by: Bruno Henriques (bmh20)
 */

/* This version assumes that all stars formed on a particular timestep
 * go into a single time history bin with star formation time equal to
 * the mid point of the bin. All times are code units (Mpc/Km/s/h) and
 * are the time till present.
 *
 * Number of time bins needed:
 * #define SFH_NMERGE 3
 * #define SFH_NBIN 19
 * At worst there are SFH_NMERGE-1 bins of each size.
 * So set SFH_NBIN to the smallest integer greater than
 * (SFH_NMERGE-1)log_2(MAXSNAPS*STEPS/(SFH_NMERGE-1)+1)
 * MAXSNAPS=61, STEPS=20, NSH_NMERGE=3 --> SFH_NBIN=19 (actually uses 19)
 * MAXSNAPS=61, STEPS=20, NSH_NMERGE=5 --> SFH_NBIN=34 (actually uses 33)
 *
 * Usage:
 *
 * On init.c:
 *   create_sfh_bins();
 * Generates the reference structure for storing the star formation histories in
 * logarithmic bins (for each snapshot/time step combination). In the code galaxy
 * structures are adjusted with respect to this structure at each step.
 * double SFH_t[MAXSNAPS][STEPS][SFH_NBIN]; //Time to present (i.e. z=0 ?) at the low-z edge of the bin (code units)
 * double SFH_dt[MAXSNAPS][STEPS][SFH_NBIN]; //Time width of the bin (code units)
 * int SFH_Nbins[MAXSNAPS][STEPS][SFH_NBIN]; //Number of bins merged in each bin (only useful for the merging algorithm)
 * int SFH_ibin[MAXSNAPS][STEPS]; //Last active bin
 *
 * On initialising galaxy p:
 * 	  sfh_initialise(p);
 * Whenever new stars are created:
 *    sfh_update_bins(p, time);
 *    Gal[p].sfh_DiskMass[Gal[p].sfh_ibin]+=added_mass;
 * (could make the above a function call but it hardly seems worth
 * it as would need a separate call for each property to be updated).
 * To merge galaxy p1 into galaxy p:
 *    sfh_merge(p,p1);
 * For debugging purposes:
 *    sfh_print(p);
 * will print the time bin structure for galaxy p.
 *
 * The following variables are all defined as part of the galaxy structures:
 * int sfh_ibin
 *    Index of highest bin are currently in use
 * double sfh_age
 *    Time in code units of last call to sph_update_bins.
 *    (Not strictly required, but useful in L-Galaxies to prevent having to
 *    pass time explicitly to the save subroutine.)
 * int sfh_dt[SFH_NBIN]
 *    Width of each time bin in code units
 * int sfh_t[SFH_NBIN]
 *    Time at low-z edge of bin in code units
 * float sfh_time[SFH_NBIN]
 *    time till output from the middle of the bin in years.
 *    Used only in save.c
 * float sfh_DiskMass[SFH_NBIN]
 *    Disk mass in bin in standard mass units.
 * float sfh_MetalsDiskMass[SFH_NBIN];
 *    Metals locked up in stars in the disk ditto.
 */

#include "allvars.h"
#include "proto.h"

void sfh_initialise(int p)
{
  /* Initialises the sfh-variables for a galaxy */
  int i, ii;
#ifdef H2_AND_RINGS
  int jj;
#endif

  for (i=0;i<SFH_NBIN;i++){
    Gal[p].sfh_dt[i]=0.;
    Gal[p].sfh_t[i]=0.;
    Gal[p].sfh_Nbins[i]=0;
    Gal[p].sfh_DiskMass[i]=0.;
#ifdef H2_AND_RINGS
    for(jj=0; jj<RNUM; jj++)
      {
	Gal[p].sfh_DiskMassRings[jj][i]=0.;
#ifdef RINGS_IN_BULGES
	Gal[p].sfh_BulgeMassRings[jj][i]=0.;
#endif
      }
#endif
    Gal[p].sfh_BulgeMass[i]=0.;
    Gal[p].sfh_ICM[i]=0.;
    for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
      {
	Gal[p].sfh_MetalsDiskMass[i][ii] = 0.;
	Gal[p].sfh_MetalsBulgeMass[i][ii] = 0.;
	Gal[p].sfh_MetalsICM[i][ii] = 0.;
      }
#ifdef INDIVIDUAL_ELEMENTS
    int kk;
    for(kk=0;kk<NUM_ELEMENTS;kk++)
      {
	Gal[p].sfh_DiskMass_elements[i][kk]=0.;
	Gal[p].sfh_BulgeMass_elements[i][kk]=0.;
	Gal[p].sfh_ICM_elements[i][kk]=0.;
      }
#endif
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
    Gal[p].sfh_MassFromInSitu[i]=0.;
    Gal[p].sfh_MassFromMergers[i]=0.;
    Gal[p].sfh_MassFromBursts[i]=0.;
#endif
#ifdef TRACK_BURST
    Gal[p].sfh_BurstMass[i]=0.;
#endif
  }

  /* Create first bin */
  Gal[p].sfh_ibin=0;

  /* Age is used for comparing galaxies during mergers, 
   * so needs to have a value set in case a merger happens before stars 
   * form (which can happen). */
  Gal[p].sfh_age=0.;
}

void sfh_merge(int p, int p1)
{
  /* Merge galaxy p1 into galaxy p */
  int i, ii;
#ifdef H2_AND_RINGS
  int jj;
#endif
  /* Perform minimal test that the two galaxies have the same time structure */
  if (Gal[p1].sfh_ibin != Gal[p].sfh_ibin) {
    printf("sfh_merge: trying to merge galaxies with different sfh bins\n");
    sfh_print(p);
    sfh_print(p1);
    exit(1);
  }

  /* The zero-ing of galaxy p1 here is not strictly necessary as galaxy p1 should
   * cease to exist after merging, but helps to make mass conservation explicit. */
  for(i=0;i<=Gal[p].sfh_ibin;i++) {
    Gal[p].sfh_DiskMass[i]+=Gal[p1].sfh_DiskMass[i];
#ifdef H2_AND_RINGS
    for(jj=0; jj<RNUM; jj++)
      {
	Gal[p].sfh_DiskMassRings[jj][i]+=Gal[p1].sfh_DiskMassRings[jj][i];
#ifdef RINGS_IN_BULGES
	Gal[p].sfh_BulgeMassRings[jj][i]+=Gal[p1].sfh_BulgeMassRings[jj][i];
#endif
      }
#endif
    Gal[p].sfh_BulgeMass[i]+=Gal[p1].sfh_BulgeMass[i];
    Gal[p].sfh_ICM[i]+=Gal[p1].sfh_ICM[i];

    Gal[p1].sfh_DiskMass[i]=0.;
#ifdef H2_AND_RINGS
    for(jj=0; jj<RNUM; jj++)
      {
	Gal[p1].sfh_DiskMassRings[jj][i]=0.;
#ifdef RINGS_IN_BULGES
	Gal[p1].sfh_BulgeMassRings[jj][i]=0.;
#endif
      }
#endif
    Gal[p1].sfh_BulgeMass[i]=0.;
    Gal[p1].sfh_ICM[i]=0.;

    for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
      {
	 Gal[p].sfh_MetalsDiskMass[i][ii] += Gal[p1].sfh_MetalsDiskMass[i][ii];
	 Gal[p].sfh_MetalsBulgeMass[i][ii] += Gal[p1].sfh_MetalsBulgeMass[i][ii];
	 Gal[p].sfh_MetalsICM[i][ii] += Gal[p1].sfh_MetalsICM[i][ii];
	 Gal[p1].sfh_MetalsDiskMass[i][ii] = 0.;
	 Gal[p1].sfh_MetalsBulgeMass[i][ii] = 0.;
	 Gal[p1].sfh_MetalsICM[i][ii] = 0.;
      }

#ifdef INDIVIDUAL_ELEMENTS
    int kk;
    for(kk=0;kk<NUM_ELEMENTS;kk++)
      {
	Gal[p].sfh_DiskMass_elements[i][kk] += Gal[p1].sfh_DiskMass_elements[i][kk];
	Gal[p].sfh_BulgeMass_elements[i][kk] += Gal[p1].sfh_BulgeMass_elements[i][kk];
	Gal[p].sfh_ICM_elements[i][kk] += Gal[p1].sfh_ICM_elements[i][kk];
	Gal[p1].sfh_DiskMass_elements[i][kk]=0.;
	Gal[p1].sfh_BulgeMass_elements[i][kk]=0.;
	Gal[p1].sfh_ICM_elements[i][kk]=0.;
      }
#endif
#ifdef TRACK_BURST
    Gal[p].sfh_BurstMass[i]+=Gal[p1].sfh_BurstMass[i];
    Gal[p1].sfh_BurstMass[i]=0.;
#endif
  }
  /* Again, not strictly necessary, but safe. */
  Gal[p1].sfh_ibin=0;
  Gal[p1].sfh_age=0.;

}

void sfh_print(int p) {
  /* Prints out populated sfh_structure.
   * Does sum of Disk + Bulge only. */
  int i;

  printf("For galaxy %d:\n",p);
  printf("sfh_ibin=%d\n",Gal[p].sfh_ibin);
  printf("sfh_age=%f\n",Gal[p].sfh_age);
  printf("  i    dt   t      Stars      Metals\n");
  for(i=0;i<SFH_NBIN;i++)
    if (Gal[p].sfh_dt[i]!=0) {
      printf("%5d %5e %5e %12f\n",i,Gal[p].sfh_dt[i],Gal[p].sfh_t[i],(Gal[p].sfh_DiskMass[i]+Gal[p].sfh_BulgeMass[i]));
      printf(".type2 [Msun] = %.2f\n",(Gal[p].sfh_MetalsDiskMass[i][0]+Gal[p].sfh_MetalsBulgeMass[i][0])*1.0e10/Hubble_h);
      printf(".type1a [Msun]  = %.2f\n",(Gal[p].sfh_MetalsDiskMass[i][1]+Gal[p].sfh_MetalsBulgeMass[i][1])*1.0e10/Hubble_h);
      printf(".agb  [Msun]   = %.2f\n",(Gal[p].sfh_MetalsDiskMass[i][2]+Gal[p].sfh_MetalsBulgeMass[i][2])*1.0e10/Hubble_h);
#ifdef INDIVIDUAL_ELEMENTS
      int kk;
      for(kk=0;kk<NUM_ELEMENTS;kk++)
	printf("Element %d [Msun]  = %.2f\n",kk,Gal[p].sfh_DiskMass_elements[i][kk]+Gal[p].sfh_BulgeMass_elements[i][kk]);
#endif
      printf(".......................\n");
    }
}

void create_sfh_bins()
{
  double previoustime, newtime, deltaT, time;
  int snap, step, sfh_ibin, i, j, sfh_Nbins[SFH_NBIN];
  int ibin_max=0;
  double sfh_t[SFH_NBIN];

  for(snap = 0; snap < MAXSNAPS; snap++) {
 	  for(step=0;step < STEPS;step++) {
 		  for(j=0;j < SFH_NBIN;j++) {
 			  SFH_t[snap][step][j]=0;
 			  SFH_dt[snap][step][j]=0;
 			  SFH_ibin[snap][step]=0;
 			  SFH_Nbins[snap][step][j] = 0;
 		  }
 	  }
   }

  for(i=0;i<SFH_NBIN;i++) {
    sfh_Nbins[i]=0;
    sfh_t[i]=0.;
  }
  sfh_ibin=0;

	//for(snap=0;snap<(LastDarkMatterSnapShot+1)-1;snap++) {
  for(snap=0;snap<(LastDarkMatterSnapShot+1);snap++) {
    previoustime = NumToTime(snap);
    newtime = NumToTime(snap+1);
    deltaT = previoustime - newtime;

    for(step=0;step<STEPS;step++) {
      int ibin;
      int flag_merged_bins; // Boolean used to check whether have merged bins
      int dt_merge; // Size of bins that we are checking for merging
      int n_merge; // Number of bins of this size

      time = previoustime - (step + 1.0) * (deltaT / STEPS);
      ibin=sfh_ibin;

      //printf("sna=%d step=%d step time=%f time low=%f\n",
      //		snap,step,(previoustime - (step + 0.5) * (deltaT / STEPS))*UnitTime_in_years/Hubble_h/1.e9,
      //		(time)*UnitTime_in_years/Hubble_h/1.e9);
      //Add one extra bin
      if(snap==0 && step==0) {
      	sfh_t[0]=time;
      	sfh_Nbins[0]=1;
      }
      else {
      	ibin+=1;
      	if(ibin==SFH_NBIN)
      		terminate("sfh_update_bins: too many bins required\n");
      	ibin_max=max(ibin_max,ibin);
      	sfh_Nbins[ibin]=1;
      	sfh_t[ibin]=time;
      }

      /* Now merge bins where we have SFH_NMERGE bins of the same size.
       * Need to do this iteratively. */
      flag_merged_bins=1;
      while(flag_merged_bins) {
      	flag_merged_bins=0;
      	dt_merge=sfh_Nbins[0];
      	i=0;
      	// Will have checked all bins once dt_merge drops to zero
      	while(!flag_merged_bins && dt_merge>0) {
      		// Count number of bins of this size
      		n_merge=0;
      		// The follwoing 2 lines are to suppress a warning message
		int itemp;
		itemp=i;
      		for(i=itemp;sfh_Nbins[i]==dt_merge;i++) n_merge+=1;
      		/* If fewer than SFH_NMERGE bins then do nothing
      		 * (SFH_NMERGE+1 bins if dt_merge=1)
      		 * else exit loop and flag for merging */
      		if (n_merge<SFH_NMERGE || (n_merge==SFH_NMERGE && dt_merge==1)) {
      			/* In new version of the code, treat smallest bins just like any others */
      			//if (n_merge<SFH_NMERGE) {
      			dt_merge/=2;
      			n_merge=0;
      		}
      		else {
      			flag_merged_bins=1;
      			i=i-n_merge;
      		}
      	}
	
      	/* At this point, if flag_merged_bins is set then
      	 * we have to merge SFH_NMERGE bins into SFH_NMERGE-1. */
      	if(flag_merged_bins) {
      		/* Merge bins i and i+1 */
      		sfh_Nbins[i]*=2;
      		sfh_t[i]=sfh_t[i+1];
      		/* Relabel all the other bins */
      		for(i=i+1;i<ibin;i++) {
      			sfh_Nbins[i]=sfh_Nbins[i+1];
      			sfh_t[i]=sfh_t[i+1];
      		}
      		sfh_Nbins[i]=0;
      		sfh_t[i]=0.;
      		ibin=i-1;
      	}
      } // End loop over bin merging

      sfh_ibin=ibin;

      for(j=0;j<=sfh_ibin;j++)
      {
      	SFH_t[snap][step][j]=sfh_t[j]; //Time to present at the low-z edge of the bin (code units)
      	SFH_Nbins[snap][step][j]=sfh_Nbins[j];//Number of bins merged in each bin (only useful for the merging algorithm)
      	if(j==0)
      		SFH_dt[snap][step][j]=NumToTime(0)-sfh_t[j];//Time width of the bin (code units)
      	else
      		SFH_dt[snap][step][j]=sfh_t[j-1]-sfh_t[j];//Time width of the bin (code units)
      	//printf("snap=%d step=%d bin=%d time=%f time_low=%f\n",
      	//		snap,step,j,(SFH_t[snap][step][j]+SFH_dt[snap][step][j]/2.)*UnitTime_in_years/Hubble_h/1.e9,
      	//		(SFH_t[snap][step][j])*UnitTime_in_years/Hubble_h/1.e9);
      }	
      SFH_ibin[snap][step]=sfh_ibin; //Last active bin

    }//end loop on steps

  }//end loop on snaps


#ifdef PARALLEL
  if(ThisTask==0)
#endif
  printf("Max number of SFH bins used = %d\n",ibin_max+1);

}

void sfh_update_bins(int p, int snap, int step, double time)
{
  /* Adds new bins as required.
   * Then merges bins whenever you have three or more of the same size.
   * Assumes that time counts from zero at the big bang. */
  int i, j, ii; // loop index
  int SFH_ibin; //desired ibin (i.e. bin in question in for loop below)
#ifdef H2_AND_RINGS
  int jj;
#endif
#ifdef INDIVIDUAL_ELEMENTS
  int kk;
#endif

  SFH_ibin=0;

  Gal[p].sfh_age=time;

  //t=time/SFH_TIME_INTERVAL;
  //ibin=Gal[p].sfh_ibin;
  for(i=0;i<SFH_NBIN;i++) //ROB: Could this not be: for(i=0;i<Gal[p].sfh_ibin;i++) ?
    if(SFH_Nbins[snap][step][i]>0) //i.e. If bin is active...
      SFH_ibin=i; //Assign with 'bin in question'

  if (Gal[p].sfh_ibin == 0) //i.e. If highest active bin is bin 0...
  {
    for(i=0;i<=SFH_ibin;i++) {
      Gal[p].sfh_t[i]=SFH_t[snap][step][i];
      Gal[p].sfh_Nbins[i]=SFH_Nbins[snap][step][i];
    }
    Gal[p].sfh_ibin=SFH_ibin;
  }
  else //i.e. If highest active bin is > bin 0...
  {
    i=0;
    
    while(i<=SFH_ibin) //Up to 'bin in question'...
    {
      if(i<=Gal[p].sfh_ibin) //...until highest active bin is reached...
      {
	if(Gal[p].sfh_Nbins[i]!= SFH_Nbins[snap][step][i]) //...and until bin has grown to required size.
	{
	  // Merge bins i and i+1
	  Gal[p].sfh_Nbins[i]+=Gal[p].sfh_Nbins[i+1];
	  Gal[p].sfh_t[i]=Gal[p].sfh_t[i+1];
	  Gal[p].sfh_DiskMass[i]+=Gal[p].sfh_DiskMass[i+1];
#ifdef H2_AND_RINGS
	  for(jj=0; jj<RNUM; jj++)
	    {
	      Gal[p].sfh_DiskMassRings[jj][i]+=Gal[p].sfh_DiskMassRings[jj][i+1];
#ifdef RINGS_IN_BULGES
	      Gal[p].sfh_BulgeMassRings[jj][i]+=Gal[p].sfh_BulgeMassRings[jj][i+1];
#endif
	    }
#endif
	  Gal[p].sfh_BulgeMass[i]+=Gal[p].sfh_BulgeMass[i+1];
	  Gal[p].sfh_ICM[i]+=Gal[p].sfh_ICM[i+1];
	  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	    {
	      Gal[p].sfh_MetalsDiskMass[i][ii] += Gal[p].sfh_MetalsDiskMass[i+1][ii];
	      Gal[p].sfh_MetalsBulgeMass[i][ii] += Gal[p].sfh_MetalsBulgeMass[i+1][ii];
	      Gal[p].sfh_MetalsICM[i][ii] += Gal[p].sfh_MetalsICM[i+1][ii];
	    }

#ifdef INDIVIDUAL_ELEMENTS
	  for(kk=0;kk<NUM_ELEMENTS;kk++)
	    {
	      Gal[p].sfh_DiskMass_elements[i][kk] += Gal[p].sfh_DiskMass_elements[i+1][kk];
	      Gal[p].sfh_BulgeMass_elements[i][kk] += Gal[p].sfh_BulgeMass_elements[i+1][kk];
	      Gal[p].sfh_ICM_elements[i][kk] += Gal[p].sfh_ICM_elements[i+1][kk];
	    }
#endif
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
	  Gal[p].sfh_MassFromInSitu[i]+=Gal[p].sfh_MassFromInSitu[i+1];
	  Gal[p].sfh_MassFromMergers[i]+=Gal[p].sfh_MassFromMergers[i+1];
	  Gal[p].sfh_MassFromBursts[i]+=Gal[p].sfh_MassFromBursts[i+1];
#endif
#ifdef TRACK_BURST
	  Gal[p].sfh_BurstMass[i]+=Gal[p].sfh_BurstMass[i+1];
#endif
	  // Relabel all the other bins
	  for(j=i+1;j<Gal[p].sfh_ibin;j++) {
	    Gal[p].sfh_Nbins[j]=Gal[p].sfh_Nbins[j+1];
	    Gal[p].sfh_t[j]=Gal[p].sfh_t[j+1];
	    Gal[p].sfh_DiskMass[j]=Gal[p].sfh_DiskMass[j+1];
#ifdef H2_AND_RINGS
	    for(jj=0; jj<RNUM; jj++)
	      {
		Gal[p].sfh_DiskMassRings[jj][j]=Gal[p].sfh_DiskMassRings[jj][j+1];
#ifdef RINGS_IN_BULGES
		Gal[p].sfh_BulgeMassRings[jj][j]=Gal[p].sfh_BulgeMassRings[jj][j+1];
#endif
	      }
#endif
	    Gal[p].sfh_BulgeMass[j] = Gal[p].sfh_BulgeMass[j+1];
	    Gal[p].sfh_ICM[j] = Gal[p].sfh_ICM[j+1];
	    for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	      {
		 Gal[p].sfh_MetalsDiskMass[j][ii] = Gal[p].sfh_MetalsDiskMass[j+1][ii];
		 Gal[p].sfh_MetalsBulgeMass[j][ii] = Gal[p].sfh_MetalsBulgeMass[j+1][ii];
		 Gal[p].sfh_MetalsICM[j][ii] = Gal[p].sfh_MetalsICM[j+1][ii];
	      }
#ifdef INDIVIDUAL_ELEMENTS
	    for(kk=0;kk<NUM_ELEMENTS;kk++)
	      {
		Gal[p].sfh_DiskMass_elements[j][kk]=Gal[p].sfh_DiskMass_elements[j+1][kk];
		Gal[p].sfh_BulgeMass_elements[j][kk]=Gal[p].sfh_BulgeMass_elements[j+1][kk];
		Gal[p].sfh_ICM_elements[j][kk]=Gal[p].sfh_ICM_elements[j+1][kk];
	      }
#endif
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
	  Gal[p].sfh_MassFromInSitu[j]=Gal[p].sfh_MassFromInSitu[j+1];
	  Gal[p].sfh_MassFromMergers[j]=Gal[p].sfh_MassFromMergers[j+1];
	  Gal[p].sfh_MassFromBursts[j]=Gal[p].sfh_MassFromBursts[j+1];
#endif
#ifdef TRACK_BURST
	    Gal[p].sfh_BurstMass[j]=Gal[p].sfh_BurstMass[j+1];
#endif
	  }

	  //set last bin to zero
	  Gal[p].sfh_Nbins[j]=0;
	  Gal[p].sfh_t[j]=0.;
	  Gal[p].sfh_DiskMass[j]=0.;
#ifdef H2_AND_RINGS
	  for(jj=0; jj<RNUM; jj++)
	    {
	      Gal[p].sfh_DiskMassRings[jj][j]=0.;
#ifdef RINGS_IN_BULGES
	      Gal[p].sfh_BulgeMassRings[jj][j]=0.;
#endif
	    }
#endif
	  Gal[p].sfh_BulgeMass[j]=0.;
	  Gal[p].sfh_ICM[j]=0.;
	  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	    {
	      Gal[p].sfh_MetalsDiskMass[j][ii] = 0.;
	      Gal[p].sfh_MetalsBulgeMass[j][ii] = 0.;
	      Gal[p].sfh_MetalsICM[j][ii] = 0.;
	    }
#ifdef INDIVIDUAL_ELEMENTS
	  for(kk=0;kk<NUM_ELEMENTS;kk++)
	    {
	      Gal[p].sfh_DiskMass_elements[j][kk]=0.;
	      Gal[p].sfh_BulgeMass_elements[j][kk]=0.;
	      Gal[p].sfh_ICM_elements[j][kk]=0.;
	    }
#endif
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
	  Gal[p].sfh_MassFromInSitu[j]+=0.;
	  Gal[p].sfh_MassFromMergers[j]+=0.;
	  Gal[p].sfh_MassFromBursts[j]+=0.;
#endif
#ifdef TRACK_BURST
	  Gal[p].sfh_BurstMass[j]=0.;
#endif
	  Gal[p].sfh_ibin=j-1;
	  
	  /* If there are no more time bins in the galaxy to merge and
	   * the last bin still doesn't have the required size
	   * re-size it according to the reference structure SFH */
	  if(Gal[p].sfh_Nbins[i+1]==0) {
	    Gal[p].sfh_Nbins[i]=SFH_Nbins[snap][step][i];
	    Gal[p].sfh_t[i]=SFH_t[snap][step][i];
	    i+=1;
	  }
	}
	else
	  i+=1;
      }
      else {
	//no more bins available in the galaxy, fill the rest times from SFH array
	for(j=i;j<=SFH_ibin;j++) {
	  Gal[p].sfh_Nbins[j]=SFH_Nbins[snap][step][j];
	  Gal[p].sfh_t[j]=SFH_t[snap][step][j];
	  Gal[p].sfh_DiskMass[j]=0.;
#ifdef H2_AND_RINGS
	  for(jj=0; jj<RNUM; jj++)
	    {
	      Gal[p].sfh_DiskMassRings[jj][j]=0.;
#ifdef RINGS_IN_BULGES
	      Gal[p].sfh_BulgeMassRings[jj][j]=0.;
#endif
	    }
#endif
	  Gal[p].sfh_BulgeMass[j]=0.;
	  Gal[p].sfh_ICM[j]=0.;
	  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
	    {
	      Gal[p].sfh_MetalsDiskMass[j][ii] = 0.;
	      Gal[p].sfh_MetalsBulgeMass[j][ii] = 0.;
	      Gal[p].sfh_MetalsICM[j][ii] = 0.;
	    }
#ifdef INDIVIDUAL_ELEMENTS
	  for(kk=0;kk<NUM_ELEMENTS;kk++)
	    {
	      Gal[p].sfh_DiskMass_elements[j][kk]=0.;
	      Gal[p].sfh_BulgeMass_elements[j][kk]=0.;
	      Gal[p].sfh_ICM_elements[j][kk]=0.;
	    }
#endif
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
	  Gal[p].sfh_MassFromInSitu[j]=0.;
	  Gal[p].sfh_MassFromMergers[j]=0.;
	  Gal[p].sfh_MassFromBursts[j]=0.;
#endif
#ifdef TRACK_BURST
	  Gal[p].sfh_BurstMass[j]=0.;
#endif
	  Gal[p].sfh_ibin=j;
	}
	i=j;
      }
    }//end while
  }//end else

  for(i=0;i<SFH_NBIN;i++) {
    if(i==0)
      Gal[p].sfh_dt[i]=NumToTime(0)-Gal[p].sfh_t[i];
    else
      Gal[p].sfh_dt[i]=Gal[p].sfh_t[i-1]-Gal[p].sfh_t[i];
  }
}

