/* This version of save.c contains code to convert star-formation histories to
 * luminosities.
 * This should really be a stand-alone piece of code that acts on L-Galaxies
 * output files. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"

void save_galaxies(int filenr, int tree)
{
#ifndef USE_MEMORY_TO_MINIMIZE_IO
  char buf[1000];
#endif
  FILE *fd;
  int i, n;
  struct GALAXY_OUTPUT galaxy_output;
#ifdef OUTPUT_MOMAF_INPUTS
  FILE *fd_momaf;
  struct MOMAF_INPUTS momaf_inputs;
#endif

  for(n = 0; n < NOUT; n++)
    {
#ifdef UPDATETYPETWO
#ifdef USE_MEMORY_TO_MINIMIZE_IO
      fd = (FILE *) 2;
      offset_auxdata = 0;
#else
#ifdef MRII
      sprintf(buf, "%s/treedata/treeaux_sf1_%03d.%d", SimulationDir, LastSnapShotNr, filenr);
#else
      sprintf(buf, "%s/treedata/treeaux_%03d.%d", SimulationDir, LastSnapShotNr, filenr);
#endif
      if(!(fd = fopen(buf, "r")))
	{
	  printf("Can't open file `%s'\n", buf);
	  fflush(stdout);
	  exit(1);
	}
#endif
      myfseek(fd, 4 * sizeof(int) + 2 * TotSnaps * sizeof(int)
	      + 2 * TotSnaps * Ntrees * sizeof(int) + 2 * sizeof(int) * NtotHalos, SEEK_CUR);
      Nids = CountIDs_snaptree[ListOutputSnaps[n] * Ntrees + tree];
      OffsetIDs = OffsetIDs_snaptree[ListOutputSnaps[n] * Ntrees + tree];
      IdList = mymalloc(sizeof(long long) * Nids, "idlist");
      PosList = mymalloc(3 * sizeof(float) * Nids, "poslist");
      VelList = mymalloc(3 * sizeof(float) * Nids, "vellist");
      myfseek(fd, OffsetIDs * sizeof(long long), SEEK_CUR);
      myfread(IdList, sizeof(long long), Nids, fd);
      myfseek(fd,
	      (TotIds - Nids) * sizeof(long long) + OffsetIDs * (3 * sizeof(float) -
								 sizeof(long long)), SEEK_CUR);
      myfread(PosList, 3 * sizeof(float), Nids, fd);
      myfseek(fd, (TotIds - Nids) * sizeof(float) * 3, SEEK_CUR);
      myfread(VelList, 3 * sizeof(float), Nids, fd);
#ifndef USE_MEMORY_TO_MINIMIZE_IO
      fclose(fd);
#endif
#endif

#ifdef USE_MEMORY_TO_MINIMIZE_IO
      fd = (FILE *) (10 + n);
      offset_galsnapdata[n] = 0;
#ifdef OUTPUT_MOMAF_INPUTS 
      fd_momaf = (FILE *) (10000 + n);
      offset_momafdata[n] = 0;
#endif
#else
      sprintf(buf, "%s/%s_z%1.2f_%d", OutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[n]], filenr);
      if(!(fd = fopen(buf, "r+")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}
#ifdef OUTPUT_MOMAF_INPUTS 
      sprintf(buf, "%s/%s_%d.%d", OutputDir,FileNameGalaxies,filenr,ListOutputSnaps[n]);
      if(!(fd_momaf = fopen(buf, "r+")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}
#endif
#endif
      myfseek(fd, (2 + Ntrees) * sizeof(int), SEEK_CUR);
      myfseek(fd, TotGalaxies[n] * sizeof(struct GALAXY_OUTPUT), SEEK_CUR);
#ifdef OUTPUT_MOMAF_INPUTS
      myfseek(fd_momaf,sizeof(int), SEEK_CUR);
      myfseek(fd_momaf, TotGalaxies[n] * sizeof(struct MOMAF_INPUTS), SEEK_CUR);
#endif
      for(i = 0; i < NumGals; i++)
	{
	  if(HaloGal[i].SnapNum == ListOutputSnaps[n])
	    {
	      prepare_galaxy_for_output(n, filenr, tree, &HaloGal[i], &galaxy_output);
	      myfwrite(&galaxy_output, sizeof(struct GALAXY_OUTPUT), 1, fd);
#ifdef OUTPUT_MOMAF_INPUTS
	      prepare_galaxy_for_momaf(n,filenr,tree,&HaloGal[i],&momaf_inputs);
	      myfwrite(&momaf_inputs,sizeof(struct MOMAF_INPUTS), 1, fd_momaf);
#endif
	      TotGalaxies[n]++;
	      TreeNgals[n][tree]++;
	    }
	}

#ifndef USE_MEMORY_TO_MINIMIZE_IO
      fclose(fd);
#ifdef OUTPUT_MOMAF_INPUTS
      fclose(fd_momaf);
#endif
#endif

#ifdef UPDATETYPETWO
      myfree(VelList);
      myfree(PosList);
      myfree(IdList);
#endif
    }

}



void finalize_galaxy_file(int filenr)
{
	int one, size_of_struct;
	one = 1;
	size_of_struct = sizeof(struct GALAXY_OUTPUT);

#ifndef USE_MEMORY_TO_MINIMIZE_IO
  char buf[1000];
#endif
  FILE *fd;

#ifndef GALAXYTREE
  int n;

  for(n = 0; n < NOUT; n++)
    {
#ifdef USE_MEMORY_TO_MINIMIZE_IO
      fd = (FILE *) (10 + n);
      offset_galsnapdata[n] = 0;
#else
      sprintf(buf, "%s/%s_z%1.2f_%d", OutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[n]], filenr);

      if(!(fd = fopen(buf, "r+")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}
#endif

      myfwrite(&Ntrees, sizeof(int), 1, fd);
      myfwrite(&TotGalaxies[n], sizeof(int), 1, fd);
      myfwrite(TreeNgals[n], sizeof(int), Ntrees, fd);

#ifndef USE_MEMORY_TO_MINIMIZE_IO
      fclose(fd);
#else
      write_galaxy_data_snap(n, filenr);
#endif
    }
#else //GALAXYTREE

#ifdef USE_MEMORY_TO_MINIMIZE_IO
  fd = (FILE *) 4;
  offset_galaxydata = 0;
#else //USE_MEMORY_TO_MINIMIZE_IO
#ifdef NEW_IO
  fd = output_file;
  rewind(fd);
#else
  fd = open_outputtree_file(filenr, "r+");
#endif //NEW_IO
#endif //USE_MEMORY_TO_MINIMIZE_IO
  // for DB compatible output, write 1, then size of a struct, then number of galaxies
  myfwrite(&one,sizeof(int),1,fd);
  myfwrite(&size_of_struct,sizeof(int),1,fd);
  myfwrite(&TotGalCount, sizeof(int), 1, fd);

#ifndef USE_MEMORY_TO_MINIMIZE_IO
#ifndef NEW_IO
  fclose(fd);
#endif
#else
  write_all_galaxy_data(filenr);
#endif
#endif
}

#ifdef OUTPUT_MOMAF_INPUTS
void finalize_momaf_file(int filenr)
{
#ifndef USE_MEMORY_TO_MINIMIZE_IO
  char buf[1000];
#endif
  FILE *fd;
  int n;
  
  for(n = 0; n < NOUT; n++)
    {
#ifdef USE_MEMORY_TO_MINIMIZE_IO
      fd = (FILE *) (10000+n);
      offset_momafdata[n] = 0;
#else
      sprintf(buf, "%s/%s_%d.%d", OutputDir, FileNameGalaxies, filenr, ListOutputSnaps[n]);
      if(!(fd = fopen(buf, "r+")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}
#endif
      myfwrite(&TotGalaxies[n], sizeof(int), 1, fd);
#ifndef USE_MEMORY_TO_MINIMIZE_IO
      fclose(fd);
#else
      write_galaxy_for_momaf(n, filenr);
#endif
    }
}
#endif

#ifdef UPDATETYPETWO

void get_coordinates(float *pos, float *vel, long long ID, int tree, int halonr, int snapnum)
{
  int m, k, start, nids;

  start = OffsetIDs_halo[TreeFirstHalo[tree] + halonr] - OffsetIDs;
  nids = CountIDs_halo[TreeFirstHalo[tree] + halonr];

  while(nids > 0)
    {
      m = nids / 2;
      if(IdList[start + m] == ID)
	{
	  for(k = 0; k < 3; k++)
	    {
	      pos[k] = PosList[3 * (start + m) + k];
	      vel[k] = sqrt(AA[snapnum]) * VelList[3 * (start + m) + k];  /* to convert to peculiar velocity */
	    }

	  if(pos[0] == 0 && pos[1] == 0 && pos[2] == 0)
	    {
	      printf
		("This treeaux-files does not (yet) contain the coordinates\n for the desired output time!\n");
	      fflush(stdout);
	      exit(0);
	    }

	  return;
	}

      if(IdList[start + m] < ID)
	{
	  nids -= m;
	  start += m;
	}
      else
	{
	  nids = m;
	}
    }

  printf("ID not found! - What's going on?  ID=%d\n", (int) ID);
  fflush(stdout);
  exit(0);
  return;
}

#endif



double xx=1.11, yy=1.23,zz=1.54;

void prepare_galaxy_for_output(int n, int filenr, int tree, struct GALAXY *g, struct GALAXY_OUTPUT *o)
{
  int j,ll;
#ifdef GALAXYTREE
  long long big;
  double scalefac;
  int i, k, lenmax, tmpfirst, next;
#endif

    double f1, f2, fmet1, fmet2;
    double X1, X2,age;
    double time, metallicity;
    int outputbin, idx, metindex, tabindex, zindex;


#ifndef NO_PROPS_OUTPUTS
	o->InfallVmax = g->InfallVmax;

	o->InfallSnap = g->InfallSnap;
	o->HotRadius =  g->HotRadius;
	

#ifdef UPDATETYPETWO
  if(g->Type == 2)
    {
      get_coordinates(g->Pos, g->Vel, g->MostBoundID, tree, g->HaloNr, g->SnapNum);

      for(j = 0; j < 3; j++)
      	{
	  o->Pos[j]=g->MergCentralPos[j] + (-g->MergCentralPos[j] + g->Pos[j])*sqrt(g->MergTime/g->OriMergTime);
	  if (o->Pos[j] < 0 || o->Pos[j] > 500)
	    {
	      printf("wrong inpos of type 2\n");
	      exit(0);
	    }

	}
    }
#endif
  if(g->Type == 2 || (g->Type == 1 && g->MergeOn == 1))
    {
      o->OriMergTime=g->OriMergTime;
      o->MergTime = g->MergTime;
    }
  else
    {
      o->OriMergTime=0.0;
      o->MergTime = 0.0;
    }
      
  o->Type = g->Type;
 
#ifndef GALAXYTREE
  o->HaloIndex = g->HaloNr;
#endif
  o->SnapNum = g->SnapNum;
  /*new model*/
  o->CentralMvir = get_virial_mass(Halo[g->HaloNr].FirstHaloInFOFgroup,1);
  /*de lucia*/
  // o->CentralMvir = get_virial_mass(Halo[g->HaloNr].FirstHaloInFOFgroup);
  for(j = 0; j < 3; j++)
    {
      if (g->Type !=2)	
	o->Pos[j] = g->Pos[j];
      o->Vel[j] = g->Vel[j];
      o->GasSpin[j] = g->GasSpin[j];
      o->StellarSpin[j] = g->StellarSpin[j];
   }
  //o->Haloj=g->Haloj;
  o->Len = g->Len;
  o->Mvir = g->Mvir;
  o->Rvir = g->Rvir;
  o->Vvir = g->Vvir;
  o->Vmax = g->Vmax;
  /*ram pressure*/
  /*if (g->InfallGasFrac <1.e-6)*/
  /* o->retainfac=0;*/
  /*else*/
    /* o->retainfac=g->HotGas/g->InfallMass/g->InfallGasFrac;*/

  //  o->ReheatedMass = g->ReheatedMass;
  // o->MetalsReheatedMass = g->MetalsReheatedMass;

  //- o->CoolingGas = g->CoolingGas;
  // o->HotStart = g->HotStart;
  // o->HotStripping = g->HotStripping;
  // o->MetalsHotStart = g->MetalsHotStart;
  // o->MetalsHotStripping = g->MetalsHotStripping;
  // o->HotICM  =g->HotICM;
  // o->MetalsHotICM = g->MetalsHotICM;

#ifdef GALAXYTREE
  //-o->HotRadius = g->HotRadius;
#endif
  o->ColdGas = g->ColdGas;
  o->BulgeSize=g->BulgeSize;
  o->StellarMass = g->StellarMass;
  o->BulgeMass = g->BulgeMass;
  o->HotGas = g->HotGas;
  o->EjectedMass = g->EjectedMass;
  o->BlackHoleMass = g->BlackHoleMass;
  //-  o->MergeSat = g->MergeSat;

#ifdef GALAXYTREE
  o->DisruptOn = g->DisruptOn;
  o->MergeOn = g->MergeOn;
#endif
  o->MetalsColdGas = g->MetalsColdGas;
  o->MetalsStellarMass = g->MetalsStellarMass;
  o->MetalsBulgeMass = g->MetalsBulgeMass;
  o->MetalsHotGas = g->MetalsHotGas;
  o->MetalsEjectedMass = g->MetalsEjectedMass;

#ifdef STAR_FORMATION_HISTORY
  o->sfh_ibin=g->sfh_ibin;

  o->sfh_time[0]=(Age[0]-NumToTime(g->SnapNum))*UnitTime_in_years/Hubble_h; //g->sfh_age;
  o->sfh_StellarMass[0]=g->sfh_StellarMass[0];
  o->sfh_MetalsStellarMass[0]=g->sfh_MetalsStellarMass[0];
  for (j=1;j<=o->sfh_ibin;j++) {
	  o->sfh_time[j]=(Age[0]-NumToTime(g->SnapNum))*UnitTime_in_years/Hubble_h - g->sfh_t[j-1]*SFH_TIME_INTERVAL;
	  o->sfh_StellarMass[j]=g->sfh_StellarMass[j];
	  o->sfh_MetalsStellarMass[j]=g->sfh_MetalsStellarMass[j];
  }
  for (j=o->sfh_ibin+1;j<SFH_NBIN;j++) {
	  o->sfh_time[j]=0.;
	  o->sfh_StellarMass[j]=0.;
	  o->sfh_MetalsStellarMass[j]=0.;
  }
#endif
  for(j = 0; j < NMAG; j++){
  o->TestMag[j] =0;
  o->TestObsMag[j] =0;
  }
  o->TestMass=0;

  if(g->SnapNum == ListOutputSnaps[n]){
	  for(j = 0; j < NMAG-4; j++){

		  for(ll=0;ll<=o->sfh_ibin;ll++){


			  X1 =  o->sfh_StellarMass[ll]*0.1/ (Hubble_h*(1 - RecycleFraction));
			  X2 = -0.4 * M_LN10;


			  if(o->sfh_StellarMass[ll] > 0.0){


				 if(ll<o->sfh_ibin)
					  time=((o->sfh_time[ll])*Hubble_h/1.e12+(o->sfh_time[ll+1])*Hubble_h/1.e12)/2.;
				  else
					  time=((o->sfh_time[ll])*Hubble_h/1.e12/2.);

				 //printf("time1=%f\n",(o->sfh_time[ll])*Hubble_h/1.e6);
				 //if(ll==o->sfh_ibin) time=(NumToTime(g->SnapNum-1)-NumToTime(g->SnapNum))/40.0;

			    if(time < (NumToTime(g->SnapNum-1)-NumToTime(g->SnapNum))/20.0)
			    	time=(NumToTime(g->SnapNum-1)-NumToTime(g->SnapNum))/40.0;
			    else
			    	if(time < 4./2.*(NumToTime(g->SnapNum-1)-NumToTime(g->SnapNum))/20.0)
				    		time=3.* (NumToTime(g->SnapNum-1)-NumToTime(g->SnapNum))/40.0;

        //printf("Snap=%d time=%f Snap=%d time=%f\n", g->SnapNum, NumToTime(g->SnapNum), g->SnapNum-1,NumToTime(g->SnapNum-1));

				  find_interpolated_lum(time,0.,
						  o->sfh_MetalsStellarMass[ll]/o->sfh_StellarMass[ll],
						  &metindex, &tabindex, &f1, &f2, &fmet1, &fmet2);

				  o->TestMag[j] +=
						  X1 * exp(X2 * (fmet1 * (f1 * MagTableZz[j][metindex][0][tabindex] +
								  f2 * MagTableZz[j][metindex][0][tabindex + 1]) 	+
								  fmet2 * (f1 * MagTableZz[j][metindex + 1][0][tabindex] +
										  f2 * MagTableZz[j][metindex + 1][0][tabindex + 1])));

				 zindex   = 63 - ListOutputSnaps[n];
				 o->TestObsMag[j] +=
				 		X1 * exp(X2 * (fmet1 * (f1 * MagTableZz[j][metindex][zindex][tabindex] +
				 					f2 * MagTableZz[j][metindex][zindex][tabindex + 1]) +
				 			       fmet2 * (f1 * MagTableZz[j][metindex + 1][zindex][tabindex] +
				 					f2 * MagTableZz[j][metindex + 1][zindex][tabindex + 1])));

				// printf("snap=%d mass=%f Time=%f (10^6yr) \n",o->SnapNum,log10(o->sfh_StellarMass[ll]*1.e10), time*1e6);
			  }
		  }

		  o->TestMag[j]=lum_to_mag(o->TestMag[j]);
		  o->TestObsMag[j]=lum_to_mag(o->TestObsMag[j]);
		  //printf("mag=%f\n",o->TestMag[j]);

		  if(lum_to_mag(g->ObsLum[0][n])-o->TestObsMag[0]>1.0)
		  {
			  printf("Mag=%f TestMag=%f\n",lum_to_mag(g->ObsLum[0][n]),o->TestObsMag[0]);
			  for(ll=0;ll<=o->sfh_ibin;ll++){
				  printf("mass=%f\n",log10(o->sfh_StellarMass[ll]*1.e10));
			  }
		  }

		 // if(g->SnapNum == 63 && o->TestMag[0]<-22.)
			//  printf("mag=%f\n",o->TestMag[0]);
	  }

/*
 * snap=27 mass=5.915346 Time=19.284953 (10^6yr)
snap=27 mass=5.937815 Time=12.714953 (10^6yr)
snap=27 mass=5.958786 Time=6.144953 (10^6yr)
 *
 * */

	 /* for(ll=0;ll<o->sfh_ibin+1;ll++)
	  {
		  o->TestMass += o->sfh_StellarMass[ll];
		 printf("snap=%d mass=%f Time=%f (10^6yr) \n",o->SnapNum, log10(1e-10+o->sfh_StellarMass[ll]*1e10),o->sfh_time[ll]/1e6);
	  }*/
	  //}
  }


  // if(g->Type<2 && HaloIDs[g->HaloNr].HaloID >= 1003318 && HaloIDs[g->HaloNr].HaloID < 1003360)
 	//  {
 	//		  printf("type=%d HaloNr=%lld Mag=%f \n", g->Type, HaloIDs[g->HaloNr].HaloID, lum_to_mag(g->Lum[0][ListOutputSnaps[n]]));
      //  printf("Snap=%d\n",g->SnapNum);


  //if(ll<o->sfh_ibin)
			//	  time=min(NumToTime(g->SnapNum-1)-NumToTime(g->SnapNum)/20.0,((o->sfh_time[ll])*Hubble_h/1.e12+(o->sfh_time[ll+1])*Hubble_h/1.e12)/2.);
			 // else
			//	  time=min(NumToTime(g->SnapNum-1)-NumToTime(g->SnapNum)/20.0,(o->sfh_time[ll])*Hubble_h/1.e12/2.);

  // if(ll<o->sfh_ibin)mytime=(o->sfh_time[ll]*Hubble_h/1e6+o->sfh_time[ll+1]*Hubble_h/1e6)/2.;
 			 // else mytime=(o->sfh_time[ll])*Hubble_h/1e6/2.;
 			/*  if(g->Type == 0 && g->HaloNr == 1986)
 			  {
 				   printf("HaloNr=%d FinalTime=%0.8f Mass=%0.8f\n",g->HaloNr, mytime, o->sfh_StellarMass[ll]*1e4);
 				   if(mytime < 201) aa+=o->sfh_StellarMass[ll]*1e4;
 				     if(ll==o->sfh_ibin)
 				       printf(" aa=%f \n",aa);

 			  }*/

  /*  NOTE: in Msun/yr */
#ifdef UseFullSfr
  for(j = 0; j < MAXSNAPS; j++)
    {
      o->Sfr[j] = g->Sfr[j] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
      o->SfrBulge[j] = g->SfrBulge[j] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
    }
#else
#ifdef SAVE_MEMORY
  o->Sfr = g->Sfr* UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
  o->SfrBulge =  g->SfrBulge * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
#else
  o->Sfr = g->Sfr[n] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
  o->SfrBulge = g->SfrBulge[n] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
#endif
#endif
  //- o->StarMerge=g->StarMerge;
  o->XrayLum = g->XrayLum;
  o->GasDiskRadius = g->GasDiskRadius;
  o->StellarDiskRadius = g->StellarDiskRadius;
  //o->halfradius = g->halfradius;
  // o->periradius = g->periradius;
  o->CoolingRadius = g->CoolingRadius;
 o->ICM = g->ICM;
  o->MetalsICM = g->MetalsICM;

#endif
 

#ifdef OUTPUT_REST_MAGS
  /* Luminosities in various bands */
  for(j = 0; j < NMAG; j++)
    {
      o->Mag[j]      = lum_to_mag(g->Lum[j][n]);
      o->MagBulge[j] = lum_to_mag(g->LumBulge[j][n]);
      o->MagDust[j]  = lum_to_mag(g->LumDust[j][n]);
#ifdef ICL
      o->MagICL[j]   = lum_to_mag(g->ICLLum[j][n]);
#endif
    }
  if(g->StellarMass > 0.0)
    {
      o->MassWeightAge = g->MassWeightAge[n] / g->StellarMass ;
      o->MassWeightAge = o->MassWeightAge / 1000. *UnitTime_in_Megayears / Hubble_h;
    }
  else
    {
      o->MassWeightAge = 0.;
    }
#endif
#ifdef OUTPUT_OBS_MAGS
#ifdef COMPUTE_OBS_MAGS
  /* Luminosities in various bands */
  for(j = 0; j < NMAG; j++)
    {
      o->ObsMag[j]       = lum_to_mag(g->ObsLum[j][n]);
      o->ObsMagBulge[j]  = lum_to_mag(g->ObsLumBulge[j][n]);
      o->ObsMagDust[j]   = lum_to_mag(g->ObsLumDust[j][n]);
#ifdef OUTPUT_MOMAF_INPUTS
      o->dObsMag[j]      = lum_to_mag(g->dObsLum[j][n]);
      o->dObsMagBulge[j] = lum_to_mag(g->dObsLumBulge[j][n]);
      o->dObsMagDust[j]  = lum_to_mag(g->dObsLumDust[j][n]);
#endif
    }
#endif
#endif





#ifndef NO_PROPS_OUTPUTS
#ifdef GALAXYTREE
  o->GalID = g->GalID;
  o->FirstProgGal = g->FirstProgGal;
  o->NextProgGal = g->NextProgGal;
  o->LastProgGal = g->LastProgGal;
  o->MainLeafId = g->MainLeaf;
  o->TreeRootId = g->TreeRoot;
  o->HaloID = HaloIDs[g->HaloNr].HaloID;

  o->DescendantGal = g->DescendantGal;
  o->Redshift = ZZ[g->SnapNum];

// TODO check how to treat next setting of multiplier prefixes properly
// Use make/input parameters?
#ifdef MRII
  // TBD assumes < 10^9 galaxies per tree, is that always correct?
  big = (((filenr * (long long) 1000000) + tree) * (long long) 1000000000);
#else
  big = (((filenr * (long long) 1000000) + tree) * (long long) 1000000);
#endif

  o->FileTreeNr = big;

  o->GalID += big;

  if(o->FirstProgGal >= 0)
    o->FirstProgGal += big;

  if(o->LastProgGal >= 0)
    o->LastProgGal += big;
  else
    o->LastProgGal = o->GalID;

  if(o->MainLeafId >= 0)
    o->MainLeafId += big;
  else
    o->MainLeafId = o->GalID;

  if(o->TreeRootId >= 0)
    o->TreeRootId += big;
  else  // TODO this should not occur, throw an error if it does?
    o->TreeRootId = -1;

  if(o->NextProgGal >= 0)
    o->NextProgGal += big;

  if(o->DescendantGal >= 0)
    o->DescendantGal += big;
  // Set pointer ("foreign key") to central galaxy of FOF group this galaxy is in.
  // TBD next seems to work on millimil.
  // Assumes that firstgalaxy in halo is its central galaxy.
  // This is not generaly true (see comments in main.c::join_galaxies_of_progenitors).
  // Therefore an extra test is included and
  o->FOFCentralGal = HaloGal[HaloAux[Halo[g->HaloNr].FirstHaloInFOFgroup].FirstGalaxy].GalID+big;
  // test whether the FOF central galaxy actually is a type-0
  // if not, search for the type-0 galaxy
  if(HaloGal[HaloAux[Halo[g->HaloNr].FirstHaloInFOFgroup].FirstGalaxy].Type != 0)
  {
  	printf("a FOFCentralGal, galaxyId=%lld, has been assigned that is no type 0!\n",o->FOFCentralGal);
  	// TBD can we do something else ?
  	//  For example store FOFCentralGalaxy on HaloAux, to be set in update_centralgal ?
  	// That introduces memory overhead for what should be very rare occurrence.
  	// TODO check this code is correct.
  	for(i = 0; i < HaloAux[Halo[g->HaloNr].FirstHaloInFOFgroup].NGalaxies; i++)
  		if(HaloGal[HaloAux[Halo[g->HaloNr].FirstHaloInFOFgroup].FirstGalaxy+i].Type == 0)
  		{
  			o->FOFCentralGal = HaloGal[HaloAux[Halo[g->HaloNr].FirstHaloInFOFgroup].FirstGalaxy+i].GalID+big;
  			break;
  		}
  }


  scalefac = 1.0 / BoxSize;

  xx = g->Pos[0] * scalefac + 1.0;
  yy = g->Pos[1] * scalefac + 1.0;
  zz = g->Pos[2] * scalefac + 1.0;
 

  i = DOUBLE_to_HASHBITS(xx);
  j = DOUBLE_to_HASHBITS(yy);
  k = DOUBLE_to_HASHBITS(zz);


  o->PeanoKey = peano_hilbert_key(i, j, k, Hashbits);
  //printf("i %d, j %d, k %d o->PeanoKey %d \n",i,j,k,o->PeanoKey);
  // add info for Gerard
#ifdef MRII
  big = g->SnapNum* (long long) 10000000000 + Halo[g->HaloNr].FileNr* (long long) 1000000 + Halo[g->HaloNr].SubhaloIndex;
#else
  big = g->SnapNum* (long long) 1000000000000 + Halo[g->HaloNr].FileNr* (long long) 100000000 + Halo[g->HaloNr].SubhaloIndex;
#endif

  o->SubID = big;

  tmpfirst = Halo[g->HaloNr].FirstHaloInFOFgroup;
  lenmax   = 0;
  next     = tmpfirst;
  while(next != -1) 
    {
      if(Halo[next].Len > lenmax)
	{
	  lenmax   = Halo[next].Len;
	  tmpfirst = next;
	}
      next = Halo[next].NextHaloInFOFgroup;
    }  
  // TODO check next relation for both Millennium and Millennium-II !
#ifdef MRII
  big = g->SnapNum* (long long) 10000000000 + Halo[tmpfirst].FileNr* (long long) 1000000 + Halo[tmpfirst].SubhaloIndex;
#else
  big = g->SnapNum* (long long) 1000000000000 + Halo[tmpfirst].FileNr* (long long) 100000000 + Halo[tmpfirst].SubhaloIndex;
#endif
  o->MMSubID = big;

#endif
#endif

#ifndef NO_PROPS_OUTPUTS
#ifdef FIX_OUTPUT_UNITS
  fix_units_for_ouput(o);
#endif
#endif

}


/**@brief Remove h from units of galaxy properties.
 * If desired (makefile option FIX_OUTPUT_UNITS is set), the output properties of the galaxies can be scaled to
 * physical units excluding any factors of h == Hubble/100 km/s/Mpc.
 */
void fix_units_for_ouput(struct GALAXY_OUTPUT *o)
{
	int j;
	// TBD do we want to change the positions that have corresponding value in halo tables that DO include the /h?
	o->Pos[0] /= Hubble_h;
	o->Pos[1] /= Hubble_h;
	o->Pos[2] /= Hubble_h;
	// TBD do we want to change the halo properties that have corresponding value in halo tables that DO include the /h?
	o->CentralMvir /= Hubble_h;
	o->Mvir /= Hubble_h;
	o->Rvir /= Hubble_h;
	// end TBD

	o->GasSpin[0] /= Hubble_h;
	o->GasSpin[1] /= Hubble_h;
	o->GasSpin[2] /= Hubble_h;
	o->StellarSpin[0] /= Hubble_h;
	o->StellarSpin[1] /= Hubble_h;
	o->StellarSpin[2] /= Hubble_h;
	o->HotRadius /= Hubble_h;
	o->ColdGas /= Hubble_h;
	o->StellarMass /= Hubble_h;
	o->BulgeMass /= Hubble_h;
	o->HotGas /= Hubble_h;
	o->EjectedMass /= Hubble_h;
	o->BlackHoleMass /= Hubble_h;
	o->ICM /= Hubble_h;
	o->MetalsColdGas /= Hubble_h;
	o->MetalsStellarMass /= Hubble_h;
	o->MetalsBulgeMass /= Hubble_h;
	o->MetalsHotGas /= Hubble_h;
	o->MetalsEjectedMass /= Hubble_h;
	o->MetalsICM /= Hubble_h;
	o->BulgeSize /= Hubble_h;
	o->StellarDiskRadius /= Hubble_h;
	o->GasDiskRadius /= Hubble_h;
	o->CoolingRadius /= sqrt(Hubble_h); //rcool units are "wrong" in the code (/h^1/2)
#ifdef STAR_FORMATION_HISTORY
	for(j=0;j<=o->sfh_ibin;j++) {
		o->sfh_StellarMass[j] /= Hubble_h;
		o->sfh_MetalsStellarMass[j] /= Hubble_h;
	}
#endif
}



#ifdef OUTPUT_MOMAF_INPUTS
void prepare_galaxy_for_momaf(int n, int filenr, int tree, struct GALAXY *g, struct MOMAF_INPUTS *o)
{
  int j;
  float dm,dmb,dz;
#ifdef GALAXYTREE
  long long big;
#endif

#ifdef UPDATETYPETWO
  if(g->Type == 2)
    {
      get_coordinates(g->Pos, g->Vel, g->MostBoundID, tree, g->HaloNr, g->SnapNum);

      for(j = 0; j < 3; j++)
        {
          o->Pos[j]=g->MergCentralPos[j] + (-g->MergCentralPos[j] + g->Pos[j])*sqrt(g->MergTime/g->OriMergTime);
          if (o->Pos[j] < 0 || o->Pos[j] > 500)
            {
              printf("wrong inpos of type 2\n");
              exit(0);
            }

        }
      
    }
#endif

  o->SnapNum = g->SnapNum;
  for(j = 0; j < 3; j++)
    {
      if (g->Type != 2)
      o->Pos[j] = g->Pos[j];
      o->Vel[j] = g->Vel[j];
    }
#ifdef COMPUTE_OBS_MAGS
  for(j = 0; j < NMAG; j++)
    {
      o->ObsMagBulge[j]  = lum_to_mag(g->ObsLumBulge[j][n]);
      o->ObsMagDust[j]   = lum_to_mag(g->ObsLumDust[j][n]);
      dmb = lum_to_mag(g->dObsLumBulge[j][n]);
      dm  = lum_to_mag(g->dObsLumDust[j][n]);
      dz  = ZZ[g->SnapNum-1] - ZZ[g->SnapNum];
      o->dObsMagBulge[j] = (dmb - o->ObsMagBulge[j])/dz;
      o->dObsMagDust[j]  = (dm - o->ObsMagDust[j])/dz;
    }
#endif

#ifdef GALAXYTREE
  o->GalID = g->GalID;
  o->HaloID = HaloIDs[g->HaloNr].HaloID;
// TODO check how to treat next, currently ok for Millennium-II
//   big = (((filenr * (long long) 1000000) + tree) * (long long) 1000000);
  big = (((filenr * (long long) 1000000) + tree) * (long long) 1000000000);
  o->GalID += big;
#else
  o->GalID  = (long long) (2);
  o->HaloID = (long long) (2);
#endif
}
#endif


#ifdef GALAXYTREE
void save_galaxy_tree(int filenr, int tree)
{
#ifndef USE_MEMORY_TO_MINIMIZE_IO
  char buf[1000];
#endif
  FILE *fd;
  int i, n, p, num;
  struct GALAXY_OUTPUT galaxy_output;
#ifdef OUTPUT_MOMAF_INPUTS
  FILE *fd_momaf;
  struct MOMAF_INPUTS momaf_inputs;
#endif

  for(i = 0; i < NumGals; i++)
    {
      HaloGal[i].Done = 0;
      HaloGal[i].LastProgGal = -1;
      HaloGal[i].MainLeaf = -1;
      HaloGal[i].TreeRoot = -1;
    }
  for (num = LastSnapShotNr; num >= 0; num--) {
		for (i = 0; i < NumGals; i++) {
			if (HaloGal[i].SnapNum == num)
				if (HaloGal[i].Done == 0) {
					walk(i);
				}
		}
	}
  for(i = 0; i < NumGals; i++)
    {
      p = HaloGal[i].FirstProgGal;
      while(p >= 0)
	{
	  HaloGal[p].DescendantGal = i;
	  p = HaloGal[p].NextProgGal;
	}
    }
  for(i = 0; i < NumGals; i++)
    {
      if(HaloGal[i].FirstProgGal >= 0)
	HaloGal[i].FirstProgGal = HaloGal[HaloGal[i].FirstProgGal].GalID;

      if(HaloGal[i].LastProgGal >= 0)
	HaloGal[i].LastProgGal = HaloGal[HaloGal[i].LastProgGal].GalID;

      if(HaloGal[i].MainLeaf >= 0)
	HaloGal[i].MainLeaf = HaloGal[HaloGal[i].MainLeaf].GalID;

      if(HaloGal[i].TreeRoot >= 0)
	HaloGal[i].TreeRoot = HaloGal[HaloGal[i].TreeRoot].GalID;

      if(HaloGal[i].NextProgGal >= 0)
	HaloGal[i].NextProgGal = HaloGal[HaloGal[i].NextProgGal].GalID;

      if(HaloGal[i].DescendantGal >= 0)
	HaloGal[i].DescendantGal = HaloGal[HaloGal[i].DescendantGal].GalID;

    }

  for(n = 0; n < NOUT; n++)
    {
#ifdef UPDATETYPETWO
#ifdef USE_MEMORY_TO_MINIMIZE_IO
      fd = (FILE *) 2;
      offset_auxdata = 0;
#else
#ifdef NEW_IO
      fd = treeaux_file;
      rewind(fd);
#else
      fd = open_treeaux_file(filenr);
#endif //NEW_IO
#endif //USE_MEMORY_TO_MINIMIZE_IO
      myfseek(fd, 4 * sizeof(int) + 2 * TotSnaps * sizeof(int)
	      + 2 * TotSnaps * Ntrees * sizeof(int) + 2 * sizeof(int) * NtotHalos, SEEK_CUR);

      Nids = CountIDs_snaptree[ListOutputSnaps[n] * Ntrees + tree];
      OffsetIDs = OffsetIDs_snaptree[ListOutputSnaps[n] * Ntrees + tree];

      IdList = mymalloc(sizeof(long long) * Nids, "idlist");
      PosList = mymalloc(3 * sizeof(float) * Nids, "poslist");
      VelList = mymalloc(3 * sizeof(float) * Nids, "vellist");

      myfseek(fd, OffsetIDs * sizeof(long long), SEEK_CUR);

      myfread(IdList, sizeof(long long), Nids, fd);

      myfseek(fd,
	      (TotIds - Nids) * sizeof(long long) + OffsetIDs * (3 * sizeof(float) -
								 sizeof(long long)), SEEK_CUR);

      myfread(PosList, 3 * sizeof(float), Nids, fd);

      myfseek(fd, (TotIds - Nids) * sizeof(float) * 3, SEEK_CUR);

      myfread(VelList, 3 * sizeof(float), Nids, fd);
#ifndef USE_MEMORY_TO_MINIMIZE_IO
#ifndef NEW_IO
      fclose(fd);
#endif
#endif
#endif

#ifdef USE_MEMORY_TO_MINIMIZE_IO
      fd = (FILE *) 4;
      offset_galaxydata = 0;
#ifdef OUTPUT_MOMAF_INPUTS
      fd_momaf = (FILE *) (10000 + n);
      offset_momafdata[n] = 0;
#endif //OUTPUT_MOMAF_INPUTS
#else //USE_MEMORY_TO_MINIMIZE_IO
#ifdef NEW_IO
      fd = output_file;
#else //NEW_IO
      fd = open_outputtree_file(filenr,"r+");
#endif //NEW_IO
#ifdef OUTPUT_MOMAF_INPUTS
      sprintf(buf, "%s/%s_%d.%d", OutputDir, FileNameGalaxies, filenr,ListOutputSnaps[n]);
      if(!(fd_momaf = fopen(buf, "r+")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}
#endif
#endif

// ~ old code with header a simple integer for number of galaxies ~
      // jump over an int (used to store the number of galaxies at the end) and
      // and over all the galaxies written so far
//      myfseek(fd, sizeof(int) + TotGalCount * sizeof(struct GALAXY_OUTPUT), SEEK_CUR);

      // for DB compatible output, pad the first line with the size of a struct.
#ifndef NEW_IO
      myfseek(fd, (1+TotGalCount) * sizeof(struct GALAXY_OUTPUT), SEEK_CUR);
#endif

#ifdef OUTPUT_MOMAF_INPUTS
      myfseek(fd_momaf,sizeof(int),SEEK_CUR);
      myfseek(fd_momaf, TotGalaxies[n]*sizeof(struct MOMAF_INPUTS), SEEK_CUR);
#endif
      for(i = 0; i < NumGals; i++)
	{
	  if(HaloGal[i].SnapNum == ListOutputSnaps[n])
	    {
	      prepare_galaxy_for_output(n, filenr, tree, &HaloGal[i], &galaxy_output);
	      myfwrite(&galaxy_output, sizeof(struct GALAXY_OUTPUT), 1, fd);
#ifdef OUTPUT_MOMAF_INPUTS
	      prepare_galaxy_for_momaf(n,filenr,tree,&HaloGal[i],&momaf_inputs);
	      myfwrite(&momaf_inputs,sizeof(struct MOMAF_INPUTS), 1, fd_momaf);
	      TotGalaxies[n]++;
#endif
	      TotGalCount++;
	    }
	}

#ifndef USE_MEMORY_TO_MINIMIZE_IO
#ifndef NEW_IO
      fclose(fd);
#endif
#ifdef OUTPUT_MOMAF_INPUTS
      fclose(fd_momaf);
#endif
#endif

#ifdef UPDATETYPETWO
      myfree(VelList);
      myfree(PosList);
      myfree(IdList);
#endif
    }
  printf("TotGalCount= %d   at tree %d\n", TotGalCount, tree);
}



int walk(int nr) {
	int last;

	last = nr;

	if (HaloGal[nr].Done == 0) {
		HaloGal[nr].Done = 1;
		HaloGal[nr].GalID = GalCount++;
		// set treeroot.
		if (HaloGal[nr].TreeRoot == -1)
			HaloGal[nr].TreeRoot = nr;

		if (HaloGal[nr].FirstProgGal >= 0) {
			HaloGal[HaloGal[nr].FirstProgGal].TreeRoot = HaloGal[nr].TreeRoot;
			last = walk(HaloGal[nr].FirstProgGal);
			HaloGal[nr].MainLeaf = HaloGal[HaloGal[nr].FirstProgGal].MainLeaf;
		} else
			HaloGal[nr].MainLeaf = nr;

		HaloGal[nr].LastProgGal = last;

		if (HaloGal[nr].NextProgGal >= 0) {
			HaloGal[HaloGal[nr].NextProgGal].TreeRoot = HaloGal[nr].TreeRoot;
			last = walk(HaloGal[nr].NextProgGal);
		}
	}

	return last;
}

#endif
