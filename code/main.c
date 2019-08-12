#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "allvars.h"
#include "proto.h"

#ifdef MCMC
#include "mcmc_vars.h"
#include "mcmc_proto.h"
#endif

#ifdef ALL_SKY_LIGHTCONE
#include "lightcone.h"
#endif


/**@file main.c
 * @brief Controlling function of L-Galaxies plus SAM Construct Galaxies,
 *        Join Galaxies of progenitors, Evolve Galaxies and check
 *        Makefile options.
 *
 * Reads the parameter file given as an argument at run time:
 * read_parameter_file(argv[1]).
 *
 * Checks for consistency between some Makefile options: check_options().
 *
 * Runs init() to initialize some stuff.
 *
 * Ifdef MCMC - it calls Senna which will start the MCMC parameter estimation,
 * calling SAM for each step of the chain to run the semi-analytic model on
 * a set of merger trees and compute a likelihood for the model with those
 * parameters with respect to a set of observations - see mcmc.c
 *
 * Otherwise, SAM is called for each of the chosen dark matter tree files and
 * for each tree, on each treefile: reads tree, constructs the galaxies, saves
 * the galaxies, frees memory for galaxies and tree.
 *
 * @param filenr  Dark Matter tree file
 * @returns       0
 * */

/**@brief Main routine of L-Galaxies*/
int main(int argc, char **argv)
{
  int filenr, *FileToProcess, *TaskToProcess, nfiles;
  char buf[1000];
  time_t start, current;
  inputFile = argv[1];


#ifdef PARALLEL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);
#else
  NTask = 1;
  ThisTask = 0;
#endif //PARALLEL


#ifdef MCMC
  time(&global_starting_time);
#endif

  if(argc > 3)
    {
      printf("\n  usage: L-Galaxies <parameterfile>\n\n");
      endrun(0);
    }

 if (ThisTask == 0)
   printf("\n\nMakefile options used:\n%s\n",COMPILETIMESETTINGS);

  /*Reads the parameter file, given as an argument at run time. */
  read_parameter_file(argv[1]);

#ifdef MR_PLUS_MRII
  //Start with MR files and later change to MRII
  LastDarkMatterSnapShot=LastDarkMatterSnapShot_MR;
  sprintf(FileWithZList, "%s", FileWithZList_MR);
  sprintf(FileWithZList_OriginalCosm, "%s", FileWithZList_OriginalCosm_MR);
#endif

  mymalloc_init();

  sprintf(FinalOutputDir, "%s", OutputDir);
#ifndef MCMC
  if(argc == 3)
	sprintf(OutputDir, "%s", argv[2]);
#else
  FirstChainNumber=0;
  if(argc == 3)
  	FirstChainNumber=atoi(argv[2]);
#endif


  //copy input parameters and makefile options to output directory
  if (ThisTask == 0)
    {
      sprintf(buf, "mkdir %s/inputs/", FinalOutputDir);
      system(buf);
      sprintf(buf, "cp %s %s/inputs/", argv[1], FinalOutputDir);
      system(buf);
      sprintf(buf, "cp %s %s/inputs/", FileWithFilterNames, FinalOutputDir);
      system(buf);
      sprintf(buf, "cp ./My_Makefile_options %s/inputs/", FinalOutputDir);
      system(buf);
      sprintf(buf, "cp ./Makefile_compilers %s/inputs/", FinalOutputDir);
      system(buf);
    }

  //time(&start);

  /* check compatibility of some Makefile Options*/
  check_options();

#ifdef COMPUTE_SPECPHOT_PROPERTIES
  //for dust_model
  mu_seed = -150;
#endif

  init();

#ifdef OUTPUT_SFH
#ifdef PARALLEL
  if(ThisTask == 0)
#endif
    write_sfh_bins();
#endif



#ifndef MCMC
  nfiles=get_nr_files_to_process(ThisTask);
  FileToProcess=mymalloc("FileToProcess", sizeof(int) * nfiles);
  TaskToProcess=mymalloc("TaskToProcess", sizeof(int) * nfiles);
  assign_files_to_tasks(FileToProcess, TaskToProcess, ThisTask, NTask, nfiles);

  int file;
  for(file = 0; file < nfiles; file++)
    {
      if(ThisTask==TaskToProcess[file])
	filenr=FileToProcess[file];
      else
	continue;
#else //MCMC
  /* In MCMC mode only one file is loaded into memory
   * and the sampling for all the steps is done on it */    
      sprintf(SimulationDir, "%s/MergerTrees_%d/", SimulationDir, 0);
      for(filenr = MCMCTreeSampleFile; filenr <= MCMCTreeSampleFile; filenr++)
	{
#endif //MCMC
	  time(&start);

#ifdef PARALLEL
#ifdef MCMC
          //first load trees with task 0 and broadcast, than put a small delay on each task
          load_tree_table(filenr);

	  do
            time(&current);
          while(difftime(current, start) < 1.0 * ThisTask);
#else
          //first put a small delay on each task so they don't all read all the != files at the same time
          //time_t current;
          do
            time(&current);
          while(difftime(current, start) < 1.0 * ThisTask);
          load_tree_table(filenr);
#endif
#else
          load_tree_table(filenr);
#endif

#ifdef MCMC
      Senna(); // run the model in MCMC MODE
#else
      SAM(filenr); // run the model in NORMAL MODE
#endif

#ifdef MCMC
      break;	//break loop on files since the MCMC is done on a single file
#else
      time(&current);
      printf("\ndone tree file %d in %ldmin and %lds\n\n", filenr, (current - start)/60, (current - start)%60);

#endif //MCMC
      free_tree_table();
      //if temporary directory given as argument
       if(argc == 3)
         {
#ifdef GALAXYTREE
    	   sprintf(buf, "mv %s/%s_galtree_%d %s", OutputDir,FileNameGalaxies, filenr, FinalOutputDir);
#else
    	   sprintf(buf, "mv %s/%s_z*_%d %s", OutputDir,FileNameGalaxies, filenr, FinalOutputDir);
#endif
    	   system(buf);
         }
    }

#ifndef MCMC
  myfree(TaskToProcess);
  myfree(FileToProcess);
#endif

#ifdef PARALLEL
  MPI_Finalize();
#endif
  return 0;

}

/**@brief SAM() loops on trees and calls construct_galaxies.*/
#ifdef MCMC
double SAM(int filenr)
#else
void SAM(int filenr)
#endif
{
  int treenr, halonr;

#ifdef MCMC
  int ii;
  MCMC_GAL = mymalloc("MCMC_Gal", sizeof(struct MCMC_GALAXY) * MCMCAllocFactor);
  for(ii=0;ii<NOUT;ii++)
  	  TotMCMCGals[ii] = 0;

  if(Sample_Cosmological_Parameters==1)
  {
  	 reset_cosmology();
#ifdef HALOMODEL
  	 initialize_halomodel();
#endif
  }

#ifdef MR_PLUS_MRII
  change_dark_matter_sim("MR");
#else
  if(Sample_Cosmological_Parameters==1 || CurrentMCMCStep==1)
  	read_sample_info();
  else
  {
  	int snap, ii;
  	for(snap=0;snap<NOUT;snap++)
  		for(ii=0;ii<NFofsInSample[snap];ii++)
  			MCMC_FOF[ii].NGalsInFoF[snap]=0;
  }
#endif
#endif

  //to be used when we have tables for the scaling in any cosmology
  //read_scaling_parameters();

#ifndef MCMC
#ifdef GALAXYTREE
  create_galaxy_tree_file(filenr);
#else
  create_galaxy_files(filenr);
#endif
#ifdef ALL_SKY_LIGHTCONE
 int nr;
 for ( nr = 0; nr < NCONES; nr++)
   create_galaxy_lightcone_files(filenr, nr);
#endif
#endif

#ifdef GALAXYTREE
  FILE *fdg = fopen("treengal.dat", "w");
#endif

//***************************************************************************************
//***************************************************************************************

  //for(treenr = 0; treenr < NTrees_Switch_MR_MRII/5.; treenr++)
  //for(treenr = NTrees_Switch_MR_MRII; treenr < Ntrees; treenr++)
  for(treenr = 0; treenr < Ntrees; treenr++)
  //for(treenr = 0; treenr <10+1;treenr++)
  //for(treenr = 0; treenr <100;treenr++)
  {
      //printf("doing tree %d of %d (MR trees=%d)\n", treenr, Ntrees, NTrees_Switch_MR_MRII);
      //if(treenr%1000==0)
    //printf("doing tree %d of %d\n", treenr, Ntrees);
#ifdef MR_PLUS_MRII
  	if(treenr == NTrees_Switch_MR_MRII)
  		change_dark_matter_sim("MRII");
#endif

  	load_tree(treenr);
#ifdef MCMC
#ifdef PRELOAD_TREES
      if(Sample_Cosmological_Parameters==1 || CurrentMCMCStep==1)
#endif
#endif
        scale_cosmology(TreeNHalos[treenr]);

      gsl_rng_set(random_generator, filenr * 100000 + treenr);
      NumMergers = 0;
      NHaloGal = 0;
#ifdef GALAXYTREE
      NGalTree = 0;
      IndexStored = 0;
#endif
      int snapnum;
      //LastSnapShotNr is the highest output snapshot
      /* we process the snapshots now in temporal order 
       * (as a means to reduce peak memory usage) */
      for(snapnum = 0; snapnum <= LastSnapShotNr; snapnum++) {
	  //	printf("\n\nsnap=%d\n",snapnum);
#ifdef MCMC
    	  /* read the appropriate parameter list for current snapnum
    	   * into the parameter variables to be used in construct_galaxies */
    	  read_mcmc_par(snapnum);
#ifdef HALOMODEL
	  //because we need halo masses even for FOFs
	  //with no galaxies it needs to be done here
          assign_FOF_masses(snapnum, treenr);
#endif
#else
	  //used to allow parameter values to vary with redshift
    	  //re_set_parameters(snapnum);
#endif
    	  for(halonr = 0; halonr < TreeNHalos[treenr]; halonr++)
	      if(HaloAux[halonr].DoneFlag == 0 && Halo[halonr].SnapNum == snapnum)
		  construct_galaxies(filenr, treenr, halonr);
      }

      /* output remaining galaxies as needed */
      while(NHaloGal)
	  output_galaxy(treenr, 0);


#ifndef MCMC
#ifdef GALAXYTREE
      save_galaxy_tree_finalize(filenr, treenr);
#ifndef PARALLEL
      if((treenr/100)*100==treenr) printf("treenr=%d  TotGalCount=%d\n", treenr, TotGalCount);
#endif
      fflush(stdout);
      fprintf(fdg, "%d\n", NGalTree);
#endif
#else//ifdef MCMC
#ifdef PRELOAD_TREES
      if(Sample_Cosmological_Parameters==1)
    	  un_scale_cosmology(TreeNHalos[treenr]);
#endif
#endif
      free_galaxies_and_tree();
  }//loop on trees

#ifdef MCMC
  double lhood = get_likelihood();

#ifdef MR_PLUS_MRII
  free(MCMC_FOF);
#else
  if(Sample_Cosmological_Parameters==1 ||  CurrentMCMCStep==ChainLength)
  	free(MCMC_FOF);
#endif

#ifdef HALOMODEL
  if (Sample_Cosmological_Parameters==1) {
    gsl_spline_free(FofSpline);
    gsl_interp_accel_free(FofAcc);
    gsl_spline_free(SigmaSpline);
    gsl_interp_accel_free(SigmaAcc);
    gsl_spline_free(ellipSpline);
    gsl_interp_accel_free(ellipAcc);
    gsl_spline_free(PowSpline);
  } //if
#endif

  myfree(MCMC_GAL);
  return lhood;

#else //MCMC

#ifdef GALAXYTREE
  close_galaxy_tree_file();
#else
  close_galaxy_files();
#endif
#ifdef ALL_SKY_LIGHTCONE
  for (nr = 0; nr < NCONES; nr++)
    close_galaxy_lightcone_files(nr);
#endif

#ifdef DEBUG
  fclose(FdGalDebug);
#endif

  return;
#endif
}


/**@brief  construct_galaxies() recursively runs the semi-analytic model.
  *        For each halo it checks if its main progenitor has been done, then
  *        if all the halos in the FOF of its main progenitor have been
  *        done and then runs the SAM in the current halo. This means that
  *        for the first time its called it will walk up the tree into the
  *        haloes in the earliest snapshot.
  *
  *        When it finds a halo that needs to be done it calls
  *        join_galaxies_of_progenitors and evolve_galaxies. */
void construct_galaxies(int filenr, int treenr, int halonr)
{
  static int halosdone = 0;
  int prog;     // progenitor halo id
  int fofhalo;  // current halo id 
  int ngal;     // galaxy counter
  int FOF_centralgal;  // galaxy at the centre of this FOF group
  int p;        // temporary counter

  HaloAux[halonr].DoneFlag = 1;
  halosdone++;

  prog = Halo[halonr].FirstProgenitor;

  while(prog >= 0) //If halo has a progenitor
    {
      if(HaloAux[prog].DoneFlag == 0) //If progenitor hasn't been done yet
		construct_galaxies(filenr, treenr, prog);
      prog = Halo[prog].NextProgenitor;	//Jump to next halo in progenitors FOF
    }

  //Now check for the progenitors of all the halos in the current FOF group
  fofhalo = Halo[halonr].FirstHaloInFOFgroup;	//Starting at the first halo in current FOF
  if(HaloAux[fofhalo].HaloFlag == 0)	//If it hasn't been done
    {
      HaloAux[fofhalo].HaloFlag = 1;	//mark as to do
      while(fofhalo >= 0)	//go through all the halos in current FOF
        {
    	  prog = Halo[fofhalo].FirstProgenitor;
    	  while(prog >= 0)	//build its progenitors
    	    {
    		  if(HaloAux[prog].DoneFlag == 0)
    			construct_galaxies(filenr, treenr, prog);
    		  prog = Halo[prog].NextProgenitor;
    	    }

    	  fofhalo = Halo[fofhalo].NextHaloInFOFgroup;	//Jump to next halo in FOF
        }
    }

  /* At this point, the galaxies for all progenitors of this halo have been
   * properly constructed. Also, the galaxies of the progenitors of all other 
   * halos in the same FOF-group have been constructed as well. We can hence go
   * ahead and construct all galaxies for the subhalos in this FOF halo, and
   * evolve them in time. */
  fofhalo = Halo[halonr].FirstHaloInFOFgroup;
  FOF_centralgal=-1;
  if(HaloAux[fofhalo].HaloFlag == 1) {	//If it is marked as an halo to do.
    ngal = 0;                           // Start count of galaxies in this FoF group
    HaloAux[fofhalo].HaloFlag = 2;      // Mark halo as done.
     /* For all the halos in this FOF group, find all their progenitor
	 galaxies and construct a list Gals[] of their current
	 properties. */
    while(fofhalo >= 0) {
	ngal = join_galaxies_of_progenitors(fofhalo, ngal, &FOF_centralgal);
	fofhalo = Halo[fofhalo].NextHaloInFOFgroup;
    }
    // For the very rare occasions when CentralGal==-1, reset to FOF_centralgal
    for (p=0;p<ngal;p++)
        if (Gal[p].CentralGal == -1) {
	  //printf("For galaxy %d in halo %d, CentralGal==-1: resetting\n",p,halonr);
	    Gal[p].CentralGal = FOF_centralgal;
	    Gal[p].MergerCentralGal=FOF_centralgal;
	    int j;
	    for(j = 0; j < 3; j++) Gal[p].MergCentralPos[j] = Gal[FOF_centralgal].Pos[j];  
	}
    /* Check that have identified type 0 galaxy correctly.
     * This incidentally verifies that we have exactly 1 type 0. */
    if (FOF_centralgal < 0) terminate("FOF_centralgal < 0\n");
    if (FOF_centralgal >= ngal) terminate("FOF_centralgal >= ngal\n");
    if (Gal[FOF_centralgal].Type != 0) terminate("Gal[FOF_centralgal].Type != 0\n");
    for (p=0;p<ngal;p++)
	if (Gal[p].Type == 0 && p != FOF_centralgal) terminate("Type 0 not FOF_centralgal\n");
    /* Update FOFCentralgal for each galaxy.  This could be done within
     * join_galaxies_of_progenitors provided that the type 0 galaxy is in the
     * FirstHaloInFOFgroup.  Safer to do it here in case that should change. */
    for (p=0;p<ngal;p++) Gal[p].FOFCentralGal=FOF_centralgal;

#ifdef SAVE_FOFHALO
    for (p=0;p<ngal;p++) Gal[p].FileNr=filenr;
    for (p=0;p<ngal;p++) Gal[p].TreeNr=treenr;
#endif

    /*Evolve the Galaxies -> SAM! */
    evolve_galaxies(Halo[halonr].FirstHaloInFOFgroup, ngal, treenr);

    for (p =0;p<ngal;p++) mass_checks(p,"main.c",__LINE__);
  }
}


/**@brief join_galaxies_of_progenitors() updates the properties of the
 *        galaxy from the dark matter halo properties and deals with
 *        merging clocks. This routine is called by construct_galaxies
 *        for every halo in the FOF being constructed. When there is no
 *        galaxy in the Halo of FirstProgenitor, the first_occupied
 *        pointer is changed to a subhalo which have the maximum mass.
 *
 *        For a central galaxy it just updates its properties. For
 *        satellites it needs to know its most massive (or only progenitor)
 *        to keep track of the merging clock. It also finds the central
 *        galaxies into which galaxies should merge. If MERGE01=1, type 1's
 *        can merge if their baryonic mass is bigger than the dark matter
 *        mass and type 2's can merge into them. Once the type 1's merge
 *        into a type 0 all its satellites will have the merging clock
 *        into the type 0 reset .
 * */
int join_galaxies_of_progenitors(int halonr, int ngalstart, int *FOF_centralgal)
{
    int ngal, prog, i, j, first_occupied, lenmax, centralgal, mostmassive;

    /* When there is no galaxy in the Halo of FirstProgenitor, the first_occupied
     * pointer is changed to a subhalo which have the maximum mass (This should
     * only happen in the case that the leaf on the firstprogenitor branch occurs
     * as a subhalo, in that case no galaxy would be assigned to it). */
    lenmax = 0;
    first_occupied = Halo[halonr].FirstProgenitor;
    prog = Halo[halonr].FirstProgenitor;
    if(prog >= 0 && HaloAux[prog].NGalaxies == 0) {
        while(prog >= 0) {
	    int currentgal;
	    for(i = 0, currentgal = HaloAux[prog].FirstGalaxy; i < HaloAux[prog].NGalaxies; i++) {
	        if(HaloGal[currentgal].Type == 0 || HaloGal[currentgal].Type == 1) {
		    if(Halo[prog].Len > lenmax) {
		        lenmax = Halo[prog].Len;
			first_occupied = prog;	//define the new first_occupied
		    }
		}
		currentgal = HaloGal[currentgal].NextGalaxy;
	    }
	    prog = Halo[prog].NextProgenitor;
	}
    }

    /* TODO is it true that IF the progenitor halos are sorted by Len
     * (apart maybe from firstprogenitor), then the first_occupied's FirstGalaxy will also be
     * the central galaxy of the new halo?*/
    /* TODO Would it be better ALWAYS to set the most massive galaxy as the first-occupied,
     * even if the FirstProgenitor halo has a galaxy in it? */
   
    // Loop through all the progenitors and get the halo mass and ID of the most massive progenitor.
    lenmax = 0;
    prog = Halo[halonr].FirstProgenitor;
    mostmassive = Halo[halonr].FirstProgenitor;
    while(prog >= 0) {
	if(Halo[prog].Len > lenmax) {
	    lenmax = Halo[prog].Len;
	    mostmassive = prog;
	}
	prog = Halo[prog].NextProgenitor;
    }

    ngal = ngalstart;
    prog = Halo[halonr].FirstProgenitor;
    while(prog >= 0) {
	int currentgal;
	for(i = 0, currentgal = HaloAux[prog].FirstGalaxy; i < HaloAux[prog].NGalaxies; i++) {
	    if(ngal >= MaxGal) {
		AllocValue_MaxGal *= ALLOC_INCREASE_FACTOR;
		MaxGal = AllocValue_MaxGal;
		if(MaxGal<ngal+1) MaxGal=ngal+1;  // How can this ever be true?
		Gal = myrealloc_movable(Gal, sizeof(struct GALAXY) * MaxGal);
	    }
	    // Copy galaxy properties from progenitor, except for those that need initialising.
	    Gal[ngal] = HaloGal[currentgal];
	    Gal[ngal].HaloNr = halonr;
	    Gal[ngal].CoolingRadius = 0.0;
	    Gal[ngal].CoolingGas = 0.0;
	    Gal[ngal].PrimordialAccretionRate = 0.0;
	    Gal[ngal].CoolingRate = 0.0;
	    Gal[ngal].CoolingRate_beforeAGN = 0.0;
	    Gal[ngal].Sfr = 0.0;
#ifdef H2_AND_RINGS  
	  for(j=0;j<RNUM;j++)
	    Gal[ngal].SfrRings[j]=0.0;
#endif
	    Gal[ngal].SfrBulge = 0.0;
	    Gal[ngal].QuasarAccretionRate=0.0;
	    Gal[ngal].RadioAccretionRate=0.0;
#ifdef GALAXYTREE
	    Gal[ngal].FirstProgGal = HaloGal[currentgal].GalTreeIndex;	/* CHECK */
#endif
	    // To fail this check means that we copy in a failed galaxy:
	    mass_checks(ngal,"main.c",__LINE__);
	    /* Update Properties of this galaxy with physical properties of halo.
	     * This deals with the central galaxies of subhalos. */
	    if(Gal[ngal].Type == 0 || Gal[ngal].Type == 1) {
		if(prog == first_occupied) {           // i.e. main progenitor of this halo.
#ifdef HALOPROPERTIES
		    Gal[ngal].HaloM_Mean200 = Halo[halonr].M_Mean200;
		    Gal[ngal].HaloM_Crit200 = Halo[halonr].M_Crit200;
		    Gal[ngal].HaloM_TopHat = Halo[halonr].M_TopHat;
		    Gal[ngal].HaloVelDisp = Halo[halonr].VelDisp;
		    Gal[ngal].HaloVmax = Halo[halonr].Vmax;
#endif
		    Gal[ngal].MostBoundID = Halo[halonr].MostBoundID;
		    for(j = 0; j < 3; j++) {
			Gal[ngal].Pos[j] = Halo[halonr].Pos[j];
			Gal[ngal].Vel[j] = Halo[halonr].Vel[j];
#ifdef HALOPROPERTIES
			Gal[ngal].HaloPos[j] = Halo[halonr].Pos[j];
			Gal[ngal].HaloVel[j] = Halo[halonr].Vel[j];
#endif
		    }
		    Gal[ngal].Len = Halo[halonr].Len;
		    /* TODO This could be place where to set this (FOF-central-)subhalo's
		     * FOFCentralGal property in case that is different from FirstGalaxy. */
		    mass_checks(ngal,"main.c",__LINE__);
		    if(halonr == Halo[halonr].FirstHaloInFOFgroup)
		      {
			*FOF_centralgal=ngal;
			update_centralgal(ngal, halonr);  // Sets Type=0 and updates properties based on halo

#ifdef TRACK_SPLASHBACKS
			//if galaxy was a satellite in a previous snap and now becomes a central, it is a SplashBack
			if(HaloGal[currentgal].Type==1 || HaloGal[currentgal].Type==2)
			  {
			    Gal[ngal].flagSplashBack=1;
			    Gal[ngal].TimeSinceSplashBack=0.;
			  }
#endif
		      }
		    else
			update_type_1(ngal, halonr, prog);  // Sets Type=1 and updates merger timescale

		    Gal[ngal].Vmax = Halo[halonr].Vmax;

		    if(DiskRadiusModel == 1 || DiskRadiusModel == 2)
		      {
			Gal[ngal].ColdGasRadius = get_initial_disk_radius(halonr, ngal);
			Gal[ngal].DiskRadius = Gal[ngal].ColdGasRadius;
		      }
		    else
		      {
			// After changing Vmax disk sizes may need resetting
			Gal[ngal].ColdGasRadius=get_gas_disk_radius(ngal);
			Gal[ngal].DiskRadius = get_stellar_disk_radius(ngal);
		      }

		}
		else // If did not come from main progenitor then set to Type 2    
		    update_type_2(ngal, halonr, prog, mostmassive);
	    }
	    mass_checks(ngal,"main.c",__LINE__);
	    if(Gal[ngal].Type < 0 || Gal[ngal].Type > 2) terminate("Unknown galaxy type\n");
	    ngal++;
	    currentgal = HaloGal[currentgal].NextGalaxy;
	}
	prog = Halo[prog].NextProgenitor;
    }

    /* If there are no progenitors with galaxies, a new galaxy is created.
     * However, if it's a subhalo, no galaxy is placed, since it would stay
     * at zero luminosity. */
    if(ngal == 0) {
	if(Halo[halonr].FirstHaloInFOFgroup == halonr) {
	    init_galaxy(ngal, halonr);
	    *FOF_centralgal=ngal;
	    ngal++;
	}
    }

    // Find and set the central galaxies for each galaxy in this halo
    if (ngal > ngalstart) {
        // Find the central galaxy of this halo
        for (i = ngalstart, centralgal = -1; i < ngal; i++) {
	    if (Gal[i].Type == 0 || Gal[i].Type == 1) {
	        if(centralgal != -1) terminate("Subhalo hosts more than one Type 0/1\n");
		centralgal = i;
	    }
	}
	//if(centralgal == -1)
	//    printf("Warning: for halo %d, No central galaxy found\n",halonr);
	// Set galaxy centres.
	for (i = ngalstart; i < ngal; i++) {
	    Gal[i].CentralGal=centralgal;
	    Gal[i].MergerCentralGal=centralgal;  // Until good reason to change.
	    /* Need to remember the following only for GUO10 as galaxy IDs get confused
	     * by copying everything back into HaloGal.  This doesn't seem to happen
	     * for any other version of the model. */
	    if (centralgal != -1) {
		for(j = 0; j < 3; j++) Gal[i].MergCentralPos[j] = Gal[centralgal].Pos[j];  
	    }
	}
	for (i = ngalstart; i<ngal; i++) mass_checks(i,"main.c",__LINE__);
    }
    
    report_memory_usage(&HighMark, "join_galaxies");

    return ngal;
}


/**@brief evolve_galaxies() deals with most of the SA recipes. This is
 *       where most of the physical recipes are called, including:
 *       infall_recipe() (gets the fraction of primordial gas that infalled),
 *       add_infall_to_hot() (adds the infalled gas to the hot phase),
 *       reincorporate_gas() (reincorporates gas ejected by SN),
 *       cooling_recipe() (gets the amount of gas that cooled - takes into
 *       account AGN feedback), cool_gas_onto_galaxy() (adds the gas that
 *       cooled into the cold phase, starformation_and_feedback() (normal
 *       SF and SN feedback), deal_with_galaxy_merger() (adds components,
 *       grows black hole and deals with SF burst), disruption() (total and
 *       instantaneous disruption of type 2 satellites) and dust() (if galaxy
 *       is in an output time, dust extinction is computed).
 *
 *       All these calculations are done in time steps of 1/STEPS the time
 *       between each snapshot (STEPS=20).
 *       
 *       * TODO - Note that because galaxies can skip a snapshot, they should
 * be evolved forward in time to the same snap as all the other galaxies,
 * treating them as if they were isolated type 0 galaxies.
 * Currently this is not done.
 * In versions of the code up to Aug 2011 each galaxy was given its own
 * timestep; however this is also incorrect as galaxies then have different
 * times during periods of interaction.  This showed up when using 
 * star-formation histories.
 */
 
/* Note: halonr is here the FOF-background subhalo (i.e. main halo) */
void evolve_galaxies(int halonr, int ngal, int treenr) {
    int p, q, nstep, centralgal, FOF_centralgal, merger_centralgal, currenthalo, prevgal, start, i;
    double infallingGas, deltaT;
    double time, previoustime, newtime;
    double AGNaccreted, t_Edd;
#ifdef STAR_FORMATION_HISTORY
    double age_in_years;
#endif
#ifdef HT09_DISRUPTION
    double CentralRadius, CentralMass, SatelliteRadius, SatelliteMass;
#endif

    // Eddington time in code units
    // code units are UnitTime_in_s/Hubble_h
    t_Edd=1.42e16*Hubble_h/UnitTime_in_s;

    //previoustime = NumToTime(Gal[0].SnapNum);
    previoustime = NumToTime(Halo[halonr].SnapNum-1);
    newtime = NumToTime(Halo[halonr].SnapNum);

    /* Time between snapshots */
    deltaT = previoustime - newtime;

  /* If Gal[0].Centralgal and Gal[0].FOFCentralGal are always the same then we don't need the latter.
     * That will be true whenever (as now) the FirstHaloInFOFGroup always holds the Type 0.
     * Have put in a trap to check that, even though the code should work fine as currently written. */
    FOF_centralgal = Gal[0].FOFCentralGal;
    if (Gal[FOF_centralgal].Type != 0) terminate("(Gal[FOF_centralgal].Type != 0\n");
    if (Gal[0].CentralGal != FOF_centralgal) terminate("(Gal[0].CentralGal != FOF_centralgal\n");

#ifdef DEBUG
  debug_galaxy(FOF_centralgal, "main.c",__LINE__);
#endif

  // This proves that FOF_centralgal is the galaxy at the centre of the FOF group.
    if(Gal[FOF_centralgal].Type != 0 || Gal[FOF_centralgal].HaloNr != halonr)
	terminate("Something wrong here ..... \n");
    
#ifdef SAVE_FOFHALO
   for (p =0;p<ngal;p++)
       Gal[p].FoFHaloNr=halonr;
#endif

    for (p =0;p<ngal;p++) mass_checks(p,"main.c",__LINE__);

  /* Update all galaxies to same star-formation history time-bins.
   * Needed in case some galaxy has skipped a snapshot. */
#ifdef STAR_FORMATION_HISTORY
    age_in_years=(Age[0]-previoustime)*UnitTime_in_years/Hubble_h; //ROB: age_in_years is in units of "real years"!
    nstep=0;
    for (p=0; p<ngal; p++) sfh_update_bins(p,Halo[halonr].SnapNum-1,nstep,age_in_years);
#endif

    /* Handle the transfer of mass between satellites and central galaxies */
    deal_with_satellites(ngal);

    /* Delete inconsequential galaxies */
    for (p =0;p<ngal;p++)
	if (Gal[p].Type ==2 && Gal[p].ColdGas+Gal[p].DiskMass+Gal[p].BulgeMass <1.e-8)
	    Gal[p].Type = 3;
	else
	    mass_checks(p,"main.c",__LINE__);
   
    /* Calculate how much hot gas needs to be accreted to give the correct baryon fraction
     * in the main halo. This is the universal fraction, less any reduction due to reionization. */
    infallingGas = infall_recipe(ngal);
    Gal[FOF_centralgal].PrimordialAccretionRate=infallingGas/deltaT;
    
    /* All the physics are computed in a number of intervals between snapshots
     * equal to STEPS */
    for (nstep = 0; nstep < STEPS; nstep++) {
	/* time to present of the current step */
	time = previoustime - (nstep + 0.5) * (deltaT / STEPS);

	/* Update all galaxies to the star-formation history time-bins of current step*/
#ifdef STAR_FORMATION_HISTORY
	age_in_years=(Age[0]-time)*UnitTime_in_years/Hubble_h;
	for (p=0; p<ngal; p++)
	    sfh_update_bins(p,Halo[halonr].SnapNum-1,nstep,age_in_years);
#endif

#ifdef DEBUG
	debug_galaxy(FOF_centralgal, "main.c",__LINE__);
#endif

      /* Infall onto central galaxy only, if required to make up a baryon deficit */
#ifndef GUO10
#ifndef GUO13
#ifndef EXCESS_MASS
	// If we have NOT defined one of the flags, then only allow positive infall
	if (infallingGas > 0.)
#endif
#endif
#endif
	add_infall_to_hot(infallingGas / STEPS);

	mass_checks(FOF_centralgal,"main.c",__LINE__);

	for (p = 0; p < ngal; p++) {
	    /* don't treat galaxies that have already merged */
	    if(Gal[p].Type == 3) continue;
	    mass_checks(p,"main.c",__LINE__);
	    
	    if (Gal[p].Type == 0 || Gal[p].Type == 1)
	      {
		if((ReIncorporationModel == 2 && Gal[p].Type==0) || ReIncorporationModel < 2)
		  if(Gal[p].EjectedMass>0.)
		    reincorporate_gas(p, deltaT / STEPS);
		// determine cooling gas given halo properties and add it to the cold phase
		mass_checks(p,"main.c",__LINE__);
		compute_cooling(p, deltaT / STEPS, ngal);
	      }
	}

	/* This must be separated as now satellite AGN can heat central
	 * galaxies.  Therefore the AGN from all satellites must be computed, in
	 * a loop inside this function, before gas is cooled into central
	 * galaxies (only suppress cooling, the gas is not actually heated) */
	if(AGNRadioModeModel <4)
	    do_AGN_heating(deltaT / STEPS, ngal,FOF_centralgal);

	for (p = 0; p < ngal; p++) {
	    cool_gas_onto_galaxy(p, deltaT / STEPS);
	    mass_checks(p,"main.c",__LINE__);
#ifdef H2_AND_RINGS
	    if(Gal[p].ColdGas>0.)
	      gas_inflow(p, deltaT / STEPS);
#endif
	    mass_checks(p,"main.c",__LINE__);
	    starformation(p, FOF_centralgal, time, deltaT / STEPS, nstep);
	    mass_checks(p,"main.c",__LINE__);
#ifdef DEBUG
	  debug_galaxy(p, "main.c",__LINE__);
#endif
	}

	// Merger detection
	for(p = 0; p < ngal; p++) {
#ifdef MERGE01
	    if(Gal[p].Type == 2 || (Gal[p].Type == 1 && Gal[p].MergeOn == 1)) { /* satellite galaxy */
#else
	    if(Gal[p].Type == 2) {
#endif
#ifndef HT09_DISRUPTION
		Gal[p].MergTime -= deltaT / STEPS;
		if(Gal[p].MergTime < 0.0) {
#else
		Gal[p].MergRadius -= min(get_deltar(p, deltaT/STEPS, FOF_centralgal),Gal[p].MergRadius);
		disruption_code (p, time, FOF_centralgal);
		if(Gal[p].MergRadius<(Gal[FOF_centralgal].DiskRadius+Gal[FOF_centralgal].BulgeSize +
				      Gal[p].DiskRadius+Gal[p].BulgeSize) || Gal[p].BulgeMass+Gal[p].DiskMass == 0 ) {
#endif
		    NumMergers++;

#ifdef MERGE01
		    if (Gal[p].Type == 1)
			// Reset everything orbiting this galaxy to have a merger centre of the FOF_centralgal
 			for(q = 0; q < ngal; q++) 
			    if (Gal[q].MergerCentralGal == p) 
				Gal[q].MergerCentralGal = FOF_centralgal;
#endif
		    merger_centralgal = Gal[p].MergerCentralGal;

		    mass_checks(p,"main.c",__LINE__);
		    mass_checks(merger_centralgal,"main.c",__LINE__);
		    mass_checks(FOF_centralgal,"main.c",__LINE__);

		    deal_with_galaxy_merger(p, time, deltaT, nstep);

		    mass_checks(p,"main.c",__LINE__);
		    mass_checks(merger_centralgal,"main.c",__LINE__);
		    mass_checks(FOF_centralgal,"main.c",__LINE__);
		}
	    }
#ifdef DEBUG
	  debug_galaxy(p, "main.c",__LINE__);
#endif
	}//loop on all galaxies to detect mergers

	/* Cool gas onto AGN */
	if (BlackHoleGrowth == 1)
	    for (p = 0; p < ngal; p++) {
		AGNaccreted=min(Gal[p].BlackHoleGas, Gal[p].BlackHoleMass*BlackHoleAccretionRate*deltaT/(STEPS*t_Edd));
		if (AGNaccreted > 0.) {
		    Gal[p].BlackHoleMass += AGNaccreted;
		    Gal[p].BlackHoleGas -= AGNaccreted;
		    // Instantaneous accretion rate.  This will get overwritten on each mini-step but that's OK
		    //Gal[p].QuasarAccretionRate = AGNaccreted*STEPS/deltaT;
		}
	    }

#ifdef DETAILED_METALS_AND_MASS_RETURN
	//DELAYED ENRICHMENT AND MASS RETURN + FEEDBACK: No fixed yield or recycling fraction anymore. FB synced with enrichment
	for (p = 0; p < ngal; p++)
	  update_yields_and_return_mass(p, FOF_centralgal, deltaT/STEPS, nstep);
#endif

#ifdef ALL_SKY_LIGHTCONE
	int nr, istep, ix, iy, iz;
	istep = Halo[halonr].SnapNum*STEPS + nstep;
	Gal[p].SnapNum = Halo[halonr].SnapNum;
	for (p = 0; p < ngal; p++)
	  for (nr = 0; nr < NCONES; nr++)
	    for (ix = 0; ix < NREPLICA; ix++)
	      for (iy = 0; iy < NREPLICA; iy++)
		for (iz = 0; iz < NREPLICA; iz++)
		  inside_lightcone(p, istep, nr, ix, iy, iz);
#endif


    }/* end move forward in interval STEPS */

    /* Location and disruption of type 2 galaxies. Type 1 galaxies are not
     * disrupted since usually bayonic component is more compact than dark
     * matter.*/
    for(p = 0; p < ngal; p++) {
	if(Gal[p].Type == 2) {
#ifndef UPDATETYPETWO
	    int jj;
	    float tmppos;
	    for(jj = 0; jj < 3; jj++) {
		tmppos = wrap(Gal[p].DistanceToCentralGal[jj],BoxSize);
		//some fudge factor to get convergent results
		//from when using mostboundpart (Bruno Henriques 07.2015)
#if defined(GUO10) || defined(GUO13) || defined(HENRIQUES13)
	      tmppos *=  sqrt(Gal[p].MergTime/Gal[p].OriMergTime);
#else
	      //some fudge factor to get convergent results
	      //from when using mostboundpart (Bruno Henriques 07.2015)
	      tmppos *=  2.*sqrt(Gal[p].MergTime/Gal[p].OriMergTime);
#endif
	      Gal[p].Pos[jj] = Gal[p].MergCentralPos[jj] + tmppos;

	      if(Gal[p].Pos[jj] < 0)
		Gal[p].Pos[jj] = BoxSize + Gal[p].Pos[jj];
	      if(Gal[p].Pos[jj] > BoxSize)
		Gal[p].Pos[jj] = Gal[p].Pos[jj] - BoxSize;
	    }
#endif
#ifdef DISRUPTION
	    disrupt(p, Gal[p].CentralGal);
#endif

        }
    }

    for (p =0;p<ngal;p++)
      if(Gal[p].DiskMass+Gal[p].BulgeMass>0.)
	{
	  int do_ColdGas=0, do_DiskMass=1, do_BulgeMass=1;
	  Gal[p].StellarHalfMassRadius=half_mass_radius(p,do_ColdGas,do_DiskMass,do_BulgeMass);
	}

    for (p =0;p<ngal;p++) mass_checks(p,"main.c",__LINE__);
  


#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef  POST_PROCESS_MAGS
    int n;
    /* If this is an output snapshot apply the dust model to each galaxy */
    for(n = 0; n < NOUT; n++) {
	if(Halo[halonr].SnapNum == ListOutputSnaps[n]) {
	    for(p = 0; p < ngal; p++) dust_model(p, n, halonr);
	    break;
        }
    }
#endif  //POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

    /* now save the galaxies of all the progenitors (and free the associated storage) */
    int prog = Halo[halonr].FirstProgenitor;

    while(prog >= 0) {
	int currentgal;
	for(i = 0, currentgal = HaloAux[prog].FirstGalaxy; i < HaloAux[prog].NGalaxies; i++) {
	    int nextgal = HaloGal[currentgal].NextGalaxy;
	    /* this will write this galaxy to an output file and free the storage associate with it */
	    output_galaxy(treenr, HaloGal[currentgal].HeapIndex);
	    currentgal = nextgal;
        }
	prog = Halo[prog].NextProgenitor;
    }

    for(p = 0, prevgal = -1, currenthalo = -1, centralgal = -1, start = NGalTree; p < ngal; p++) {
	if(Gal[p].HaloNr != currenthalo) {
	    currenthalo = Gal[p].HaloNr;
	    HaloAux[currenthalo].FirstGalaxy = -1;
	    HaloAux[currenthalo].NGalaxies = 0;
	}
	
	mass_checks(p,"main.c",__LINE__);

	if(Gal[p].Type != 3) {
	    if(NHaloGal >= MaxHaloGal) {
		int oldmax = MaxHaloGal;
		AllocValue_MaxHaloGal *= ALLOC_INCREASE_FACTOR;
		MaxHaloGal = AllocValue_MaxHaloGal;
		if(MaxHaloGal<NHaloGal+1) MaxHaloGal=NHaloGal+1;
		HaloGal = myrealloc_movable(HaloGal, sizeof(struct GALAXY) * MaxHaloGal);
		HaloGalHeap = myrealloc_movable(HaloGalHeap, sizeof(int) * MaxHaloGal);
		for(i = oldmax; i < MaxHaloGal; i++) HaloGalHeap[i] = i;
	    }
	    Gal[p].SnapNum = Halo[currenthalo].SnapNum;
#ifndef GUO10
#ifdef UPDATETYPETWO
	    update_type_two_coordinate_and_velocity(treenr, p, Gal[0].CentralGal);
#endif
#endif
	    /* when galaxies are outputed, the slot is filled with the
	     * last galaxy in the heap. New galaxies always take the last spot */
	    int nextgal = HaloGalHeap[NHaloGal];
	    HaloGal[nextgal] = Gal[p];
	    HaloGal[nextgal].HeapIndex = NHaloGal;
	    if(HaloAux[currenthalo].FirstGalaxy < 0) HaloAux[currenthalo].FirstGalaxy = nextgal;
	    if(prevgal >= 0) HaloGal[prevgal].NextGalaxy = nextgal;
	    prevgal = nextgal;
	    HaloAux[currenthalo].NGalaxies++;
	    NHaloGal++;
#ifdef GALAXYTREE
	    if(NGalTree >= MaxGalTree) {
		AllocValue_MaxGalTree *= ALLOC_INCREASE_FACTOR;
		MaxGalTree = AllocValue_MaxGalTree;
		if(MaxGalTree<NGalTree+1) MaxGalTree=NGalTree+1;
		GalTree = myrealloc_movable(GalTree, sizeof(struct galaxy_tree_data) * MaxGalTree);
	    }
	    HaloGal[nextgal].GalTreeIndex = NGalTree;

	    memset(&GalTree[NGalTree], 0, sizeof(struct galaxy_tree_data));
	    GalTree[NGalTree].HaloGalIndex = nextgal;
	    GalTree[NGalTree].SnapNum = Halo[currenthalo].SnapNum;
	    GalTree[NGalTree].NextProgGal = -1;
	    GalTree[NGalTree].DescendantGal = -1;

	    GalTree[NGalTree].FirstProgGal = Gal[p].FirstProgGal;
	    if(Gal[p].Type == 0)
		centralgal = NGalTree;
	    NGalTree++;
#endif
	}
    }

#ifdef GALAXYTREE
  for(p = start; p < NGalTree; p++)
    {
      if(centralgal < 0)
    	  terminate("centralgal < 0");
      GalTree[p].FOFCentralGal = centralgal;
    }
#endif

  report_memory_usage(&HighMark, "evolve_galaxies");
}

  /**
   * @brief Check whether makefile options are compatible.
   */ 
void check_options()
  {
#ifdef OUTPUT_OBS_MAGS
#ifndef COMPUTE_OBS_MAGS
    printf("> Error : option OUTPUT_OBS MAGS requires option COMPUTE_OBS_MAGS \n");
    exit(32);
#endif
#endif

/*#ifdef OVERWRITE_OUTPUT
#ifdef PARALLEL
    terminate("PARALLEL cannot be run with OVERWRITE_OUTPUT\n");
#endif
#endif*/
  
#ifdef OUTPUT_MOMAF_INPUTS
#ifndef POST_PROCESS_MAGS
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("> Warning : option OUTPUT_MOMAF_INPUTS should really only be used with POST_PROCESS_MAGS \n");
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
#endif
#endif


#ifdef OUTPUT_MOMAF_INPUTS
#ifndef COMPUTE_OBS_MAGS 
    printf("> Error : option OUTPUT_MOMAF_INPUTS requires option COMPUTE_OBS_MAGS \n");
    exit(32);
#endif
#endif
  
#ifdef KITZBICHLER
#ifndef OUTPUT_MOMAF_INPUTS
    printf("> Error : option KITZBICHLER requires option OUTPUT_MOMAF_INPUTS \n");
    exit(32);
#endif
#ifndef POST_PROCESS_MAGS
    printf("> Error : option KITZBICHLER requires option POST_PROCESS_MAGS \n");
    exit(32);
#endif
#endif

#ifdef OUTPUT_L_CONE_INPUTS
#ifndef OUTPUT_OBS_MAGS
    printf("> Error : option OUTPUT_L_CONE_INPUTS requires option OUTPUT_OBS_MAGS \n");
    exit(32);
#endif
#endif

#ifdef OUTPUT_L_CONE_INPUTS
#ifndef GALAXYTREE
    printf("> Error : option OUTPUT_L_CONE_INPUTS requires option GALAXYTREE \n");
    exit(32);
#endif
#endif

#ifndef DISRUPTION
#ifdef ICL
  terminate("> Warning : DISRUPTION off then ICL makes no sense \n");
#endif
#endif

#ifndef LOADIDS
#ifdef GALAXYTREE
  terminate("> Warning : GALAXYTREE requires LOADIDS \n");
#endif
#endif

#ifndef STAR_FORMATION_HISTORY
#ifdef POST_PROCESS_MAGS
  terminate("> Warning : POST_PROCESS_MAGS  requires STAR_FORMATION_HISTORY \n");
#endif
#endif

#ifdef MCMC
#ifndef LOADIDS
  terminate("> Warning : MCMC  requires LOADIDS \n");
#endif

#ifdef HALOMODEL
#ifdef MR_PLUS_MRII
  terminate("> Warning : HALOMODEL doesn't work yet with MR_PLUS_MRII\n");
#endif
#endif
#endif

#ifdef DEBUG_PRINT
#ifdef DEBUG_READ_AND_CHECK
  terminate("DEBUG_PRINT & DEBUG_READ_AND_CHECK cannot run together");
#endif
#endif

#ifdef PHOTTABLES_PRECOMPUTED
#ifdef SPEC_PHOTABLES_ON_THE_FLY
  terminate("> Warning : PHOTTABLES_PRECOMPUTED cannot run with SPEC_PHOTABLES_ON_THE_FLY\n");
#endif
#endif

#ifdef LIGHT_OUTPUT
#ifdef POST_PROCESS_MAGS
  terminate("> Warning : LIGHT_OUTPUT cannot run with POST_PROCESS_MAGS \n");
#endif
#ifdef OUTPUT_MOMAF_INPUTS
  terminate("> Warning : LIGHT_OUTPUT cannot run with OUTPUT_MOMAF_INPUTS \n");
#endif
#endif //LIGHT_OUTPUT

/* Description of the code to appear in the first page of the documentation. */

}


void output_galaxy(int treenr, int heap_index)
{
  int gal_index = HaloGalHeap[heap_index];

  if(heap_index >= NHaloGal)
    terminate("heap_index >= NHaloGal");

  if(HaloGal[gal_index].HeapIndex != heap_index)	// consistency check
    terminate("HeapIndex != heap_index");

#ifdef GUO10
#ifdef UPDATETYPETWO
  update_type_two_coordinate_and_velocity(treenr, gal_index, HaloGal[0].CentralGal);
#endif
#endif

#ifdef GALAXYTREE
  GalTree[HaloGal[gal_index].GalTreeIndex].IndexStored = IndexStored++;
  /* Note: to use a minimum output mass for the GALAXYTREE option then one would
   * have to go through the code and check that all the pointers are handled
   * correctly.  Also, the code currently attempts to read back in all the
   * galaxies in save_galaxy_tree_finalize() and so crashes. */
  //if((HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)>=MinGalOutputMass)
    save_galaxy_tree_append(gal_index);
#else
#ifdef MCMC
  /* if MCMC galaxies are saved into a structure to then be compared
   * with observations. The output will come in form of a file with
   * parameters, not galaxy properties */
  if((HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)>=MinGalOutputMass)
    save_galaxy_for_mcmc(gal_index);
#else //Normal snap output
  int n;
  for(n = 0; n < NOUT; n++)
    if(ListOutputSnaps[n] == HaloGal[gal_index].SnapNum)
      if((HaloGal[gal_index].DiskMass+HaloGal[gal_index].BulgeMass)>=MinGalOutputMass)
	save_galaxy_append(treenr, gal_index, n);
#endif
#endif

  /* fill the gap in the heap with the galaxy in the last occupied slot */

  int last = NHaloGal - 1;
  int last_gal_index = HaloGalHeap[last];
  HaloGalHeap[last] = gal_index;
  HaloGalHeap[heap_index] = last_gal_index;

  /* make sure that the back-pointer of the last galaxy is updated */
  HaloGal[last_gal_index].HeapIndex = heap_index;

  NHaloGal--;
}
