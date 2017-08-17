#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"

#include "mcmc_vars.h"
#include "mcmc_proto.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define RNMX (1.0-EPS)
#define EPS 1.2e-7


/**@file mcmc.c
 * @brief This is the main mcmc file, it reads in observational tests,
 *        starts a chain, proposes new sets of parameters, calls the
 *        SAM, gets the likelihood for galaxy properties obtained with
 *        the current parameters, compares the likelihood with the
 *        previous run and decides on whether or not to accept the
 *        proposed parameters.
 *
 * To run you need to chose input_mcmc_wmp7.par and select the
 * desired sample of halos. Also choose the correct desired_output_snap
 * according to the redshift at which you want to run the constraints.
 * Turn on MCMC on the Makefile.
 *
 * Senna is the controlling routine. First is calls read_observations,
 * to read in all the observational data sets into the struct MCMC_Obs.
 * Then, an initial SAM is run, to obtain a first likelihood value.
 * This will be the starting point of the chain, given by the initial
 * values of p[] (this should be some non 0 likelihood region to avoid
 * a long burn in). The likelihood for each set of parameters is computed
 * in mcmc_likelihood.c (get_likelihood()) at the end of each SAM run.
 *
 * After this first step, a chain is started with mcmc() being called
 * for each step to control all the processes. Namely, propose_new_parameters(),
 * SAM() and then at each step the likelihood for the given set of
 * parameters compared with previous and accepted or not. Two options
 * available according to the value of MCMCMode in input.par. If MCMCMode
 * =0, normal MCMC is done. If MCMCMode = 1, the new set of parameters
 * is only accepted if the likelihood is higher than for the previous.
 * This makes the chain going up in likelihood very quickly and its
 * useful to find a region of non 0 likelihood when the parameter is
 * still unknown.
 *
 * The name of the parameters sampled can be seen in the beginning of
 * function SAM() where the values of the internal variables of the code
 * containing the parameters are changed according to the new set of
 * parameters proposed. At each step p[] will contain the previously
 * accepted set of parameters and prop[] will contain the newly proposed
 * set of parameters (given by propose_new_parameters()). If desired,
 * less parameters then the default number, given by MCMCNpar can be sampled.
 * This can be done in propose_new_parameters() by never assigning a new
 * value to a given parameter. */
#ifdef MCMC
void Senna()
{
  int ii, jj, snap, IndividualAcceptRate = 0, TotAcceptRate = 0;
  FILE *fmcmc;
  char buf[1000];
  char DirNum;
  time_t local_initial, final;
  double lhood2, qratio, AcceptanceProbability, ran;
  int AcceptanceLogic;

  if(ThisTask==0)
    {
      printf("\n\n\n");
      printf("**********************************************************\n");
      printf("*                                                        *\n");
      printf("*                   Starting Senna                       *\n");
      printf("*                                                        *\n");
      printf("*             MCMC parameter estimation                  *\n");
      printf("*  Applied to a Semi-Analytic Model of Galaxy Formation  *\n");
      printf("*                                                        *\n");
      printf("**********************************************************\n\n");
    }

  MCMCseed = -((ThisTask+FirstChainNumber) * 100 + 25);
  MaxLikelihood=0.;

#ifdef HALOMODEL //to compute correlation function for MCMC
  if (Sample_Cosmological_Parameters==0)
    initialize_halomodel();
  printf("halo model initialized\n");
#endif

  //This file will have the values for the parameters
  //It will be the output from the mcmc sampling used in the analysis
  //it is also where the parameters are read from, using previous chains
  getcwd(buf, sizeof(buf));
  DirNum=buf[strlen(buf)-1];
  sprintf(buf, "%s/senna_g%c_%d.txt", OutputDir, DirNum, ThisTask+FirstChainNumber);
  //open to read and write. if file doesn't exists open just to write and read from backup
  if((fmcmc = fopen(buf, "r+")) == NULL)
    if((fmcmc = fopen(buf, "w")) == NULL)
      {
	char sbuf[1000];
	sprintf(sbuf, "can't open file `%s'\n", buf);
	terminate(sbuf);
      }

  //read (from previous output) and initialize mcmc parameters, also lhood1
  initialize_mcmc_par_and_lhood (fmcmc);

  //Read observational constraints (SMF, LFs, BHBM relation, etc)
  read_observations();

  //open files that will contain the comparison with observations at each step (used to compute the likelihoods)
  open_files_with_comparison_to_observations();

  for(ii=1;ii<ChainLength+1;ii++)
    {
      CurrentMCMCStep=ii;
      GlobalMCMCStep+=1;
      time(&local_initial);
      printf("\n\n\nMCMC %d STARTED on Task %d\n\n\n", ii, ThisTask+FirstChainNumber);

      IndividualAcceptRate = 0;

      //get a new set of parameters and return qratio - the prior
      qratio = propose_new_parameters();

      //runs the SAM with the new parameters and gives the likelyhood for them

      lhood2=SAM(MCMCSampleFile);

      printf("LIKELY1=%0.5e\nLIKELY2=%0.5e\n",lhood1,lhood2);

      if(isnan(lhood1))	lhood1=0.0;
      if(isnan(lhood2))	lhood2=0.0;

      /* By default qratio = 1, meaning we assume a flat prior. Therefore, the acceptance
       * probability is just given by the ratio of the likelihoods from two steps. */
      if(1.0 > (qratio * (lhood2 / lhood1)))
	AcceptanceProbability = qratio * (lhood2 / lhood1);
      else
	AcceptanceProbability = 1.0;

      //generate random number to assess probability
      ran = ran3(&MCMCseed);

      //accepts the proposed parameters with probability=acceptance rate
      if(MCMCMode == 0)
	AcceptanceLogic = (ran < AcceptanceProbability);
      //only accepts new parameters if likelihood increases
      if(MCMCMode == 1)
	AcceptanceLogic = (lhood2 > lhood1);

      //AcceptanceLogic=1;
      //if new set of parameters is accepted change current set of parameter and lhood
      if(AcceptanceLogic)
	{
	  for(jj=0;jj<MCMCNpar;++jj)
	    for(snap=0;snap<NOUT;snap++)
	      MCMC_PAR[jj].Value[snap] = MCMC_PAR[jj].PropValue[snap];

	  lhood1=lhood2;
	  IndividualAcceptRate=1;

	  //if current step is the maximum likelihood so far, copy the comparison with
	  //observations into a special directory in output
	  if(lhood1>MaxLikelihood)
	    {
	      MaxLikelihood=lhood1;
	      create_bestfit_files();
	    }
	}

      //print the values to file and to screen (to screen only if IndividualAcceptRate=1)
      print_parameters(IndividualAcceptRate, fmcmc);

      TotAcceptRate = TotAcceptRate + IndividualAcceptRate;
      printf("\nCurrent acceptance rate of this chain=%0.2f%%\n", ((float) TotAcceptRate / ii) * 100);

      time(&final);
      printf("Task %d chain %d (%d) took %lds\n", ThisTask+FirstChainNumber, CurrentMCMCStep, GlobalMCMCStep, final - local_initial);
      printf("Global Time Elapsed %lds, %f hours\n", final - global_starting_time, (final - global_starting_time)/3600.);

#ifdef PARALLEL
      if(ThisTask==0)
	if((final - global_starting_time)/3600. > MachineTimeOut)
	  {
	    sprintf(buf, "%s %s %s%d.bash", JobSubmitCommand, JobSubmitPipe, JobSubmitFile, FirstChainNumber);
	    printf("resubmit command: %s\n",buf);
	    system(buf);
	    fflush(stdout);
	    terminate("\n\nMachine TimeOut reached, restarting\n\n");
	  }

#endif
    }// ChainLength (MCMC steps)

  fclose(fmcmc);

  /* For each set of parameters writes the semi-analytic
   * predictions used for the likelihood calculation*/
  //write_comparison_to_obs();
  close_files_with_comparison_to_observations();

  myfree(MCMC_Obs);

#ifdef HALOMODEL //to compute correlation function for MCMC
  if (Sample_Cosmological_Parameters==0) {
      gsl_spline_free(FofSpline);
      gsl_interp_accel_free(FofAcc);
      gsl_spline_free(SigmaSpline);
      gsl_interp_accel_free(SigmaAcc);
      gsl_spline_free(ellipSpline);
      gsl_interp_accel_free(ellipAcc);
      gsl_spline_free(PowSpline);
  } //if
#endif

  printf("\nFinal acceptance rate of this chain=%f%%\n", ((float) TotAcceptRate / ChainLength) * 100);
  printf("\n\nMCMC OVER\n\n");

}


////////
//MCMC//
////////





/*@brief Function to print parameters whenever
 *       a step is accepted in the MCMC*/
void print_parameters (int AcceptanceLogic, FILE *fmcmc)
{
  int i, snap, chainweight=1;

  fprintf(fmcmc,"%d %0.8g",chainweight, -(log10(lhood1)));

  //print parameters into output file
  for(i=0;i<MCMCNpar;++i)
    {
      if(MCMC_PAR[i].Sampling_Switch==1)
	{
	  if(strcmp(MCMC_PAR[i].Type,"Physical")==0)
	    {
	      if(Sample_Physical_Parameters==1)
		{
		  if(Time_Dependant_PhysPar==1)
		    for(snap=0;snap<NOUT;snap++)
		      fprintf(fmcmc," %0.6f", log10(MCMC_PAR[i].Value[snap]));
		  else
		    fprintf(fmcmc," %0.6f", log10(MCMC_PAR[i].Value[0]));
		}
	    }
	  else if (strcmp(MCMC_PAR[i].Type,"Cosmological")==0)
	    {
	      if(Sample_Cosmological_Parameters==1)
		fprintf(fmcmc, " %0.6f", log10(MCMC_PAR[i].Value[0]));
	    }
	}
    }
  fprintf(fmcmc,"\n");

  if(AcceptanceLogic==2)
    {

      fprintf(fmcmc,"%d %0.8g",chainweight, lhood1);

      //print parameters into output file
      for(i=0;i<MCMCNpar;++i)
	{
	  if(MCMC_PAR[i].Sampling_Switch==1)
	    {
	      if(strcmp(MCMC_PAR[i].Type,"Physical")==0)
		{
		  if(Sample_Physical_Parameters==1)
		    {
		      if(Time_Dependant_PhysPar==1)
			for(snap=0;snap<NOUT;snap++)
			  fprintf(fmcmc," %0.2g", MCMC_PAR[i].Value[snap]);
		      else
			fprintf(fmcmc," %0.2g", MCMC_PAR[i].Value[0]);
		    }
		}
	      else if (strcmp(MCMC_PAR[i].Type,"Cosmological")==0)
		{
		  if(Sample_Cosmological_Parameters==1)
		    fprintf(fmcmc, " %0.2g", MCMC_PAR[i].Value[0]);
		}
	    }
	}
      fprintf(fmcmc,"\n");
    }

  fflush(fmcmc);
  fflush(stdout);


  //print to screen
  if(AcceptanceLogic==1)
    {

      printf("\n******************************************************\n");
      printf("Accepted!!!\n");
      for(i=0;i<MCMCNpar;++i)
	if(MCMC_PAR[i].Sampling_Switch==1)
	  {
	    if (strcmp(MCMC_PAR[i].Type,"Physical")==0)
	      {
		if(Sample_Physical_Parameters==1)
		  {
		    if(Time_Dependant_PhysPar==1)
		      for(snap=0;snap<NOUT;snap++)
			printf("%0.2g ",MCMC_PAR[i].PropValue[snap]);
		    else
		      printf("%0.2g ",MCMC_PAR[i].PropValue[0]);
		  }
	      }
	    else if (strcmp(MCMC_PAR[i].Type,"Cosmological")==0)
	      {
		if(Sample_Cosmological_Parameters==1)
		  printf("%0.2g ",MCMC_PAR[i].PropValue[0]);
	      }
	  }
      printf("\n");


      //print log of parameters to screen
      printf("%d %0.8g ",chainweight, -(log10(lhood1)));
      for(i=0;i<MCMCNpar;++i)
	if(MCMC_PAR[i].Sampling_Switch==1)
	  {
	    if (strcmp(MCMC_PAR[i].Type,"Physical")==0)
	      {
		if(Sample_Physical_Parameters==1)
		  {
		    if(Time_Dependant_PhysPar==1)
		      for(snap=0;snap<NOUT;snap++)
			printf("%0.6f ",log10(MCMC_PAR[i].PropValue[snap]));
		    else
		      printf("%0.6f ",log10(MCMC_PAR[i].PropValue[0]));
		  }
	      }
	    else if (strcmp(MCMC_PAR[i].Type,"Cosmological")==0)
	      {
		if(Sample_Cosmological_Parameters==1)
		  printf("%0.4f ",log10(MCMC_PAR[i].PropValue[0]));
	      }
	  }


      printf("\n******************************************************\n\n\n\n");
    } //end if(AcceptanceLogic==1) - print ot screen

}


void create_bestfit_files()
{
  char buf[1000], sbuf[1000];
  int ii, ObsNr, snap;
  FILE *fa;

  sprintf(buf, "%s/bestfit/task%d_bestfit.txt",OutputDir, ThisTask);
  //if best fit directory doesn't exist, create it
  if(!(fa = fopen(buf, "w")))
    {
      sprintf(buf, "mkdir %s/bestfit/", OutputDir);
      system(buf);

      sprintf(buf, "%s/bestfit/task%d_bestfit.txt",OutputDir, ThisTask);
      if(!(fa = fopen(buf, "w")))
	{
	  sprintf(sbuf, "can't open file `%s'\n", buf);
	  terminate(sbuf);
	}
      print_parameters(2, fa);
    }
  else
    {
      print_parameters(2, fa);
    }
  //move comparison to bestfit directory
  sprintf(buf, "mv %s/task%d_mcmc_plus_obs_* %s/bestfit/", OutputDir, ThisTask,OutputDir);
  system(buf);

  //move the current best fit from all tasks into final_bestfit
  double tmp_MaxLike=0., MaxLike=0.;
  int dummy, MaxLikeTask=ThisTask;
  //check if another task has already printed a higher likelihood
  for(ii=0;ii<NTask;ii++)
    {
      sprintf(buf, "%s/bestfit/task%d_bestfit.txt",OutputDir, ii);
      if((fa = fopen(buf, "r")))
	{
	  fscanf(fa,"%d %lg",&dummy,&tmp_MaxLike);
	  if(pow(10,-1.*tmp_MaxLike)>MaxLike)
	    {
	      MaxLike=pow(10,-1.*tmp_MaxLike);
	      MaxLikeTask=ii;
	    }
	  fclose(fa);
	}
    }

  //if MaxLikeTask==ThisTask thistask is the highest likelihood, copy to final_bestfit
  if(MaxLikeTask==ThisTask)
    {
      sprintf(buf, "cp %s/bestfit/task%d_bestfit.txt %s/bestfit/final_bestfit.txt",OutputDir, ThisTask, OutputDir, ThisTask);
      system(buf);
      for(ObsNr=0;ObsNr<MCMCNConstraints;ObsNr++)
	for(snap=0;snap<NOUT;snap++)
	  if(MCMC_Obs[ObsNr].ObsTest_Switch_z[snap]==1)
	    {
	      sprintf(buf, "cp %s/bestfit/task%d_mcmc_plus_obs_%s_z%1.2f.txt %s/bestfit/final_mcmc_plus_obs_%s_z%1.2f.txt",
		      OutputDir,ThisTask,MCMC_Obs[ObsNr].Name,(double)((int)((MCMCConstraintsZZ[snap]*10)+0.5)/10.),
		      OutputDir,MCMC_Obs[ObsNr].Name,(double)((int)((MCMCConstraintsZZ[snap]*10)+0.5)/10.));
	      system(buf);
	    }
    }

}



/* initialize MCMC_PAR.Value and MCMC_PAR.PropValue with the same values
 *
 * the number of parameters, limits and switches (whitch to sample) are
 * read from MCMCParameterPriorsAndSwitches while the actual values are read from
 * previous output or MCMCStartingParFile
 *  */
void initialize_mcmc_par_and_lhood (FILE *fmcmc)
{
  int i, jj, snap, dumb_weight, EoF=0;
  double aux_p;
  char buf[1000];
  FILE *fa;

  sprintf(buf, "%s", MCMCParameterPriorsAndSwitches);
  if(!(fa = fopen(buf, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "can't open file `%s'\n", buf);
      terminate(sbuf);
    }

  fgets(buf, 300, fa);
  fscanf(fa,"%d\n",&MCMCNpar);

  //initialize structure to contain parameter names, values, priors and other properties
  MCMC_PAR = mymalloc("MCMC_PAR", sizeof(struct MCMC_PAR) * MCMCNpar);

  //read names and switches
  fgets(buf, 300, fa);
  for(i=0;i<MCMCNpar;i++)
    {
      fscanf(fa,"%s %lg %lg %lg %s %d\n",MCMC_PAR[i].Name, &MCMC_PAR[i].PropValue[0],
	     &MCMC_PAR[i].PriorMin, &MCMC_PAR[i].PriorMax, MCMC_PAR[i].Type, &MCMC_PAR[i].Sampling_Switch);
    }

  fclose(fa);  //done reading from MCMCParameterPriorsAndSwitches


  //read actual values from previous outputs and check if are inside priors
  jj=0;
  GlobalMCMCStep=0;
  do
    {
      //read lhood1
      if(fscanf(fmcmc,"%d %lg ", &dumb_weight, &lhood1)==EOF)
	EoF=1;//if there is no new line to read, fscanf gives error and EoF changes to 1
      else
	{
	  jj++;
	  //read parameter values
	  for(i=0;i<MCMCNpar;++i)
	    if(MCMC_PAR[i].Sampling_Switch==1)
	      fscanf(fmcmc,"%lg",&MCMC_PAR[i].Value[0]);
	  //printf("par[%d]=%f\n",i, MCMC_PAR[i].Value[0]);
	  GlobalMCMCStep+=1;
	}
    }
  while (EoF==0);

  //if there is something in that file,this is a re-start so:
  if(jj>0)
    MCMC_Initial_Par_Displacement=0.;

  //if there is nothing in that file read from backup file 00
  if(jj==0)
    {
      sprintf(buf, "%s", MCMCStartingParFile);
      if((fa = fopen(buf, "r")) == NULL)
	{
	  char sbuf[1000];
	  sprintf(sbuf, "can't open file `%s'\n", buf);
	  terminate(sbuf);
	}
      fscanf(fa,"%d %lg ", &dumb_weight, &lhood1);
      for(i=0;i<MCMCNpar;++i)
	if(MCMC_PAR[i].Sampling_Switch==1)
	  fscanf(fa,"%lg",&MCMC_PAR[i].Value[0]);
      fclose(fa);
    }

  //convert from log
  lhood1=pow(10,-lhood1);
  for(i=0;i<MCMCNpar;++i)
    if(MCMC_PAR[i].Sampling_Switch==1)
      MCMC_PAR[i].Value[0]=pow(10,MCMC_PAR[i].Value[0]);

  //for(i=0;i<MCMCNpar;++i)
  //		if(MCMC_PAR[i].Sampling_Switch==1)
  //			printf("par[%d]=%f\n",i, MCMC_PAR[i].Value[0]);

#ifdef PARALLEL
  //If PARALLEL only print initial parameter values for TASK 0
  if(ThisTask == 0)
    {
#endif
      printf("Initial Parameter Values:\n");

      for(i=0;i<MCMCNpar;i++)
	{
	  if(MCMC_PAR[i].Sampling_Switch==1)
	    {
	      printf("%s",MCMC_PAR[i].Name);
	      for(jj=0;jj<28-string_length(MCMC_PAR[i].Name);jj++)
		printf(" ");
	      printf("= %0.6f\n",MCMC_PAR[i].Value[0]);
	    }

	  if(MCMC_PAR[i].Sampling_Switch==1)
	    if(MCMC_PAR[i].Value[0]<MCMC_PAR[i].PriorMin ||
		MCMC_PAR[i].Value[0]>MCMC_PAR[i].PriorMax)
	      {
		printf("value=%0.6f priormin=%0.4f priormax=%0.4f\n",MCMC_PAR[i].Value[0],MCMC_PAR[i].PriorMin,MCMC_PAR[i].PriorMax);
		char sbuf[1000];
		sprintf(sbuf, "parameter '%s' outside prior range \n", MCMC_PAR[i].Name);
		terminate(sbuf);
	      }
	}
      printf("\n");
#ifdef PARALLEL
    }
#endif

  //if MCMC_Initial_Par_Displacement>0. introduce displacement in parameter values
  for(i=0;i<MCMCNpar;++i)
    {
      if (strcmp(MCMC_PAR[i].Type,"Physical")==0)
	for(snap=0;snap<NOUT;snap++)
	  {
	    if(Sample_Physical_Parameters==1)
	      {
		if(Time_Dependant_PhysPar==1 || snap==0)
		  {
		    if(MCMC_PAR[i].Sampling_Switch==1)
		      {
			aux_p=MCMC_PAR[i].Value[0];
			do
			  MCMC_PAR[i].Value[snap] = aux_p * exp(MCMC_Initial_Par_Displacement * gassdev(&MCMCseed));
			while(MCMC_PAR[i].Value[snap] < MCMC_PAR[i].PriorMin
			    || MCMC_PAR[i].Value[snap] > MCMC_PAR[i].PriorMax);
		      }
		    else
		      MCMC_PAR[i].Value[snap] = MCMC_PAR[i].Value[0];

		    MCMC_PAR[i].PropValue[snap] = MCMC_PAR[i].Value[snap];
		  }
		else //if(Time_Dependant_PhysPar==0 || snap>0)
		  {
		    MCMC_PAR[i].Value[snap] = MCMC_PAR[i].Value[0];
		    MCMC_PAR[i].PropValue[snap] = MCMC_PAR[i].Value[0];
		  }
	      }
	    else //if(Sample_Physical_Parameters==0)
	      {
		MCMC_PAR[i].Value[snap] = MCMC_PAR[i].Value[0];
		MCMC_PAR[i].PropValue[snap] = MCMC_PAR[i].Value[0];
	      }
	  }
      else if (strcmp(MCMC_PAR[i].Type,"Cosmological")==0)
	{
	  if(Sample_Cosmological_Parameters==1)
	    {
	      if(MCMC_PAR[i].Sampling_Switch==1)
		{
		  aux_p=MCMC_PAR[i].Value[0];
		  do
		    MCMC_PAR[i].Value[0] = aux_p * exp(MCMC_Initial_Par_Displacement * gassdev(&MCMCseed));
		  while(MCMC_PAR[i].Value[0] < MCMC_PAR[i].PriorMin
		      || MCMC_PAR[i].Value[0] > MCMC_PAR[i].PriorMax);
		}
	      MCMC_PAR[i].PropValue[0] = MCMC_PAR[i].Value[0];
	    }
	  else
	    MCMC_PAR[i].PropValue[0] = MCMC_PAR[i].Value[0];
	  //when running cosmology check to see if it was ok to just use snap=0 by default
	}
    }


  //zlist_0012_0047.txt
  /* p[0][6]=6.6;
      p[1][6]=5.9;
      p[2][6]=3.7;
      p[3][6]=7.4;
      p[0][9]=0.99;
      p[1][9]=0.95;
      p[2][9]=0.6;
      p[3][9]=0.93;*/

  //LOAD INTITIAL VALUES FOR ALL PARAMETERS FROM A DIFFERENT FILE
  if(Time_Dependant_PhysPar==1)
    {
      sprintf(buf, "./input/MCMC_inputs/mcmc_allz_par.txt");
      if(!(fa = fopen(buf, "r")))
	{
	  char sbuf[1000];
	  sprintf(sbuf, "can't open file `%s'\n", buf);
	  terminate(sbuf);
	}

      for(i=0;i<MCMCNpar;i++)
	for(snap=0;snap<NOUT;snap++)
	  if(i<1 || (i>2 && i<13))
	    {
	      fscanf(fa,"%lg",&MCMC_PAR[i].Value[snap]);
	      aux_p=MCMC_PAR[i].Value[snap];
	      do
		MCMC_PAR[i].Value[snap] = aux_p * exp(MCMC_Initial_Par_Displacement * gassdev(&MCMCseed));
	      while(MCMC_PAR[i].Value[snap] < MCMC_PAR[i].PriorMin || MCMC_PAR[i].Value[snap] > MCMC_PAR[i].PriorMax);
	      MCMC_PAR[i].PropValue[snap] = MCMC_PAR[i].Value[snap];
	    }
      fclose(fa);
    }


  //PROP = P in the first step

  /*for(i=0;i<MCMCNpar;++i)
  {
   	for(snap=0;snap<NOUT;snap++)
   	{
   		//MCMC_PAR[i].PropValue[snap] = MCMC_PAR[i].Value[snap];
   		printf("p[%d][%d]=%0.3lg  ", i,snap, MCMC_PAR[i].PropValue[snap]);
   	}
   printf("\n");
  }*/

  for(i=0;i<MCMCNpar;i++)
    if(MCMC_PAR[i].Sampling_Switch==1)
      {
#ifdef H2_AND_RINGS
	if (strcmp(MCMC_PAR[i].Name,"SfrColdCrit")==0)
	  terminate("MCMC parameter not free for current Makefile options: SfrColdCrit is not a free parameter if H2_AND_RINGS");
#endif
#ifdef DETAILED_METALS_AND_MASS_RETURN
	if (strcmp(MCMC_PAR[i].Name,"Yield")==0)
	  terminate("MCMC parameter not free for current Makefile options: Yield is not a free parameter if DETAILED_METALS_AND_MASS_RETURN");
#endif
      }
}






///////////
//PROPOSE//
///////////


/*@brief Function to propose new parameters given by
 *       a random normal distribution gassdev*/

double propose_new_parameters ()
{
  double qratio;
  int i, snap;

  for(i=0;i<MCMCNpar;++i)
    if(MCMC_PAR[i].Sampling_Switch==1)
      {
	if (strcmp(MCMC_PAR[i].Type,"Physical")==0)
	  {
	    if(Sample_Physical_Parameters==0)
	      for(snap=0;snap<NOUT;snap++)
		MCMC_PAR[i].PropValue[snap]= MCMC_PAR[i].Value[0];
	    else
	      for(snap=0;snap<NOUT;snap++)
		{
		  if(Time_Dependant_PhysPar==0 && snap>0)
		    MCMC_PAR[i].PropValue[snap] = MCMC_PAR[i].PropValue[0];
		  else
		    {
		      do
			MCMC_PAR[i].PropValue[snap] = MCMC_PAR[i].Value[snap] * exp(MCMC_LogStep_Size * gassdev(&MCMCseed));
		      while(MCMC_PAR[i].PropValue[snap] < MCMC_PAR[i].PriorMin
			  || MCMC_PAR[i].PropValue[snap] > MCMC_PAR[i].PriorMax);
		    }
		}
	  }
	else if (strcmp(MCMC_PAR[i].Type,"Cosmological")==0)
	  if(Sample_Cosmological_Parameters==0)
	    MCMC_PAR[i].PropValue[0]= MCMC_PAR[i].Value[0];
      }

  if(Sample_Cosmological_Parameters==1)
    do
      {
	for(i=0;i<MCMCNpar;++i)
	  if(MCMC_PAR[i].Sampling_Switch==1)
	    if (strcmp(MCMC_PAR[i].Type,"Cosmological")==0)
	      do
		MCMC_PAR[i].PropValue[0] = MCMC_PAR[i].Value[0] * exp(MCMC_LogStep_Size * gassdev(&MCMCseed));
	      while(MCMC_PAR[i].PropValue[0] < MCMC_PAR[i].PriorMin
		  || MCMC_PAR[i].PropValue[0] > MCMC_PAR[i].PriorMax);

	reset_cosmology();

      }while(fabs(ZZ[ListOutputSnaps[0]]-ListOutputRedshifts[0]) > 0.1);


  //to make just some parameters the same at all z

  /* if(Time_Dependant_PhysPar==1)
   for(i=0;i<MCMCNpar;++i)
  {
  	for(snap=0;snap<NOUT;snap++)
  	{
  		//if(i<3 || (i>3 && i<6) || (i>6 && i<9) )
	    //if(i<6 || (i>6 && i<9) || i>9)
  		if(i!=8 && i!=11)
  		//if(i!=11)
  		 {
  		 //IF THIS IS UNCOMMENTED ALL THE INTITIAL VALUES SHOULD BE THE SAME except for this parameters
  		  //COMMENT READING FROM FILE IN NEXT ROUTINE
  			 MCMC_PAR[i].PropValue[snap]=MCMC_PAR[i].PropValue[0];
  		 }
  	}
  }*/


  qratio = 1;
  //qratio= prop[0]/p[0]*prop[1]/p[1]*prop[2]/p[2]*prop[3]/p[3]*prop[4]/p[4];
  return qratio;
}



void read_mcmc_par (int snapnum)
{
  int snap, i;

  //betwenn snapnum[i] and snapnum[i+1] parameters have the values of snap[i+1]
  for(snap=0;snap<NOUT;snap++)
    if(snapnum < ListOutputSnaps[NOUT-snap-1]+1)
      break;
  snap=NOUT-snap-1;

  if(snap<0)
    snap=0;

  for(i=0;i<MCMCNpar;i++)
    if(MCMC_PAR[i].Sampling_Switch==1)
      {
	if(strcmp(MCMC_PAR[i].Name,"SfrEfficiency")==0)
	  SfrEfficiency = MCMC_PAR[i].PropValue[snap];
	if(strcmp(MCMC_PAR[i].Name,"SfrColdCrit")==0)
	  SfrColdCrit = MCMC_PAR[i].PropValue[snap];
	else if(strcmp(MCMC_PAR[i].Name,"SfrBurstEfficiency")==0)
	  SfrBurstEfficiency = MCMC_PAR[i].PropValue[snap];
	else if(strcmp(MCMC_PAR[i].Name,"SfrBurstSlope")==0)
	  SfrBurstSlope = MCMC_PAR[i].PropValue[snap];

	else if(strcmp(MCMC_PAR[i].Name,"AgnEfficiency")==0)
	  AgnEfficiency = MCMC_PAR[i].PropValue[snap];
	else if(strcmp(MCMC_PAR[i].Name,"BlackHoleGrowthRate")==0)
	  BlackHoleGrowthRate = MCMC_PAR[i].PropValue[snap];
	else if(strcmp(MCMC_PAR[i].Name,"BlackHoleDisruptGrowthRate")==0)
	  BlackHoleDisruptGrowthRate = MCMC_PAR[i].PropValue[snap];
	else if(strcmp(MCMC_PAR[i].Name,"BlackHoleCutoffVelocity")==0)
	  BlackHoleCutoffVelocity = MCMC_PAR[i].PropValue[snap];

	else if(strcmp(MCMC_PAR[i].Name,"FeedbackReheatingEpsilon")==0)
	  FeedbackReheatingEpsilon = MCMC_PAR[i].PropValue[snap];
	else if(strcmp(MCMC_PAR[i].Name,"ReheatPreVelocity")==0)
	  ReheatPreVelocity = MCMC_PAR[i].PropValue[snap];
	else if(strcmp(MCMC_PAR[i].Name,"ReheatSlope")==0)
	  ReheatSlope = MCMC_PAR[i].PropValue[snap];

	else if(strcmp(MCMC_PAR[i].Name,"FeedbackEjectionEfficiency")==0)
	  FeedbackEjectionEfficiency = MCMC_PAR[i].PropValue[snap];
	else if(strcmp(MCMC_PAR[i].Name,"EjectPreVelocity")==0)
	  EjectPreVelocity = MCMC_PAR[i].PropValue[snap];
	else if(strcmp(MCMC_PAR[i].Name,"EjectSlope")==0)
	  EjectSlope = MCMC_PAR[i].PropValue[snap];

	else if(strcmp(MCMC_PAR[i].Name,"ReIncorporationFactor")==0)
	  ReIncorporationFactor	= MCMC_PAR[i].PropValue[snap];
	else if(strcmp(MCMC_PAR[i].Name,"ReincZpower")==0)
	  ReincZpower = MCMC_PAR[i].PropValue[snap];
	else if(strcmp(MCMC_PAR[i].Name,"ReincVelocitypower")==0)
	  ReincVelocitypower = MCMC_PAR[i].PropValue[snap];

	else if(strcmp(MCMC_PAR[i].Name,"FracZtoHot")==0)
	  FracZtoHot = MCMC_PAR[i].PropValue[snap];

	else if(strcmp(MCMC_PAR[i].Name,"Yield")==0)
	  Yield = MCMC_PAR[i].PropValue[snap];

	else if(strcmp(MCMC_PAR[i].Name,"ThreshMajorMerger")==0)
	  ThreshMajorMerger = MCMC_PAR[i].PropValue[snap];

	else if(strcmp(MCMC_PAR[i].Name,"MergerTimeMultiplier")==0)
	  MergerTimeMultiplier = MCMC_PAR[i].PropValue[snap];

	else if(strcmp(MCMC_PAR[i].Name,"RamPressureStrip_CutOffMass")==0)
	  RamPressureStrip_CutOffMass = MCMC_PAR[i].PropValue[snap];

	else if(strcmp(MCMC_PAR[i].Name,"RamPressureRadiusThreshold")==0)
	  RamPressureRadiusThreshold = MCMC_PAR[i].PropValue[snap];

	else if(strcmp(MCMC_PAR[i].Name,"Reionization_z0")==0)
	  Reionization_z0 = MCMC_PAR[i].PropValue[snap];
	else if(strcmp(MCMC_PAR[i].Name,"Reionization_zr")==0)
	  {
	    Reionization_zr = MCMC_PAR[i].PropValue[snap];
	    if(Reionization_zr>(Reionization_z0-0.5))
	      {
		Reionization_zr=Reionization_z0-0.5;
		MCMC_PAR[i].PropValue[snap]=Reionization_z0-0.5;
	      }
	  }

	else if(strcmp(MCMC_PAR[i].Name,"GasInflowVel")==0)
	  GasInflowVel = MCMC_PAR[i].PropValue[snap];



	//printf("EjectSlope=%g\n",EjectSlope);
      }
}



/*@brief Read in the IDs and Weights of the FOFs groups that
 *       constitute the sample for which galaxies will be
 *       compared with observations
 *
 *       There is a single file with merger trees that needs to be created before
 *       running the code. The file contains the trees corresponding to the IDs
 *       read here. It contains both MR and MRII trees. The MR trees are stored first.
 *       NOTE THAT THERE IS MORE THAN 1 FOF GROUP PER TREE. Another file is
 *       created with the number of trees in MR. After the code (in main()) has
 *       looped over that number of trees the necessary variables are changed to MRII.*/
void read_sample_info (void)
{
  int DumbTreeNrColector, DumbFileNrColector, MaxFoFNr;
  int i, snap;
  FILE *fa;
  char buf[1000];

  MaxFoFNr=0;


#ifdef MR_PLUS_MRII
  if(Switch_MR_MRII==1)
    {
      sprintf(buf, "%s/%ssample_allz_nh_Switch_MR_MRII_%d.dat", MCMCSampleDir, MCMCSampleFilePrefix, MCMCSampleFile);
      if(!(fa = fopen(buf, "r")))
	{
	  char sbuf[1000];
	  sprintf(sbuf, "can't open file `%s'\n", buf);
	  terminate(sbuf);
	}
      fscanf(fa, "%d \n", &NTrees_Switch_MR_MRII);
      fclose(fa);
    }
#endif

  for(snap=0;snap<NOUT;snap++)
    {
      sprintf(buf, "%s/%ssample_allz_nh_%d%d.dat", MCMCSampleDir, MCMCSampleFilePrefix,
	      MCMCSampleFile, ListOutputSnaps[snap]);
      if(!(fa = fopen(buf, "r")))
	{
	  char sbuf[1000];
	  sprintf(sbuf, "can't open file `%s'\n", buf);
	  terminate(sbuf);
	}

      fscanf(fa, "%d \n", &NFofsInSample[snap]);

      if(MaxFoFNr<NFofsInSample[snap])
	MaxFoFNr=NFofsInSample[snap];

      fclose(fa);
    }


  //structure as the size of highest number of halos in a given snapshot
  //for all the other snapshot it will not be full
  MCMC_FOF = malloc(sizeof(struct MCMC_FOF_struct) * MaxFoFNr);

  for(snap=0;snap<NOUT;snap++)
    {

      sprintf(buf, "%s/%ssample_allz_nh_%d%d.dat", MCMCSampleDir, MCMCSampleFilePrefix,
	      MCMCSampleFile, ListOutputSnaps[snap]);
      if(!(fa = fopen(buf, "r")))
	{
	  char sbuf[1000];
	  sprintf(sbuf, "can't open file `%s'\n", buf);
	  terminate(sbuf);
	}

      fscanf(fa, "%d \n", &NFofsInSample[snap]);

      for(i=0;i<NFofsInSample[snap];i++)
	{
	  fscanf(fa, "%lld %d %d %lg\n", &MCMC_FOF[i].FoFID[snap], &DumbTreeNrColector, &DumbFileNrColector, &MCMC_FOF[i].Weight[snap]);
	  MCMC_FOF[i].Weight[snap]/=BoxSize*BoxSize*BoxSize;
	  MCMC_FOF[i].NGalsInFoF[snap]=0;
	  MCMC_FOF[i].IndexOfCentralGal[snap]=-1;
	}

      fclose(fa);
    }
}



/*@brief Read in the arrays of observational data. They will be compared
 *       with the outputs from the SAM On the get_likelihood routine*/
void read_observations (void)
{
  int i, j, kk, snap, number_of_tests, n_snaps;
  //FILE *f[MCMCNConstraints];
  FILE *fa;
  float BinValueColector;
  char buf[1000],	sbuf[1000], aux_testName[1000], aux_testType[1000];

  //allocate structure to contain observational data
  MCMC_Obs = mymalloc("MCMC_Obs", sizeof(struct MCMC_OBSCONSTRAINTS) * MCMCNConstraints);


  /*Read file MCMCObsConstraints.txt that contains a header with
   * number_of_tests, number_of_chi_tests, number_of_binom_tests,
   * n_snaps and the redshifts that will be used to constrain the MCMC
   * Then there is a list of all the observational constraints, with
   * the type of test to be used and a switch controlling which redshifts
   * will be used for each constraint */
  sprintf(buf, "%s", MCMCObsConstraints);
  if(!(fa = fopen(buf, "r")))
    {
      sprintf(sbuf, "can't open file `%s'\n", buf);
      terminate(sbuf);
    }


  fgets(buf, 500, fa);
  fscanf(fa,"%d\n",&number_of_tests);

  if(number_of_tests!=MCMCNConstraints)
    {
      sprintf(sbuf, "check MCMCNConstraints & number_of_tests in mcmc_var.h and %s\n",MCMCObsConstraints);
      terminate(sbuf);
    }

  fgets(buf, 500, fa);
  fgets(buf, 500, fa);
  fscanf(fa,"%d\n",&n_snaps);
  if(n_snaps > NOUT)
    {
      sprintf(sbuf, "n_snaps > NOUT\n");
      terminate(sbuf);
    }



  //Check if all the required snaps are in the current ListOutputSnaps
  for(i=0;i<n_snaps;i++)
    {
      fscanf(fa,"%lg\n",&MCMCConstraintsZZ[i]);
      kk=0;
      for(j=0;j<NOUT;j++)
	{
	  if(MCMCConstraintsZZ[i] >= (double)((int)((ZZ[ListOutputSnaps[j]]*10)+0.5)/10.)-0.1 &&
	      MCMCConstraintsZZ[i] <= (double)((int)((ZZ[ListOutputSnaps[j]]*10)+0.5)/10.)+0.1)
	    kk+=1;
	}
      if(kk==0)
	{
	  sprintf(sbuf, "redshift %0.2f required for MCMC not in outputlist \n",MCMCConstraintsZZ[i]);
	  terminate(sbuf);
	}
    }

  fgets(buf, 500, fa);

  //Scan test names, types and redshift switches
  for(i=0;i<MCMCNConstraints;i++)
    {
      fscanf(fa,"%s %s",MCMC_Obs[i].Name, MCMC_Obs[i].TestType);
      for(j=0;j<NOUT;j++)
	fscanf(fa,"%d",&MCMC_Obs[i].ObsTest_Switch_z[j]);
    }
  fclose(fa);




  //now read weights for different constraints
  sprintf(buf, "%s", MCMCWeightsObsConstraints);
  if(!(fa = fopen(buf, "r")))
    {
      sprintf(sbuf, "can't open file `%s'\n", buf);
      terminate(sbuf);
    }

  fgets(buf, 500, fa);
  fscanf(fa,"%d\n",&number_of_tests);

  if(number_of_tests!=MCMCNConstraints)
    {
      sprintf(sbuf, "check MCMCNConstraints & number_of_tests in mcmc_var.h and %s\n",MCMCWeightsObsConstraints);
      terminate(sbuf);
    }
  fgets(buf, 500, fa);
  for(i=0;i<MCMCNConstraints;i++)
    {
      fscanf(fa,"%s %s",aux_testName, aux_testType);
      for(j=0;j<NOUT;j++)
	fscanf(fa,"%lf",&MCMC_Obs[i].ObsTest_Weight_z[j]);
    }
  fclose(fa);



  //now read the observations
  for(snap=0;snap<NOUT;snap++)
    {
      //the round of ZZ[] is to assure that independently of the cosmology used you still
      //round to z=1.0 or 2.0,etc...
      for(i=0;i<MCMCNConstraints;++i)
	{
	  if(MCMC_Obs[i].ObsTest_Switch_z[snap]==1)
	    {

	      sprintf(buf, "%s/%s_z%1.2f.txt",ObsConstraintsDir,MCMC_Obs[i].Name,
		      (double)((int)((MCMCConstraintsZZ[snap]*10)+0.5)/10.) );
	      if((fa=fopen(buf,"r"))==NULL)
		{
		  sprintf(sbuf, "can't open file `%s'\n", buf);
		  terminate(sbuf);
		}

	      /* This values will give the size of the observational arrays
	       * They will be used in get_likelihood. */
	      fscanf(fa, "%d", &Nbins[snap][i]);
	      if(Nbins[snap][i]>MCMCMaxObsBins)
		{
		  sprintf(sbuf, "NBins for Obs[%d] Snap[%d] > MCMCMaxObsBins", i,snap);
		  terminate(sbuf);
		}

	      //Read observational data
	      for(j = 0; j < Nbins[snap][i]; j++)
		{
		  //Chi_Sq and Maximum Likelihood TESTS
		  if(strcmp(MCMC_Obs[i].TestType,"chi_sq")==0 || strcmp(MCMC_Obs[i].TestType,"maxlike")==0)
		    fscanf(fa, "%lg %lg %lg %lg", &MCMC_Obs[i].Bin_low[snap][j], &MCMC_Obs[i].Bin_high[snap][j],
			   &MCMC_Obs[i].Obs[snap][j], &MCMC_Obs[i].Error[snap][j]);
		  //Binomial TESTS
		  else if(strcmp(MCMC_Obs[i].TestType,"binomial")==0)
		    fscanf(fa, "%f %lg %lg", &BinValueColector, &MCMC_Obs[i].ObsUp[snap][j], &MCMC_Obs[i].ObsDown[snap][j]);
		}
	      fclose(fa);

	    }//end if(MCMC_Obs[i].ObsTest_Switch_z[snap]==1)
	}//end loop on tests
    }//end loop on snaps

}


/* A different file is written for each observational constraint and for each redshift
 * Each file as the following structure: *
 * int ChainLength+1
 * int Nbins[snap][constraint]
 * phi[0],phi[1],phi[Nbins[snap][constraint]],likelihood(of this constraint),likelihood(of this step)
 *
 * Plus one file to contain the total likelihood at each step*/
void open_files_with_comparison_to_observations()
{
  int constraint, snap;
  char buf[1000], sbuf[1000];

  sprintf(buf, "%sMCMC_LIKELIHOOD_%d.txt", OutputDir, ThisTask+FirstChainNumber);
  if((FILE_MCMC_LIKELIHOOD = fopen(buf, "w")) == NULL)
    {
      char sbuf[1000];
      sprintf(sbuf, "can't open file `%s'\n", buf);
      terminate(sbuf);
    }

  for(constraint = 0; constraint < MCMCNConstraints; constraint++)
    for(snap=0;snap<NOUT;snap++)
      {
	if(MCMC_Obs[constraint].ObsTest_Switch_z[snap]==1)
	  {
	    sprintf(buf, "%sPredictionsPerStep_%s_z%1.2f_%d.txt", OutputDir, MCMC_Obs[constraint].Name,
		    (double)((int)((MCMCConstraintsZZ[snap]*10)+0.5)/10.), ThisTask+FirstChainNumber);
	    if((FILE_MCMC_PredictionsPerStep[snap][constraint] = fopen(buf, "w")) == NULL)
	      {
		char sbuf[1000];
		sprintf(sbuf, "can't open file `%s'\n", buf);
		terminate(sbuf);
	      }
	    fprintf(FILE_MCMC_PredictionsPerStep[snap][constraint], " %d\n", ChainLength+1);
	    fprintf(FILE_MCMC_PredictionsPerStep[snap][constraint], " %d\n", Nbins[snap][constraint]);
	  }
      }
}

void close_files_with_comparison_to_observations()
{
  int constraint, snap;

  fclose(FILE_MCMC_LIKELIHOOD);

  for(constraint = 0; constraint < MCMCNConstraints; constraint++)
    for(snap=0;snap<NOUT;snap++)
      if(MCMC_Obs[constraint].ObsTest_Switch_z[snap]==1)
	fclose(FILE_MCMC_PredictionsPerStep[snap][constraint]);

}


#ifdef MR_PLUS_MRII
void change_dark_matter_sim(char SimName[])
{

  if (strcmp(SimName,"MR")==0)
    {
      Switch_MR_MRII=1;
      sprintf(FileWithZList, "%s", FileWithZList_MR);
      PartMass=PartMass_MR;
      BoxSize=BoxSize_MR;

      sprintf(FileWithZList_OriginalCosm, "%s", FileWithZList_OriginalCosm_MR);
      PartMass_OriginalCosm=PartMass_OriginalCosm_MR;
      BoxSize_OriginalCosm=BoxSize_OriginalCosm_MR;

      LastDarkMatterSnapShot=LastDarkMatterSnapShot_MR;
    }
  else if (strcmp(SimName,"MRII")==0)
    {
      Switch_MR_MRII=2;
      sprintf(FileWithZList, "%s", FileWithZList_MRII);
      PartMass=PartMass_MRII;
      BoxSize=BoxSize_MRII;

      sprintf(FileWithZList_OriginalCosm, "%s", FileWithZList_OriginalCosm_MRII);
      PartMass_OriginalCosm=PartMass_OriginalCosm_MRII;
      BoxSize_OriginalCosm=BoxSize_OriginalCosm_MRII;

      LastDarkMatterSnapShot=LastDarkMatterSnapShot_MRII;
    }


  //do part of init() again.
  //Can't do the all functions because of metallicity arrays in read_cooling_functions()
  read_zlist();
  read_zlist_original_cosm();
  read_output_snaps();

  //CREATE ARRAYS OF SFH TIME STRUCTURE:
#ifdef  STAR_FORMATION_HISTORY
  create_sfh_bins();
#endif

#ifdef COMPUTE_SPECPHOT_PROPERTIES
  //read in photometric tables
#ifdef PHOTTABLES_PRECOMPUTED
  setup_LumTables_precomputed(SimName);
#endif
#ifdef SPEC_PHOTABLES_ON_THE_FLY
  setup_Spec_LumTables_onthefly();
#endif
#endif

  //READ IN THE YIELD TABLES, AND FORM NORMALISED YIELD ARRAYS:
#ifdef DETAILED_METALS_AND_MASS_RETURN
  //read_yield_tables();
  //#ifdef INDIVIDUAL_ELEMENTS
  //  SNe_rates();
  //#endif
  init_integrated_yields();
  integrate_yields();
#endif


  /*After all the millennium trees have been done read sample for MRII trees
   * from treenr=NTrees_MR to NTrees_MR+NTrees_MRII=Ntrees */
  if (strcmp(SimName,"MR")==0)
    {
      sprintf(MCMCSampleFilePrefix,"%s",MCMCSampleFilePrefix_MR);
      MCMCSampleFile=MCMCSampleFile_MR;
    }
  else if (strcmp(SimName,"MRII")==0)
    {
      sprintf(MCMCSampleFilePrefix,"%s",MCMCSampleFilePrefix_MRII);
      MCMCSampleFile=MCMCSampleFile_MRII;
    }

  if (strcmp(SimName,"MRII")==0)
    {
      free(MCMC_FOF);
    }

  read_sample_info();


}
#endif

#ifdef HALOMODEL
void assign_FOF_masses(snapnum, treenr)
{
  int fof,snap, halonr;

  for(halonr = 0; halonr < TreeNHalos[treenr]; halonr++)
    if(HaloAux[halonr].DoneFlag == 0 && Halo[halonr].SnapNum == snapnum)
      {
	for(snap=0;snap<NOUT;snap++)
	  {
	    if(snapnum==ListOutputSnaps[snap])
	      {
		for(fof=0;fof<NFofsInSample[snap]; fof++)
		  if(HaloIDs[halonr].FirstHaloInFOFgroup == MCMC_FOF[fof].FoFID[snap])
		    {
		      MCMC_FOF[fof].M_Crit200[snap] = log10(Halo[halonr].M_Crit200*1.e10);
		      MCMC_FOF[fof].M_Mean200[snap] = log10(Halo[halonr].M_Mean200*1.e10);
#ifdef MCRIT
MCMC_FOF[fof].M_Mean200[snap] = log10(Halo[halonr].M_Crit200*1.e10);
#endif
		    }
	      }
	  }
      }

}
#endif



//////////
//GASDEV//
//////////

//Gives a random normal deviate using ran3 (ran1 NR)


double gassdev(long *idum)
{
  static int iset = 0;
  static double gset;
  double fac, r, v1, v2;

  if(iset == 0)
    {
      do
	{
	  v1 = 2.0 * ran3(idum) - 1.0;
	  v2 = 2.0 * ran3(idum) - 1.0;
	  r = v1 * v1 + v2 * v2;
	}
      while(r >= 1.0 || r == 0.0);
      fac = sqrt(-2.0 * log(r) / r);

      //Box Muller deviates to get two normal deviates
      gset = v1 * fac;
      iset = 1;
      return v2 * fac;
    }

  else
    {
      iset = 0;
      return gset;
    }
}






////////
//RAN3//
////////

double ran3(long *idum)
{
  int j;
  long k;
  static long iy = 0;
  static long iv[NTAB];
  double temp;

  if(*idum <= 0 || !iy)
    {
      if(-(*idum) < 1)
	*idum = 1;
      else
	*idum = -(*idum);
      for(j = NTAB + 7; j >= 0; j--)
	{
	  k = (*idum) / IQ;
	  *idum = IA * (*idum - k * IQ) - IR * k;
	  if(*idum < 0)
	    *idum += IM;
	  if(j < NTAB)
	    iv[j] = *idum;
	}

      iy = iv[0];
    }

  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if(*idum < 0)
    *idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = *idum;

  if((temp = AM * iy) > RNMX)
    return RNMX;
  else
    return temp;
}
#endif //MCMC



#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef RNMX
#undef EPS
