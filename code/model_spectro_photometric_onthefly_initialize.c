#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <stddef.h>

#include "allvars.h"
#include "proto.h"

/* void read_z_list()
 *
 * void read_vega_spectra(double *LambdaVega, double *FluxVega)
 * void read_filters()
 * void read_InputSSP_spectra(int zz) */


/* Setup table with redshifts to calculate observed frame magnitudes.
 * The list is different from the one on the code (ZZ[]) since no negative
 * times are allowed, which can happen for scaled cosmologies. */
void setup_RedshiftTab()
{
  int snap;

  for(snap=0;snap<(LastDarkMatterSnapShot+1);snap++)
    {
      RedshiftTab[snap]=ZZ[snap];
      //for the photometry calculation do not allow negative times (it can happen for scaled cosmologies)
      if(ZZ[snap]<0.0)
	RedshiftTab[snap]=0.;
    }
}


void read_vega_spectra(double *LambdaVega, double *FluxVega)
{
  int i;
  char buf[1000];
  FILE *fa;

  sprintf(buf, "%s/FullSEDs/VEGA_A0V_KUR_BB.SED",SpecPhotDir);
  if((fa=fopen(buf,"r"))==NULL)
    {
      char sbuf[1000];
      sprintf(sbuf, "Can't open file %s\n", buf);
      terminate(sbuf);
    }

  for (i=0;i<NLambdaVega;i++)
    {
      fscanf(fa,"%lf %lf" , &LambdaVega[i], &FluxVega[i]);
      //convert to ergs.cm^-2.s^-1.AA-1
      FluxVega[i]*=2e-17;
      FluxVega[i]=FluxVega[i]*LambdaVega[i]*LambdaVega[i]/(C*1e8);
    }
  fclose(fa);
}

void read_filters(double LambdaFilter[NMAG][MAX_NLambdaFilter], double FluxFilter[NMAG][MAX_NLambdaFilter])
{
  int j, bandn, NFilters;
  FILE *fa, *fb;
  char buf[1000], buf2[1000], FilterFile[1000], FilterName[1000];

  sprintf(buf, "%s",FileWithFilterNames);
  if((fa=fopen(buf,"r"))==NULL)
    {
      char sbuf[1000];
      sprintf(sbuf, "Can't open file %s\n", buf);
      terminate(sbuf);
    }

  fscanf(fa,"%d" ,&NFilters);
  if (NFilters != NMAG) {printf("NFilters not equal to  NMAG, line %d of read_filters.c!!! ",__LINE__);exit(0);}

  for(bandn=0;bandn<NMAG;bandn++)
    {
      fscanf(fa,"%s %f %s" ,FilterFile, &FilterLambda[bandn],FilterName);
      sprintf(buf2, "%s/Filters/%s",SpecPhotDir,FilterFile);

      if((fb=fopen(buf2,"r"))==NULL)
	{
	  char sbuf[1000];
	  sprintf(sbuf, "Can't open file %s\n", buf2);
	  terminate(sbuf);
  	}

      fscanf(fb,"%d" ,&NLambdaFilter[bandn]);
      if(NLambdaFilter[bandn]>MAX_NLambdaFilter)
	{
	  char sbuf[1000];
	  sprintf(sbuf, "NLambdaFilter[%d]>MAX_NLambdaFilter \n", bandn);
	  terminate(sbuf);
  	}

      for(j=0;j<NLambdaFilter[bandn];j++)
	fscanf(fb,"%lf %lf" ,&LambdaFilter[bandn][j], &FluxFilter[bandn][j]);

      fclose(fb);
  }

  fclose(fa);

}

void read_MetalTab()
{
  int i,dumb_ssp_nmetallicites;
  FILE *fa;
  char buf[1000];
#ifdef M05
  char *SSP = {"M05"};
#endif
#ifdef CB07
  char *SSP = {"CB07"};
#endif
#ifdef BC03
  char *SSP = {"BC03"};
#endif

  sprintf(buf, "%s/FullSEDs/%s_%s_Metallicity_list.dat", SpecPhotDir, SSP, SpecPhotIMF);
  if(!(fa = fopen(buf, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "file `%s' not found.\n", buf);
      terminate(sbuf);
    }

  fscanf(fa, "%d", &dumb_ssp_nmetallicites);
  if(dumb_ssp_nmetallicites != SSP_NMETALLICITES)
    {
      terminate("nmetallicites on file not equal to SSP_NMETALLICITES");
    }

  for(i=0;i<SSP_NMETALLICITES;i++)
    {
      fscanf(fa, "%f", &SSP_logMetalTab[i]);
      SSP_logMetalTab[i]=log10(SSP_logMetalTab[i]);
    }

  fclose(fa);
}

/* Reads in the SSP full spectra for the current metallicity
 * LambdaInputSSP[NAGE][NLambdaInputSSP] &  FluxInputSSP[NAGE][NLambdaInputSSP];
 * original units of FluxInputSSP[i] are (erg.s^-1.AA^-1);
 * Flux*Lambda^2/Clight*1e8 converts it to (erg.s^-1.Hz^-1) */
void read_InputSSP_spectra(double LambdaInputSSP[SSP_NAGES][SSP_NLambda], double FluxInputSSP[SSP_NAGES][SSP_NLambda], int MetalLoop)
{
  double Dumb1, age;
  int i, ageloop;
  FILE *fa;
  char buf1[1000];
#ifdef M05
  char *SSP = {"M05"};
#endif
#ifdef CB07
  char *SSP = {"CB07"};
#endif
#ifdef BC03
  char *SSP = {"BC03"};
#endif

  sprintf(buf1, "%s/FullSEDs/%s_%s_FullSED_m%0.4f.dat",SpecPhotDir, SSP, SpecPhotIMF, pow(10,SSP_logMetalTab[MetalLoop]));

  if((fa=fopen(buf1,"r"))==NULL)
    {
      char sbuf[1000];
      sprintf(sbuf, "Can't open file %s\n", buf1);
      terminate(sbuf);
    }

  for(ageloop=0;ageloop<SSP_NAGES;ageloop++)
    for(i=0;i<SSP_NLambda;i++)
      {
    	LambdaInputSSP[ageloop][i]=0.0;
    	FluxInputSSP[ageloop][i]=0.0;
      }

  for(ageloop=0;ageloop<SSP_NAGES;ageloop++)
    {
      for (i=0;i<SSP_NLambda;i++)
	{
	  fscanf(fa,"%lf %lf %lf %lf\n" , &age, &Dumb1,&LambdaInputSSP[ageloop][i], &FluxInputSSP[ageloop][i]);
	  FluxInputSSP[ageloop][i]=1e11*FluxInputSSP[ageloop][i]*LambdaInputSSP[ageloop][i]*LambdaInputSSP[ageloop][i]/(C*1.e8);
  	}
      if(MetalLoop==0) //only read age table once
	if(age>0.)
	  SSP_logAgeTab[ageloop]=log10(age / 1.0e6 / UnitTime_in_Megayears * Hubble_h);
	else
	  SSP_logAgeTab[ageloop]=0.;
    }

  fclose(fa);
}
