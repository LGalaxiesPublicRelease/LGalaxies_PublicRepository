#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ph.h"

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

int init_interp_profmean(int in_nprof, float *in_profmean,
                          float *profradius, int in_nrad, float *maxradius, 
                          float in_rscale, float in_pmscale);
float interp_prof(float radius);

/********************************************************************/
IDL_LONG idl_interp_profmean(int      argc,
                             void *   argv[])
{
	IDL_LONG nprof, nrad; 
  float *profmean, *profradius, *maxradius, rscale, pmscale,radius,*value;
	
	IDL_LONG i;
	IDL_LONG retval=1;

	/* 0. allocate pointers from IDL */
	i=0;
	nprof=*((int *)argv[i]); i++;
	profmean=((float *)argv[i]); i++;
	profradius=((float *)argv[i]); i++;
	nrad=*((int *)argv[i]); i++;
	pmscale=*(float *)argv[i]; i++;
	rscale=*(float *)argv[i]; i++;
	maxradius=((float *)argv[i]); i++;
  radius=*(float *)argv[i]; i++;
	value=((float *)argv[i]); i++;
	
	/* 1. run the fitting routine */
	retval=(IDL_LONG) init_interp_profmean(nprof,profmean,profradius, 
                                         nrad,maxradius,rscale,pmscale);
  *value=interp_prof(radius) ;
    
	/* 2. free memory and leave */
	free_memory();
	return retval;
}

/***************************************************************************/

