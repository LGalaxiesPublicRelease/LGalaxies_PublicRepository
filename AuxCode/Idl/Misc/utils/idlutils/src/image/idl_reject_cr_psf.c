#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ph.h"

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

int reject_cr_psf(float *image, float *image_ivar, int xnpix, int ynpix,
                  float nsig, float *psfvals, float cfudge, float c2fudge,
                  int *rejected, int *ignoremask);

/********************************************************************/
IDL_LONG idl_reject_cr_psf(int      argc,
                           void *   argv[])
{
	IDL_LONG i;
	IDL_LONG retval=1;

  float *image, *image_ivar, nsig, *psfvals, cfudge, c2fudge;
  IDL_LONG *rejected, *ignoremask;
  IDL_LONG xnpix, ynpix;

	/* 0. allocate pointers from IDL */
	i=0;
	image=((float *)argv[i]); i++;
	image_ivar=((float *)argv[i]); i++;
	xnpix=*((IDL_LONG *)argv[i]); i++;
	ynpix=*((IDL_LONG *)argv[i]); i++;
	nsig=*(float *)argv[i]; i++;
	psfvals=(float *)argv[i]; i++;
	cfudge=*((float *)argv[i]); i++;
  c2fudge=*(float *)argv[i]; i++;
	rejected=((IDL_LONG *)argv[i]); i++;
	ignoremask=((IDL_LONG *)argv[i]); i++;
	
	/* 1. run the fitting routine */
	retval=(IDL_LONG) reject_cr_psf(image,image_ivar,(int) xnpix,(int) ynpix,
                                  nsig,psfvals,cfudge,c2fudge,(int *)rejected, 
                                  (int *)ignoremask);
  
	/* 2. free memory and leave */
	free_memory();
	return retval;
}

/***************************************************************************/

