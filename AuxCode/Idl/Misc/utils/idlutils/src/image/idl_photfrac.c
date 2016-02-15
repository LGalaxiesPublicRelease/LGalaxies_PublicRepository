#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ph.h"

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

/********************************************************************/
IDL_LONG idl_photfrac (int      argc,
                       void *   argv[])
{
	IDL_LONG xnpix,ynpix,xcen,ycen;
	float *frac,radius;
	
	IDL_LONG i;
	IDL_LONG retval=1;

	/* 0. allocate pointers from IDL */
	i=0;
	xnpix=*((int *)argv[i]); i++;
	ynpix=*((int *)argv[i]); i++;
	radius=*((float *)argv[i]); i++;
	frac=(float *)argv[i]; i++;
	xcen=*((int *)argv[i]); i++;
	ycen=*((int *)argv[i]); i++;
	
	/* 1. run the fitting routine */
	retval=(IDL_LONG) photfrac(xnpix,ynpix,radius,frac,xcen,ycen);
	
	/* 2. free memory and leave */
	free_memory();
	return retval;
}

/***************************************************************************/

