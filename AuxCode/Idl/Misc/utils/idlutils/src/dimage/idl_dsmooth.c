#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

int dsmooth(float *image, int nx, int ny, float sigma, float *smooth);

/********************************************************************/
IDL_LONG idl_dsmooth (int      argc,
                      void *   argv[])
{
	IDL_LONG nx,ny;
	float *image, *smooth, sigma;
	
	IDL_LONG i;
	IDL_LONG retval=1;

	/* 0. allocate pointers from IDL */
	i=0;
	image=((float *)argv[i]); i++;
	nx=*((int *)argv[i]); i++;
	ny=*((int *)argv[i]); i++;
  sigma=*((float *)argv[i]); i++;
  smooth=((float *)argv[i]); i++;
	
	/* 1. run the fitting routine */
	retval=(IDL_LONG) dsmooth(image, nx, ny, sigma, smooth);
	
	/* 2. free memory and leave */
	free_memory();
	return retval;
}

/***************************************************************************/

