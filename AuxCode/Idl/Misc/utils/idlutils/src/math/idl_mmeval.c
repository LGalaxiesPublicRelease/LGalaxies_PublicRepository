#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h" 

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

/********************************************************************/
IDL_LONG idl_mmeval (int argc, 
                       void *argv[])
{
  IDL_LONG nx,ny,nz,*rowstart,*nxrow,*x;
  float *aa, *bb, *val;
    
	IDL_LONG i,j,k,l;
  float comp;
	IDL_LONG retval=1;

	/* 0. allocate pointers from IDL */
	i=0;
	val=((float *)argv[i]); i++;
	bb=((float *)argv[i]); i++;
	aa=((float *)argv[i]); i++;
	nx=*((IDL_LONG *)argv[i]); i++;
	ny=*((IDL_LONG *)argv[i]); i++;
	nz=*((IDL_LONG *)argv[i]); i++;
	x=((IDL_LONG *)argv[i]); i++;
	rowstart=((IDL_LONG *)argv[i]); i++;
	nxrow=((IDL_LONG *)argv[i]); i++;

  for(j=0L; j<ny; j++) {
    for(l=0L; l<nxrow[j]; l++) {
      i=x[rowstart[j]+l];
      comp=0.;
      for(k=0;k<nx;k++) 
        comp+=bb[k+i*nx]*aa[k+j*nx];
      val[rowstart[j]+l]=comp;
    }
  }

	return retval;
}

/***************************************************************************/

