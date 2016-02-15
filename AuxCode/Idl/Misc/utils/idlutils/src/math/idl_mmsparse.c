#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h" 

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

int mmcomp(float *b, int nxb, int nyb, float *a, int *ax, int nxa, 
           int xstart, int i, float *comp);

/********************************************************************/
IDL_LONG idl_mmsparse (int argc, 
                       void *argv[])
{
  IDL_LONG nx,ny,nz,*rowstart,*nxrow,*x;
  float *cc, *bb, *val;
    
	IDL_LONG i,j,k,l;
  float comp;
	IDL_LONG retval=1;

	/* 0. allocate pointers from IDL */
	i=0;
	cc=((float *)argv[i]); i++;
	bb=((float *)argv[i]); i++;
	nx=*((IDL_LONG *)argv[i]); i++;
	ny=*((IDL_LONG *)argv[i]); i++;
	nz=*((IDL_LONG *)argv[i]); i++;
	val=((float *)argv[i]); i++;
	x=((IDL_LONG *)argv[i]); i++;
	rowstart=((IDL_LONG *)argv[i]); i++;
	nxrow=((IDL_LONG *)argv[i]); i++;

  for(i=0L; i<nz; i++) {
    for(j=0L; j<ny; j++) {
      comp=0.;
      for(l=0; l<nxrow[j]; l++) {
        k=x[rowstart[j]+l];
        comp+=val[rowstart[j]+l]*bb[k+i*nx];
      }
      cc[i+j*nz]=comp;
    }
  }
	
	return retval;
}

/***************************************************************************/

