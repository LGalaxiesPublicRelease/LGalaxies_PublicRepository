#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <polygon.h>
#include "export.h"

IDL_LONG gvert(polygon *poly, int vcirc, double *tol, int nvmax, int per, int nve, 
          int *nv, vec ve[/*nvmax * nve*/], double angle[/*nvmax*/], 
          int ipv[/*nvmax*/], int gp[/*poly->np*/], int *nev, int *nev0, 
          int ev[/*nvmax*/]);
IDL_LONG vmidc(polygon *poly, int nv, int nve, 
         vec ve[/*nv * nve*/], int ipv[/*nv*/], int ev[/*nv*/], int *nvm, 
         vec **vm_p);

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

#define DEG2RAD .01745329251994

/********************************************************************/
IDL_LONG idl_vmidc
  (int      argc,
   void *   argv[])
{
  double *outvm;

	polygon *poly;
  IDL_LONG nvm;
	double *vm_p;
	double    tol;
	IDL_LONG vcirc;
	IDL_LONG nv;
	IDL_LONG nve;
	double *ve_p;
	double *angle_p;
	IDL_LONG *ipv_p;
	IDL_LONG *gp_p;
	IDL_LONG nev;
	IDL_LONG nev0;
	IDL_LONG *ev_p;
	IDL_LONG i;
	IDL_LONG retval=1;

  poly=NULL;
  vm_p=NULL;
  ve_p=NULL;
  angle_p=NULL;
  ipv_p=NULL;
  gp_p=NULL;
  ev_p=NULL;
  outvm=NULL;
  nvm=0;

	poly=(polygon *) malloc(sizeof(polygon));
  nve=2;

	/* 0. allocate pointers from IDL */  
	i=0;
	poly->rp = (vec *) argv[i]; i++;
	poly->cm = (double *)argv[i]; i++;
	poly->np = *((IDL_LONG *)argv[i]); i++;
	poly->npmax = poly->np;
  outvm =(double *)argv[i]; i++;
  nv=poly->np*3*(poly->np*3+1)+4;

	/* 1. we know sizes, so we can do the real thing */
  ve_p=(double *) malloc(3*nv*nve*sizeof(double));
  angle_p=(double *) malloc(nv*sizeof(double));
  ipv_p=(IDL_LONG *) malloc(nv*sizeof(IDL_LONG));
  gp_p=(IDL_LONG *) malloc(poly->np*sizeof(IDL_LONG));
  ev_p=(IDL_LONG *) malloc(nv*sizeof(IDL_LONG));
	retval=gvert(poly, vcirc, &tol, nv, 0, nve, (int *) &nv, (vec *) ve_p, 
               angle_p, (int *) ipv_p, (int *) gp_p, (int *) &nev, 
               (int *) &nev0, (int *) ev_p);
#if 0
  printf("%d\n",retval);
  for(i=0;i<nve;i++)
    printf("%e %e %e\n",ve_p[i*3+0],ve_p[i*3+1],ve_p[i*3+2]);
  fflush(stdout);
#endif
  retval=vmidc(poly, nv, nve, (vec *) ve_p, (int *) ipv_p, (int *) ev_p, 
              (int *) &nvm, (vec **) &vm_p);
  
#if 0
  printf("%d\n",retval);
  for(i=0;i<nvm;i++)
    printf("%e %e %e\n",vm_p[i*3+0],vm_p[i*3+1],vm_p[i*3+2]);
  fflush(stdout);
#endif
  for(i=0;i<3;i++)
    outvm[i]=vm_p[0*3+i];
	
	FREEVEC(poly);
	FREEVEC(vm_p);
	FREEVEC(ve_p);
	FREEVEC(angle_p);
	FREEVEC(ipv_p);
	FREEVEC(gp_p);
	FREEVEC(ev_p);
	free_memory();
	return retval;
}

/******************************************************************************/

