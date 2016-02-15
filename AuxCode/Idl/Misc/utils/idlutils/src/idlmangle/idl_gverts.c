#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <polygon.h>
#include "export.h"

static int ve_int;
static int nve_int;
#if 0
static double *ve_p_int=NULL;
static double *angle_p_int=NULL;
static int *ipv_p_int=NULL;
static int *ev_p_int=NULL;
#endif

IDL_LONG gvert(polygon *poly, int vcirc, double *tol, int nvmax, int per, int nve, 
          int *nv, vec ve[/*nvmax * nve*/], double angle[/*nvmax*/], 
          int ipv[/*nvmax*/], int gp[/*poly->np*/], int *nev, int *nev0, 
          int ev[/*nvmax*/]);

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
#if 0
	FREEVEC(ve_p_int);
	FREEVEC(angle_p_int);
	FREEVEC(ipv_p_int);
	FREEVEC(ev_p_int);
#endif
}

#define DEG2RAD .01745329251994

/********************************************************************/
IDL_LONG idl_gverts
  (int      argc,
   void *   argv[])
{
	polygon *poly;
	double    tol;
	IDL_LONG vcirc;
	IDL_LONG *nv;
	IDL_LONG nve;
	vec *ve_p;
	double *angle_p;
	IDL_LONG *ipv_p;
	IDL_LONG *gp_p;
	IDL_LONG *nev;
	IDL_LONG *nev0;
	IDL_LONG *ev_p;
	
	IDL_LONG i;
	IDL_LONG retval=1;

	poly=(polygon *) malloc(sizeof(polygon));

	/* 0. allocate pointers from IDL */
	i=0;
	poly->rp = (vec *) argv[i]; i++;
	poly->cm = (double *)argv[i]; i++;
	poly->np = *((IDL_LONG *)argv[i]); i++;
	poly->npmax = poly->np;
	tol = *((double *)argv[i]); i++;
	vcirc = *((IDL_LONG *)argv[i]); i++;
	nv = ((IDL_LONG *)argv[i]); i++;
	nve_int =nve = *((IDL_LONG *)argv[i]); i++;
	ve_p = ((vec *)argv[i]); i++;
	angle_p = ((double *)argv[i]); i++;
	ipv_p = ((IDL_LONG *)argv[i]); i++;
	gp_p = ((IDL_LONG *)argv[i]); i++;
	nev = ((IDL_LONG *)argv[i]); i++;
	nev0 = ((IDL_LONG *)argv[i]); i++;
	ev_p = ((IDL_LONG *)argv[i]); i++;
	
	/* 1. we know sizes, so we can do the real thing */
	retval=gvert(poly, vcirc, &tol, *nv, 0, nve, (int *) nv, ve_p, angle_p, 
               (int *) ipv_p, (int *) gp_p, (int *) nev, (int *) nev0, 
               (int *) ev_p);
	
	FREEVEC(poly);
	free_memory();
	return retval;
}

/******************************************************************************/

