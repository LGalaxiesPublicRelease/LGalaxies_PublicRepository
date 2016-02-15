#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <polygon.h>
#include "export.h"

IDL_LONG garea(polygon *poly, double *tol, IDL_LONG verb, double *area); 

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

#define DEG2RAD .01745329251994

/********************************************************************/
IDL_LONG idl_garea
  (int      argc,
   void *   argv[])
{
	polygon *poly;
	double    tol;
	IDL_LONG  verbose;
	double    *area;
	
	IDL_LONG i;
	IDL_LONG retval=1;

	poly=(polygon *) malloc(sizeof(polygon));

	/* 0. allocate pointers from IDL */
	i=0;
	poly->rp = (vec *)argv[i]; i++;
	poly->cm = (double *)argv[i]; i++;
	poly->np = *((IDL_LONG *)argv[i]); i++;
	poly->npmax = poly->np;
	tol = *((double *)argv[i]); i++;
	verbose = *((IDL_LONG *)argv[i]); i++;
	area = (double *)argv[i]; i++;
	
	/* 1. call garea */
	retval=garea(poly, &tol, verbose, area);
	FREEVEC(poly);
	return retval;
}

/******************************************************************************/

