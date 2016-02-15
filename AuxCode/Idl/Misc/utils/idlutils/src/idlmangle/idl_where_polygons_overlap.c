#include <stdio.h>
#include <math.h>
#include <manglefn.h>
#include <stdlib.h>
#include <polygon.h>
#include "export.h"

polygon *new_poly(int npmax);

/*---------------------------------------------------------------------------
  Suppress obviously superfluous caps from polygon,
  by setting cm = 2 for second of two coincident caps.
  In addition, detect null caps (those which contain nothing),
  and complementary caps (those which exclude each other,
  in which case replace the polygon with a single null cap.

  Two caps are considered coincident if their axes (rp) and latitudes (cm)
  are EXACTLY equal.  Caps with axes pointing in opposite directions are
  not detected.

  Two caps are considered complementary if their axes (rp) are EXACTLY equal,
  and their latitudes (cm) are EXACTLY opposing (cm of one is -cm of the other).
  Axes pointing in opposite directions are not detected.

  All this makes garea, gspher et al happy, provided that near coincident caps,
  including those with axes pointing in opposite directions, have been
  modified by `snap' to coincide exactly, with coaligned axes.

  Input: poly is a pointer to a polygon.
  Output: poly with obviously superfluous caps suppressed;
  the number and order of caps remains unchanged
  UNLESS polygon is replaced by null polygon.
  Return value: 0 if nothing changed;
  1 if one or more caps were suppressed;
  2 if nothing changed, and polygon is null;
  3 if polygon was changed to null polygon.

  MODIFIED FROM trim_poly
*/
int trim_poly_ez(polygon *poly)
{
  int ip, iret, jp;
  double tol;

  /* initialize return value to no change */
  iret = 0;
  tol=1.e-10;

  /* check for cap which excludes everything */
  for (jp = 0; jp < poly->np; jp++) {
    if (poly->cm[jp] == 0. || poly->cm[jp] <= -2.) {
	    if (poly->np == 1		/* polygon is already single null cap */
          && poly->rp[0][0] == 0.
          && poly->rp[0][1] == 0.
          && poly->rp[0][2] == 1.
          && poly->cm[0] == 0.) {
        return(2);
	    } else {			/* change polygon to single null cap */
        poly->rp[0][0] = 0.;
        poly->rp[0][1] = 0.;
        poly->rp[0][2] = 1.;
        poly->cm[0] = 0.;
        poly->np = 1;
        return(3);
	    }
    }
  }

  /* for each cap jp, check for coincident caps */
  for (jp = 0; jp < poly->np; jp++) {
    /* don't check superfluous cap */
    if (poly->cm[jp] >= 2.) continue;
    for (ip = jp+1; ip < poly->np; ip++) {
	    /* don't check superfluous cap */
	    if (poly->cm[ip] >= 2.) continue;
	    /* cap axes coincide */
	    if (fabs(poly->rp[ip][0]-poly->rp[jp][0])<tol
          && fabs(poly->rp[ip][1]-poly->rp[jp][1])<tol
          && fabs(poly->rp[ip][2]-poly->rp[jp][2])<tol) {
        /* suppress coincident cap ip */
        if (fabs(poly->cm[ip]-poly->cm[jp])<tol) {
          poly->cm[ip] = 2.;
          iret = 1;
        } else if (fabs(poly->cm[ip]+poly->cm[jp])<tol) {
          /* complementary cap means polygon is null */
          poly->rp[0][0] = 0.;
          poly->rp[0][1] = 0.;
          poly->rp[0][2] = 1.;
          poly->cm[0] = 0.;
          poly->np = 1;
          return(3);
        }
	    }
    }
  }

  return(iret);
}

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

#define DEG2RAD .01745329251994

/********************************************************************/
IDL_LONG idl_where_polygons_overlap
(int      argc,
 void *   argv[])
{
	polygon *poly;
  IDL_LONG ncaps,maxncaps, nmatch, *matchncaps, *ismatch;
  double *x, *cm, *xmatch, *cmmatch, *areamatch;
	
	IDL_LONG i,j,k;
	IDL_LONG retval=1;
	double    tol;
	IDL_LONG  verbose;

  tol=0.;
  verbose=0;

	/* 0. allocate pointers from IDL */
	i=0;
  x=(double *) argv[i]; i++;
  cm=(double *) argv[i]; i++;
  ncaps=*((IDL_LONG *) argv[i]); i++;
  xmatch=(double *) argv[i]; i++;
  cmmatch=(double *) argv[i]; i++;
  maxncaps=*((IDL_LONG *) argv[i]); i++;
  nmatch=*((IDL_LONG *) argv[i]); i++;
  matchncaps=((IDL_LONG *) argv[i]); i++;
  ismatch=((IDL_LONG *) argv[i]); i++;
  areamatch=((double *) argv[i]); i++;

	poly=new_poly(ncaps+maxncaps);

  for(i=0;i<nmatch;i++) {
    poly->np=ncaps+matchncaps[i];
    for(j=0;j<ncaps;j++) {
      poly->cm[j]=cm[j];
      for(k=0;k<3;k++) 
        poly->rp[j][k]=x[j*3+k];
    }
    for(j=0;j<matchncaps[i];j++) {
      poly->cm[ncaps+j]=cmmatch[i*maxncaps+j];
      for(k=0;k<3;k++) 
        poly->rp[ncaps+j][k]=xmatch[i*maxncaps*3+j*3+k];
    }
    trim_poly_ez(poly);
    retval=garea(poly, &tol, verbose, &(areamatch[i]));
    if(areamatch[i]>0.) 
      ismatch[i]=1; 
    else 
      ismatch[i]=0; 
  }
  
  free_poly(poly);
	return retval;
}

/******************************************************************************/

