#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ph.h"

/*
 * photfrac.c
 *
 * Yield fraction of flux which each pixel should 
 * contribute to an aperture of a given size (only
 * return one quadrant).
 *
 * Mike Blanton
 * 8/2003 */

#define PI 3.14159265358979

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
void p_cisi(float x, float *ci, float *si);

static float p_fi, p_fj;
static float p_radius;

float photfrac_func(float x) {
  float sincarg,y_int1,y_int2,dum_int,val;

  sincarg=PI*(p_fj+sqrt(p_radius*p_radius-x*x));
  p_cisi(sincarg, &dum_int, &y_int1);
  sincarg=PI*(p_fj-sqrt(p_radius*p_radius-x*x));
  p_cisi(sincarg, &dum_int, &y_int2);
  sincarg=PI*(p_fi-x);
  if(sincarg==0.) 
    val=(y_int1-y_int2)/PI;
  else 
    val=sin(sincarg)/(sincarg*PI)*(y_int1-y_int2);
  
  return(val);
} /* end photfrac_func */

int photfrac(int xnpix, 
             int ynpix,
             float radius, 
             float *frac, 
             long xcen, 
             long ycen)
{
  int i,j;

  p_radius=radius;

  for(i=0;i<xnpix;i++) {
    p_fi=(float) (i-xcen);
    for(j=0;j<ynpix;j++) {
      p_fj=(float) (j-ycen);
      frac[j*xnpix+i]= 
        p_qromo(photfrac_func, -radius, radius, p_midpnt);
    }
  } /* end for i,j */

	return(1);
} /* end photfrac */
