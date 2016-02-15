
#ifndef __INCevilmath_h
#define __INCevilmath_h

#include "export.h"

#define min(a,b) ( ((a) < (b)) ? (a) : (b) )
#define max(a,b) ( ((a) > (b)) ? (a) : (b) )

void recenter_fweight
  (IDL_LONG    nx,
   float    *  imrow,
   float    *  imerr,
   float       radius,
   float       xinit,
   float    *  xcen,
   float    *  xerr);

int fixedGauss2(float x[], float y[], float ymod[], float invvar[], float *xcen,
        int ndat, int nTrace, int nPoly, float *sigma,
        float a[], int ia[], int ma, float **covar, float *chisq,
        int *xmin, int *xmax, float **abig);
float **matrix_nr(int nrow, int ncol);
void chebyshevFunc(float x, float *coeff, int nCoeff);
void Profile(float *x, int ndat, float *y, float xcen, int xmin,
                int xmax, float sigma);


#endif /* __INCevilmath_h */

