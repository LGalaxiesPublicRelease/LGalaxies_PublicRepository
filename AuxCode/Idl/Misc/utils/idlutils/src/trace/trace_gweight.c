#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "evilmath.h"

IDL_LONG trace_gweight
  (int      argc,
   void *   argv[])
{
   IDL_LONG    nx;
   IDL_LONG    ny;
   float    ** fimage;
   float    ** invvar;
   float       sigma;
   IDL_LONG    ncen;
   float    *  xcen;
   IDL_LONG *  ycen;
   float    *  xerr;

   long        iy, i;
   long        icen;
   IDL_LONG    retval = 1;
   long        lower, upper;
   float       profile;
   float       denom, sum, frac, fact;
   float       mean, weight,meanvar;
   float       xtemp;
   long        ytemp;
   long        bad;


   /* Allocate pointers from IDL */
   nx = *((IDL_LONG *)argv[0]);
   ny = *((IDL_LONG *)argv[1]);
   fimage = (float **)malloc(ny * sizeof(float *)); /* build pointers only */
   invvar = (float **)malloc(ny * sizeof(float *)); /* build pointers only */
   for (iy=0; iy < ny; iy++) {
      fimage[iy] = (float *)argv[2] + iy*nx;
      invvar[iy] = (float *)argv[3] + iy*nx;
   }
   sigma = *((float *)argv[4]);
   ncen = *((IDL_LONG *)argv[5]);
   xcen = ((float *)argv[6]);
   ycen = ((IDL_LONG *)argv[7]);
   xerr = ((float *)argv[8]);

   denom = 1.0/sqrt(2.0*3.1416*sigma*sigma);

   /* Loop through each center value */
   for (icen=0; icen < ncen; icen ++) {
      xtemp = xcen[icen];
      ytemp = ycen[icen];
      xerr[icen] = 999.0;	/* error out peaks too close to the edge. */
      lower =  xtemp - 3.0*sigma;
      upper =  (long)(xtemp + 3.0*sigma) + 1;
      if (lower >= 0 && upper < nx) {

	   mean = 0.0;
           weight = 0.0;
           fact = 0.0;
	   meanvar = 0.0;
           bad = 0;
	   for (i = lower; i <= upper; i++) {

/*  Integrate Gaussian Profile by summing up 5 parts  */
	      if (invvar[ytemp][i] > 0.0) {
	        for (sum=0.0,frac=-0.4;frac<0.5;frac += 0.2) {
                     fact = (xtemp + frac - i)/sigma;
                     sum += exp(-0.5 * fact * fact);
                 }

/*  Now add to top and bottom  */

                 profile = 0.2 * sum * denom;
	         mean += profile * fimage[ytemp][i] * i;
                 weight += profile * fimage[ytemp][i];
	         meanvar += profile*profile * (i-xtemp) * (i-xtemp) / invvar[ytemp][i] ;
	      } else bad = 1;

	   }

	   if (weight > 0.0 && bad != 1) {
               xcen[icen] = mean/weight;
               xerr[icen] = sqrt(meanvar)/weight;
           } else {
               xerr[icen] = 999.0;
           }
      }
   }

   /* Free temporary memory */
   free(fimage);
   free(invvar);

   return retval;
}

