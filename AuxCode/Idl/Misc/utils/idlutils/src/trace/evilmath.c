
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "export.h"
#include "evilmath.h"

#define TINY 1.0e-20
 
/* Recenter one position in a row of data using a flux-weighted center */
void recenter_fweight
  (IDL_LONG    nx,
   float    *  imrow,
   float    *  invvar,
   float       radius,
   float       xinit,
   float    *  xcen,
   float    *  xerr)
{
   IDL_LONG    ii;
   IDL_LONG    ix1;
   IDL_LONG    ix2;
   IDL_LONG    npix;
   int         qbad;
   float       x1;
   float       x2;
   float       sumxw;
   float       sumsxsx;
   float       sumw;
   float       xdiff;
   float    *  convol;

   /* Only recenter if the guess value is within the bounds of the image */
   if (xinit > radius - 0.5 && xinit < nx - 0.5 - radius) {

      /* Determine which pixel numbers over which to sum */
      x1 = xinit - radius + 0.5;
      x2 = xinit + radius + 0.5;
      ix1 = floor(x1);
      ix2 = floor(x2);

      /* If either end of the convol is out of bounds in the image, then
         shrink both sides of the summing convol until it is in bounds.
       */
      if (x1 < 0.0) {
         x2 += x1;
         x1 = 0.0;
      }
      if (x2 > nx) {
         x1 += x2 - nx;
         x2 = nx;
         ix2 = nx - 1;
      }
      npix = ix2 - ix1 + 1;

      convol = (float *) malloc(sizeof(float) * npix);
      sumw = 0.0;
      sumxw = 0.0;
      qbad = 0;
      for (ii=0; ii < npix; ii++) {
         /* Determine the weight of a boxcar window function (convol) for this
          * pixel.  Note that the values of "convol" will sum to 2*radius
          * unless the edge of the image is reached.
          */
         if (ii == 0) {
            convol[ii] = 1.0 + ix1 - x1;
         } else if (ii == npix-1) {
            convol[ii] = x2 - ix2;
         } else {
            convol[ii] = 1.0;
         }
         sumw = sumw + convol[ii] * imrow[ix1+ii];

         xdiff = ix1 + ii - xinit;
         sumxw += convol[ii] * imrow[ix1+ii] * xdiff;

         if (invvar[ix1+ii] <= 0.0) qbad = 1;
      }

      if (sumw > 0.0 && qbad == 0) {
         *xcen = sumxw / sumw + xinit;

         /* Compute the error in the flux-weighted center */
         sumsxsx = 0.0;
         for (ii=0; ii < npix; ii++) {
            xdiff = ix1 + ii - *xcen;
            sumsxsx += xdiff * xdiff * convol[ii] * convol[ii]
             / (sumw * sumw * invvar[ix1+ii]);
         }
         *xerr = sqrt(sumsxsx);
      } else {
         *xcen = xinit;
         *xerr = 999.0;
      }

      free(convol);
   } else {
     *xcen = xinit;
     *xerr = 999.0;
   }

   /* Trap a centroid going out of bounds from the fitting region. */
   /* This will only happen if there are negative flux values. */
   if (*xcen > xinit + radius + 0.5 || *xcen < xinit - radius - 0.5) {
     *xcen = xinit;
     *xerr = 999.0;
   }
}

