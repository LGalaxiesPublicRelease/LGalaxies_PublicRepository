#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "nr.h"

/******************************************************************************/
IDL_LONG ccorrelate
  (int      argc,
   void *   argv[])
{
   IDL_LONG    nx;
   float    *  xvec;
   float    *  xweight;
   IDL_LONG    ny;
   float    *  yvec;
   float    *  yweight;
   IDL_LONG    nlag;
   IDL_LONG *  lags;
   float    *  result;

   IDL_LONG    ilag;
   IDL_LONG    startx;
   IDL_LONG    starty;
   IDL_LONG    ii;
   IDL_LONG    ncomp;
   IDL_LONG    retval = 1;
   float    *  wvec;
   double      wx;
   double      wy;
   double      sumw;
   double      sumwx;
   double      sumwy;
   double      sumwxwx;
   double      sumwywy;
   double      meanx;
   double      meany;
   float       res;

   /* Allocate pointers from IDL */
   nx = *((IDL_LONG *)argv[0]);
   xvec = (float *)argv[1];
   xweight = (float *)argv[2];
   ny = *((IDL_LONG *)argv[3]);
   yvec = (float *)argv[4];
   yweight = (float *)argv[5];
   nlag = *((IDL_LONG *)argv[6]);
   lags = (IDL_LONG *)argv[7];
   result = (float *)argv[8];

   /* Allocate memory for temporary vector */
   wvec = malloc((nx > ny ? nx : ny)  * sizeof(float));

   for (ilag=0; ilag < nlag; ilag++) {

      if (lags[ilag] < 0) {
         startx = -lags[ilag];
         starty = 0;
      } else {
         startx = 0;
         starty = lags[ilag];
      }

      ncomp = (nx - startx - 1) < (ny - starty - 1) ?
       (nx - startx - 1) : (ny - starty - 1);

      for (ii=0; ii < ncomp; ii++)
       wvec[ii] = xweight[startx+ii] * yweight[starty+ii];

      sumw = 0.0;
      sumwx = 0.0;
      sumwy = 0.0;
      sumwxwx = 0.0;
      sumwywy = 0.0;
      for (ii=0; ii < ncomp; ii++) {
         sumw += wvec[ii];
         wx = wvec[ii] * xvec[startx+ii];
         wy = wvec[ii] * yvec[starty+ii];
         sumwx += wx;
         sumwy += wy;
         sumwxwx += wx * wx;
         sumwywy += wy * wy;
      }

      meanx = sumwx / sumw;
      meany = sumwy / sumw;

      res = 0.0;
      for (ii=0; ii < ncomp; ii++) {
         res += wvec[ii] * wvec[ii] * (xvec[startx+ii] - meanx)
          * (yvec[starty+ii] - meany);
      }
      result[ilag] = res / sqrt( (sumwxwx - sumwx * sumwx / sumw)
                          * (sumwywy - sumwy * sumwy / sumw) );
   }

   /* Free temporary memory */
   free(wvec);

   return retval;
}

/******************************************************************************/

