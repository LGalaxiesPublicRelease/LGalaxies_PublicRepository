#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "nr.h"

float vector_avsigclip_mask
  (IDL_LONG    nData,
   float    *  pData,
   char     *  pMaskIn,
   char     *  pMaskOut,
   float       sigrejlo,
   float       sigrejhi,
   IDL_LONG    maxiter);
void vector_mean_and_disp_mask
  (IDL_LONG    nData,
   float    *  pData,
   char     *  pMask,
   float    *  pMean,
   float    *  pDispersion);
float vector_mean_mask
  (IDL_LONG    nData,
   float    *  pData,
   char     *  pMask);

/******************************************************************************/
IDL_LONG arravsigclip_mask
  (int         argc,
   void    *   argv[])
{
   IDL_LONG    ndim;
   IDL_LONG *  dimvec;
   float    *  array;
   IDL_LONG    dim;
   float       sigrejlo;
   float       sigrejhi;
   IDL_LONG    maxiter;
   float    *  avearr;
   char     *  maskin;
   char     *  maskout;

   IDL_LONG    i;
   IDL_LONG    nlo;
   IDL_LONG    nmid;
   IDL_LONG    nhi;
   IDL_LONG    ilo;
   IDL_LONG    imid;
   IDL_LONG    ihi;
   IDL_LONG    i1;
   IDL_LONG    indx;
   float    *  tempvec;
   char     *  tempm1;
   char     *  tempm2;
   IDL_LONG    retval = 1;

   /* Allocate pointers from IDL */
   ndim = *((IDL_LONG *)argv[0]);
   dimvec = (IDL_LONG *)argv[1];
   array = (float *)argv[2];
   dim = *((IDL_LONG *)argv[3]);
   sigrejlo = *((float *)argv[4]);
   sigrejhi = *((float *)argv[5]);
   maxiter = *((IDL_LONG *)argv[6]);
   avearr = (float *)argv[7];
   maskin = (char *)argv[8];
   maskout = (char *)argv[9];

   nlo = 1;
   for (i=0; i < dim-1; i++) nlo *= dimvec[i];
   nhi = 1;
   for (i=dim; i < ndim; i++) nhi *= dimvec[i];
   nmid = dimvec[dim-1];

   /* Allocate memory for temporary vector */
   tempvec = malloc(nmid * sizeof(float));
   tempm1 = malloc(nmid * sizeof(char));
   tempm2 = malloc(nmid * sizeof(char));

   /* Loop through all pixels in array */
   for (ilo=0; ilo < nlo; ilo++) {
      for (ihi=0; ihi < nhi; ihi++) {
         i1 = ilo + ihi * nmid * nlo;

         /* Construct the vector of values from which to compute */
         for (imid=0; imid < nmid; imid++) {
            indx = i1 + imid * nlo;
            tempvec[imid] = array[indx];
            tempm1[imid] = maskin[indx];
            tempm2[imid] = maskout[indx];
         }

         /* Compute a single average value with sigma clipping */
         avearr[ilo + ihi * nlo] =
          vector_avsigclip_mask(nmid, tempvec, tempm1, tempm2,
           sigrejlo, sigrejhi, maxiter);

         /* Copy output mask values */
         for (imid=0; imid < nmid; imid++) {
            indx = i1 + imid * nlo;
            maskout[indx] = tempm2[imid];
         }
      }
   }

   /* Free temporary memory */
   free(tempvec);
   free(tempm1);
   free(tempm2);

   return retval;
}

/******************************************************************************/
float vector_avsigclip_mask
  (IDL_LONG   nData,
   float    * pData,
   char     * pMaskIn,
   char     * pMaskOut,
   float      sigrejlo,
   float      sigrejhi,
   IDL_LONG   maxiter)
{
   IDL_LONG   i;
   IDL_LONG   nmask;
   IDL_LONG   nbad;
   IDL_LONG   iiter;
   float      mval;
   float      mdisp;

   /* Copy the input mask into the output mask */
   nmask = 0;
   for (i=0; i<nData; i++) {
      pMaskOut[i] = pMaskIn[i];
      if (pMaskIn[i] > 0) nmask++;
   }

   /* First compute the mean and dispersion */
   if (maxiter > 0) {
      vector_mean_and_disp_mask(nData, pData, pMaskOut, &mval, &mdisp);
   } else {
      mval = vector_mean_mask(nData, pData, pMaskOut);
   }

   /* Iterate with rejection */
   for (iiter=0; iiter < maxiter; iiter++) {

      /* Copy the input mask into the output mask */
      for (i=0; i<nData; i++) pMaskOut[i] = pMaskIn[i];

      /* Reject outliers, but only if the dispersion is greater than 0.
       * Otherwise, we could get down to only 1 good data point, with
       * a dispersion of zero, and then reject that point!
       */
      nbad = nmask;
      for (i=0; i < nData; i++) {
         if (pMaskIn[i] == 0 && mdisp > 0 &&
          (pData[i] <= mval - sigrejlo*mdisp
           || pData[i] >= mval + sigrejhi*mdisp) ) {
            pMaskOut[i] = 1;
            nbad++;
         }
      }

      /* If all points are now rejected, then un-reject the points just
       * masked out, such that we only mask those points specified by pMaskIn.
       * This is no guarantee that any points survive, since pMaskIn might
       * already mask all points.
       */
#if 0
      if (nbad == nData) {
         for (i=0; i < nData; i++) {
            if (pMaskIn[i] == 0 && pMaskOut[i] == 1) {
               pMaskOut[i] = 0;
               nbad--;
            }
         }
      }
#endif

      if (nbad == nData) {
         iiter = maxiter;
      } else {
         vector_mean_and_disp_mask(nData, pData, pMaskOut, &mval, &mdisp);
      }
   }

   return mval;
}

/******************************************************************************/
/* Find the mean and dispersion of the nData elements.
 * Return a dispersion of zero if there are fewer than two data elements.
 */
void vector_mean_and_disp_mask
  (IDL_LONG    nData,
   float    *  pData,
   char     *  pMask,
   float    *  pMean,
   float    *  pDispersion)
{
   IDL_LONG    i;
   IDL_LONG    ngood;
   float       vmean;
   float       vdisp;
   float       vtemp;

   ngood = 0;
   vmean = 0.0;
   for (i=0; i<nData; i++) {
      if (pMask[i] == 0) {
         vmean += pData[i];
         ngood++;
      }
   }

   if (ngood > 0) {
      vmean /= ngood;
      vdisp = 0.0;
      if (ngood > 1) {
         for (i=0; i<nData; i++) {
            if (pMask[i] == 0) {
               vtemp = pData[i] - vmean;
               vdisp += vtemp * vtemp;
            }
         }
         vdisp = sqrt(vdisp / (ngood-1));
      } else {
         vdisp = 0.0;
      }
   } else {
      /* Compute the mean and dispersion for all points
       * (even though all are masked)
       */
      for (i=0; i<nData; i++) vmean += pData[i];
      vmean /= nData;
      vdisp = 0.0;
      for (i=0; i<nData; i++) {
         vtemp = pData[i] - vmean;
         vdisp += vtemp * vtemp;
      }
      if (nData > 1) vdisp = sqrt(vdisp / (nData-1));
      else vdisp = 0.0;
   }

   *pMean = vmean;
   *pDispersion = vdisp;
}

/******************************************************************************/
/* Find the mean value of the nData elements of a floating point array pData[].
 */
float vector_mean_mask
  (IDL_LONG    nData,
   float    *  pData,
   char     *  pMask)
{
   float vmean;
   IDL_LONG i;
   IDL_LONG ngood;

   vmean = 0.0;
   ngood = 0;
   for (i=0; i<nData; i++) {
      if (pMask[i] == 0) {
         vmean += pData[i];
         ngood++;
      }
   }
   if (ngood > 0) vmean /= ngood;

   return vmean;
}

/******************************************************************************/
