#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "nr.h"

float vector_avsigclip
  (IDL_LONG    nData,
   float    *  pData,
   float       sigrejlo,
   float       sigrejhi,
   IDL_LONG    maxiter);
void vector_mean_and_dispersion
  (IDL_LONG    nData,
   float    *  pData,
   float    *  pMean,
   float    *  pDispersion);
float vector_mean
  (IDL_LONG    nData,
   float    *  pData);
float vector_sum
  (IDL_LONG    nData,
   float    *  pData);

/******************************************************************************/
IDL_LONG arravsigclip
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

   nlo = 1;
   for (i=0; i < dim-1; i++) nlo *= dimvec[i];
   nhi = 1;
   for (i=dim; i < ndim; i++) nhi *= dimvec[i];
   nmid = dimvec[dim-1];

   /* Allocate memory for temporary vector */
   tempvec = malloc(nmid * sizeof(float));

   /* Loop through all pixels in array */
   for (ilo=0; ilo < nlo; ilo++) {
      for (ihi=0; ihi < nhi; ihi++) {
         i1 = ilo + ihi * nmid * nlo;
         /* Construct the vector of values from which to compute */
         for (imid=0; imid < nmid; imid++) {
            indx = i1 + imid * nlo;
            tempvec[imid] = array[indx];
         }
         /* Compute a single average value with sigma clipping */
         avearr[ilo + ihi * nlo] =
          vector_avsigclip(nmid, tempvec, sigrejlo, sigrejhi, maxiter);
      }
   }

   /* Free temporary memory */
   free(tempvec);

   return retval;
}

/******************************************************************************/
float vector_avsigclip
  (IDL_LONG    nData,
   float    *  pData,
   float       sigrejlo,
   float       sigrejhi,
   IDL_LONG    maxiter)
{
   IDL_LONG   i;
   IDL_LONG   nGood;
   IDL_LONG   iiter;
   float      mval;
   float      mdisp;
   float    * pGood;

   /* First compute the mean and dispersion */
   if (maxiter > 0) {
      /* Allocate memory for temporary vector */
      pGood = malloc(nData * sizeof(float));

      vector_mean_and_dispersion(nData, pData, &mval, &mdisp);
   } else {
      mval = vector_mean(nData, pData);
   }

   /* Iterate with rejection */
   for (iiter=0; iiter < maxiter; iiter++) {
      /* Copy all values within rejection limits into temporary vector */
      nGood = 0;
      for (i=0; i < nData; i++) {
         if (pData[i] > mval - sigrejlo*mdisp
          && pData[i] < mval + sigrejhi*mdisp)
          pGood[nGood++] = pData[i];
      }
      if (nGood == 0) {
         iiter = maxiter;
      } else {
         vector_mean_and_dispersion(nGood, pGood, &mval, &mdisp);
      }
   }

   /* Free memory */
   if (maxiter > 0) free(pGood);

   return mval;
}

/******************************************************************************/
/* Find the mean and dispersion of the nData elements.
 * Return a dispersion of zero if there are fewer than two data elements.
 */
void vector_mean_and_dispersion
  (IDL_LONG    nData,
   float    *  pData,
   float    *  pMean,
   float    *  pDispersion)
{
   IDL_LONG    i;
   float       vmean;
   float       vdisp;
   float       vtemp;

   vmean = vector_mean(nData, pData);
   vdisp = 0.0;
   if (nData > 1) {
      for (i=0; i<nData; i++) {
         vtemp = pData[i] - vmean;
         vdisp += vtemp * vtemp;
      }
      vdisp = sqrt(vdisp/(nData-1));
   }

   *pMean = vmean;
   *pDispersion = vdisp;
}

/******************************************************************************/
/* Find the mean value of the nData elements of a floating point array pData[].
 */
float vector_mean
  (IDL_LONG    nData,
   float    *  pData)
{
   float       vmean;

   vmean = vector_sum(nData, pData);
   vmean /=  nData;

   return vmean;
}

/******************************************************************************/
/* Find the summed value of the nData elements of a floating point array pData[].
 */
float vector_sum
  (IDL_LONG    nData,
   float    *  pData)
{
   IDL_LONG    i;
   float       vsum;

   vsum = 0.0;
   for (i=0; i<nData; i++) vsum += pData[i];

   return vsum;
}

