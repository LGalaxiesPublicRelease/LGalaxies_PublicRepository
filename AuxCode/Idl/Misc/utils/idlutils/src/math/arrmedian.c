#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "nr.h"

float vector_median
  (IDL_LONG   nData,
   float *    pData);

/******************************************************************************/
IDL_LONG arrmedian
  (int      argc,
   void *   argv[])
{
   IDL_LONG    ndim;
   IDL_LONG *  dimvec;
   float    *  array;
   IDL_LONG    dim;
   float    *  medarr;

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
   medarr = (float *)argv[4];

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
         /* Construct the vector of values from which to compute a median */
         for (imid=0; imid < nmid; imid++) {
            indx = i1 + imid * nlo;
            tempvec[imid] = array[indx];
         }
         /* Compute a single median value */
         medarr[ilo + ihi * nlo] = vector_median(nmid, tempvec);
      }
   }

   /* Free temporary memory */
   free(tempvec);

   return retval;
}

/******************************************************************************/
/* Find the median of the "nData" elements of a floating point array pData[].
 * The data vector is returned unchanged.
 */
float vector_median
  (IDL_LONG   nData,
   float *    pData)
{
   unsigned long nlong;
   unsigned long klong;
   float    retval;

   /* Numerical Recipes routines expect everything to be 1-indexed.
    * Pass the array pData as 1-indexed, and the element counter klong
    * as a 1-indexed element number.
    */
   nlong = nData;
   klong = (nlong+1)/2;
   retval = selip(klong, nlong, pData-1);
   if (nData % 2 == 0) {
      retval = 0.5 * (retval + selip(klong+1, nlong, pData-1));
   }

   return retval;
}

