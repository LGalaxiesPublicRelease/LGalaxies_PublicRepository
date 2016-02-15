#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "evilmath.h"

IDL_LONG trace_fweight
  (int      argc,
   void *   argv[])
{
   IDL_LONG    nx;
   IDL_LONG    ny;
   float    ** fimage;
   float    ** invvar;
   float       radius;
   IDL_LONG    ncen;
   float    *  xcen;
   IDL_LONG *  ycen;
   float    *  xerr;

   long        iy;
   long        icen;
   IDL_LONG    retval = 1;

   /* Allocate pointers from IDL */
   nx = *((IDL_LONG *)argv[0]);
   ny = *((IDL_LONG *)argv[1]);
   fimage = (float **)malloc(ny * sizeof(float *)); /* build pointers only */
   invvar = (float **)malloc(ny * sizeof(float *)); /* build pointers only */
   for (iy=0; iy < ny; iy++) {
      fimage[iy] = (float *)argv[2] + iy*nx;
      invvar[iy] = (float *)argv[3] + iy*nx;
   }
   radius = *((float *)argv[4]);
   ncen = *((IDL_LONG *)argv[5]);
   xcen = ((float *)argv[6]);
   ycen = ((IDL_LONG *)argv[7]);
   xerr = ((float *)argv[8]);

   /* Loop through each center value */
   for (icen=0; icen < ncen; icen ++) {
      recenter_fweight(nx, fimage[ycen[icen]], invvar[ycen[icen]], radius,
       xcen[icen], &xcen[icen], &xerr[icen]);
   }

   /* Free temporary memory */
   free(fimage);
   free(invvar);

   return retval;
}

