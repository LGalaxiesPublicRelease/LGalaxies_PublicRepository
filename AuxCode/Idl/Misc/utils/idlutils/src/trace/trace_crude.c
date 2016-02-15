#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "evilmath.h"

IDL_LONG trace_crude
  (int      argc,
   void *   argv[])
{
   IDL_LONG    nx;
   IDL_LONG    ny;
   float    ** fimage;
   float    ** invvar;
   float       radius;
   IDL_LONG    ntrace;
   float    *  xstart;
   IDL_LONG *  ystart;
   float    ** xset;
   float    ** xerr;
   float    *  maxerr;
   float    *  maxshift;
   float    *  maxshift0;

   IDL_LONG    iy;
   IDL_LONG    itrace;
   IDL_LONG    retval = 1;
   float       xinit;
   float       xshift;
   float       xfit;
   float       xfiterr;

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
   ntrace = *((IDL_LONG *)argv[5]);
   xstart = ((float *)argv[6]);
   ystart = ((IDL_LONG *)argv[7]);
   xset = (float **)malloc(ntrace * sizeof(float *)); /* build pointers only */
   xerr = (float **)malloc(ntrace * sizeof(float *)); /* build pointers only */
   for (itrace=0; itrace < ntrace; itrace++) {
      xset[itrace] = (float *)argv[8] + itrace*ny;
      xerr[itrace] = (float *)argv[9] + itrace*ny;
   }
   maxerr = ((float *)argv[10]);
   maxshift = ((float *)argv[11]);
   maxshift0 = ((float *)argv[12]);

   /* Loop through each trace */
   for (itrace=0; itrace < ntrace; itrace ++) {

      /* RECENTER INITIAL ROW */
      iy = ystart[itrace];
      xinit = xstart[itrace];
      recenter_fweight(nx, fimage[iy], invvar[iy], radius, xinit,
       &xfit, &xfiterr);
      if (xfiterr < *maxerr && xfiterr != 0) {
         xshift = xfit - xinit;
         xshift = max(xshift, -*maxshift0);
         xshift = min(xshift, *maxshift0);
      } else {
         xshift = 0.0;
      }
      xset[itrace][iy] = xinit + xshift;
      xerr[itrace][iy] = xfiterr;

      /* LOOP FROM INITIAL (COL,ROW) NUMBER TO LARGER ROW NUMBERS */
      for (iy=ystart[itrace]+1; iy < ny; iy++) {
         xinit = xset[itrace][iy-1];
         recenter_fweight(nx, fimage[iy], invvar[iy], radius, xinit,
          &xfit, &xfiterr);
         if (xfiterr < *maxerr && xfiterr != 0) {
            xshift = xfit - xinit;
            xshift = max(xshift, -*maxshift);
            xshift = min(xshift, *maxshift);
         } else {
            xshift = 0.0;
         }
         xset[itrace][iy] = xinit + xshift;
         xerr[itrace][iy] = xfiterr;
      }

      /* LOOP FROM INITIAL (COL,ROW) NUMBER TO SMALLER ROW NUMBERS */
      for (iy=ystart[itrace]-1; iy >= 0; iy--) {
         xinit = xset[itrace][iy+1];
         recenter_fweight(nx, fimage[iy], invvar[iy], radius, xinit,
          &xfit, &xfiterr);
         if (xfiterr < *maxerr && xfiterr != 0) {
            xshift = xfit - xinit;
            xshift = max(xshift, -*maxshift);
            xshift = min(xshift, *maxshift);
         } else {
            xshift = 0.0;
         }
         xset[itrace][iy] = xinit + xshift;
         xerr[itrace][iy] = xfiterr;
      }

   }

   /* Free temporary memory */
   free(fimage);
   free(invvar);
   free(xset);
   free(xerr);

   return retval;
}

