/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Minimum and maximum values of cm = 1-cos(th) between polygon
  and a unit vector rp.

  This is a c interface to fortran subroutine gcmlim.

   Input: poly is a polygon.
	  *tol = angle within which to merge multiple intersections.
	  rp = unit vector.
  Output: minimum and maximum values of cm = 1-cos(th).
  Return value:  0 if ok;
		-1 if failed to allocate memory.
*/
int gcmlim(polygon *poly, double *tol, vec rp, double *cmmin, double *cmmax)
{
    /* work arrays */
    int *iord;
    double *phi;

    /* allocate memory for work arrays */
    iord = (int *) malloc(sizeof(int) * poly->np * 2);
    if (!iord) {
	fprintf(stderr, "gcmlim: failed to allocate memory for %d ints\n", poly->np * 2);
	return(-1);
    }
    phi = (double *) malloc(sizeof(double) * poly->np * 2);
    if (!phi) {
	fprintf(stderr, "gcmlim: failed to allocate memory for %d doubles\n", poly->np * 2);
	return(-1);
    }

    /* fortran routine */
    gcmlim_(poly->rp, poly->cm, &poly->np, rp, cmmin, cmmax, tol, phi, iord);

    /* free work arrays */
    free(iord);
    free(phi);

    return(0);
}
