#include <stdio.h>
#include <stdlib.h>
#include "polygon.h"

/* local functions */
int gvphi();

/* external functions */
extern void gvphi_();

/*------------------------------------------------------------------------------
  Point at the centre of that part of a boundary segment which
  (a) lies inside a polygon, and
  (b) contains or is closest to a specified unit vector.

  This is a c interface to fortran subroutine gvphi.

   Input: poly is a polygon.
	  rp[3], cm define a boundary.
	  vi[3] is a unit vector.
	  *tol = angle within which to merge multiple intersections.
  Output: angle = length of that part of boundary segment which
		  (a) lies inside the polygon, and
		  (b) contains, or is closest to, unit vector vi.
		= 0. if the boundary lies entirely outside the polygon.
	  v[3] = unit vector at centre of said boundary segment.
  Return value: 0 if ok;
		-1 if failed to allocate memory.
*/
int gvphi(poly, rp, cm, vi, tol, angle, v)
polygon *poly;
double rp[3], cm, vi[3];
double *tol;
double *angle;
double v[3];
{
    /* work arrays */
    int *iord;
    double *phi;

    /* allocate memory for work arrays */
    iord = (int *) malloc(sizeof(int) * poly->np * 2);
    if (!iord) {
	fprintf(stderr, "gvphi: failed to allocate memory for %d integers\n", poly->np * 2);
	return(-1);
    }
    phi = (double *) malloc(sizeof(double) * poly->np * 2);
    if (!phi) {
	fprintf(stderr, "gvphi: failed to allocate memory for %d doubles\n", poly->np * 2);
	return(-1);
    }

    /* fortran routine */
    gvphi_(angle, v, poly->rp, poly->cm, &poly->np, rp, &cm, vi, tol, phi, iord);

    /* free work arrays */
    free(iord);
    free(phi);

    return(0);
}
