#include <stdio.h>
#include <stdlib.h>
#include "logical.h"
#include "polygon.h"

/* local functions */
int garea();

/* external functions */
extern void garea_();

/*------------------------------------------------------------------------------
  Area of polygon.

  This is a c interface to fortran subroutine garea.

   Input: poly is a polygon.
	  *verb = 0 to suppress messages from garea, even fatal ones;
		  1 to allow messages from garea.
  Input/Output: *tol = angle within which to merge multiple intersections.
  Output: *area = area of polygon.
  Return value: 0 if ok;
		1 if fatal intersection of boundaries;
		-1 if failed to allocate memory.
*/
int garea(poly, tol, verb, area)
int *verb;
double *tol, *area;
polygon *poly;
{
    logical ldegen;
    /* work arrays */
    int *iord;
    double *phi;

    /* allocate memory for work arrays */
    iord = (int *) malloc(sizeof(int) * poly->np * 2);
    if (!iord) {
	fprintf(stderr, "garea: failed to allocate memory for %d integers\n", poly->np * 2);
	return(-1);
    }
    phi = (double *) malloc(sizeof(double) * poly->np * 2);
    if (!phi) {
	fprintf(stderr, "garea: failed to allocate memory for %d doubles\n", poly->np * 2);
  free(iord);
	return(-1);
    }

    /* fortran routine */
    garea_(area, poly->rp, poly->cm, &poly->np, tol, verb, phi, iord, &ldegen);

    /* free work arrays */
    free(iord);
    free(phi);

    /* fatal intersection of boundaries */
    if (ldegen) return(1);

    return(0);
}
