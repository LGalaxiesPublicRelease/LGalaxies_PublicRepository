#include <stdlib.h>
#include <stdio.h>
#include "harmonics.h"
#include "logical.h"
#include "polygon.h"

/* local functions */
int gsphr(), gsphra();

/* external functions */
extern void gspher_(), gsphera_();

/*------------------------------------------------------------------------------
  Spherical harmonics of polygon.

  This is a simplified c interface to fortran subroutine gspher.
  It returns the spherical harmonics, and does not worry about bound and vert.

   Input: poly is a polygon.
	  lmax = maximum harmonic number.
  Input/Output: *tol = angle within which to merge multiple intersections.
  Output: w = array containing spherical harmonics of polygon.
  Return value: 0 if ok;
		1 if fatal intersection of boundaries;
		-1 if could not allocate temporary memory.
*/
int gsphr(poly, lmax, tol, w)
polygon *poly;
int lmax;
double *tol;
#ifdef GCC
double w[IM * NW];
#else
double w[];
#endif
{
    logical ldegen;
    int i, ibv, im, iphi, lmax1, npc, nw;
    double area, bound[2], vert[2];
    /* work arrays */
    int *iord;
    double *v, *phw;

    /* allocate memory for work arrays */
    iord = (int *) malloc(sizeof(int) * poly->np * 2);
    if (!iord) {
	fprintf(stderr, "gsphr: failed to allocate memory for %d integers\n", poly->np * 2);
	return(-1);
    }
    phw = (double *) malloc(sizeof(double) * poly->np * 2);
    if (!phw) {
	fprintf(stderr, "gsphr: failed to allocate memory for %d doubles\n", poly->np * 2);
	return(-1);
    }
    v = (double *) malloc(sizeof(double) * (lmax + 1));
    if (!v) {
	fprintf(stderr, "gsphr: failed to allocate memory for %d doubles\n", lmax + 1);
	return(-1);
    }

    /* parameters */
    lmax1 = lmax + 1;
    im = IM;
    nw = NW;
    npc = 0;
    ibv = 0;
    iphi = 0;

    /* set harmonics to zero */
    for (i = 0; i < IM * NW; i++) w[i] = 0.;

    /* the fortran routine */
    gspher_(&area, bound, vert, w, &lmax1, &im, &nw, poly->rp, poly->cm, &poly->np, &npc, &ibv, &iphi, tol, phw, iord, v, &ldegen);

    /* gspher_ may set area = 0 but not w = 0, because of numerical roundoff */
    if (area == 0.) {
	for (i = 0; i < IM * NW; i++) w[i] = 0.;
    }

    /* free work arrays */
    free(iord);
    free(phw);
    free(v);

    /* fatal intersection of boundaries */
    if (ldegen) return(1);

    return(0);
}

/*------------------------------------------------------------------------------
  Accelerated computation of spherical harmonics of rectangle.

  This is a simplified c interface to fortran subroutine gsphera.
  It returns the spherical harmonics, and does not worry about bound and vert.

  The acceleration involves some overhead, and works only if two or more
  rectangles with the same elmin & elmax are computed in succession.
  The overhead means that the accelerated computation is actually slightly
  slower for just a single rectangle.

   Input: lmax = maximum harmonic number.
  Output: w = array containing spherical harmonics of polygon.
  Return value: 0 if ok;
		-1 if could not allocate temporary memory.
*/
int gsphra(azmin, azmax, elmin, elmax, lmax, w)
double azmin, azmax, elmin, elmax;
int lmax;
#ifdef GCC
double w[IM * NW];
#else
double w[];
#endif
{
    /* array used for acceleration */
    static double *dw = 0x0;

    int i, ibv, im, lmax1, nw;
    double area, bound[2], vert[2];
    /* work array */
    double *v;

    /* allocate memory for work arrays */
    v = (double *) malloc(sizeof(double) * (lmax + 1));
    if (!v) {
	fprintf(stderr, "gsphra: failed to allocate memory for %d doubles\n", lmax + 1);
	return(-1);
    }

    /* parameters */
    lmax1 = lmax + 1;
    im = IM;
    nw = NW;
    ibv = 0;

    /* set harmonics to zero */
    for (i = 0; i < IM * NW; i++) w[i] = 0.;

    /* dw contains array that is pre-computed, then used by all rects with same elmin, elmax */
    if (!dw) {
	dw = (double *) malloc(sizeof(double) * NW);
	if (!dw) {
	    fprintf(stderr, "gsphra: failed to allocate memory for %d doubles\n", NW);
	    return(-1);
	}
    }

    /* fortran routine */
    gsphera_(&area, bound, vert, w, &lmax1, &im, &nw, &ibv, &azmin, &azmax, &elmin, &elmax, v, dw);

    /* free work array */
    free(v);

    return(0);
}
