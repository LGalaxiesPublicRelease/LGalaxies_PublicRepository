#include <stdio.h>
#include <stdlib.h>
#include "logical.h"
#include "polygon.h"

/* number of extra vertices to allocate, to allow for expansion */
#define DNV		4

/* local functions */
int gvlims();
int gvlim();

/* external functions */
extern void gvlim_();

/*------------------------------------------------------------------------------
  Points on polygon nearest to and farthest from unit direction vi,
  and minimum and maximum values of cm=1-cos(theta)
  between polygon and unit vector vi.

  This is a wrapper around gvlim,
  that calls gvlim until the arrays are large enough.

   Input: poly is a polygon.
	  vcirc = 1 to return vertices and midpoints also for bounding circles
		  which have no intersections;
		= 0 not so.
	  vi[3] = unit vector.
  Input/Output: *tol = angle within which to merge multiple intersections.
  Output: *nv = number of vertices.
	  *vmin_p, vmax_p = pointers to arrays vmin[nv][3], vmax[nv][3]
		  giving nearest and farthest points from unit vector vi
		  on each edge; memory for the arrays is allocated.
	  *cmmin_p, *cmmax_p = pointer to arrays cmmin[nv], cmmax[nv]
		  giving minimum and maximum cm=1-cos(theta)
		  between each edge and unit vector vi;
		  memory for the arrays is allocated.
	  *nev = number of connected sequences of vertices.
	  *nev0 = number of bounding circles which have no intersections.
	  *ev_p = pointer to array ev[nv] of end indices;
		  memory for the array is allocated.
  Return value: 0 if ok;
		1 if fatal degenerate intersection of boundaries;
		-1 if failed to allocate memory.
*/
int gvlims(poly, vcirc, vi, tol, nv, vmin_p, vmax_p, cmmin_p, cmmax_p, nev, nev0, ev_p)
polygon *poly;
int vcirc;
double vi[3];
double *tol;
int *nv, *nev, *nev0;
double **vmin_p, **vmax_p, **cmmin_p, **cmmax_p;
int **ev_p;
{
    static int nvmax = 0;
    static int *ev = 0x0;
    static double *vmin = 0x0, *vmax = 0x0, *cmmin = 0x0, *cmmax = 0x0;

    int ier;

    /* first guess at number of vertices */
    *nv = poly->np;
    
    /* keep trying till the arrays are big enough */
    do {

	/* make sure that allocated arrays contain enough space */
	if (!vmin || !vmax || !cmmin || !cmmax || !ev || *nv > nvmax) {
	    if (vmin) free(vmin);
	    if (vmax) free(vmax);
	    if (cmmin) free(cmmin);
	    if (cmmax) free(cmmax);
	    if (ev) free(ev);
	    vmin = (double *) malloc(sizeof(double) * (*nv + DNV) * 3);
	    if (!vmin) {
		fprintf(stderr, "gvlims: failed to allocate memory for %d vectors\n", *nv + DNV);
		return(-1);
	    }
	    vmax = (double *) malloc(sizeof(double) * (*nv + DNV) * 3);
	    if (!vmax) {
		fprintf(stderr, "gvlims: failed to allocate memory for %d vectors\n", *nv + DNV);
		return(-1);
	    }
	    cmmin = (double *) malloc(sizeof(double) * (*nv + DNV));
	    if (!cmmin) {
		fprintf(stderr, "gvlims: failed to allocate memory for %d doubles\n", *nv + DNV);
		return(-1);
	    }
	    cmmax = (double *) malloc(sizeof(double) * (*nv + DNV));
	    if (!cmmax) {
		fprintf(stderr, "gvlims: failed to allocate memory for %d doubles\n", *nv + DNV);
		return(-1);
	    }
	    ev = (int *) malloc(sizeof(int) * (*nv + DNV));
	    if (!ev) {
		fprintf(stderr, "gvlims: failed to allocate memory for %d integers\n", *nv + DNV);
		return(-1);
	    }
	    nvmax = *nv + DNV;
	}

	/* compute vertices of polygon */
	ier = gvlim(poly, vcirc, vi, nvmax, tol, nv, vmin, vmax, cmmin, cmmax, nev, nev0, ev);
	if (ier) return(ier);

    } while (*nv > nvmax);

    /* point arguments at allocated arrays */
    *vmin_p = vmin;
    *vmax_p = vmax;
    *cmmin_p = cmmin;
    *cmmax_p = cmmax;
    *ev_p = ev;

    return(0);
}

/*------------------------------------------------------------------------------
  Points on polygon nearest to and farthest from unit direction vi,
  and minimum and maximum values of cm=1-cos(theta)
  between polygon and unit vector vi.

  This is a c interface to fortran subroutine gvlim.

   Input: poly is a polygon.
	  vcirc = 1 to return vertices and midpoints also for bounding circles
		  which have no intersections;
		= 0 not so.
	  vi[3] = unit vector.
	  nvmax = dimension of v[nvmax][3] and ev[nmax].
  Input/Output: *tol = angle within which to merge multiple intersections.
  Output: *nv = number of vertices.
	  vmin[nv][3], vmax[nv][3] = arrays giving nearest and farthest
		    points from unit vector vi on each edge.
	  cmmin[nv], cmmax[nv] = minimum and maximum cm=1-cos(theta)
		    between each edge and unit vector vi.
	  *nev = number of connected sequences of vertices.
	  *nev0 = number of bounding circles which have no intersections.
	  ev[nev] = end indices of each connected sequence of vertices.
  Return value: 0 if ok;
		1 if fatal degenerate intersection of boundaries;
		-1 if failed to allocate memory.
*/
int gvlim(poly, vcirc, vi, nvmax, tol, nv, vmin, vmax, cmmin, cmmax, nev, nev0, ev)
polygon *poly;
double vi[3];
int nvmax;
double *tol;
int *nv, *nev, *nev0;
#ifdef GCC
double vmin[nvmax][3], vmax[nvmax][3];
double cmmin[nvmax], cmmax[nvmax];
int ev[nvmax];
#else
double vmin[][3], vmax[][3];
double cmmin[], cmmax[];
int ev[];
#endif
{
    logical ldegen;
    /* work arrays */
    int *iord, *iwk;
    double *phi, *wk;

    /* allocate memory for work arrays */
    iord = (int *) malloc(sizeof(int) * poly->np * 2);
    if (!iord) {
	fprintf(stderr, "gvlim: failed to allocate memory for %d integers\n", poly->np * 2);
	return(-1);
    }
    phi = (double *) malloc(sizeof(double) * poly->np * 2);
    if (!phi) {
	fprintf(stderr, "gvlim: failed to allocate memory for %d doubles\n", poly->np * 2);
	return(-1);
    }
    iwk = (int *) malloc(sizeof(int) * nvmax * 3);
    if (!iwk) {
	fprintf(stderr, "gvlim: failed to allocate memory for %d integers\n", nvmax * 3);
	return(-1);
    }
    wk = (double *) malloc(sizeof(double) * nvmax);
    if (!wk) {
	fprintf(stderr, "gvlim: failed to allocate memory for %d doubles\n", nvmax);
	return(-1);
    }

    /* fortran routine */
    gvlim_(vmin, vmax, cmmin, cmmax, ev, &nvmax, nv, nev, nev0, poly->rp, poly->cm, &poly->np, vi, &vcirc, tol, phi, iord, wk, iwk, &ldegen);

    /* free work arrays */
    free(iord);
    free(phi);
    free(iwk);
    free(wk);

    /* fatal intersection of boundaries */
    if (ldegen) return(1);

    return(0);
}
