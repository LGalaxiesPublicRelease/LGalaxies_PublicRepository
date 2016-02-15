#include <stdio.h>
#include <stdlib.h>
#include "logical.h"
#include "polygon.h"

/* number of extra vertices to allocate, to allow for expansion */
#define DNV		4

/* local functions */
int gverts();
int gvert();

/* external functions */
extern void gvert_();

/*------------------------------------------------------------------------------
  Points on edges of polygon.

  This is a wrapper around gvert,
  that calls gvert until the arrays are large enough.

   Input: poly is a polygon.
	  vcirc = 1 to return vertices and midpoints also for bounding circles
		  which have no intersections;
		= 0 not so.
	  nve = desired number of points on each edge, including vertex.
  Input/Output: *tol = angle within which to merge multiple intersections.
  Output: *nv = number of vertices.
	  *ve_p = pointer to array ve[nv][nve][3] of points on edges;
		  memory for the array is allocated.
	  *angle_p = pointer to array angle[nv] of lengths of edges;
		  memory for the array is allocated.
	  *ipv_p = pointer to array ipv[nv] containing index of i'th edge;
		   that is, ve[i] lie on circle number ipv;
		   memory for the array is allocated.
	  *nev = number of connected sequences of vertices.
	  *nev0 = number of bounding circles which have no intersections.
	  *ev_p = pointer to array ev[nv] of end indices;
		  memory for the array is allocated.
  Return value: 0 if ok;
		1 if fatal degenerate intersection of boundaries;
		-1 if could not allocate memory.
*/
int gverts(poly, vcirc, tol, nv, nve, ve_p, angle_p, ipv_p, nev, nev0, ev_p)
polygon *poly;
int vcirc;
double *tol;
int nve;
int *nv, *nev, *nev0;
double **ve_p, **angle_p;
int **ipv_p, **ev_p;
{
    static int nvmax = 0;
    static int *ipv = 0x0, *ev = 0x0;
    static double *v = 0x0, *ve = 0x0, *angle = 0x0;

    int ier;

    /* first guess at number of vertices */
    *nv = poly->np;
    
    /* keep trying till the arrays are big enough */
    do {

	/* make sure that allocated arrays contain enough space */
	if (!ve || !angle || !ipv || !ev || *nv > nvmax) {
	    if (v) free(v);
	    if (ve) free(ve);
	    if (angle) free(angle);
	    if (ipv) free(ipv);
	    if (ev) free(ev);
	    ve = (double *) malloc(sizeof(double) * (*nv + DNV) * nve * 3);
	    if (!ve) {
		fprintf(stderr, "failed to allocate memory for %d x %d vectors\n", *nv + DNV, nve);
		return(-1);
	    }
	    angle = (double *) malloc(sizeof(double) * (*nv + DNV));
	    if (!angle) {
		fprintf(stderr, "failed to allocate memory for %d doubles\n", *nv + DNV);
		return(-1);
	    }
	    ipv = (int *) malloc(sizeof(int) * (*nv + DNV));
	    if (!ipv) {
		fprintf(stderr, "failed to allocate memory for %d integers\n", *nv + DNV);
		return(-1);
	    }
	    ev = (int *) malloc(sizeof(int) * (*nv + DNV));
	    if (!ev) {
		fprintf(stderr, "failed to allocate memory for %d integers\n", *nv + DNV);
		return(-1);
	    }
	    nvmax = *nv + DNV;
	}

	/* compute vertices of polygon */
	ier = gvert(poly, nvmax, vcirc, tol, nv, nve, ve, angle, ipv, nev, nev0, ev);
	if (ier) return(1);

    } while (*nv > nvmax);

    /* point ve_p, angle_p, ipv_p, and ev_p at ve, angle, ipv, and ev */
    *ve_p = ve;
    *angle_p = angle;
    *ipv_p = ipv;
    *ev_p = ev;
		

    return(0);
}

/*------------------------------------------------------------------------------
  Vertices of polygon.

  This is a c interface to fortran subroutine gvert.

   Input: poly is a polygon.
	  nvmax = dimension of v[nvmax][3] and ev[nmax].
	  vcirc = 1 to return vertices and midpoints also for bounding circles
		  which have no intersections;
		= 0 not so.
	  nve = desired number of points on each edge, including vertex.
  Input/Output: *tol = angle within which to merge multiple intersections.
  Output: *nv = number of vertices.
	  ve[nv][nve][3] = points on edges of polygon.
	  angle[nv] = angular lengths of edges of polygon.
	  ipv[nv] = indices of vertices/edges;
		    that is, vertex v[i] and edge points ve[i]
		    lie on boundary number ipv.
	  *nev = number of connected sequences of vertices.
	  *nev0 = number of bounding circles which have no intersections.
	  ev[nev] = end indices of each connected sequence of vertices.
  Return value: 0 if ok;
		1 if fatal intersection of boundaries;
		-1 if failed to allocate memory.
*/
int gvert(poly, nvmax, vcirc, tol, nv, nve, ve, angle, ipv, nev, nev0, ev)
polygon *poly;
int nvmax;
int vcirc;
int nve;
double *tol;
int *nv, *nev, *nev0;
#ifdef GCC
double ve[nvmax][3];
double angle[nvmax];
int ipv[nvmax];
int ev[nvmax];
#else
double ve[][3];
double angle[];
int ipv[];
int ev[];
#endif
{
    int iv;
    logical ldegen;
    /* work arrays */
    int *iord, *iwk;
    double *phi, *wk;

    /* allocate memory for work arrays */
    iord = (int *) malloc(sizeof(int) * poly->np * 2);
    if (!iord) {
	fprintf(stderr, "gvert: failed to allocate memory for %d integers\n", poly->np * 2);
	return(-1);
    }
    phi = (double *) malloc(sizeof(double) * poly->np * 2);
    if (!phi) {
	fprintf(stderr, "gvert: failed to allocate memory for %d doubles\n", poly->np * 2);
	return(-1);
    }
    iwk = (int *) malloc(sizeof(int) * nvmax * 3);
    if (!iwk) {
	fprintf(stderr, "gvert: failed to allocate memory for %d integers\n", nvmax * 3);
	return(-1);
    }
    wk = (double *) malloc(sizeof(double) * nvmax);
    if (!wk) {
	fprintf(stderr, "gvert: failed to allocate memory for %d doubles\n", nvmax);
	return(-1);
    }

    /* fortran routine */
    gvert_(ve, angle, ipv, ev, &nvmax, nv, &nve, nev, nev0, poly->rp, poly->cm, &poly->np, &vcirc, tol, phi, iord, wk, iwk, &ldegen);

    /* convert vertex/edge indices from fortran to c convention */
    for (iv = 0; iv < *nv && iv < nvmax; iv++) ipv[iv]--;

    /* free work arrays */
    free(iord);
    free(phi);
    free(iwk);
    free(wk);

    /* fatal intersection of boundaries */
    if (ldegen) return(1);

    return(0);
}
