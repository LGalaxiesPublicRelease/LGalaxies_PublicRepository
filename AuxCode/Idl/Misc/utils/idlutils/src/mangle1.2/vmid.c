#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "polygon.h"

/* initial angular tolerance within which to merge multiple intersections */
extern double mtol;

/* number of extra vertices to allocate, to allow for expansion */
#define DNV		4

/* local functions */
int vmid();

/* external declarations */
extern int gvlims(), gvphi();
extern void rp_to_gc();

/*------------------------------------------------------------------------------
  Point(s) somewhere in the middle of a polygon.
  The number of points is equal to the number of connected sequences of
  vertices, as computed by gverts, except that if there are zero connected
  sequences, then there is a single midpoint.  Zero connected sequences
  happens when the option to return midpoints for non-intersecting circles is
  turned off in gverts, and each boundary of the polygon is a disjoint circle.

  If the connected sequences of vertices delineate holes in the polygon,
  then the number of connected sequences of vertices may exceed
  the number of connected parts of the polygon.  In this case vmid supplies
  more midpoints than there are connected parts of the polygon.

  The intention is to aim for a point (or points) squarely in the middle of
  the polygon, but the algorithm, though mildly paranoid, is not cast iron.

  Designed to be called after gverts().

   Input: poly = pointer to polygon.
	  nv = number of vertices of polygon, as output by gverts().
	  nve = number of points per edge, as input to gverts();
		must be even and >= 2.
	  ve = points on edges of polygon, as output by gverts().
	  ev = end indices of connected sequences, as output by gverts().
  Output: *nvm = number of middle points in *vm
	       = max(nev, 1).
	  *vm = pointer to unit vectors giving middle points;
		the required memory is allocated.
		If the algorithm failed to find a point inside the polygon,
		then the vector for that point is {0.,0.,0.}.
  Return value: 0 if ok (even if no inside point was found);
		-1 if failed to allocate memory.
*/
int vmid(poly, nv, nve, ve, ev, nvm, vm_p)
polygon *poly;
int nv, nve;
#ifdef GCC
double ve[nv * nve][3];
int ev[nv];
#else
double ve[][3];
int ev[];
#endif
int *nvm;
double **vm_p;
{
    static int do_vcirc = 1;
    static int nvmmax = 0;
    static double *vm = 0x0;

    int i, ier, iev, iv, ivl, ivm, ivu, neva, nev0, nva;
    double angle, cm, cmbest, rp[3], s, tol, vc[3], vp[3];
    int *eva;
    double *vmin, *vmax, *cmmin, *cmmax;

    /* make sure nve is even and >= 2 */
    if (nve < 2 || nve % 2 != 0) {
	fprintf(stderr, "vmid: nve = %d must be even and >= 2\n", nve);
	fprintf(stderr, "stop\n");
	exit(1);
    }

    /* number of midpoints */
    for (ivm = 0; ev[ivm] != nv; ivm++);
    *nvm = ivm + 1;

    /* make sure that vm contains enough space */
    if (!vm || *nvm > nvmmax) {
	if (vm) free(vm);
	vm = malloc(sizeof(double) * (*nvm + DNV) * 3);
	if (!vm) {
	    fprintf(stderr, "failed to allocate memory for %d vectors\n", *nvm + DNV);
	    return(-1);
	}
	nvmmax = *nvm + DNV;
    }

    /* polygon is whole sphere: choose point to be north pole */
    if (poly->np == 0) {
	vm[0] = 0.;
	vm[1] = 0.;
	vm[2] = 1.;

    /* polygon has 1 cap: choose point at centre of circle */
    } else if (poly->np == 1) {
	s = (poly->cm[0] >= 0.)? 1. : -1.;
	vm[0] = s * poly->rp_(0, 0);
	vm[1] = s * poly->rp_(1, 0);
	vm[2] = s * poly->rp_(2, 0);

    /* polygon has >= 2 caps */
    } else {

	/* index of central point of edge */
	iev = nve / 2;

	/* no vertices implies each boundary is a disjoint circle */
	if (nv == 0) {
	    /* set `edge points' equal to centres of circles */
	    for (iv = 0; iv < poly->np; iv++) {	/* gverts allocated enough memory for this */
		s = (poly->cm[iv] >= 0.)? 1. : -1.;
		for (i = 0; i < 3; i++) {
		    ve[iv * nve + iev][i] = s * poly->rp_(i, iv);
		}
	    }
	}

	/* do each connected sequence of vertices */
	for (ivm = 0; ivm < *nvm; ivm++) {

	    /* initialize the midpoint to zero */
	    for (i = 0; i < 3; i++) vm[i + 3 * ivm] = 0.;

	    /* vertices in connected sequence */
	    if (nv > 0) {
		ivl = (ivm == 0)? 0 : ev[ivm - 1];
		ivu = ev[ivm];
	    } else {		/* nv = 0 */
		ivl = 0;
		ivu = poly->np;	/* gverts allocated enough memory for this */
	    }

	    /* central vector is summed vector of points at centres of edges */
	    for (i = 0; i < 3; i++) vc[i] = 0.;
	    for (iv = ivl; iv < ivu; iv++) {
		for (i = 0; i < 3; i++) {
		    vc[i] += ve[iv * nve + iev][i];
		}
	    }
	    /* length of central vector */
	    s = sqrt(vc[0]*vc[0] + vc[1]*vc[1] + vc[2]*vc[2]);
	    /* s = 0 should probably never happen, but just in case ... */
	    if (s == 0.) {
		vc[0] = 0.;
		vc[1] = 0.;
		vc[2] = 1.;
	    /* normalize central vector to unit length */
	    } else {
		for (i = 0; i < 3; i++) {
		    vc[i] /= s;
		}
	    }

	    /* try several possible central points, and choose the best */
	    cmbest = -1.;
	    for (iv = ivl; iv < ivu; iv++) {
		/* great circle joining centre of edge to central vector vc */
		rp_to_gc(ve[iv * nve + iev], vc, rp, &cm);
		/* central point vp of great circle inside polygon */
		tol = mtol;
		ier = gvphi(poly, rp, cm, ve[iv * nve + iev], &tol, &angle, vp);
		if (ier == -1) return(-1);
		/* distance of central point to edges */
		ier = gvlims(poly, do_vcirc, vp, &tol, &nva, &vmin, &vmax, &cmmin, &cmmax, &neva, &nev0, &eva);
		if (ier == -1) return(-1);
		if (ier) continue;
		/* closest distance of central point to edges */
		cm = 2.;
		for (i = ivl; i < ivu; i++) {
		    if (cmmin[i] < cm) cm = cmmin[i];
		}
		/* best bet so far is the point farthest from any edge */
		if (cm > cmbest) {
		    cmbest = cm;
		    for (i = 0; i < 3; i++) {
			vm[i + 3 * ivm] = vp[i];
		    }
		}
	    }

	}

    }

    /* point vm_p at vm array */
    *vm_p = vm;

    return(0);
}

/*------------------------------------------------------------------------------
  Point(s) at the barycentre of the centres of the edges of a polygon.
  The point(s) may lie either inside or outside the polygon.
  The number of points is equal to the number of connected sequences of
  vertices, as computed by gverts, except that if there are zero connected
  sequences, then there is a single point.  Zero connected sequences
  happens when the option to return midpoints for non-intersecting circles is
  turned off in gverts, and each boundary of the polygon is a disjoint circle.

  Designed to be called after gverts().

   Input: poly = pointer to polygon.
	  nv = number of vertices of polygon, as output by gverts().
	  nve = number of points per edge, as input to gverts();
	  ve = points at centres of edges of polygon, as output by gverts().
	  ev = end indices of connected sequences, as output by gverts().
  Output: *nvm = number of middle points in *vm
	       = max(nev, 1).
	  *vm = pointer to unit vectors giving barycentres;
		the required memory is allocated.
		If the algorithm failed to find a point inside the polygon,
		then the vector for that point is {0.,0.,0.}.
  Return value: 0 if ok (even if no inside point was found);
		-1 if failed to allocate memory.
*/
int vmidc(poly, nv, nve, ve, ev, nvm, vm_p)
polygon *poly;
int nv, nve;
#ifdef GCC
double ve[nv * nve][3];
int ev[nv];
#else
double ve[][3];
int ev[];
#endif
int *nvm;
double **vm_p;
{
    static int nvmmax = 0;
    static double *vm = 0x0;

    int i, iev, iv, ivl, ivm, ivu;
    double s;
    double *vc;

    /* number of midpoints */
    for (ivm = 0; ev[ivm] != nv; ivm++);
    *nvm = ivm + 1;

    /* make sure that vm contains enough space */
    if (!vm || *nvm > nvmmax) {
	if (vm) free(vm);
	vm = malloc(sizeof(double) * (*nvm + DNV) * 3);
	if (!vm) {
	    fprintf(stderr, "failed to allocate memory for %d vectors\n", *nvm + DNV);
	    return(-1);
	}
	nvmmax = *nvm + DNV;
    }

    /* polygon is whole sphere: choose point to be north pole */
    if (poly->np == 0) {
	vm[0] = 0.;
	vm[1] = 0.;
	vm[2] = 1.;

    /* polygon has 1 cap: choose point at centre of circle */
    } else if (poly->np == 1) {
	vm[0] = poly->rp_(0, 0);
	vm[1] = poly->rp_(1, 0);
	vm[2] = poly->rp_(2, 0);

    /* no vertices implies each boundary is a disjoint circle */
    } else if (nv == 0) {
	/* choose point at centre of 1st circle */
	vm[0] = poly->rp_(0, 0);
	vm[1] = poly->rp_(1, 0);
	vm[2] = poly->rp_(2, 0);

    /* polygon has >= 2 caps */
    } else {

	/* index of central point of edge */
	iev = nve / 2;

	/* do each connected sequence of vertices */
	for (ivm = 0; ivm < *nvm; ivm++) {

	    /* point vc at appropriate element of vm */
	    vc = &vm[3 * ivm];

	    /* vertices in connected sequence */
	    ivl = (ivm == 0)? 0 : ev[ivm - 1];
	    ivu = ev[ivm];

	    /* central vector is summed vector of points at centres of edges */
	    for (i = 0; i < 3; i++) vc[i] = 0.;
	    for (iv = ivl; iv < ivu; iv++) {
		for (i = 0; i < 3; i++) {
		    vc[i] += ve[iv * nve + iev][i];
		}
	    }
	    /* length of central vector */
	    s = sqrt(vc[0]*vc[0] + vc[1]*vc[1] + vc[2]*vc[2]);
	    /* s = 0 should probably never happen, but just in case ... */
	    if (s == 0.) {
		vc[0] = 0.;
		vc[1] = 0.;
		vc[2] = 1.;
	    /* normalize central vector to unit length */
	    } else {
		for (i = 0; i < 3; i++) {
		    vc[i] /= s;
		}
	    }

	}

    }

    /* point vm_p at vm array */
    *vm_p = vm;

    return(0);
}
