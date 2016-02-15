#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "polygon.h"
#include "vertices.h"

#include "pi.h"
#define TWOPI		(2. * PI)
#define HALFTWOPI	PI
#define PIBYTWO		(PI / 2.)

/* local functions */
int poly_to_rect();
double cmij();
polygon *vert_to_poly(), *edge_to_poly();
vertices *rps_to_vert();
void rp_to_vert(), azel_to_gc(), rp_to_gc(), edge_to_rpcm(), rp_to_rpcm(),
	circ_to_rpcm(), rpcm_to_circ(), az_to_rpcm(), el_to_rpcm();

/* external functions */ 
extern polygon *new_poly();
extern vertices *new_vert();
extern void free_vert();
extern void msg(char *, ...);

/*------------------------------------------------------------------------------
  Convert vertices structure to polygon.
  Assume that vertices are joined by great circles.

  Return value: pointer to polygon,
		or null if error occurred.
*/
polygon *vert_to_poly(vert)
vertices *vert;
{
    int iv;
    polygon *poly = 0x0;

    /* make polygon with nv boundaries */
    poly = new_poly(vert->nv);
    if (!poly) return(0x0);

    /* convert each pair of adjacent vertices to a great circle */
    poly->np = vert->nv;
    for (iv = 0; iv < vert->nv; iv++) {
	azel_to_gc(&vert->v[iv], &vert->v[(iv+1) % vert->nv],
	    &poly->rp_(0, iv), &poly->cm[iv]);
    }

    return(poly);
}

/*------------------------------------------------------------------------------
  Convert edges, stored in vertices structure, to polygon.

  Input: vert = pointer to vertices structure.
	 nve = number of points per edge.
	 ev = cumulative number of points on each connected boundary.
  Return value: pointer to polygon,
		or null if error occurred.
*/
polygon *edge_to_poly(vert, nve, ev)
vertices *vert;
int *ev;
{
    int dev, i, evo, iedge, iv, iv0, iv1, iv2, jv, nedge;
    polygon *poly = 0x0;

    /* number of edges */
    nedge = vert->nv / nve;

    /* make polygon with nedge boundaries */
    poly = new_poly(nedge);
    if (!poly) return(0x0);
    poly->np = nedge;

    evo = 0;
    iv = 0;
    iedge = 0;
    /* do each connected boundary */
    for (jv = 0; evo < vert->nv; jv++) {
	dev = ev[jv] - evo;
	/* points on connected boundary */
	while (iv < ev[jv]) {
	    i = iv - evo;
	    iv0 = evo + i;
	    iv1 = evo + (i + 1) % dev;
	    iv2 = evo + (i + 2) % dev;
	    if (nve == 1) {
		/* convert pair of adjacent vertices to a great circle */
		azel_to_gc(&vert->v[iv0], &vert->v[iv1],
		    &poly->rp_(0, iedge), &poly->cm[iedge]);
	    } else {
		/* convert triple of vertices to circle */
		edge_to_rpcm(&vert->v[iv0], &vert->v[iv1], &vert->v[iv2],
		    &poly->rp_(0, iedge), &poly->cm[iedge]);
	    }
	    iv += nve;
	    iedge++;
	}
	evo = ev[jv];
    }

    return(poly);
}

/*------------------------------------------------------------------------------
  Convert unit vectors to vertices structure.

   Input: rp = array of nv unit vectors.
	  nv = number of unit vectors.
  Return value: pointer to a vertices structure,
		or null if an error occurred.
*/
vertices *rps_to_vert(nv, rp)
int nv;
#ifdef GCC
double rp[nv][3];
#else
double rp[][3];
#endif
{
/* number of extra vertices to allocate, to allow for expansion */
#define DNV             4
    static vertices *vert = 0x0;

    int iv;

    /* ensure that vert contains enough space */
    if (!vert || nv > vert->nvmax) {
	free_vert(vert);
	vert = new_vert(nv + DNV);
	if (!vert) {
	    fprintf(stderr, "failed to allocate memory for vertices structure of %d vertices\n", nv + DNV);
	    return(0x0);
	}
    }

    /* convert vectors to vertices */
    vert->nv = nv;
    for (iv = 0; iv < nv; iv++) {
	rp_to_vert(rp[iv], &vert->v[iv]);
    }

    /* phase vertices to the first */
    for (iv = 1; iv < nv; iv++) {
	vert->v[iv].az -= rint((vert->v[iv].az - vert->v[0].az) / TWOPI) * TWOPI;
    }

    return(vert);
}

/*------------------------------------------------------------------------------
  Convert unit vector to vertex.
*/
void rp_to_vert(rp, v)
double rp[3];
azel *v;
{
    v->az = atan2(rp[1], rp[0]);
    v->el = atan2(rp[2], sqrt(rp[0] * rp[0] + rp[1] * rp[1]));
}

/*------------------------------------------------------------------------------
  Determine rp, cm for great circle passing through two az-el vertices.
  The great circle goes right-handedly from v0 to v1.
*/
void azel_to_gc(v0, v1, rp, cm)
azel *v0, *v1;
double rp[3], *cm;
{
    azel *v;
    int iv;
    double rpv[2][3];

    for (iv = 0; iv < 2; iv++) {
	v = (iv == 0)? v0: v1;
	rpv[iv][0] = cos(v->el) * cos(v->az);
	rpv[iv][1] = cos(v->el) * sin(v->az);
	rpv[iv][2] = sin(v->el);
    }

    rp_to_gc(rpv[0], rpv[1], rp, cm);
}

/*------------------------------------------------------------------------------
  Determine rp, cm for great circle passing through two unit vectors.
  The great circle goes right-handedly from rp0 to rp1.
  If the two unit vectors coincide, suppress the boundary,
  on the assumption that coincident points are redundant
  (for example, a mask-maker may specify a triangle with 4 vertices,
  with 2 vertices being coincident).
*/
void rp_to_gc(rp0, rp1, rp, cm)
double rp0[3], rp1[3];
double rp[3], *cm;
{
    int i;
    double rpa;

    /* cofactors */
    rp[0] = rp0[1]*rp1[2] - rp1[1]*rp0[2];
    rp[1] = rp1[0]*rp0[2] - rp0[0]*rp1[2];
    rp[2] = rp0[0]*rp1[1] - rp1[0]*rp0[1];

    rpa = sqrt(rp[0]*rp[0] + rp[1]*rp[1] + rp[2]*rp[2]);

    /* indeterminate solution (rp0 and rp1 are same or antipodeal points) */
    if (rpa == 0.) {
    /* suppress boundary, on assumption that coincident points are redundant */
	*cm = 2.;
    /* normalize rp to 1 */
    } else {
	for (i = 0; i < 3; i++) rp[i] /= rpa;
	/* cm = 1 for great circle */
	*cm = 1.;
    }
}

/*------------------------------------------------------------------------------
  Determine rp, cm for circle passing through three az-el points.
  The circle goes right-handedly from v0 to v1 to v2.
*/
void edge_to_rpcm(v0, v1, v2, rp, cm)
azel *v0, *v1, *v2;
double rp[3], *cm;
{
    azel *v;
    int iv;
    double rpv[3][3];

    for (iv = 0; iv < 3; iv++) {
	switch (iv) {
	case 0:	v = v0;	break;
	case 1:	v = v1;	break;
	case 2:	v = v2;	break;
	}
	rpv[iv][0] = cos(v->el) * cos(v->az);
	rpv[iv][1] = cos(v->el) * sin(v->az);
	rpv[iv][2] = sin(v->el);
    }

    rp_to_rpcm(rpv[0], rpv[1], rpv[2], rp, cm);
}

/*------------------------------------------------------------------------------
  Determine rp, cm for circle passing through three unit vectors.
  The circle goes right-handedly from rp0 to rp1 to rp2.
  If two of the unit vectors coincide, join with a great circle.
  If three of the unit vectors coincide, suppress the boundary.
*/
void rp_to_rpcm(rp0, rp1, rp2, rp, cm)
double rp0[3], rp1[3], rp2[3];
double rp[3], *cm;
{
    int coincide, i, j;
    double det, rpa;
    double *rpi, *rpj;

    /* check whether any two of the three unit vectors coincide */
    coincide = 0;
    for (j = 0; j < 3; j++) {
	switch (j) {
	case 0:	rpj = rp0;	break;
	case 1:	rpj = rp1;	break;
	case 2:	rpj = rp2;	break;
	}
	for (i = 0; i < j; i++) {
	    switch (i) {
	    case 0:	rpi = rp0;	break;
	    case 1:	rpi = rp1;	break;
	    case 2:	rpi = rp2;	break;
	    }
	    /* vector coincide */
	    if (rpi[0] == rpj[0] && rpi[1] == rpj[1] && rpi[2] == rpj[2]) {
		coincide = 1;
		/* set rpj equal to the third vector */
		if (i == 0 && j == 1) {
		    rpj = rp2;
		} else if (i == 0 && j == 2) {
		    rpj = rp1;
		} else if (i == 1 && j == 2) {
		    rpj = rp0;
		}
	    }
	    if (coincide) break;
	}
	if (coincide) break;
    }

    /* if any two vectors coincide, join with great circle */
    if (coincide) {
	rp_to_gc(rpi, rpj, rp, cm);

    /* non-coincident vectors */
    } else {
	/* cofactors, arranged to reduce roundoff */
	rp[0] = rp0[1] * (rp1[2] - rp2[2])
	      + rp1[1] * (rp2[2] - rp0[2])
	      + rp2[1] * (rp0[2] - rp1[2]);
	rp[1] = rp0[0] * (rp2[2] - rp1[2])
	      + rp1[0] * (rp0[2] - rp2[2])
	      + rp2[0] * (rp1[2] - rp0[2]);
	rp[2] = (rp0[0] * rp1[1] - rp1[0] * rp0[1])
	      + (rp1[0] * rp2[1] - rp2[0] * rp1[1])
	      + (rp2[0] * rp0[1] - rp0[0] * rp2[1]);

	/* |rp| */
	rpa = sqrt(rp[0]*rp[0] + rp[1]*rp[1] + rp[2]*rp[2]);

	/* indeterminate solution (2 of 3 unit vectors coincide: shouldn't happen) */
	if (rpa == 0.) {
	    /* suppress boundary, not knowing what to do with it */
	    *cm = 2.;
	/* normal solution */
	} else {
	    /* determinant */
	    det = (rp0[0] * rp1[1] - rp1[0] * rp0[1]) * rp2[2]
		+ (rp1[0] * rp2[1] - rp2[0] * rp1[1]) * rp0[2]
		+ (rp2[0] * rp0[1] - rp0[0] * rp2[1]) * rp1[2];
	    /* cos(th) = det / |rp| */
	    *cm = 1. - det / rpa;
	    /* normalize rp to 1 */
	    for (i = 0; i < 3; i++) rp[i] /= rpa;
	}
    }
}

/*------------------------------------------------------------------------------
  Convert circle to rp, cm.

   Input: angle = (azimuth, elevation, radius) in radians.
  Output: rp, cm as used by garea, gspher et al.
*/
void circ_to_rpcm(angle, rp, cm)
double angle[3];
double rp[3], *cm;
{
    double s;

    /* Cartesian coordinates of azimuth, elevation */
    rp[0] = cos(angle[1]) * cos(angle[0]);
    rp[1] = cos(angle[1]) * sin(angle[0]);
    rp[2] = sin(angle[1]);
    /* 1 - cos(radius) = 2 sin^2(radius/2) */
    s = sin(angle[2] / 2.);
    *cm = s * s * 2.;
    *cm = (angle[2] >= 0.)? *cm : -*cm;
}

/*------------------------------------------------------------------------------
  Convert rp, cm to circle.

   Input: rp, cm as used by garea, gspher et al.
  Output: angle = (azimuth, elevation, radius) in radians.
*/
void rpcm_to_circ(rp, cm, angle)
double rp[3], *cm;
double angle[3];
{
    double s;

    angle[0] = atan2(rp[1], rp[0]);
    angle[1] = atan2(rp[2], sqrt(rp[0] * rp[0] + rp[1] * rp[1]));
    s = sqrt(fabs(*cm) / 2.);
    if (s > 1.) s = 1.;
    angle[2] = 2. * asin(s);
    angle[2] = (*cm >= 0.)? angle[2] : -angle[2];
}

/*------------------------------------------------------------------------------
  Convert line of constant azimuth to rp, cm.

   Input: az = azimuth in radians.
	  m = 0 for minimum elevation;
	    = 1 for maximum elevation.
  Output: rp, cm as used by garea, gspher et al.
*/
void az_to_rpcm(az, m, rp, cm)
double az;
int m;
double rp[3], *cm;
{
    /* axis along equator */
    rp[0] = - sin(az);
    rp[1] = cos(az);
    rp[2] = 0.;
    /* 1 - cos(th) = 1 for great circle */
    *cm = 1.;
    /* min, max */
    *cm = (m == 0)? *cm : - *cm;
}

/*------------------------------------------------------------------------------
  Convert line of constant elevation to rp, cm.

   Input: el = elevation in radians.
	  m = 0 for minimum elevation;
	    = 1 for maximum elevation.
  Output: rp, cm as used by garea, gspher et al.
*/
void el_to_rpcm(el, m, rp, cm)
double el;
int m;
double rp[3], *cm;
{
    /* north pole */
    rp[0] = 0.;
    rp[1] = 0.;
    rp[2] = 1.;
    /* 1 - cos(th) = 1 - sin(el) */
    *cm = 1 - sin(el);
    /* min, max */
    *cm = (m == 0)? *cm : - *cm;
}

/*------------------------------------------------------------------------------
   1-cos(theta_ij) = 2 sin^2(theta_ij/2) between two unit vectors.
*/
double cmij(rpi, rpj)
double rpi[3], rpj[3];
{
    double cm, dx, dy, dz;

    dx = rpi[0] - rpj[0];
    dy = rpi[1] - rpj[1];
    dz = rpi[2] - rpj[2];
    cm = (dx*dx + dy*dy + dz*dz) / 2.;

    return(cm);
}

/*------------------------------------------------------------------------------
  Determine whether a polygon is a rectangle,
  and if so return <azmin> <azmax> <elmin> <elmax>.

   Input: poly = pointer to polygon.
  Output: *azmin, *azmax, *elmin, *elmax;
	  if *azmin >= *azmax, or *elmin >= *elmax, the rectangle is empty;
	  if *azmin = *azmax + 2*pi, there is no azimuthal constraint.
  Return value: 0 if polygon is not a rectangle,
		1 if polygon is a rectangle,
		2 if polygon is a rectangle with superfluous boundaries.
*/
int poly_to_rect(poly, azmin, azmax, elmin, elmax)
polygon *poly;
double *azmin, *azmax, *elmin, *elmax;
{
    int iaz, ielmin, ielmax, ip;
    double az, el;

    *azmin = -TWOPI;
    *azmax = TWOPI;
    *elmin = -PIBYTWO;
    *elmax = PIBYTWO;
    iaz = 0;
    ielmin = 0;
    ielmax = 0;
    for (ip = 0; ip < poly->np; ip++) {
	/* null polygon */
	if (poly->cm[ip] == 0. || poly->cm[ip] <= -2.) return(0);
	/* skip superfluous cap */
	if (poly->cm[ip] >= 2.) continue;
	/* line of constant azimuth */
	if (poly->rp_(2, ip) == 0.
	    && (poly->cm[ip] == 1. || poly->cm[ip] == -1.)) {
	    /* az is in (-pi,pi] */
	    az = atan2(-poly->rp_(0, ip), poly->rp_(1, ip));
	    if (poly->cm[ip] == -1.) az += (az <= 0.)? PI : -PI;
	    /* shift azmin, azmax to phase of az, az+pi */
	    if (az >= *azmax) {
		while (az >= *azmax) {
		    *azmin += TWOPI;
		    *azmax += TWOPI;
		}
	    } else if (az + PI <= *azmin) {
		while (az + PI <= *azmin) {
		    *azmin -= TWOPI;
		    *azmax -= TWOPI;
		}
	    }
	    /* revise azmin, azmax */
	    if (az > *azmin) *azmin = az;
	    if (az + PI < *azmax) *azmax = az + PI;
	    iaz++;
	/* line of constant elevation */
	} else if (poly->rp_(0, ip) == 0. && poly->rp_(1, ip) == 0.) {
	    el = asin(1. - fabs(poly->cm[ip]));
	    if (poly->rp_(2, ip) < 0.) el = -el;
	    if ((poly->rp_(2, ip) > 0. && poly->cm[ip] >= 0.)
		|| (poly->rp_(2, ip) < 0. && poly->cm[ip] < 0.)) {
		if (el > *elmin) *elmin = el;
		ielmin++;
	    } else {
		if (el < *elmax) *elmax = el;
		ielmax++;
	    }
	} else {
	    return(0);
	}
    }

    /* azimuthal constraints include nothing */
    if (*azmin > *azmax) {
	*azmin = *azmax = 0.;
    /* azimuthal constraints include everything */
    } else if (*azmin + TWOPI <= *azmax) {
	*azmin = -PI;
	*azmax = PI;
    /* phase azmin to (-pi,pi] */
    } else if (*azmin > PI) {
        *azmin -= TWOPI;
        *azmax -= TWOPI;
    }

    /* should not happen */
    if (iaz > 2 || ielmin > 1 || ielmax > 1) {
	msg("poly_to_rect: polygon contains");
	if (iaz > 2) msg(" %d az boundaries", iaz);
	if (ielmin > 1) msg(" %d elmin boundaries", ielmin);
	if (ielmax > 1) msg(" %d elmax boundaries", ielmax);
	msg("\n");
	return(2);
    }

    return(1);
}
