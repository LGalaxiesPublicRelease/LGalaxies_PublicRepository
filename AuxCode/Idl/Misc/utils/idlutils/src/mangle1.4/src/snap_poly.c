/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Make almost coincident caps of 2 polygons coincide.
  Caps of poly2 are adjusted to equal those of poly1.

  Input:  poly1, poly2 = pointers to polygon structures.
	  axtol = angle in radians [actually 2 sin(angle/2)]:
		  if angle twixt polar axes of caps <= axtol,
		  then make axis of poly2 cap
		  exactly parallel to axis of poly1 cap.
	  btol =  angle in radians:
		  if two axes of caps of poly1 and poly2 are parallel,
		  and if angle between latitudes of caps <= btol,
		  then make latitude of poly2 cap
		  exactly equal to latitude of poly1 cap.
  Output: adjusted caps of poly2 (i.e. poly2->rp, poly2->cm).
  Return value: number of caps adjusted.
*/
int snap_poly(polygon *poly1, polygon *poly2, double axtol, double btol)
{
    int adjusted, ip, ip1, ip2, nadj, sp;
    double cm, dl, drp, dx, dy, dz;

    nadj = 0;
    for (ip1 = 0; ip1 < poly1->np; ip1++) {	/* for each cap of poly1 ... */
	/* superfluous cap */
	if (poly1->cm[ip1] == 0. || fabs(poly1->cm[ip1]) >= 2.) continue;
	for (ip2 = 0; ip2 < poly2->np; ip2++) {	/* ... and each cap of poly2 */
	    /* superfluous cap */
	    if (poly2->cm[ip2] == 0. || fabs(poly2->cm[ip2]) >= 2.) continue;
	    for (ip = 0; ip < 2; ip++) {	/* check rp2 = +- rp1 */
		adjusted = 0;
		sp = (ip == 0)? 1 : -1;
		/* [2 sin(alpha/2)]^2, where alpha is angle twixt axes */
		dx = poly2->rp[ip2][0] - sp * poly1->rp[ip1][0];
		dy = poly2->rp[ip2][1] - sp * poly1->rp[ip1][1];
		dz = poly2->rp[ip2][2] - sp * poly1->rp[ip1][2];
		drp = sqrt(dx * dx + dy * dy + dz * dz);
		if (drp <= axtol) {		/* axes are nearly parallel */
		    if (!(poly2->rp[ip2][0] == sp * poly1->rp[ip1][0]
		      && poly2->rp[ip2][1] == sp * poly1->rp[ip1][1]
		      && poly2->rp[ip2][2] == sp * poly1->rp[ip1][2])) {
			/* make axis of poly2 cap exactly parallel to poly1
			   (made exactly equal below if caps nearly coincide) */
			poly2->rp[ip2][0] = sp * poly1->rp[ip1][0];
			poly2->rp[ip2][1] = sp * poly1->rp[ip1][1];
			poly2->rp[ip2][2] = sp * poly1->rp[ip1][2];
			adjusted = 1;
		    }
		    /* angle between latitudes of caps */
		    if (sp == 1) {		/* axes are aligned */
			dl = 2. * (asin(sqrt(fabs(poly2->cm[ip2]) / 2.))
			    - asin(sqrt(fabs(poly1->cm[ip1]) / 2.)));
		    } else {			/* axes are anti-aligned */
			dl = 2. * (asin(sqrt((2. - fabs(poly2->cm[ip2])) / 2.))
			    - asin(sqrt(fabs(poly1->cm[ip1]) / 2.)));
		    }
		    if (fabs(dl) <= btol) {	/* caps nearly coincide */
			if (sp == -1) {
			    /* reflect axis of poly2 cap */
			    poly2->rp[ip2][0] = - poly2->rp[ip2][0];
			    poly2->rp[ip2][1] = - poly2->rp[ip2][1];
			    poly2->rp[ip2][2] = - poly2->rp[ip2][2];
			    adjusted = 1;
			}
			/* set latitude of poly2 cap equal to poly1 */
			cm = (poly2->cm[ip2] >= 0.)?
			    sp * fabs(poly1->cm[ip1]):
			    - sp * fabs(poly1->cm[ip1]);
			if (poly2->cm[ip2] != cm) {
			    poly2->cm[ip2] = cm;
			    adjusted = 1;
			}
		    }
		    if (adjusted) nadj++;
		}
	    }
	}
    }
    return(nadj);
}
/*------------------------------------------------------------------------------
  Snap edge of poly2 to cap boundary of poly1.
  Caps of poly2 are adjusted to equal those of poly1.

  Input:  poly1, poly2 = pointers to polygon structures.
	  thtol = edge tolerance in radians.
	  ytol = edge to length tolerance;
		 if the two vertices and centre point of an edge of poly2 are
		 all closer to a boundary of poly1 than the lesser of
		 (1) thtol, and
		 (2) ytol times the length of the edge,
		 and if in addition at least one of the three points lies
		 inside poly1 (sans said boundary),
		 then make boundary of the poly2 cap equal to that of poly1.
	  mtol = initial tolerance angle for multiple intersections in radians.
  Output: adjusted caps of poly2 (i.e. poly2->rp, poly2->cm).
  Return value: number of caps adjusted,
		or -1 if error occurred.
*/
int snap_polyth(polygon *poly1, polygon *poly2, double thtol, double ytol, double mtol)
{
    const int per = 0;
    const int nve = 2;

    int adjusted, do_vcirc, i, ier, in, ip1, ip2, iv, ivp, nadj, nev, nev0, nv;
    int *ipv, *gp, *ev;
    double cm, cm1, dth, dthmax, sp, tol;
    double *angle;
    vec *v, *ve;

    /* vertices and centres of edges of poly2 */
    do_vcirc = 0;
    tol = mtol;
    ier = gverts(poly2, do_vcirc, &tol, per, nve, &nv, &ve, &angle, &ipv, &gp, &nev, &nev0, &ev);
    if (ier != 0) return(-1);

    /* convert angle of each edge to scalar length angle * sin(theta) */
    for (iv = 0; iv < nv; iv++) {
	ip2 = ipv[iv];
	cm = fabs(poly2->cm[ip2]);
	angle[iv] = angle[iv] * sqrt(cm * (2. - cm));
    }

    nadj = 0;
    /* for each edge of poly2 ... */
    for (iv = 0; iv < nv; iv++) {
	ivp = (iv + 1) % nv;
	ip2 = ipv[iv];

	/* ... and each axis of poly1 */
	for (ip1 = 0; ip1 < poly1->np; ip1++) {
	    adjusted = 0;

	    /* distance from edge of poly2 to cap of poly1 */
	    cm1 = poly1->cm[ip1];
            poly1->cm[ip1] = 2.;        /* suppress cap of poly1 */
	    in = 0;
	    dthmax = 0.;
	    for (i = 0; i < 3; i++) {
		/* vertex, centre point, vertex of edge of poly2 */
		v = &ve[(iv * nve + i) % (nv * nve)];
		in |= gptin(poly1, *v);	/* in if any one point is in */
		cm = cmij(*v, poly1->rp[ip1]);
		dth = 2. * (sqrt(cm/2.) - sqrt(fabs(cm1/2.)));
		dth = fabs(dth);	/* angle from point to cap of poly1 */
		if (dth > dthmax) dthmax = dth;
	    }
            poly1->cm[ip1] = cm1;       /* restore cap of poly1 */

	    /* three points of poly2 edge are all close to boundary of poly1 */
	    if (in && dthmax <= thtol && dthmax <= ytol * angle[iv]) {
		sp = poly1->rp[ip1][0] * poly2->rp[ip2][0] + poly1->rp[ip1][1] * poly2->rp[ip2][1] + poly1->rp[ip1][2] * poly2->rp[ip2][2];
		sp = (sp >= 0.)? 1. : -1.;
		if (!(poly2->rp[ip2][0] == poly1->rp[ip1][0]
		  && poly2->rp[ip2][1] == poly1->rp[ip1][1]
		  && poly2->rp[ip2][2] == poly1->rp[ip1][2])) {
		    /* make axis of poly2 cap exactly equal to that of poly1 */
		    poly2->rp[ip2][0] = poly1->rp[ip1][0];
		    poly2->rp[ip2][1] = poly1->rp[ip1][1];
		    poly2->rp[ip2][2] = poly1->rp[ip1][2];
		    adjusted = 1;
		}
		/* set latitude of poly2 cap equal to that of poly1 */
		cm = (poly2->cm[ip2] >= 0.)?
		    sp * fabs(poly1->cm[ip1]):
		    - sp * fabs(poly1->cm[ip1]);
		if (poly2->cm[ip2] != cm) {
		    poly2->cm[ip2] = cm;
		    adjusted = 1;
		}
		if (adjusted) nadj++;
	    }

	}

    }

    /* trim adjusted polygon */
    if (nadj > 0) trim_poly(poly2);

    return(nadj);
}
