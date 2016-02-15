#include <math.h>
#include <stdio.h>
#include "polygon.h"

/* initial angular tolerance within which to merge multiple intersections */
extern double mtol;

/* external functions */
extern int room_poly();
extern int snap_poly();
extern int trim_poly();
extern int vmidc();
extern int gptin();
extern int gverts(), gvlims();
extern void copy_poly(), poly_poly();
extern void gaxisi_();

/*------------------------------------------------------------------------------
  Partition disconnected polygon into connected polygons.
  Partitioning is done by lassooing each connected boundary of the polygon
  with an extra cap.  A lassoo is discarded if it lies inside the polygon,
  so the procedure works also for polygons that are not simply-connected.

  The value of the flag overwrite_original determines whether the
  original polygon in poly is overwritten or not.
  if poly is overwritten, then it is a copy of the polygon in polys[npoly].
  
	Input: *poly is a polygon.
	npolys = maximum number of polygons available in polys array.
	all_oneborder = 1 to lassoo all one-border polygons;
	0 to lassoo only those one-border polygons
	with more caps than vertices.
	In either case, all multi-border polygons are lassooed.
	adjust_lassoo = how to tighten lassoo:
	= 1 for balkanize,
	2 for ransack.
	overwrite_original = 2 to overwrite original polygon poly in all cases,
	whether or not lassoo succeeds;
	1 to overwrite poly only if lassoo succeeds;
	0 never to overwrite original poly.
  Output: *poly and polys[i], i = 0 to npoly-1, are the connected parts of poly.
	*npoly_try = attempted number of extra polygons in polys.
  Return value: npoly = number of connected polygons in polys,
	or -1 if error occurred.
	If npoly < *npoly_try, then the algorithm failed
	to lassoo a polygon.
*/
int partition_poly(poly, polys, npolys, all_oneborder, adjust_lassoo, 
									 overwrite_original, npoly_try)
		 int npolys, all_oneborder, adjust_lassoo;
		 polygon **poly;
#ifdef GCC
		 polygon *polys[npolys];
#else
		 polygon *polys[];
#endif
		 int *npoly_try;
{
	/* number of extra caps to allocate to polygon, to allow for expansion */
#define DNP		4
	static int do_vcirc = 1;
	static int nve = 2;
	static int itmax = 30;
	static polygon *extracap = 0x0;

	int i, ier, iev, in, ip, it, iv, ivm, ivthat, ivthis, nev, nev0, np, npoly, 
		nv, nvm;
	int *ipv, *ev;
	double *ve, *vm, *angle;
	double *vmax, *vmin, *cmmin, *cmmax;
	double cmthat, cmthis, cth, s, sth, th, thm, ththat, ththis, tol, v[3], 
		x[3], y[3];

	*npoly_try = 0;

	/* vertices and centres of edges of polygon */
	tol = mtol;
	ier = gverts(*poly, do_vcirc, &tol, &nv, nve, &ve, &angle, &ipv, &nev, 
							 &nev0, &ev);
	if (ier) return(-1);

	/* no boundaries */
	if (nev == 0) return(0);

	/* polygon has 1 connected border, and not too many caps */
	if (!all_oneborder && nev == 1 && (*poly)->np <= nv + 1) return(0);

	/* barycentres of connected boundaries of poly */
	ier = vmidc(*poly, nv, nve, ve, ev, &nvm, &vm);
	if (ier == -1) return(-1);

	/* number of polygons to try to split into */
	*npoly_try = nvm;

	/* attempt to partition polygon around each barycentre vm[ivm] */
	npoly = 0;
	for (ivm = 0; ivm < nvm; ivm++) {

		/* repeat until find isolating boundary */
		it = 0;
		do {

	    /* points on each edge nearest to and farthest from vm[ivm] */
	    tol = mtol;
	    ier = gvlims(*poly, do_vcirc, &vm[3 * ivm], &tol, &nv, &vmin, &vmax, 
									 &cmmin, &cmmax, &nev, &nev0, &ev);
	    if (ier == -1) return(-1);
	    if (ier) break;

	    /* distance that encloses edges of this part of polygon */
	    cmthis = 0.;
	    for (iv = (ivm == 0)? 0 : ev[ivm - 1]; iv < ev[ivm]; iv++) {
				if (cmmax[iv] > cmthis) {
					ivthis = iv;
					cmthis = cmmax[iv];
				}
	    }
	    /* distance that excludes edges of other parts of polygon */
	    cmthat = 2.;
	    for (iev = 0; iev < nev; iev++) {
				if (iev == ivm) continue;
				for (iv = (iev == 0)? 0 : ev[iev - 1]; iv < ev[iev]; iv++) {
					if (cmmin[iv] < cmthat) {
						ivthat = iv;
						cmthat = cmmin[iv];
					}
				}
	    }
	    ththis = 2. * asin(sqrt(cmthis / 2.));
	    ththat = 2. * asin(sqrt(cmthat / 2.));

	    /* printf("%21.15lg %21.15lg %21.15lg %4d %21.15lg %21.15lg\n", vm[3 * ivm], vm[1 + 3 * ivm], vm[2 + 3 * ivm], ivm, ththis, ththat); */

	    /* found boundary that isolates this part of polygon */
	    if (cmthis < cmthat) {

				/* Cartesian axes with z-axis to barycentre */
				gaxisi_(&vm[3 * ivm], x, y);
				/* opening angle of isolating boundary */
				th = (ththis + ththat) / 2.;
				cth = cos(th);
				sth = sin(th);
				/* point on isolating boundary */
				for (i = 0; i < 3; i++) {
					v[i] = cth * vm[i + 3 * ivm] + sth * x[i];
				}
				/* does (point on) isolating boundary lie inside or outside poly? */
				in = gptin(*poly, v);

				/* don't need isolating boundary that lies inside poly */
				if (in) {
					/* decrement number of polygons to try for */
					(*npoly_try)--;

					/* make new polygon with extra boundary */
				} else {
					/* not enough polygons for a new one */
					if (npoly >= npolys) return(*npoly_try + 1);

					/* put isolating boundary into new polygon */
					np = 1;
					ier = room_poly(&extracap, np, DNP, 0);
					if (ier == -1) {
						fprintf(stderr, "partition_poly: failed to allocate memory for polygon of %d caps\n", np + DNP);
						return(-1);
					}
					for (i = 0; i < 3; i++) {
						extracap->rp_(i, 0) = vm[i + 3 * ivm];
					}
					th = (ththis + ththat) / 2.;
					switch (adjust_lassoo) {
						/* for balkanize: tiny angles give garea problems */
					case 1:	thm = ththis + .001;	break;
						/* for ransack: want tight lassoo */
					case 2:	thm = ththis * 1.05;	break;
					}
					if (th > thm) th = thm;
					s = sin(th / 2.);
					extracap->cm[0] = 2. * s * s;
					extracap->np = 1;

					/* snap new boundary to poly */
					snap_poly(*poly, extracap);

					/* make sure new polygon contains enough space */
					np = (*poly)->np + 1;
					ier = room_poly(&polys[npoly], np, 0, 0);
					if (ier == -1) {
						fprintf(stderr, "partition_poly: failed to allocate memory for polygon of %d caps\n", np);
						return(-1);
					}

					/* combination of poly with new boundary */
					poly_poly(*poly, extracap, polys[npoly]);

					/* increment number of polygons */
					npoly++;

				}

				/* failed to find isolating boundary */
	    } else {
				/* vector from that to this */
				for (i = 0; i < 3; i++) {
					v[i] = vmax[i + 3 * ivthis] - vmin[i + 3 * ivthat];
				}
				s = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
				/* translate centre point along vector */
				/* the 0.02 puts centre just beyond equal distance */
				for (i = 0; i < 3; i++) {
					v[i] = vm[i + 3 * ivm] + ((cmthis - cmthat) / s + 0.02 * (it + 1)) 
						* v[i];
				}
				s = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
				for (i = 0; i < 3; i++) {
					vm[i + 3 * ivm] = v[i] / s;
				}

	    }

		} while (cmthis >= cmthat && nvm > 1 && it++ < itmax);

	}

	/* go with original polygon */
	if (npoly == 0 && *npoly_try == 1) {
		(*npoly_try)--;
		return(npoly);
	}

	/* if one or more polygons failed, add complement of all new caps to poly */
	if (npoly < *npoly_try) {
		/* not enough polygons for a new one */
		if (npoly >= npolys) return(*npoly_try + 1);

		/* put complement of all new caps into new polygon */
		np = npoly;
		ier = room_poly(&extracap, np, DNP, 0);
		if (ier == -1) {
	    fprintf(stderr, "partition_poly: failed to allocate memory for polygon of %d caps\n", np + DNP);
	    return(-1);
		}
		for (np = 0; np < npoly; np++) {
	    ip = polys[np]->np - 1;
	    for (i = 0; i < 3; i++) {
				extracap->rp_(i, np) = polys[np]->rp_(i, ip);
	    }
	    extracap->cm[np] = - polys[np]->cm[ip];
		}
		extracap->np = npoly;

		/* snap new caps to poly */
		snap_poly(*poly, extracap);

		/* make sure new polygon contains enough space */
		np = (*poly)->np + npoly;
		ier = room_poly(&polys[npoly], np, 0, 0);
		if (ier == -1) {
	    fprintf(stderr, "partition_poly: failed to allocate memory for polygon of %d caps\n", np);
	    return(-1);
		}

		/* poly with complement of new caps from other polygons */
		poly_poly(*poly, extracap, polys[npoly]);

		/* increment number of polygons */
		npoly++;
	}

	/* trim new polygons */
	for (ip = 0; ip < npoly; ip++) {
		trim_poly(polys[ip]);
	}

	/* copy final polygon part into poly */
	if (npoly > 0
			&& (overwrite_original == 2
					|| (overwrite_original == 1 && npoly == *npoly_try))) {
		/* make sure poly contains enough space */
		np = polys[npoly - 1]->np;
		ier = room_poly(poly, np, 0, 0);
		if (ier == -1) {
	    fprintf(stderr, "partition_poly: failed to allocate memory for polygon of %d caps\n", np);
	    return(-1);
		}
		copy_poly(polys[npoly - 1], *poly);

		/* decrement number of polygons */
		npoly--;
		(*npoly_try)--;
	}

	return(npoly);
}
