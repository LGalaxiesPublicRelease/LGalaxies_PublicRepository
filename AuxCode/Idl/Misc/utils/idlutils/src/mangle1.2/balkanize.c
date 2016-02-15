#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <values.h>
#include "format.h"
#include "polygon.h"
#include "scale.h"
#include "defaults.h"

/* number of extra caps to allocate to polygon, to allow for expansion */
#define DNP		4

/* getopt options */
static char *optstr = "dqm:s:e:v:p:i:o:";

/* local functions */
int split_poly(), fragment_poly(), balkanize();
void usage(), parse_args();

/* external functions */
extern int partition_poly();
extern int prune_poly(), trim_poly();
extern int rdmask(), wrmask();
extern int room_poly();
extern int garea();
extern void advise_fmt();
extern void copy_poly(), poly_poly(), poly_polyn();
extern void msg(char *, ...);
extern void scale();
extern void memmsg();

polygon *polys[NPOLYSMAX];

/*------------------------------------------------------------------------------
  Main program.
*/
int main(argc, argv)
		 int argc;
		 char *argv[];
{
	int ifile, nfiles, npoly, npolys;

	/* default output format */
	fmt.out = keywords[POLYGON];
	/* default is to renumber output polygons with new id numbers */
	fmt.newid = 'n';

	/* parse arguments */
	parse_args(argc, argv);

	/* at least one input and output filename required as arguments */
	if (argc - optind < 2) {
		if (optind > 1 || argc - optind == 1) {
	    fprintf(stderr, "%s requires at least 2 arguments: infile and outfile\n", argv[0]);
	    usage();
	    exit(1);
		} else {
	    usage();
	    exit(0);
		}
	}

	msg("---------------- balkanize ----------------\n");

	/* tolerance angle for multiple intersections */
	if (mtol != 0.) {
		scale (&mtol, munit, 's');
		munit = 's';
		msg("multiple intersections closer than %g%c will be treated as coincident\n", mtol, munit);
		scale (&mtol, munit, 'r');
		munit = 'r';
	}

	/* advise data format */
	advise_fmt(&fmt);

	/* read polygons */
	npoly = 0;
	nfiles = argc - 1 - optind;
	for (ifile = optind; ifile < optind + nfiles; ifile++) {
		npolys = rdmask(argv[ifile], &fmt, &polys[npoly], NPOLYSMAX - npoly);
		if (npolys == -1) exit(1);
		npoly += npolys;
	}
	if (nfiles >= 2) {
		msg("total of %d polygons read\n", npoly);
	}
	if (npoly == 0) {
		msg("stop\n");
		exit(0);
	}

	/* balkanize polygons */
	npolys = balkanize(polys, npoly, &polys[npoly], NPOLYSMAX - npoly);
	if (npolys == -1) exit(1);

	/* write polygons */
	ifile = argc - 1;
	npolys = wrmask(argv[ifile], &fmt, &polys[npoly], npolys);
	if (npolys == -1) exit(1);
	/* memmsg(); */

	exit(0);
}

/*------------------------------------------------------------------------------
*/
void usage()
{
	printf("usage:\n");
	printf("balkanize [-d] [-q] [-s<n>] [-e<n>] [-vo|-vn] [-p<n>] [-i<f>[<n>][u]] [-o<f>[u]] infile1 [infile2] ... outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  If poly1 overlaps poly2, split poly1 into two parts.

   Input: *poly1, poly2 are 2 polygons.
  Output: *poly1 and *poly3 are 2 split polygons of poly1, if poly1 is split,
		with *poly1 the part outside poly2,
		and *poly3 the part intersecting poly2;
		*poly1 and *poly3 remain untouched if poly1 is not split.
  Return value: -1 = error occurred;
		0 = poly1 and poly2 have zero intersection;
		1 = poly2 contains poly1;
		2 = poly1 split into 2.
*/
int split_poly(poly1, poly2, poly3)
		 polygon **poly1, *poly2, **poly3;
{
	static polygon *poly = 0x0;

	int i, ier, ip, itrim, np, verb;
	double area, area_tot, cm, tol;

	/* poly2 is whole sphere, therefore contains poly1 */
	if (poly2->np == 0) return(1);

	/* make sure poly contains enough space for intersection */
	np = (*poly1)->np + poly2->np;
	ier = room_poly(&poly, np, DNP, 0);
	if (ier == -1) goto out_of_memory;

	/* intersection of poly1 and poly2 */
	poly_poly(*poly1, poly2, poly);

	/* suppress coincident boundaries, to make garea happy */
	itrim = trim_poly(poly);

	/* intersection of poly1 and poly2 is null polygon */
	if (itrim >= 2) return(0);

	/* area of intersection */
	tol = mtol;
	verb = 1;
	ier = garea(poly, &tol, &verb, &area_tot);
	if (ier) goto error;

	/* poly1 and poly2 have zero intersection */
	if (area_tot == 0.) return(0);

	/* find boundary of poly2 which intersects poly1 */
	verb = 0;
	np = (*poly1)->np;
	for (ip = 0; ip < poly2->np; ip++) {
		cm = poly->cm[np + ip];
		poly->cm[np + ip] = 2.;		/* suppress boundary to be tested */
		ier = garea(poly, &tol, &verb, &area);	/* area of intersection sans boundary */
		poly->cm[np + ip] = cm;		/* restore tested boundary */
		if (area != area_tot) {		/* boundary intersects poly1 */

	    /* make sure poly3 contains enough space */
	    np = (*poly1)->np + 1;
	    ier = room_poly(poly3, np, DNP, 0);
	    if (ier == -1) goto out_of_memory;

	    /* poly3 is intersection of poly1 and ip'th cap of poly2 */
	    poly_polyn(*poly1, poly2, ip, *poly3);

	    /* make sure poly1 contains enough space */
	    np = (*poly1)->np + 1;
	    ier = room_poly(poly1, np, DNP, 1);
	    if (ier == -1) goto out_of_memory;

	    /* poly1 is intersection of poly1 and complement of ip'th cap of poly2 */
	    np = (*poly1)->np;
	    for (i = 0; i < 3; i++) {
				(*poly1)->rp_(i, np) = poly2->rp_(i, ip);
	    }
	    (*poly1)->cm[np] = - poly2->cm[ip];
	    (*poly1)->np++;

	    /* prune poly1 */
	    itrim = prune_poly(*poly1);
	    if (itrim == -1) goto error;
	    /* poly1 may be null because of numerical roundoff */
	    if (itrim >= 2) {
				/* restore an original copy of poly1 from poly3 */
				copy_poly(*poly3, *poly1);
				/* cut the extra cap */
				(*poly1)->np--;
				/* skip to the next cap in the for loop */
				continue;
	    }

	    /* prune poly3 */
	    itrim = prune_poly(*poly3);
	    if (itrim == -1) goto error;
	    /* poly3 is unaccountably null */
	    if (itrim >= 2) goto null_poly3;

	    /* poly1 successfully split into poly1 and poly3 */
	    return(2);
		}
	}

	/* poly2 contains poly1 */
	return(1);

	/* ---------------- error returns ---------------- */
 error:
	return(-1);

 out_of_memory:
	fprintf(stderr, "split_poly: failed to allocate memory for polygon of %d caps\n", np + DNP);
	return(-1);

 null_poly3:
	fprintf(stderr, "split_poly: poly3 is null, which should not happen\n");
	return(-1);
}

/*------------------------------------------------------------------------------
  Fragment poly1 into several disjoint polygons,
  each of which is either wholly outside or wholly inside poly2.

   Input: *poly1, poly2 are 2 polygons.
	  discard = 0 to retains all parts of poly1;
	  	  = 1 to discard intersection of poly1 with poly2.
	  npolys = maximum number of polygons available in polys array.
  Output: *poly1 and polys[npoly] are disjoint polygons of poly1;
		all but the last polygon lie outside poly2;
		if discard = 0:
		    if poly1 intersects poly2, then the last polygon,
		    polys[npoly - 1] (or *poly1 if npoly = 0),
		    is the intersection of poly1 and poly2;
		if discard = 1:
		    if poly1 intersects poly2, then the last+1 polygon,
		    polys[npoly],
		    is the discarded intersection of poly1 and poly2;
		    if poly1 lies entirely inside poly2 (so npoly = 0),
		    then *poly1 is set to null.
  Return value: npoly = number of disjoint polygons, excluding poly1,
			or -1 if error occurred in split_poly().
*/
int fragment_poly(poly1, poly2, discard, polys, npolys)
		 int discard, npolys;
		 polygon **poly1, *poly2;
#ifdef GCC
		 polygon *polys[npolys];
#else
		 polygon *polys[];
#endif
{
	int npoly, nsplit;
	polygon **poly;

	/* iteratively subdivide polygons of poly1 */
	npoly = 0;
	poly = poly1;
	while (1) {
		/* check space is available */
		if (npoly >= npolys) return(npoly + 1);
		/* split */
		nsplit = split_poly(poly, poly2, &polys[npoly]);
		/* error */
		if (nsplit == -1) return(-1);
		/* done */
		if (nsplit == 0 || nsplit == 1) {
	    if (nsplit == 1 && discard) {
				if (npoly > 0) {
					npoly--;
				} else {
					*poly = 0x0;
				}
	    }
	    return(npoly);
		}
		poly = &polys[npoly++];
	}

}

/*------------------------------------------------------------------------------
  Balkanize overlapping polygons into many disjoint connected polygons.

   Input: poly = array of pointers to polygons.
	  npoly = number of polygons.
	  npolys = maximum number of output polygons.
  Output: polys = array of pointers to polygons.
  Return value: number of disjoint connected polygons,
		or -1 if error occurred.
*/
int balkanize(poly, npoly, polys, npolys)
		 int npoly, npolys;
#ifdef GCC
		 polygon *poly[npoly], *polys[npolys];
#else
		 polygon *poly[], *polys[];
#endif
{
	/* partition_polly should not lassoo all one-border polygons */
#define ALL_ONEBORDER		0
	/* how partition_poly should tighten lassoo */
#define ADJUST_LASSOO		1
	/* partition_poly should overwrite all original polygons */
#define OVERWRITE_ORIGINAL	2
#define WARNMAX			8
	int discard, dm, dn, dnp, failed, i, ier, inull, ip, iprune, j, k, m, n, np, npoly_try;

    /* start by pruning all input polygons */
	np = 0;
	inull = 0;
	for (i = 0; i < npoly; i++) {
		iprune = prune_poly(poly[i]);
		/* error */
		if (iprune == -1) return(-1);
		/* zero area polygon */
		if (iprune >= 2) {
	    if (WARNMAX > 0 && inull == 0) msg("warning from balkanize: following polygons have zero area & are being discarded:\n");
	    if (inull < WARNMAX) {
				msg(" %d", (fmt.newid == 'o')? poly[i]->id : i);
	    } else if (inull == WARNMAX) {
				msg(" ... more\n");
	    }
	    inull++;
		} else {
	    np++;
		}
	}
	if (WARNMAX > 0 && inull > 0 && inull <= WARNMAX) msg("\n");
	if (inull > 0) {
		msg("balkanize: %d polygons with zero area are being discarded;\n", inull);
	}

	/* number of polygons */
	msg("balkanizing %d polygons ...\n", np);

	/* nullify all output polygons */
	for (i = npoly; i < npolys; i++) {
		polys[i] = 0x0;
	}

	/*
      m = starting index of current set of fragments of i'th polygon
      dm = number of current set of fragments of i'th polygon
      n = starting index of new subset of fragments of i'th polygon
      dn = number of new subset of fragments of i'th polygon
    */

	msg("balkanize stage 1 (fragment into non-overlapping polygons):\n");
	n = 0;
	dnp = 0;
	ip = 0;
	/* fragment each polygon in turn */
	for (i = 0; i < npoly; i++) {
		/* skip null polygons */
		if (poly[i]->np > 0 && poly[i]->cm[0] == 0.) continue;
		/* update indices */
		m = n;
		dm = 1;
		n = m + dm;
		/* make sure output polygon has enough room */
		ier = room_poly(&polys[m], poly[i]->np, DNP, 0);
		if (ier == -1) {
	    fprintf(stderr, "balkanize: failed to allocate memory for polygon of %d caps\n", poly[i]->np + DNP);
	    return(-1);
		}
		/* copy polygon i into output polygon */
		copy_poly(poly[i], polys[m]);

		/* fragment successively against other polygons */
		for (j = 0; j < npoly; j++) {
	    /* skip self, or null polygons */
	    if (j == i || (poly[j]->np > 0 && poly[j]->cm[0] == 0.)) continue;
	    /* keep only one copy of the intersection of i & j */
	    /* intersection inherits weight of polygon being fragmented,
	       so keeping later polygon ensures intersection inherits
	       weight of later polygon */
	    if (i < j) {
				discard = 1;
	    } else {
				discard = 0;
	    }
	    /* fragment each part of i'th polygon */
	    for (k = m; k < m + dm; k++) {
				/* skip null polygons */
				if (!polys[k] || (polys[k]->np > 0 && polys[k]->cm[0] == 0.)) continue;
				/* fragment */
				dn = fragment_poly(&polys[k], poly[j], discard, &polys[n], npolys - n);
				/* error */
				if (dn == -1) {
					fprintf(stderr, "UHOH at polygon %d; continuing\n", (fmt.newid == 'o')? polys[i]->id : ip);
					continue;
					/* return(-1); */
				}
				/* increment index of next subset of fragments */
				n += dn;
				/* increment polygon count */
				np += dn;
				dnp += dn;
				if (!polys[k]) {
					np--;
					dnp--;
				}
				/* check whether exceeded maximum number of polygons */
				if (n > npolys) {
					fprintf(stderr, "balkanize: total number of polygons exceeded maximum %d\n", npoly + npolys);
					fprintf(stderr, "if you need more space, enlarge NPOLYSMAX in defaults.h, and recompile\n");
					return(-1);
				}
	    }
	    dm = n - m;
		}
		ip++;
	}
	msg("added %d polygons to make %d\n", dnp, np);

    /* partition disconnected polygons into connected parts  */
	msg("balkanize stage 2 (partition disconnected polygons into connected parts):\n");
	m = n;
	dnp = 0;
	ip = 0;
	failed = 0;
	for (i = 0; i < m; i++) {
		/* skip null polygons */
		if (!polys[i] || (polys[i]->np > 0 && polys[i]->cm[0] == 0.)) continue;
		/* partition disconnected polygons */
		dn = partition_poly(&polys[i], &polys[n], npolys - n, ALL_ONEBORDER, ADJUST_LASSOO, OVERWRITE_ORIGINAL, &npoly_try);
		/* error */
		if (dn == -1) {
	    fprintf(stderr, "UHOH at polygon %d; continuing\n", (fmt.newid == 'o')? polys[i]->id : ip);
	    continue;
	    /* return(-1); */
		}
		/* failed to partition polygon into desired number of parts */
		if (dn < npoly_try) {
	    fprintf(stderr, "failed to split polygon %d into %d parts\n", (fmt.newid == 'o')? polys[i]->id : ip, npoly_try - dn + 1);
	    failed++;
		}
		/* increment index of next subset of fragments */
		n += dn;
		/* increment polygon count */
		np += dn;
		dnp += dn;
		/* check whether exceeded maximum number of polygons */
		if (n > npolys) {
	    fprintf(stderr, "balkanize: total number of polygons exceeded maximum %d\n", npoly + npolys);
	    fprintf(stderr, "if you need more space, enlarge NPOLYSMAX in defaults.h, and recompile\n");
	    return(-1);
		}
		ip++;
	}
	msg("added %d polygons to make %d\n", dnp, np);
	if (failed > 0) {
		msg("balkanize: failed to split %d polygons into desired number of connected parts\n", failed);
		msg(".............................................................................\n");
		msg("Failure to split polygon probably means:\n");
		msg("either (1) you forgot to run snap on all your input polygon files;\n");
		msg("    or (2) the polygon is too small for the numerics to cope with;\n");
		msg("    or (3) you have a weird-shaped polygon.\n");
		msg("You may ignore this warning message if the weights of polygons in the input\n");
		msg("polygon file(s) are already correct, and you do not want to reweight them.\n");
		msg("Similarly, you may ignore this warning message if you do want to reweight the\n");
		msg("polygons, but the weights of the different parts of each unsplit polygon are\n");
		msg("the same.  If you want to reweight the different parts of an unsplit polygon\n");
		msg("with different weights, then you will need to split that polygon by hand.\n");
		msg("Note that balkanize may think that a polygon is disconnected when in fact it\n");
		msg("is just not simply-connected; in this case it is neither necessary nor\n");
		msg("possible to split the polygon, and again you may ignore this message.\n");
		msg("Whatever the case, the output file of balkanized polygons constitutes\n");
		msg("a valid mask of non-overlapping polygons, which is safe to use.\n");
		msg(".............................................................................\n");
	}

	/* final prune */
	j = 0;
	inull = 0;
	for (i = 0; i < n; i++) {
		if (!polys[i]) continue;
		iprune = prune_poly(polys[i]);
		if (iprune == -1) {
	    fprintf(stderr, "failed to prune polygon %d; continuing ...\n", (fmt.newid == 'o')? polys[i]->id : j);
	    /* return(-1); */
		}
		if (iprune >= 2) {
	    inull++;
		} else {
	    polys[j] = polys[i];
	    j++;
		}
	}
	if (inull > 0) msg("balkanize: %d balkanized polygons have zero area, and are being discarded\n", inull);
	n = j;

	/* assign new polygon id numbers in place of inherited ids */
	if (fmt.newid == 'n') {
		for (i = 0; i < n; i++) {
	    polys[i]->id = i;
		}
	}

	if (n != -1) msg("balkanize: balkans contain %d polygons\n", n);

	return(n);
}
