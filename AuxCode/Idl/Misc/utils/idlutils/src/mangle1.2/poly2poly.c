#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <values.h>
#include "format.h"
#include "polygon.h"
#include "scale.h"
#include "defaults.h"

/* getopt options */
static char *optstr = "dqm:j:k:ns:e:v:p:i:o:";

/* local functions */
int intersect_poly();
void usage(), parse_args();

/* external functions */
extern int rdmask(), wrmask();
extern int prune_poly();
extern polygon *new_poly();
extern void advise_fmt();
extern void msg(char *, ...);
extern void scale();
extern void free_poly();
extern void poly_poly();

polygon *polys[NPOLYSMAX];

/*------------------------------------------------------------------------------
  Main program.
*/
int main(argc, argv)
int argc;
char *argv [];
{
    int ifile, ipoly, nfiles, npoly, npolys;

    /* default output format */
    fmt.out = keywords[POLYGON];

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

    msg("---------------- poly2poly ----------------\n");

    /* tolerance angle for multiple intersections */
    if (mtol != 0.) {
	scale (&mtol, munit, 's');
	munit = 's';
	msg("multiple intersections closer than %g%c will be treated as coincident\n", mtol, munit);
	scale (&mtol, munit, 'r');
	munit = 'r';
    }

    /* weight limits */
    if (is_weight_min && is_weight_max) {
	/* min <= max */
	if (weight_min <= weight_max) {
	    msg("will keep only polygons with weights inside [%g, %g]\n", weight_min, weight_max);
	/* min > max */
	} else {
	    msg("will keep only polygons with weights >= %g or <= %g\n", weight_min, weight_max);
	    msg("         (only polygons with weights outside (%g, %g))\n", weight_max, weight_min);
	}
    } else if (is_weight_min) {
	msg("will keep only polygons with weights >= %g\n", weight_min);
    } else if (is_weight_max) {
	msg("will keep only polygons with weights <= %g\n", weight_max);
    }
    /* area limits */
    if (is_area_min && is_area_max) {
	/* min <= max */
	if (area_min < area_max) {
	    msg("will keep only polygons with areas inside [%g, %g]\n", area_min, area_max);
	/* min > max */
	} else {
	    msg("will keep only polygons with areas >= %g or <= %g\n", area_min, area_max);
	    msg("         (only polygons with areas outside (%g, %g))\n", area_max, area_min);
	}
    } else if (is_area_min) {
	msg("will keep only polygons with areas >= %g\n", area_min);
    } else if (is_area_max) {
	msg("will keep only polygons with areas <= %g\n", area_max);
    }

    /* advise data format */
    advise_fmt(&fmt);

    /* read polygons */
    npoly = 0;
    nfiles = argc - 1 - optind;
    for (ifile = optind; ifile < optind + nfiles; ifile++) {
	npolys = rdmask(argv[ifile], &fmt, &polys[npoly], NPOLYSMAX - npoly);
	if (npolys == -1) exit(1);
	/* intersect polygons of infile1 with those of subsequent infiles */
	if (ifile > optind && intersect) {
	    npoly = intersect_poly(polys, npoly, &polys[npoly], npolys);
	    if (npoly == -1) exit(1);
	/* increment number of polygons */
	} else {
	    npoly += npolys;
	}
    }
    if (nfiles >= 2 && !intersect) {
        msg("total of %d polygons read\n", npoly);
    }
    if (npoly == 0) {
	msg("stop\n");
	exit(0);
    }

    /* apply new id numbers to output polygons */
    if (fmt.newid == 'n') {
	for (ipoly = 0; ipoly < npoly; ipoly++) {
	    polys[ipoly]->id = ipoly;
	}
    }

    /* write polygons */
    ifile = argc - 1;
    npoly = wrmask(argv[ifile], &fmt, polys, npoly);
    if (npoly == -1) exit(1);

    exit(0);
}

/*------------------------------------------------------------------------------
*/
void usage()
{
    printf("usage:\n");
    printf("poly2poly [-d] [-q] [-j[<min>][,<max>]] [-k[min][,<max>]] [-n] [-s<n>] [-e<n>] [-vo|-vn] [-p<n>] [-i<f>[<n>][u]] [-o<f>[u]] infile1 [infile2] ... outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Intersect polygons of poly1 with any polygon(s) of poly2 having the
  same id number.
*/
int intersect_poly(poly1, npoly1, poly2, npoly2)
int npoly1, npoly2;
#ifdef GCC
polygon *poly1[npoly1], *poly2[npoly2];
#else
polygon *poly1[], *poly2[];
#endif
{
    int inull, iprune, i, j, np;
    polygon *poly;

    /* intersect each poly1 with any poly2 having same id number */
    for (i = 0; i < npoly1; i++) {
	for (j = 0; j < npoly2; j++) {
	    if (poly1[i]->id == poly2[j]->id) {

		/* make sure poly contains enough space for intersection */
		np = poly1[i]->np + poly2[j]->np;
		poly = new_poly(np);
		if (!poly) goto out_of_memory;

		/* intersection of poly1 and poly2 */
		poly_poly(poly1[i], poly2[j], poly);

		/* free poly1 */
		free_poly(poly1[i]);

		/* point poly1 at intersection */
		poly1[i] = poly;
	    }
	}
    }

    /* free poly2 polygons */
    for (j = 0; j < npoly2; j++) {
        free_poly(poly2[j]);
	poly2[j] = 0x0;
    }

    /* prune poly1 polygons */
    j = 0;
    inull = 0;
    for (i = 0; i < npoly1; i++) {
        iprune = prune_poly(poly1[i]);
        if (iprune == -1) return(-1);
	if (iprune >= 2) {
	    free_poly(poly1[i]);
	    inull++;
	} else {
	    poly1[j] = poly1[i];
	    j++;
	}
    }
    if (inull > 0) msg("%d intersected polygons have zero area, and are being discarded\n", inull);
    npoly1 = j;

    return(npoly1);

    /* error returns */
    out_of_memory:
    fprintf(stderr, "intersect_poly: failed to allocate memory for polygon of %d caps\n", np);
    return(-1);
}
