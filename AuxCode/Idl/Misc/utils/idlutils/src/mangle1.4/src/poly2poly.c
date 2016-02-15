/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "manglefn.h"
#include "defaults.h"

polygon *polys[NPOLYSMAX];

/* getopt options */
const char *optstr = "dqm:j:k:ns:e:v:p:i:o:";

/* local functions */
void	usage(void);
#ifdef	GCC
int	intersect_poly(int npoly1, polygon *[npoly1], int npoly2, polygon *[npoly2], double);
#else
int	intersect_poly(int npoly1, polygon *[/*npoly1*/], int npoly2, polygon *[/*npoly2*/], double);
#endif

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
    int ifile, ipoly, nfiles, npoly, npolys;

    /* default output format */
    fmt.out = keywords[POLYGON];

    /* parse arguments */
    parse_args(argc, argv);

    /* at least one input and output filename required as arguments */
    if (argc - optind < 2) {
	if (optind > 1 || argc - optind == 1) {
	    fprintf(stderr, "%s requires at least 2 arguments: polygon_infile and polygon_outfile\n", argv[0]);
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
	npolys = rdmask(argv[ifile], &fmt, NPOLYSMAX - npoly, &polys[npoly]);
	if (npolys == -1) exit(1);
	/* intersect polygons of infile1 with those of subsequent infiles */
	if (ifile > optind && intersect) {
	    npoly = intersect_poly(npoly, polys, npolys, &polys[npoly], mtol);
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
	msg("STOP\n");
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
    npoly = wrmask(argv[ifile], &fmt, npoly, polys,0);
    if (npoly == -1) exit(1);

    return(0);
}

/*------------------------------------------------------------------------------
*/
void usage(void)
{
    printf("usage:\n");
    printf("poly2poly [-d] [-q] [-m<a>[u]] [-j[<min>][,<max>]] [-k[min][,<max>]] [-n] [-s<n>] [-e<n>] [-vo|-vn] [-p[+|-][<n>]] [-i<f>[<n>][u]] [-o<f>[u]] polygon_infile1 [polygon_infile2 ...] polygon_outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Intersect polygons of poly1 with any polygon(s) of poly2 having the
  same id number.

  This subroutine implements the -n option of poly2poly.
*/
int intersect_poly(int npoly1, polygon *poly1[/*npoly1*/], int npoly2, polygon *poly2[/*npoly2*/], double mtol)
{
    int ier, inull, iprune, i, j, np;

    /* intersect each poly1 with any poly2 having same id number */
    for (i = 0; i < npoly1; i++) {
	for (j = 0; j < npoly2; j++) {
	    if (poly1[i]->id == poly2[j]->id) {

		/* make sure poly1 contains enough space for intersection */
		np = poly1[i]->np + poly2[j]->np;
		ier = room_poly(&poly1[i], np, 0, 1);
		if (ier == -1) goto out_of_memory;

		/* intersection of poly1 and poly2 */
		poly_poly(poly1[i], poly2[j], poly1[i]);

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
        iprune = prune_poly(poly1[i], mtol);
        if (iprune == -1) {
	    fprintf(stderr, "intersect_poly: failed to prune polygon %d; continuing ...\n", (fmt.newid == 'o')? poly1[i]->id : j);
	}
	if (iprune >= 2) {
	    free_poly(poly1[i]);
	    poly1[i] = 0x0;
	    inull++;
	} else {
	    poly1[j] = poly1[i];
	    j++;
	}
    }
    if (inull > 0) msg("%d intersected polygons have zero area, and are being discarded\n", inull);
    npoly1 = j;

    return(npoly1);

    /* ---------------- error returns ---------------- */
    out_of_memory:
    fprintf(stderr, "intersect_poly: failed to allocate memory for polygon of %d caps\n", np);
    return(-1);
}
