#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <values.h>
#include "format.h"
#include "polygon.h"
#include "scale.h"
#include "vertices.h"
#include "defaults.h"

/* getopt options */
static char *optstr = "dqz:m:s:e:v:p:i:o:";

/* local functions */
int weight();
void usage(), parse_args();

/* external functions */
extern int rdmask(), wrmask();
extern int vmid();
extern int gverts();
extern double weight_fn();
extern void advise_fmt();
extern void msg(char *, ...);
extern void rp_to_vert();
extern void scale(), scale_azel();

/*------------------------------------------------------------------------------
  Main program.
*/
int main(argc, argv)
int argc;
char *argv [];
{
    int ifile, nfiles, npoly, npolys;
    polygon *polys[NPOLYSMAX];

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

    /* survey must have been specified */
    if (!survey) {
	fprintf(stderr, "%s requires -z<survey> option to specify name of survey\n", argv[0]);
	fprintf(stderr, "please look in weight_fn.c for the names of known surveys\n");
	exit(1);
    }

    msg("---------------- weight ----------------\n");

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

    /* weight polygons */
    npoly = weight(polys, npoly, survey);

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
    printf("weight [-d] [-q] -z<survey> [-s<n>] [-e<n>] [-vo|-vn] [-p<n>] [-i<f>[<n>][u]] [-o<f>[u]] infile1 [infile2] ... outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Weight polygons.

   Input: poly = array of pointers to polygons.
	  npoly = pointer to number of polygons.
	  survey = name of survey.
  Output: polys = array of pointers to polygons;
  Return value: number of polygons weighted,
		or -1 if error occurred.
*/
int weight(poly, npoly, survey)
int npoly;
#ifdef GCC
polygon *poly[npoly];
#else
polygon *poly[];
#endif
char *survey;
{
    static int nve = 2;

    int do_vcirc, i, imid, ipoly, iverts, ivm, nev, nev0, nomid, nv, nvm, nzero;
    int *ipv, *ev;
    double tol;
    double *ve, *vm, *angle;
    azel v;

    nomid = 0;
    nzero = 0;
    for (ipoly = 0; ipoly < npoly; ipoly++) {
	/* vertices of polygon */
	do_vcirc = 0;
	tol = mtol;
	iverts = gverts(poly[ipoly], do_vcirc, &tol, &nv, nve, &ve, &angle, &ipv, &nev, &nev0, &ev);
	if (iverts != 0) return(-1);
	/* point somewhere in the middle of the polygon */
	imid = vmid(poly[ipoly], nv, nve, ve, ev, &nvm, &vm);
	if (imid == -1) return(-1);
	/* check found point inside the polygon */
	imid = 0;
	for (ivm = 0; ivm < nvm; ivm++) {
	    if (vm[3 * ivm] != 0. || vm[1 + 3 * ivm] != 0. || vm[2 + 3 * ivm] != 0.) {
		imid = 1;
		if (ivm > 0) for (i = 0; i < 3; i++) vm[i] = vm[i + 3 * ivm];
		break;
	    }
	}
	/* found point */
        if (imid == 1) {
	    /* convert unit vector to az, el */
	    rp_to_vert(vm, &v);
	    /* scale angles from radians to degrees */
	    scale_azel(&v, 'r', 'd');
	    /* weight at that point */
	    poly[ipoly]->weight = weight_fn(v.az, v.el, survey);
	    if (poly[ipoly]->weight == 0.) nzero++;
	/* failed to find point */
	} else {
	    if (nomid == 0) msg("weight: failed to find interior point for the following polygons:\n");
	    msg(" %d", (fmt.newid == 'o')? ipoly : poly[ipoly]->id);
	    nomid++;
	}
    }
    if (nomid > 0) msg("\n");
    if (nomid > 0) {
	msg("weight: failed to find interior point for %d polygons\n", nomid);
        msg("FAILURE TO FIND INTERIOR POINT PROBABLY MEANS YOU HAVE A WEIRD-SHAPED POLYGON.\n");
	msg("PLEASE FILL IN THE CORRECT WEIGHT BY HAND.\n");
    }

    /* assign new polygon id numbers in place of inherited ids */
    if (fmt.newid == 'n') {
	for (ipoly = 0; ipoly < npoly; ipoly++) {
	    poly[ipoly]->id = ipoly;
	}
    }

    /* warn about zero weights */
    if (nzero > 0) msg("weight: %d polygons have zero weight\n", nzero);

    return(npoly);
}
