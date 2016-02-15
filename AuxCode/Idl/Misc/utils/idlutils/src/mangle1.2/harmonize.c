#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <values.h>
#include "format.h"
#include "harmonics.h"
#include "polygon.h"
#include "scale.h"
#include "defaults.h"

/* advise how many polygons done if lmax >= this */
#define LMAX_ADVICE		250

/* getopt options */
static char *optstr = "dql:m:s:e:i:";

/* local functions */
int harmonize();
void usage(), parse_args();

/* external functions */
extern int rdmask();
extern int poly_to_rect();
extern int wrspher();
extern int gsphr();
extern int gsphra();
extern void advise_fmt();
extern void findtop();
extern void msg(char *fmt, ...);
extern void scale();

polygon *polys[NPOLYSMAX];

/*------------------------------------------------------------------------------
  Main program.
*/
int main(argc, argv)
int argc;
char *argv [];
{
    int ifile, nfiles, npoly, npolys, nws;
    double area, *w;

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

    msg("---------------- harmonize ----------------\n");

    /* advise harmonic number */
    msg("maximum harmonic number %d\n", lmax);

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

    /* spherical harmonics of region */
    npoly = harmonize(polys, npoly, lmax, &w);
    if (npoly == -1) exit(1);

    /* advise area */
    area = w[0] * 2. * sqrt(PI);
    msg("area of (weighted) region is %.15lg str\n", area);

    /* write polygons */
    ifile = argc - 1;
    nws = wrspher(argv[ifile], lmax, w);
    if (nws == -1) exit(1);

    exit(0);
}

/*------------------------------------------------------------------------------
*/
void usage()
{
    printf("usage:\n");
    printf("harmonize [-d] [-q] [-l<n>] [-s<n>] [-e<n>] [-i<f>[<n>][u]] infile1 [infile2] ... outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Spherical harmonics of sum of weighted polygons.

  Return value: number of polygons for which spherical harmonics were computed,
		or -1 if error occurred.
*/
int harmonize(poly, npoly, lmax, w_p)
int npoly, lmax;
#ifdef GCC
polygon *poly[npoly];
#else
polygon *poly[];
#endif
double **w_p;
{

    /* array contains cumulative spherical harmonics of entire region */
    static double *w;

    int accelerate, i, ier, ipoly, iq, ir, isrect, iw, naccelerate, ndone, ner, nrect;
    double azmin, azmax, elmin, elmax, azmn, azmx, elmn, elmx, tol;
    /* work array contains harmonics of single polygon */
    double *dw;
    /* work arrays to deal with possible acceleration */
    int *iord, *ir_to_ip;
    double *elord;

    /* allocate array containing cumulative spherical harmonics of entire region */
    w = (double *) malloc(sizeof(double) * IM * NW);
    if (!w) {
	fprintf(stderr, "harmonize: failed to allocate memory for %d harmonics\n", IM * NW);
	return(-1);
    }

    /* work arrays */
    dw = (double *) malloc(sizeof(double) * IM * NW);
    if (!dw) {
	fprintf(stderr, "harmonize: failed to allocate memory for %d harmonics\n", IM * NW);
	return(-1);
    }
    iord = (int *) malloc(sizeof(int) * npoly);
    if (!iord) {
	fprintf(stderr, "harmonize: failed to allocate memory for %d integers\n", npoly);
	return(-1);
    }
    ir_to_ip = (int *) malloc(sizeof(int) * npoly);
    if (!ir_to_ip) {
	fprintf(stderr, "harmonize: failed to allocate memory for %d integers\n", npoly);
	return(-1);
    }
    elord = (double *) malloc(sizeof(double) * npoly);
    if (!elord) {
	fprintf(stderr, "harmonize: failed to allocate memory for %d doubles\n", npoly);
	return(-1);
    }

    /* zero harmonics of region */
    for (iw = 0; iw < IM * NW; iw++) {
	w[iw] = 0;
    }

    /* determine which polygons are rectangles, for which acceleration may be possible */
    nrect = 0;
    ir = npoly;
    for (ipoly = 0; ipoly < npoly; ipoly++) {
	isrect = poly_to_rect(poly[ipoly], &azmin, &azmax, &elmin, &elmax);
	if (isrect && poly[ipoly]->weight != 0.) {
	    ir_to_ip[nrect] = ipoly;
	    elord[nrect] = elmin * 1.e8 + elmax;
	    nrect++;
	} else {
	    ir--;
	    ir_to_ip[ir] = ipoly;
	}
    }
    msg("harmonize: %d polygons are rectangles, for which acceleration may be possible\n",
	nrect);

    /* order rectangles by elmin, elmax */
    findtop(elord, nrect, iord, nrect);
    /* convert from fortran to c convention */
    for (ir = 0; ir < nrect; ir++) iord[ir]--;

    /* do each polygon */
    ndone = 0;
    naccelerate = 0;
    ner = 0;
    if (lmax >= LMAX_ADVICE) msg("doing polygon number (of %d):\n", npoly);
    for (i = 0; i < npoly; i++) {
	if (lmax >= LMAX_ADVICE) msg(" %d", i);
	accelerate = 0;
	/* rectangle */
	if (i < nrect) {
	    ir = iord[i];
	    ipoly = ir_to_ip[ir];
	    poly_to_rect(poly[ipoly], &azmin, &azmax, &elmin, &elmax);
	    /* does previous rectangle have same elevation limits? */
	    if (i > 0) {
		iq = iord[i - 1];
		iq = ir_to_ip[iq];
		poly_to_rect(poly[iq], &azmn, &azmx, &elmn, &elmx);
		/* if so, use acceleration */
		if (elmn == elmin && elmx == elmax) accelerate = 1;
	    }
	    /* if not, does next rectangle have same elevation limits? */
	    if (!accelerate && i + 1 < nrect) {
		iq = iord[i + 1];
		iq = ir_to_ip[iq];
		poly_to_rect(poly[iq], &azmn, &azmx, &elmn, &elmx);
		/* if so, worth accelerating */
		if (elmn == elmin && elmx == elmax) accelerate = 1;
	    }
	    /* accelerated computation */
	    if (accelerate) {
		ier = gsphra(azmin, azmax, elmin, elmax, lmax, dw);
		if (ier == -1) return(-1);
	    /* standard computation */
	    } else {
		tol = mtol;
		ier = gsphr(poly[ipoly], lmax, &tol, dw);
		if (ier == -1) return(-1);
	    }
	/* non-rectangle */
	} else {
	    ipoly = ir_to_ip[i];
	    /* zero weight polygon requires no computation */
	    if (poly[ipoly]->weight == 0.) {
		ndone++;
		continue;
	    } else {
		tol = mtol;
		ier = gsphr(poly[ipoly], lmax, &tol, dw);
		if (ier == -1) return(-1);
	    }
	}
	/* computation failed */
	if (ier) {
	    ner++;
	    if (lmax >= LMAX_ADVICE) msg("\n");
	    fprintf(stderr, "computation failed for polygon %d; discard it\n", i);
	/* success */
	} else {
	    naccelerate += accelerate;
	    ndone++;
	    /* increment harmonics of region */
	    for (iw = 0; iw < IM * NW; iw++) {
		w[iw] += dw[iw] * poly[ipoly]->weight;
	    }
	}
    }
    if (lmax >= LMAX_ADVICE) msg("\n");

    /* number of computations that were accelerated */
    msg("harmonize: computation was accelerated for %d rectangles\n", naccelerate);
    /* advise */
    if (ner > 0) {
	msg("discarded %d polygons for which computations failed\n", ner);
    }
    msg("spherical harmonics of %d weighted polygons accumulated\n", ndone);

    /* point w_p at spherical harmonics */
    *w_p = w;

    /* free work arrays */
    free(dw);
    free(iord);
    free(ir_to_ip);
    free(elord);

    return(npoly);
}
