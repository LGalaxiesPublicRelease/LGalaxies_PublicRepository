/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "manglefn.h"
#include "defaults.h"

/* getopt options */
const char *optstr = "dql:m:s:e:i:";

/* local functions */
void	usage(void);

polygon *polys[NPOLYSMAX];

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
    int ifile, nfiles, npoly, npolys, nws;
    double area;
    harmonic *w;

    /* parse arguments */
    parse_args(argc, argv);

    /* at least one input and output filename required as arguments */
    if (argc - optind < 2) {
	if (optind > 1 || argc - optind == 1) {
	    fprintf(stderr, "%s requires at least 2 arguments: polygon_infile, and Wlm_outfile\n", argv[0]);
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
	scale(&mtol, munit, 's');
	munit = 's';
	msg("multiple intersections closer than %g%c will be treated as coincident\n", mtol, munit);
	scale(&mtol, munit, 'r');
	munit = 'r';
    }

    /* advise data format */
    advise_fmt(&fmt);

    /* read polygons */
    npoly = 0;
    nfiles = argc - 1 - optind;
    for (ifile = optind; ifile < optind + nfiles; ifile++) {
	npolys = rdmask(argv[ifile], &fmt, NPOLYSMAX - npoly, &polys[npoly]);
	if (npolys == -1) exit(1);
	npoly += npolys;
    }
    if (nfiles >= 2) {
	msg("total of %d polygons read\n", npoly);
    }
    if (npoly == 0) {
	msg("STOP\n");
	exit(0);
    }

    /* allocate array containing spherical harmonics of complete mask */
    w = (harmonic *) malloc(sizeof(harmonic) * NW);
    if (!w) {
        fprintf(stderr, "harmonize: failed to allocate memory for %d harmonics\n", NW);
        exit(1);
    }

    /* spherical harmonics of region */
    npoly = harmonize_polys(npoly, polys, mtol, lmax, w);
    if (npoly == -1) exit(1);

    /* advise area */
    area = w[0][0] * 2. * sqrt(PI);
    msg("area of (weighted) region is %.15lg str\n", area);

    /* write polygons */
    ifile = argc - 1;
    nws = wrspher(argv[ifile], lmax, w);
    if (nws == -1) exit(1);

    return(0);
}

/*------------------------------------------------------------------------------
*/
void usage(void)
{
    printf("usage:\n");
    printf("harmonize [-d] [-q] [-l<lmax>] [-m<a>[u]] [-s<n>] [-e<n>] [-i<f>[<n>][u]] polygon_infile1 [polygon_infile2 ...] Wlm_outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"
