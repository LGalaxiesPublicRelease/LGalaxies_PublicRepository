#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <values.h>
#include "format.h"
#include "harmonics.h"
#include "inputfile.h"
#include "scale.h"
#include "defaults.h"

/* redefine default maximum harmonic */
#undef LMAX
#define LMAX		MAXINT

inputfile file = {
    '\0',	/* input filename */
    0x0,	/* input file stream */
    '\0',	/* line buffer */
    64,		/* size of line buffer (will expand as necessary) */
    0,		/* line number */
    0		/* maximum number of characters to read (0 = no limit) */
};

/* getopt options */
static char *optstr = "dqw:l:g:x:u:p:";

/* local functions */
int map();
void usage(), parse_args();

/* external functions */
extern int rdangle(), rdline();
extern int rdspher();
extern double wrho();
extern void msg(char *fmt, ...);
extern void scale();
extern void wrangle();

/*------------------------------------------------------------------------------
  Main program.
*/
int main(argc, argv)
int argc;
char *argv [];
{
    int nmap, nws;
    double **w_p;

    /* lmax will be read from Wlm_filename, or set by command line */
    lmax = LMAX;

    /* parse arguments */
    parse_args(argc, argv);

    /* one input and output filename required as arguments */
    if (argc - optind != 2) {
	if (optind > 1 || argc - optind == 1) {
	    fprintf(stderr, "%s requires 2 arguments: infile and outfile\n", argv[0]);
	    usage();
	    exit(1);
	} else {
	    usage();
	    exit(0);
	}
    }

    /* file containing harmonics must have been specified */
    if (!Wlm_filename) {
	fprintf(stderr, "%s requires -w<Wlm_filename> option to specify file containing harmonics\n", argv[0]);
	exit(1);
    }

    msg("---------------- map ----------------\n");

    /* advise */
    if (lsmooth == 0.) {
	msg("no smoothing\n");
    } else {
	msg("smoothing harmonic number lsmooth = %lg\n", lsmooth);
	if (esmooth != 2.) {
	    msg("smoothing exponent = %lg\n", esmooth);
	}
    }

    /* read harmonics */
    nws = rdspher(Wlm_filename, &lmax, &w_p);
    if (nws == -1) exit(1);

    /* map */
    nmap = map(argv[argc - 2], argv[argc - 1], &fmt, w_p, lmax, lsmooth, esmooth);
    if (nmap == -1) exit(1);

    exit(0);
}

/*------------------------------------------------------------------------------
*/
void usage()
{
    printf("usage:\n");
    printf("map [-d] [-q] -w<Wlmfile> [-l<lmax>] [-g<lsmooth>] [-u<inunit>[,<outunit>]] [-p<n>] infile outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Map.  Implemented as interpretive read/write, to permit interactive behaviour.

   Input: in_filename = name of file to read from;
			"" or "-" means read from standard input.
	  out_filename = name of file to write to;
			"" or "-" means write to standard output.
	  fmt = pointer to format structure.
	  w = array containing spherical harmonics.
	  lmax = maximum harmonic number.
	  lsmooth = smoothing harmonic number (0. = no smooth).
	  esmooth = smoothing exponent (2. = gaussian).
  Return value: number of items written,
		or -1 if error occurred.
*/
int map(in_filename, out_filename, fmt, w, lmax, lsmooth, esmooth)
char *in_filename, *out_filename;
format *fmt;
int lmax;
double lsmooth, esmooth;
#ifdef GCC
double w[IM * NW];
#else
double w[];
#endif
{
/* precision of map values written to file */
#define PRECISION	8
#define AZEL_STR_LEN	32
    static char input[] = "input", output[] = "output";
    char *word, *next;
    char az_str[AZEL_STR_LEN], el_str[AZEL_STR_LEN];
    int ird, len, mmax, nmap, width;
    double az, el, rho;
    char *out_fn;
    FILE *outfile;

    /* open in_filename for reading */
    if (!in_filename || strcmp(in_filename, "-") == 0) {
	file.file = stdin;
	file.name = input;
    } else {
	file.file = fopen(in_filename, "r");
	if (!file.file) {
	    fprintf(stderr, "cannot open %s for reading\n", in_filename);
	    return(-1);
	}
	file.name = in_filename;
    }
    file.line_number = 0;

    /* open out_filename for writing */
    if (!out_filename || strcmp(out_filename, "-") == 0) {
	outfile = stdout;
	out_fn = output;
    } else {
	outfile = fopen(out_filename, "w");
	if (!outfile) {
	    fprintf(stderr, "cannot open %s for writing\n", out_filename);
	    return(-1);
	}
	out_fn = out_filename;
    }

    /* advise angular units */
    msg("will take units of input az, el angles in %s to be ", file.name);
    switch (fmt->inunit) {
#include "angunit.h"
    }
    msg("\n");
    if (fmt->outunit != fmt->inunit) {
	msg("units of output az, el angles will be ");
	switch (fmt->outunit) {
#include "angunit.h"
	}
	msg("\n");
    }

    /* setting mmax = lmax ensures complete set of azimuthal harmonics */
    mmax = lmax;

    /* width of map value */
    width = PRECISION + 6;

    /* write header */
    az = 0.;
    wrangle(az, fmt->outunit, fmt->outprecision, az_str, AZEL_STR_LEN);
    len = strlen(az_str);
    if (fmt->outunit == 'h') {
	sprintf(az_str, "az(hms)");
	sprintf(el_str, "el(dms)");
    } else {
	sprintf(az_str, "az(%c)", fmt->outunit);
	sprintf(el_str, "el(%c)", fmt->outunit);
    }
    fprintf(outfile, "%*s %*s %*s\n", len, az_str, len, el_str, width - 4, "wrho");

    /* interpretive read/write loop */
    nmap = 0;
    while (1) {
	/* read line */
	ird = rdline(&file);
	/* serious error */
	if (ird == -1) return(-1);
	/* EOF */
	if (ird == 0) break;

	/* read <az> */
	word = file.line;
	ird = rdangle(word, &next, fmt->inunit, &az);
	/* skip header */
	if (ird != 1 && nmap == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* read <el> */
	word = next;
	ird = rdangle(word, &next, fmt->inunit, &el);
	/* skip header */
	if (ird != 1 && nmap == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* convert az and el from input units to radians */
	scale(&az, fmt->inunit, 'r');
	scale(&el, (fmt->inunit == 'h')? 'd' : fmt->inunit, 'r');

	/*
	  The entire of map.c is an interface to the next line of code.
	  Bizarre, huh?
	*/
        /* compute the value of the window function at this point */
	rho = wrho(az, el, w, lmax, mmax, lsmooth, esmooth);

	/* convert az and el from radians to output units */
	scale(&az, 'r', fmt->outunit);
	scale(&el, 'r', (fmt->outunit == 'h')? 'd' : fmt->outunit);

	/* write result */
	wrangle(az, fmt->outunit, fmt->outprecision, az_str, AZEL_STR_LEN);
	wrangle(el, fmt->outunit, fmt->outprecision, el_str, AZEL_STR_LEN);
	fprintf(outfile, "%s %s %- #*.*lg\n", az_str, el_str, width, PRECISION, rho);
	fflush(outfile);

        /* increment counter of results */
	nmap++;
    }

    if (outfile != stdout) {
	msg("map: %d values written to %s\n", nmap, out_fn);
    }

    return(nmap);
}
