#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <values.h>
#include "format.h"
#include "inputfile.h"
#include "polygon.h"
#include "scale.h"
#include "defaults.h"

/* getopt options */
static char *optstr = "dqf:u:p:";

/* declared in rdmask() */
extern inputfile file;

/* local functions */
int rotate();
int parse_fopt();
void usage(), parse_args();
void rotate_azel();

/* external functions */
extern int rdangle(), rdline();
extern double places();
extern void msg(char *fmt, ...);
extern void scale();
extern void wrangle();
extern void azel_();
extern void fframe_();

/*------------------------------------------------------------------------------
  Main program.
*/
int main(argc, argv)
int argc;
char *argv [];
{
    int itr, np;

    /* parse arguments */
    parse_args(argc, argv);
    /* parse option <fopt> to -f<fopt> */
    if (fopt) itr = parse_fopt();

    /* one input and one output filename required as arguments */
    if (argc - optind != 2) {
	if (optind > 1 || argc - optind >= 1) {
	    fprintf(stderr, "%s requires 2 arguments: infile and outfile\n", argv[0]);
	    usage();
	    exit(1);
	} else {
	    usage();
	    exit(0);
	}
    }

    msg("---------------- rotate ----------------\n");

    /* rotate */
    np = rotate(argv[argc - 2], argv[argc - 1], &fmt);
    if (np == -1) exit(1);

    exit(0);
}

/*------------------------------------------------------------------------------
*/
void usage()
{
    printf("usage:\n");
    printf("rotate [-d] [-q] [-f[<inframe>[,<outframe>]|<azn>,<eln>,<azp>[u]]] [-u<inunit>[,<outunit>]] [-p<n>] infile outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
*/
#include "parse_fopt.c"

/*------------------------------------------------------------------------------
  Rotate az, el positions from one frame to another.
  The az, el positions are read from in_filename,
  and rotated az, el positions are written to out_filename.
  Implemented as interpretive read/write, to permit interactive behaviour.

   Input: in_filename = name of file to read from;
			"" or "-" means read from standard input.
	  out_filename = name of file to write to;
			"" or "-" means write to standard output.
	  fmt = pointer to format structure.
  Return value: number of lines written,
		or -1 if error occurred.
*/
int rotate(in_filename, out_filename, fmt)
char *in_filename, *out_filename;
format *fmt;
{
#define AZEL_STR_LEN	32
    static char input[] = "input", output[] = "output";
    char *word, *next;
    char az_str[AZEL_STR_LEN], el_str[AZEL_STR_LEN];
    int ird, len, np;
    double azi, azf, eli, elf;
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

    /* advise custom transformation */
    if (fmt->outframe == -1) {
	/* multiple transformation */
	if (strchr(fopt, ':')) {
	    msg(" -f%s is equivalent to\n", fopt);
	    msg(" -f%.16lg,%.16lg,%.16lg%c\n",
		places(fmt->azn, 14), places(fmt->eln, 14), places(fmt->azp, 14), fmt->trunit);
	/* single transformation */
	} else {
	    msg("rotate -f%.16lg,%.16lg,%.16lg%c\n",
		fmt->azn, fmt->eln, fmt->azp, fmt->trunit);
	}
    /* advise standard transformation */
    } else if (fmt->inframe != fmt->outframe) {
	msg("rotate from %s to %s\n",
	    frames[fmt->inframe], frames[fmt->outframe]);
	msg(" -f%s,%s is equivalent to\n",
	    frames[fmt->inframe], frames[fmt->outframe]);
        msg(" -f%.16lg,%.16lg,%.16lg%c\n",
	    places(fmt->azn, 14), places(fmt->eln, 14), places(fmt->azp, 14), fmt->trunit);
    }

    /* angular units */
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

    /* write header */
    azf = 0.;
    wrangle(azf, fmt->outunit, fmt->outprecision, az_str, AZEL_STR_LEN);
    len = strlen(az_str);
    if (fmt->outunit == 'h') {
	sprintf(az_str, "az(hms)");
	sprintf(el_str, "el(dms)");
    } else {
	sprintf(az_str, "az(%c)", fmt->outunit);
	sprintf(el_str, "el(%c)", fmt->outunit);
    } 
    fprintf(outfile, "%*s %*s\n", len, az_str, len, el_str);

    /* interpretive read/write loop */
    np = 0;
    while (1) {
	/* read line */
	ird = rdline(&file);
	/* serious error */
	if (ird == -1) return(-1);
	/* EOF */
	if (ird == 0) break;

	/* read <az> */
	word = file.line;
	ird = rdangle(word, &next, fmt->inunit, &azi);
	/* skip header */
	if (ird != 1 && np == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* read <el> */
	word = next;
	ird = rdangle(word, &next, fmt->inunit, &eli);
	/* skip header */
	if (ird != 1 && np == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* convert az and el from input units to degrees */
	scale(&azi, fmt->inunit, 'd');
	scale(&eli, (fmt->inunit == 'h')? 'd' : fmt->inunit, 'd');

	/* rotate az and el */
	rotate_azel(fmt, &azi, &eli, &azf, &elf);

	/* convert az and el from degrees to output units */
	scale(&azf, 'd', fmt->outunit);
	scale(&elf, 'd', (fmt->outunit == 'h')? 'd' : fmt->outunit);

	/* write result */
	wrangle(azf, fmt->outunit, fmt->outprecision, az_str, AZEL_STR_LEN);
	wrangle(elf, fmt->outunit, fmt->outprecision, el_str, AZEL_STR_LEN);
	fprintf(outfile, "%s %s\n", az_str, el_str);
	fflush(outfile);

        /* increment counters of results */
	np++;
    }

    if (outfile != stdout) {
	msg("rotate: %d positions written to %s\n", np, out_fn);
    }

    return(np);
}

/*------------------------------------------------------------------------------
   Interface to fortran routines that actually do the rotation.
*/
void rotate_azel(fmt, azi, eli, azf, elf)
format *fmt;
double *azi, *eli, *azf, *elf;
{
    /* custom */
    if (fmt->outframe == -1) {
	azel_(azi, eli, &fmt->azn, &fmt->eln, &fmt->azp, azf, elf);

    /* built-ins */
    } else {
	fframe_(&fmt->inframe, azi, eli, &fmt->outframe, azf, elf);

    }
}
