/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "inputfile.h"
#include "manglefn.h"
#include "defaults.h"

polygon *poly[NPOLYSMAX];

/* getopt options */
const char *optstr = "dqu:p:";

/* declared in rdmask */
extern inputfile file;

/* local functions */
void	usage(void);
#ifdef	GCC
int	poly_ids(char *, char *, format *, int npoly, polygon *[npoly]);
#else
int	poly_ids(char *, char *, format *, int npoly, polygon *[/*npoly*/]);
#endif

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
    int ifile, nfiles, npoly, npolys;

    /* parse arguments */
    parse_args(argc, argv);

    /* at least two input and one output filename required as arguments */
    if (argc - optind < 3) {
	if (optind > 1 || argc - optind >= 1) {
	    fprintf(stderr, "%s requires at least 3 arguments: polygon_infile, azel_infile, and outfile\n", argv[0]);
	    usage();
	    exit(1);
	} else {
	    usage();
	    exit(0);
	}
    }

    msg("---------------- polyid ----------------\n");

    /* advise data format */
    advise_fmt(&fmt);

    /* read polygons */
    npoly = 0;
    nfiles = argc - 2 - optind;
    for (ifile = optind; ifile < optind + nfiles; ifile++) {
	npolys = rdmask(argv[ifile], &fmt, NPOLYSMAX - npoly, &poly[npoly]);
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

    /* polygon id numbers */
    npolys = poly_ids(argv[argc - 2], argv[argc - 1], &fmt, npoly, poly);
    if (npolys == -1) exit(1);

    return(0);
}

/*------------------------------------------------------------------------------
*/
void usage(void)
{
    printf("usage:\n");
    printf("polyid [-d] [-q] [-u<inunit>[,<outunit>]] [-p[+|-][<n>]] polygon_infile1 [polygon_infile2 ...] azel_infile outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Id numbers of polygons containing az, el positions.
  The az, el positions are read from in_filename,
  and the results are written to out_filename.
  Implemented as interpretive read/write, to permit interactive behaviour.

   Input: in_filename = name of file to read from;
			"" or "-" means read from standard input.
	  out_filename = name of file to write to;
			"" or "-" means write to standard output.
	  fmt = pointer to format structure.
	  poly = array of pointers to polygons.
	  npoly = number of polygons in poly array.
  Return value: number of lines written,
		or -1 if error occurred.
*/
int poly_ids(char *in_filename, char *out_filename, format *fmt, int npoly, polygon *poly[/*npoly*/])
{
#define AZEL_STR_LEN	32
    char input[] = "input", output[] = "output";
    char *word, *next;
    char az_str[AZEL_STR_LEN], el_str[AZEL_STR_LEN];
    int i, idmax, idmin, idwidth, ird, len, nid, nids, nid0, nid2, np;
    int *id;
    azel v;
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

    /* largest width of polygon id number */
    idmin = 0;
    idmax = 0;
    for (i = 0; i < npoly; i++) {
	if (!poly[i]) continue;
	if (poly[i]->id < idmin) idmin = poly[i]->id;
	if (poly[i]->id > idmax) idmax = poly[i]->id;
    }
    idmin = ((idmin < 0)? floor(log10((double)-idmin)) + 2 : 1);
    idmax = ((idmax > 0)? floor(log10((double)idmax)) + 1 : 1);
    idwidth = ((idmin > idmax)? idmin : idmax);

    /* write header */
    v.az = 0.;
    wrangle(v.az, fmt->outunit, fmt->outprecision, AZEL_STR_LEN, az_str);
    len = strlen(az_str);
    if (fmt->outunit == 'h') {
	sprintf(az_str, "az(hms)");
	sprintf(el_str, "el(dms)");
    } else {
	sprintf(az_str, "az(%c)", fmt->outunit);
	sprintf(el_str, "el(%c)", fmt->outunit);
    }
    fprintf(outfile, "%*s %*s", len, az_str, len, el_str);
    if (npoly > 0) fprintf(outfile, " polygon_ids");
    fprintf(outfile, "\n");

    /* interpretive read/write loop */
    np = 0;
    nids = 0;
    nid0 = 0;
    nid2 = 0;
    while (1) {
	/* read line */
	ird = rdline(&file);
	/* serious error */
	if (ird == -1) return(-1);
	/* EOF */
	if (ird == 0) break;

	/* read <az> */
	word = file.line;
	ird = rdangle(word, &next, fmt->inunit, &v.az);
	/* skip header */
	if (ird != 1 && np == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* read <el> */
	word = next;
	ird = rdangle(word, &next, fmt->inunit, &v.el);
	/* skip header */
	if (ird != 1 && np == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* convert az and el from input units to radians */
	scale_azel(&v, fmt->inunit, 'r');

        /* id numbers of the polygons containing position az, el */
	nid = poly_id(npoly, poly, v.az, v.el, &id);

	/* convert az and el from radians to output units */
	scale_azel(&v, 'r', fmt->outunit);

	/* write result */
	wrangle(v.az, fmt->outunit, fmt->outprecision, AZEL_STR_LEN, az_str);
	wrangle(v.el, fmt->outunit, fmt->outprecision, AZEL_STR_LEN, el_str);
	fprintf(outfile, "%s %s", az_str, el_str);
	for (i = 0; i < nid; i++) {
	    fprintf(outfile, " %*d", idwidth, id[i]);
	}
	fprintf(outfile, "\n");
	fflush(outfile);

        /* increment counters of results */
	np++;
	nids += nid;
	if (nid == 0) {
	    nid0++;
	} else if (nid >= 2) {
	    nid2++;
	}
    }

    /* advise */
    if (nid0 > 0) msg("%d points were not inside any polygon\n", nid0);
    if (nid2 > 0) msg("%d points were inside >= 2 polygons\n", nid2);

    if (outfile != stdout) {
	msg("polyid: %d id numbers at %d positions written to %s\n", nids, np, out_fn);
    }

    return(np);
}
