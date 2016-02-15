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
static char *optstr = "dqu:p:";

/* declared in rdmask */
extern inputfile file;

/* local functions */
int poly_ids(), poly_id();
void usage(), parse_args();

/* external functions */
extern int rdmask();
extern int rdangle(), rdline();
extern int gptin();
extern void advise_fmt();
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
    int np, npoly;
    polygon *poly[NPOLYSMAX];

    /* parse arguments */
    parse_args(argc, argv);

    /* two input and one output filename required as arguments */
    if (argc - optind != 3) {
	if (optind > 1 || argc - optind >= 1) {
	    fprintf(stderr, "%s requires 3 arguments: polygon_infile, infile and outfile\n", argv[0]);
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
    npoly = rdmask(argv[argc - 3], &fmt, poly, NPOLYSMAX);
    if (npoly == -1) exit(1);
    if (npoly == 0) {
	msg("stop\n");
	exit(0);
    }

    /* polygon id numbers */
    np = poly_ids(argv[argc - 2], argv[argc - 1], &fmt, poly, npoly);
    if (np == -1) exit(1);

    exit(0);
}

/*------------------------------------------------------------------------------
*/
void usage()
{
    printf("usage:\n");
    printf("polyid [-d] [-q] [-u<inunit>[,<outunit>]] [-p<n>] polygon_infile infile outfile\n");
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
int poly_ids(in_filename, out_filename, fmt, poly, npoly)
char *in_filename, *out_filename;
format *fmt;
int npoly;
#ifdef GCC
polygon *poly[npoly];
#else
polygon *poly[];
#endif
{
#define AZEL_STR_LEN	32
    static char input[] = "input", output[] = "output";
    char *word, *next;
    char az_str[AZEL_STR_LEN], el_str[AZEL_STR_LEN];
    int i, ird, len, nid, nids, nid0, nid2, np;
    int *id;
    double az, el;
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
    fprintf(outfile, "%*s %*s polygon_ids\n", len, az_str, len, el_str);

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
	ird = rdangle(word, &next, fmt->inunit, &az);
	/* skip header */
	if (ird != 1 && np == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* read <el> */
	word = next;
	ird = rdangle(word, &next, fmt->inunit, &el);
	/* skip header */
	if (ird != 1 && np == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* convert az and el from input units to radians */
	scale(&az, fmt->inunit, 'r');
	scale(&el, (fmt->inunit == 'h')? 'd' : fmt->inunit, 'r');

        /* id numbers of the polygons containing position az, el */
	nid = poly_id(poly, npoly, az, el, &id);

	/* convert az and el from radians to output units */
	scale(&az, 'r', fmt->outunit);
	scale(&el, 'r', (fmt->outunit == 'h')? 'd' : fmt->outunit);

	/* write result */
	wrangle(az, fmt->outunit, fmt->outprecision, az_str, AZEL_STR_LEN);
	wrangle(el, fmt->outunit, fmt->outprecision, el_str, AZEL_STR_LEN);
	fprintf(outfile, "%s %s", az_str, el_str);
	for (i = 0; i < nid; i++) {
	    fprintf(outfile, " %d", id[i]);
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

/*------------------------------------------------------------------------------
  Id numbers of polygons containing position az, el.

   Input: poly = array of pointers to npoly polygons.
	  npoly = number of polygons in poly array.
	  az, el = angular position in radians.
  Output: id_p = pointer to array containing id numbers of polygons.
  Return value: number of polygons that contain az, el position.
*/
int poly_id(poly, npoly, az, el, id_p)
int npoly;
#ifdef GCC
polygon *poly[npoly];
#else
polygon *poly[];
#endif
double az, el;
int **id_p;
{
/* number of extra polygon id numbers to allocate, to allow for expansion */ 
#define DNID		16
    static int nidmax = 0;
    static int *id = 0x0;

    int ipoly, nid;
    double rp[3];

    /* unit vector corresponding to angular position az, el */
    rp[0] = cos(el) * cos(az);
    rp[1] = cos(el) * sin(az);
    rp[2] = sin(el);

    nid = 0;

    /* keep trying till the id array is big enough */
    do {
        /* make sure that allocated id array contain enough space */
        if (!id || nid > nidmax) {
            if (id) free(id);
            id = (int *) malloc(sizeof(int) * (nid + DNID));
            if (!id) {
                fprintf(stderr, "poly_id: failed to allocate memory for %d integers\n", nid + DNID);
                return(-1);
            }
	    nidmax = nid + DNID;
	}

	nid = 0;
	/* do each polygon in turn */
	for (ipoly = 0; ipoly < npoly; ipoly++) {
	    /* id number of each polygon that contains az, el position */
	    if (gptin(poly[ipoly], rp)) {
		if (nid < nidmax) id[nid] = poly[ipoly]->id;
		nid++;
	    }
	}

    } while (nid > nidmax);

    /* point id_p at id array */
    *id_p = id;

    /* number of polygons containing az, el position */
    return(nid);
}
