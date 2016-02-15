/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "manglefn.h"
#include "defaults.h"

polygon *polys[NPOLYSMAX];

/* getopt options */
const char *optstr = "dqm:s:e:v:p:i:o:";

/* local functions */
void	usage(void);
#ifdef	GCC
int	balkanize(int npoly, polygon *[npoly], int npolys, polygon *[npolys]);
#else
int	balkanize(int npoly, polygon *[/*npoly*/], int npolys, polygon *[/*npolys*/]);
#endif

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
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
	    fprintf(stderr, "%s requires at least 2 arguments: polygon_infile and polygon_outfile\n", argv[0]);
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

    /* balkanize polygons */
    npolys = balkanize(npoly, polys, NPOLYSMAX - npoly, &polys[npoly]);
    if (npolys == -1) exit(1);

    /* write polygons */
    ifile = argc - 1;
    npolys = wrmask(argv[ifile], &fmt, npolys, &polys[npoly],0);
    if (npolys == -1) exit(1);
    /* memmsg(); */

    return(0);
}

/*------------------------------------------------------------------------------
*/
void usage(void)
{
    printf("usage:\n");
    printf("balkanize [-d] [-q] [-m<a>[u]] [-s<n>] [-e<n>] [-vo|-vn] [-p[+|-][<n>]] [-i<f>[<n>][u]] [-o<f>[u]] polygon_infile1 [polygon_infile2 ...] polygon_outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Balkanize overlapping polygons into many disjoint connected polygons.

   Input: npoly = number of polygons.
	  poly = array of pointers to polygons.
	  npolys = maximum number of output polygons.
  Output: polys = array of pointers to polygons.
  Return value: number of disjoint connected polygons,
		or -1 if error occurred.
*/
int balkanize(int npoly, polygon *poly[/*npoly*/], int npolys, polygon *polys[/*npolys*/])
{
/* part_poly should lasso one-boundary polygons only if they have too many caps */
#define ALL_ONEBOUNDARY		1
/* how part_poly should tighten lasso */
#define ADJUST_LASSO		1
/* part_poly should force polygon to be split even if no part can be lassoed */
#define	FORCE_SPLIT		1
/* partition_poly should overwrite all original polygons */
#define OVERWRITE_ORIGINAL	2
#define WARNMAX			8
    int discard, dm, dn, dnp, failed, i, ier, inull, ip, iprune, j, k, m, n, np;

    /* start by pruning all input polygons */
    np = 0;
    inull = 0;
    for (i = 0; i < npoly; i++) {
	iprune = prune_poly(poly[i], mtol);
	/* error */
	if (iprune == -1) {
	    fprintf(stderr, "balkanize: initial prune failed at polygon %d\n", poly[i]->id);
	    return(-1);
	}
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
		dn = fragment_poly(&polys[k], poly[j], discard, npolys - n, &polys[n], mtol);
		/* error */
		if (dn == -1) {
		    fprintf(stderr, "balkanize: UHOH at polygon %d; continuing ...\n", (fmt.newid == 'o')? polys[i]->id : ip);
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
		    fprintf(stderr, "if you need more space, enlarge NPOLYSMAX in defines.h, and recompile\n");
		    return(-1);
		}
	    }
	    /* copy down non-null polygons */
	    dm = 0;
	    for (k = m; k < n; k++) {
		if (polys[k]) {
		    polys[m + dm] = polys[k];
		    dm++;
		}
	    }
	    /* nullify but don't free, because freeing polys[k] will free polys[m + dm] */
	    for (k = m + dm; k < n; k++) {
		polys[k] = 0x0;
	    }
	    n = m + dm;
	    if (dm == 0) break;
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
	ier = partition_poly(&polys[i], npolys - n, &polys[n], mtol, ALL_ONEBOUNDARY, ADJUST_LASSO, FORCE_SPLIT, OVERWRITE_ORIGINAL, &dn);
	/* error */
	if (ier == -1) {
	    fprintf(stderr, "balkanize: UHOH at polygon %d; continuing ...\n", (fmt.newid == 'o')? polys[i]->id : ip);
	    continue;
	    /* return(-1); */
	/* failed to partition polygon into desired number of parts */
	} else if (ier == 1) {
	    fprintf(stderr, "balkanize: failed to partition polygon %d fully; partitioned it into %d parts\n", (fmt.newid == 'o')? polys[i]->id : ip, dn + 1);
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
	    fprintf(stderr, "if you need more space, enlarge NPOLYSMAX in defines.h, and recompile\n");
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
	msg("Whatever the case, the output file of balkanized polygons constitutes\n");
	msg("a valid mask of non-overlapping polygons, which is safe to use.\n");
	msg(".............................................................................\n");
    }

    /* final prune */
    j = 0;
    inull = 0;
    for (i = 0; i < n; i++) {
	iprune = prune_poly(polys[i], mtol);
	if (iprune == -1) {
	    fprintf(stderr, "balkanize: failed to prune polygon %d; continuing ...\n", (fmt.newid == 'o')? polys[i]->id : j);
	    /* return(-1); */
	}
	if (iprune >= 2) {
	    free_poly(polys[i]);
	    polys[i] = 0x0;
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
