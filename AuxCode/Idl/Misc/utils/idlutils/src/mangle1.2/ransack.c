#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <values.h>
#include "format.h"
#include "polygon.h"
#include "vertices.h"
#include "scale.h"
#include "defaults.h"

/* getopt options */
static char *optstr = "dqm:c:r:u:p:";

/* local functions */
int ransack(), search();
void cmminf();
void usage(), parse_args();

/* external functions */
extern int rdmask();
extern int partition_poly();
extern int room_poly();
extern int copy_poly();
extern int prune_poly();
extern int garea();
extern int gptin();
extern double drandom();
extern void advise_fmt();
extern void msg(char *fmt, ...);
extern void free_poly();
extern void gaxisi_();
extern void rp_to_vert();
extern void scale(), scale_azel();
extern void wrangle();

polygon *poly[NPOLYSMAX];

/*------------------------------------------------------------------------------
  Main program.
*/
int main(argc, argv)
		 int argc;
		 char *argv [];
{
	int ifile, nfiles, np, npoly, npolys;

	/* parse arguments */
	parse_args(argc, argv);

	/* at least one input and output filename required as arguments */
	if (argc - optind < 2) {
		if (optind > 1 || argc - optind >= 1) {
			fprintf(stderr, "%s requires at least 2 arguments: infile and outfile\n", argv[0]);
	    usage();
	    exit(1);
		} else {
	    usage();
	    exit(0);
		}
	}

	msg("---------------- ransack ----------------\n");

	/* advise data format */
	advise_fmt(&fmt);

	/* tolerance angle for multiple intersections */
	if (mtol != 0.) {
		scale (&mtol, munit, 's');
		munit = 's';
		msg("multiple intersections closer than %g%c will be treated as coincident\n", mtol, munit);
		scale (&mtol, munit, 'r');
		munit = 'r';
	}

	/* warn about seed not being set */
	if (seed_read == 0) {
		msg("warning: seed was not set on command line: using default seed %d\n", seed);
	}

	/* read polygons */
	npoly = 0;
	nfiles = argc - 1 - optind;
	for (ifile = optind; ifile < optind + nfiles; ifile++) {
		npolys = rdmask(argv[ifile], &fmt, &poly[npoly], NPOLYSMAX - npoly);
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

	/* random points in polygons */
	ifile = argc - 1;
	np = ransack(argv[ifile], &fmt, poly, npoly, NPOLYSMAX);
	if (np == -1) exit(1);

	exit(0);
}

/*------------------------------------------------------------------------------
 */
void usage()
{
	printf("usage:\n");
	printf("ransack [-d] [-q] [-c<n>] [-r<n>] [-s<n>] [-e<n>] [-u<inunit>[,<outunit>]] [-p<n>] [-i<f>[<n>][u]] infile1 infile2 ... outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
 */
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Generate random az, el positions within mask defined by poly.
  The results are written to out_filename.

	Input: out_filename = name of file to write to;
	"" or "-" means write to standard output.
	fmt = pointer to format structure.
	poly = array of pointers to polygons.
	npoly = number of polygons in poly array.
	npolysmax = maximum number of polygons in poly array.
  Return value: number of random points generated,
	or -1 if error occurred.
*/
int ransack(out_filename, fmt, poly, npoly, npolysmax)
		 char *out_filename;
		 format *fmt;
		 int npoly, npolysmax;
#ifdef GCC
		 polygon *poly[npolysmax];
#else
		 polygon *poly[];
#endif
{
	/* lassoo all one-border polygons */
#define ALL_ONEBORDER		1
	/* how partition_poly should tighten lassoo */
#define ADJUST_LASSOO		2
	/* number of extra caps to allocate to polygon, to allow for expansion */
#define DNP			4
	/* length of state vector for random number generator */
#define STATELEN		256
	static char state[STATELEN], stateo[STATELEN];
#define AZEL_STR_LEN		32
	static char output[] = "output";
	char az_str[AZEL_STR_LEN], el_str[AZEL_STR_LEN];
	int dn, i, idmin, idmax, ier, in, inull, ip, ipmin, ipoly, ipolyo, iprune, irandom, lassooed, len, np, npolys, npoly_try, nwpoly, overwrite_original, tries, verb, width;
	double area, cmmin, cmmino, cmmint, cmi, phi, rp[3], rpoly, si, tol, w, wcum, x, xi[3], y, yi[3], z;
	double *wpoly;
	azel v;
	char *out_fn;
	FILE *outfile;

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
	if (fmt->outunit != fmt->inunit) {
		msg("units of output az, el angles will be ");
		switch (fmt->outunit) {
#include "angunit.h"
		}
		msg("\n");
	}

	/* initialize random number generator used by ransack() */
	initstate(seed, state, STATELEN);
	/* initialize random number generator used by ikrand() */
	initstate(seed, stateo, STATELEN);

	/* discard polygons with zero weight * area */
	np = 0;
	ier = 0;
	inull = 0;
	for (ipoly = 0; ipoly < npoly; ipoly++) {
		if (poly[ipoly]->weight == 0.) {
			inull++;
	    free_poly(poly[ipoly]);
	    poly[ipoly] = 0x0;
		} else {
			/* prune polygon */
	    iprune = prune_poly(poly[ipoly]);
	    /* error */
	    if (iprune == -1) {
				ier++;
				free_poly(poly[ipoly]);
				poly[ipoly] = 0x0;
				fprintf(stderr, "failed to prune polygon %d; discard it\n", ipoly);
				/* return(-1); */
				/* zero area polygon */
	    } else if (iprune >= 2) {
				inull++;
				free_poly(poly[ipoly]);
				poly[ipoly] = 0x0;
	    } else {
				np++;
	    }
		}
	}
	if (ier > 0) {
		msg("discarding %d unprunable polygons\n", ier);
	}
	if (inull > 0) {
		msg("discarding %d polygons with zero weight * area\n", inull);
	}

	if (np == 0) {
		fprintf(stderr, "no polygons to generate random points inside!\n");
		return(-1);
	}

	/* pre-lassoo polygons if there are many random points */
	if (nrandom >= np) {
		lassooed = 1;
		msg("lassooing %d polygons ...\n", np);
		/* partition_poly should overwrite original polygon only if lassoo succeeds */
		overwrite_original = 1;
		npolys = npoly;
		np = 0;
		for (ipoly = 0; ipoly < npoly; ipoly++) {
	    /* skip null polygon */
	    if (!poly[ipoly]) continue;
	    /* area/(2 pi) of smallest cap of polygon */
	    cmminf(poly[ipoly], &ipmin, &cmmino);
	    /* lassoo polygon */
	    dn = partition_poly(&poly[ipoly], &poly[npolys], npolysmax - npolys, ALL_ONEBORDER, ADJUST_LASSOO, overwrite_original, &npoly_try);
	    /* increment polygon count */
	    np++;
	    /* polygon was split successfully into several */
	    if (dn == npoly_try && dn > 0) {
				/* area/(2 pi) of combined smallest caps of lassooed polygon */
				cmminf(poly[ipoly], &ipmin, &cmmint);
				for (ip = 0; ip < dn; ip++) {
					cmminf(poly[npolys + ip], &ipmin, &cmmin);
					cmmint += cmmin;
				}
				/* lassooed polygons are genuine improvement */
				if (cmmint <= .9 * cmmino) {
					/* increment number of polygons */
					np += dn;
					npolys += dn;
					/* check whether exceeded maximum number of polygons */
					if (npolys > npolysmax) {
						fprintf(stderr, "ransack: total number of polygons exceeded maximum %d\n", npolysmax);
						fprintf(stderr, "if you need more space, enlarge NPOLYSMAX in defaults.h, and recompile\n");
						return(-1);
					}
				}
	    }
		}

    /* two few random points to make it worth pre-lassooing */
	} else {
		lassooed = 0;
		npolys = npoly;
	}

	/* largest width of polygon id number */
	idmin = 0;
	idmax = 0;
	for (ipoly = 0; ipoly < npoly; ipoly++) {
		if (!poly[ipoly]) continue;
		if (poly[ipoly]->id < idmin) idmin = poly[ipoly]->id;
		if (poly[ipoly]->id > idmax) idmax = poly[ipoly]->id;
	}
	sprintf(az_str, "%d", idmin);
	sprintf(el_str, "%d", idmax);
	idmin = strlen(az_str);
	idmax = strlen(el_str);
	width = (idmin > idmax)? idmin : idmax;

	/* write header */
	wrangle(0., fmt->outunit, fmt->outprecision, az_str, AZEL_STR_LEN);
	len = strlen(az_str);
	if (fmt->outunit == 'h') {
		sprintf(az_str, "az(hms)");
		sprintf(el_str, "el(dms)");
	} else {
		sprintf(az_str, "az(%c)", fmt->outunit);
		sprintf(el_str, "el(%c)", fmt->outunit);
	}
	fprintf(outfile, "%*s %*s %*s\n", len, az_str, len, el_str, width, "id");

	/* allocate memory for wpoly */
	nwpoly = npolys + DNP;
	wpoly = (double *) malloc(sizeof(double) * nwpoly);
	if (!wpoly) {
		fprintf(stderr, "ransack: failed to allocate memory for %d doubles\n", nwpoly);
		return(-1);
	}

	/* accept error messages from garea */
	/* unprunable polygons were already discarded, so garea should give no errors */
	verb = 1;

	/* cumulative area times weight of polygons */
	w = 0.;
	for (ipoly = 0; ipoly < npolys; ipoly++) {
		/* skip null polygons */
		if (poly[ipoly]) {
	    /* area of polygon */
	    tol = mtol;
	    ier = garea(poly[ipoly], &tol, &verb, &area);
	    if (ier) goto error;
	    /* accumulate weight times area */
	    w += poly[ipoly]->weight * area;
		}
		wpoly[ipoly] = w;
	}
	wcum = w;

	/* random points */
	if (strcmp(out_fn, output) != 0) {
		msg("generating %d random points in %d polygons ...\n", nrandom, np);
	}
	/* partition_poly should never overwrite original polygon */
	overwrite_original = 0;
	ipolyo = -1;
	dn = -1;
	for (irandom = 0; irandom < nrandom; irandom++) {

		/* which polygon to put random point in */
		setstate(state);
		rpoly = drandom() * wcum;
		setstate(stateo);
		ipoly = search(npolys, wpoly, rpoly);

		/* lassoo polygon, if not already lassooed */
		if (!lassooed && ipoly != ipolyo) {
	    /* area/(2 pi) of smallest cap of polygon */
	    cmminf(poly[ipoly], &ipmin, &cmmino);
	    /* lassoo polygon */
	    dn = partition_poly(&poly[ipoly], &poly[npolys], npolysmax - npolys, ALL_ONEBORDER, ADJUST_LASSOO, overwrite_original, &npoly_try);
	    /* polygon was lassooed successfully */
	    if (dn == npoly_try && dn > 0) {
				/* area/(2 pi) of combined smallest caps of lassooed polygon */
				cmmint = 0.;
				for (ip = 0; ip < dn; ip++) {
					cmminf(poly[npolys + ip], &ipmin, &cmmin);
					cmmint += cmmin;
				}
				/* lassooed polygons are genuine improvement */
				if ((dn == 1 && cmmint < cmmino) || cmmint <= .9 * cmmino) {
					/* enlarge memory for wpoly */
					if (npolys + dn > nwpoly) {
						nwpoly = npolys + dn + DNP;
						wpoly = (double *) realloc(wpoly, sizeof(double) * nwpoly);
						if (!wpoly) {
							fprintf(stderr, "ransack: failed to reallocate memory for %d doubles\n", nwpoly);
							return(-1);
						}
					}
					/* cumulative weight times area of lassooed polygons */
					w = (ipoly == 0)? 0. : wpoly[ipoly-1];
					for (ip = npolys; ip < npolys + dn; ip++) {
						/* area of polygon */
						tol = mtol;
						ier = garea(poly[ip], &tol, &verb, &area);
						if (ier) goto error;
						/* accumulate area times weight */
						w += poly[ip]->weight * area;
						wpoly[ip] = w;
					}
					/* lassooed polygons are too large to bother with */
				} else {
					dn = 0;
				}
				/* polygon was not lassooed */
	    } else {
				dn = 0;
	    }
		}
		ipolyo = ipoly;

		/* lassooed polygons */
		if (dn > 0) {
	    ipoly = npolys + search(dn, &wpoly[npolys], rpoly);
	    /* guard against roundoff */
	    if (ipoly >= npolys + dn) {
				fprintf(stderr, "%d should be < %d (i.e. %.15lg < %.15lg)\n", ipoly, npolys + dn, rpoly, wpoly[npolys + dn - 1]);
				ipoly = npolys + dn - 1;
	    }
		}

		/* smallest cap of polygon */
		cmminf(poly[ipoly], &ipmin, &cmmin);

		/* random point within polygon */
		tries = 0;
		do {
	    tries++;
	    /* random point within smallest cap */
	    setstate(state);
	    phi = TWOPI * drandom();
	    cmi = cmmin * drandom();
	    setstate(stateo);
	    /* coordinates of random point in cap frame */
	    si=sqrt(cmi * (2. - cmi));
	    x = si * cos(phi);
	    y = si * sin(phi);
	    z = 1. - cmi;
	    if (poly[ipoly]->cm[ipmin] < 0.) z = -z;
	    /* Cartesian axes with z-axis along cap axis */
	    gaxisi_(&poly[ipoly]->rp_(0, ipmin), &xi, &yi);
	    /* coordinates of random point */
	    for (i = 0; i < 3; i++) rp[i] = x * xi[i] + y * yi[i] + z * poly[ipoly]->rp_(i, ipmin);
	    /* whether random point is inside polygon */
	    in = gptin(poly[ipoly], rp);
		} while (!in);

		/* convert unit vector to az, el */
		rp_to_vert(rp, &v);
		v.az -= floor(v.az / TWOPI) * TWOPI;

		/* convert az and el from radians to output units */
		scale_azel(&v, 'r', fmt->outunit);

		/* write result */
		wrangle(v.az, fmt->outunit, fmt->outprecision, az_str, AZEL_STR_LEN);
		wrangle(v.el, fmt->outunit, fmt->outprecision, el_str, AZEL_STR_LEN);
		fprintf(outfile, "%s %s %*d\n", az_str, el_str, width, poly[ipolyo]->id);
		/* fprintf(outfile, "%s %s %d %d %d %d %lg %lg %lg %lg %d %d\n", az_str, el_str, irandom, ipolyo, ipoly, tries, wcum, rpoly / wcum, area, TWOPI * cmmin / area, ipmin, poly[ipoly]->np); */

	}

	/* advise */
	if (outfile != stdout) {
		msg("ransack: %d random positions written to %s\n", nrandom, out_fn);
	}

	return(nrandom);

	/* error returns */
 error:
	return(-1);
}

/*------------------------------------------------------------------------------
  Find position within ordered array by binary chop.
  Assumes points are uncorrelated between calls
  (otherwise it would be better to start search from previous point).

  Return value: i such that array[i-1] <= point < array[i]
	0 if point < array[0]
	n if point >= array[n-1]
*/
int search(n, array, point)
		 int n;
#ifdef GCC
		 double array[n];
#else
		 double array[];
#endif
		 double point;
{
	int i, im, ip;

	/* point below minimum */
	if (point < array[0]) {
		return(0);

    /* point above maximum */
	} else if (point >= array[n-1]) {
		return(n);

    /* point between limits */
	} else {
		im = 0;
		ip = n-1;
		/* binary chop */
		while (im + 1 < ip) {
	    i = (im + ip) / 2;
	    if (point < array[i]) {
				ip = i;
	    } else if (point >= array[i]) {
				im = i;
	    }
		};
		return(ip);
	} 
}

/*------------------------------------------------------------------------------
  Find smallest cap of polygon.
*/
void cmminf(poly, ipmin, cmmin)
		 polygon *poly;
		 int *ipmin;
		 double *cmmin;
{
	int ip;
	double cmi;

	*cmmin = 3.;
	for (ip = 0; ip < poly->np; ip++) {
		cmi = (poly->cm[ip] > 0.)? poly->cm[ip] : 2. + poly->cm[ip];
		if (cmi < *cmmin) {
			*ipmin = ip;
			*cmmin = cmi;
		}
	}
}
