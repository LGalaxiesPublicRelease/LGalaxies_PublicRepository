#include <math.h>
#include <stdio.h>
#include <string.h>
#include "format.h"
#include "pi.h"
#include "polygon.h"
#include "vertices.h"

#define TWOPI		(2. * PI)

#define WARNMAX		8
#define	AZEL_STR_LEN	32

/* suppress error messages from garea */
static int verb = 0;

/* initial angular tolerance within which to merge multiple intersections */
extern double mtol;

/* min, max weights to keep */
extern int is_weight_min, is_weight_max;
extern double weight_min, weight_max;

/* min, max areas to keep */
extern int is_area_min, is_area_max;
extern double area_min, area_max;

/* local functions */
int wrmask(), wr_circ(), wr_edge(), wr_rect(), wr_poly(), wr_Reg();
int discard_poly();

/* external functions */ 
extern int poly_to_rect();
extern int vmid();
extern int garea(), gverts();
extern void msg(char *, ...);
extern void free_poly();
extern void rpcm_to_circ(), rp_to_vert();
extern void scale(), scale_azel();
extern void wrangle();

/*-----------------------------------------------------------------------------
  Write mask data.

	Input: filename = name of file to write to;
	"" or "-" means write to standard output.
	fmt = pointer to format structure.
	polys = polygons to write.
	npolys = number of polygons to write.
  Return value: number of polygons written,
	or -1 if error occurred.
*/ 
int wrmask(filename, fmt, polys, npolys, selfdestruct)
		 char *filename;
		 format *fmt;
		 int npolys;
		 int selfdestruct;
#ifdef GCC
		 polygon *polys[npolys];
#else
		 polygon *polys[];
#endif
{
	int npoly;

	/* discard polygons with weight or area outside specified limits */
	npoly = discard_poly(polys, npolys);

	if (strcmp(fmt->out, "circle") == 0) {
		npoly = wr_circ(filename, fmt, polys, npolys, npoly);

	} else if (strcmp(fmt->out, "edges") == 0
						 || strcmp(fmt->out, "vertices") == 0) {
		npoly = wr_edge(filename, fmt, polys, npolys, npoly);

	} else if (strcmp(fmt->out, "rectangle") == 0) {
		npoly = wr_rect(filename, fmt, polys, npolys, npoly);

	} else if (strcmp(fmt->out, "polygon") == 0) {
		npoly = wr_poly(filename, polys, npolys, npoly, selfdestruct);

	} else if (strcmp(fmt->out, "Region") == 0) {
		npoly = wr_Reg(filename, polys, npolys, npoly);

	} else {
		npoly = 0;

	}

	return(npoly);
}

/*------------------------------------------------------------------------------
  Write mask data in circle format.

	Input: filename = name of file to write to;
	"" or "-" means write to standard output.
	fmt = pointer to format structure.
	polys = polygons to write.
	npolys = number of polygons.
	npolyw = number of polygons to write.
  Return value: number of polygons written,
	or -1 if error occurred.
*/ 
int wr_circ(filename, fmt, polys, npolys, npolyw)
		 char *filename;
		 format *fmt;
		 int npolys, npolyw;
#ifdef GCC
		 polygon *polys[npolys];
#else
		 polygon *polys[];
#endif
{
	char unit;
	char az_str[AZEL_STR_LEN], el_str[AZEL_STR_LEN], th_str[AZEL_STR_LEN];
	int i, ier, ip, ipoly, nbadarea, npoly;
	double area, angle[3], tol;
	FILE *file;
	char *circle_fmt = "circle %d ( %d caps, %.15lg weight, %.15lf str):\n";

	/* open filename for writing */
	if (!filename || strcmp(filename, "-") == 0) {
		file = stdout;
	} else {
		file = fopen(filename, "w");
		if (!file) {
	    fprintf(stderr, "wr_circ: cannot open %s for writing\n", filename);
	    return(-1);
		}
	}

	/* write number of polygons */
	fprintf(file, "%d polygons\n", npolyw);

	/* write angular unit */
	fprintf(file, "unit %c\n", fmt->outunitp);

	npoly = 0;
	nbadarea = 0;
	for (ipoly = 0; ipoly < npolys; ipoly++) {
		/* discard null polygons */
		if (!polys[ipoly]) continue;
		/* area of polygon */
		tol = mtol;
		ier = garea(polys[ipoly], &tol, &verb, &area);
		if (ier == -1) return(-1);
		if (ier) {
	    fprintf(stderr, "wr_circ: area of polygon %d is incorrect\n", 
							polys[ipoly]->id);
	    nbadarea++;
		}
		/* number of caps, weight, and area of polygon */
		fprintf(file, circle_fmt,
						polys[ipoly]->id, polys[ipoly]->np, polys[ipoly]->weight, area);
		/* write boundaries of polygon */
		for (ip = 0; ip < polys[ipoly]->np; ip++) {
	    rpcm_to_circ(&polys[ipoly]->rp_(0, ip), &polys[ipoly]->cm[ip], angle);
	    if (angle[0] < 0.) angle[0] += TWOPI;
	    for (i = 0; i < 3; i++) {
				unit = fmt->outunitp;
				if (i > 0 && fmt->outunitp == 'h') unit = 'd';
				scale(&angle[i], 'r', unit);
	    }
	    wrangle(angle[0], fmt->outunitp, fmt->outprecision, az_str, AZEL_STR_LEN);
	    wrangle(angle[1], fmt->outunitp, fmt->outprecision, el_str, AZEL_STR_LEN);
	    wrangle(angle[2], (fmt->outunitp == 'h')? 'd' : fmt->outunitp, fmt->outprecision, th_str, AZEL_STR_LEN);
	    fprintf(file, " %s %s %s", az_str, el_str, th_str);
		}
		fprintf(file, "\n");
		/* increment polygon count */
		npoly++;
	}

	/* warn about polygons with incorrect area */
	if (nbadarea > 0) {
		msg("%d polygons have incorrect area, but kept\n", nbadarea);
	}

	/* advise */
	msg("%d polygons written to %s\n",
			npoly, (file == stdout)? "output": filename);

	/* close file */
	if (file != stdout) fclose(file);

	return(npoly);
}

/*------------------------------------------------------------------------------
  Write mask data in edges or vertices format.

	Input: filename = name of file to write to;
	"" or "-" means write to standard output.
	fmt = pointer to format structure.
	polys = polygons to write.
	npolys = number of polygons.
	npolyw = number of polygons to write.
  Return value: number of polygons written,
	or -1 if error occurred.
*/ 
int wr_edge(filename, fmt, polys, npolys, npolyw)
		 char *filename;
		 format *fmt;
		 int npolys, npolyw;
#ifdef GCC
		 polygon *polys[npolys];
#else
		 polygon *polys[];
#endif
{
	char az_str[AZEL_STR_LEN], el_str[AZEL_STR_LEN];
	int do_vcirc, i, ier, imid, ipoly, iv, ive, ivm, jv, manybounds, nbadverts, nev, nev0, npoly, nv, nve, nvm;
	int *ipv, *ev;
	double az0, tol;
	double *ve, *vm, *angle;
	azel v;
	FILE *file;
	char *vertices_fmt = "vertices %d ( %d vertices, %.15lg weight, %s %s mid):\n";
	char *edges_fmt = "edges %d ( %d points/edge, %d edges, %.15lg weight, %s %s mid):\n";

	/* open filename for writing */
	if (!filename || strcmp(filename, "-") == 0) {
		file = stdout;
	} else {
		file = fopen(filename, "w");
		if (!file) {
	    fprintf(stderr, "wr_edge: cannot open %s for writing\n", filename);
	    return(-1);
		}
	}

	/* whether to write vertices also for circles with no intersections */
	if (strcmp(fmt->out, "vertices") == 0) {
		do_vcirc = 0;
	} else {
		do_vcirc = 1;
	}

	/* write number of polygons */
	fprintf(file, "%d polygons\n", npolyw);

	/* write angular unit */
	fprintf(file, "unit %c\n", fmt->outunitp);

	/* nve must be even to ensure gverts returns a centre point of edge */
	if (fmt->outnve % 2 == 0) { 
		nve = fmt->outnve;
	} else {
		nve = fmt->outnve * 2;
	}

	manybounds = 0;
	npoly = 0;
	nbadverts = 0;
	for (ipoly = 0; ipoly < npolys; ipoly++) {
		/* discard null polygons */
		if (!polys[ipoly]) continue;

		/* points on edges of polygon */
		tol = mtol;
		ier = gverts(polys[ipoly], do_vcirc, &tol, &nv, nve, &ve, &angle, &ipv, &nev, &nev0, &ev);
		if (ier == -1) return(-1);
		if (ier) {
	    nbadverts++;
	    continue;
		}

		/* point somewhere in the middle of the polygon */
		imid = vmid(polys[ipoly], nv, nve, ve, ev, &nvm, &vm);
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
	    rp_to_vert(vm, &v);
	    scale_azel(&v, 'r', fmt->outunitp);
		}

		/* warn about multi-boundary polygon */
		if (nev > 1) {
	    if (WARNMAX > 0 && manybounds == 0) {
				msg("the following polygons have > 1 boundary (not simply-connected)\n");
				msg("   separate boundaries will be split over separate lines:\n");
	    }
	    if (manybounds < WARNMAX) {
				msg(" %d", polys[ipoly]->id);
	    } else if (manybounds == WARNMAX) {
				msg(" ... more\n");
	    }
	    manybounds++;
		}

		/* number of edges, weight, and midpoint of polygon */
		wrangle(v.az, fmt->outunitp, fmt->outprecision, az_str, AZEL_STR_LEN);
		wrangle(v.el, fmt->outunitp, fmt->outprecision, el_str, AZEL_STR_LEN);
		if (strcmp(fmt->out, "vertices") == 0) {
	    fprintf(file, vertices_fmt,
							polys[ipoly]->id, nv, polys[ipoly]->weight, az_str, el_str);
		} else {
	    fprintf(file, edges_fmt,
							polys[ipoly]->id, fmt->outnve, nv, polys[ipoly]->weight, az_str, el_str);
		}

		/* write points, splitting separate boundaries over separate lines */
		for (iv = jv = 0; iv < nv; jv++) {
	    while (iv < ev[jv]) {
				for (ive = 0; ive < nve; ive += nve / fmt->outnve) {
					/* convert unit vector to azel vertex */
					rp_to_vert(&ve[3 * (ive + nve * iv)], &v);
					/* phase points to the first */
					if (iv == 0 && ive == 0) {
						az0 = v.az;
					} else {
						v.az -= rint((v.az - az0) / TWOPI) * TWOPI;
					}
					scale_azel(&v, 'r', fmt->outunitp);
					wrangle(v.az, fmt->outunitp, fmt->outprecision, az_str, AZEL_STR_LEN);
					wrangle(v.el, fmt->outunitp, fmt->outprecision, el_str, AZEL_STR_LEN);
					fprintf(file, " %s %s", az_str, el_str);
				}
				iv++;
	    }
	    fprintf(file, "\n");
		}
		/* increment polygon count */
		npoly++;
	}
	/* warn about multi-boundary polygon */
	if (WARNMAX > 0 && manybounds > 0 && manybounds <= WARNMAX) msg("\n");
	if (manybounds > 0) msg("%d polygons had more than one boundary (not simply-connected)\n", manybounds);

	/* warn about polygons producing fatal error */
	if (nbadverts > 0) {
		msg("%d polygons producing fatal error in gvert discarded\n");
	}

	/* advise */
	msg("%d polygons written to %s\n",
			npoly, (file == stdout)? "output": filename);

	/* close file */
	if (file != stdout) fclose(file);

	return(npoly);
}

/*------------------------------------------------------------------------------
  Write mask data in rectangle format.
  Only polygons which are rectangles are written: other polygons are discarded.

	Input: filename = name of file to write to;
	"" or "-" means write to standard output.
	fmt = pointer to format structure.
	polys = polygons to write.
	npolys = number of polygons.
	npolyw = number of polygons to write.
  Return value: number of polygons written,
	or -1 if error occurred.
*/ 
int wr_rect(filename, fmt, polys, npolys, npolyw)
		 char *filename;
		 format *fmt;
		 int npolys, npolyw;
#ifdef GCC
		 polygon *polys[npolys];
#else
		 polygon *polys[];
#endif
{
	char unit;
	char azmin_str[AZEL_STR_LEN], azmax_str[AZEL_STR_LEN], elmin_str[AZEL_STR_LEN], elmax_str[AZEL_STR_LEN];
	int ier, ipoly, isrect, nbadarea, nrect;
	double area, azmin, azmax, elmin, elmax, tol;
	FILE *file;
	char *rect_fmt = "rectangle %d ( %d caps, %.15lg weight, %.15lf str):\n";

	/* count how many polygons are rectangles */
	nrect = 0;
	for (ipoly = 0; ipoly < npolys; ipoly++) {
		if (!polys[ipoly]) continue;
		isrect = poly_to_rect(polys[ipoly], &azmin, &azmax, &elmin, &elmax);
		if (isrect) nrect++;
	}

	/* no rectangles */
	if (nrect == 0) {
		msg("there are no rectangles among the %d polygons\n", npolyw);
		return(0);
	} else if (nrect < npolyw) {
		msg("%d of %d polygons are rectangles; discarding %d non-rectangle polygons\n", nrect, npolyw, npolyw - nrect);
	} else if (nrect == npolyw) {
		msg("all %d polygons are rectangles\n", npolyw);
	}

	/* open filename for writing */
	if (!filename || strcmp(filename, "-") == 0) {
		file = stdout;
	} else {
		file = fopen(filename, "w");
		if (!file) {
	    fprintf(stderr, "wr_rect: cannot open %s for writing\n", filename);
	    return(-1);
		}
	}

	/* write number of rectangles */
	fprintf(file, "%d rectangles\n", nrect);

	/* write angular unit */
	fprintf(file, "unit %c\n", fmt->outunitp);

	nrect = 0;
	nbadarea = 0;
	/* write rectangles */
	for (ipoly = 0; ipoly < npolys; ipoly++) {
		/* discard null polygons */
		if (!polys[ipoly]) continue;

		/* is polygon a rectangle? */
		isrect = poly_to_rect(polys[ipoly], &azmin, &azmax, &elmin, &elmax);
		/* skip polygons that are not rectangles */
		if (!isrect) continue;

		/* area of polygon */
		tol = mtol;
		ier = garea(polys[ipoly], &tol, &verb, &area);
		if (ier == -1) return(-1);
		if (ier) {
	    fprintf(stderr, "wr_rect: area of polygon %d is incorrect\n", polys[ipoly]->id);
	    nbadarea++;
		}

		/* number of caps, weight, and area of polygon */
		fprintf(file, rect_fmt,
						polys[ipoly]->id, polys[ipoly]->np, polys[ipoly]->weight, area);

		/* scale angles to desired units */
		unit = fmt->outunitp;
		scale(&azmin, 'r', unit);
		scale(&azmax, 'r', unit);
		if (fmt->outunitp == 'h') unit = 'd';
		scale(&elmin, 'r', unit);
		scale(&elmax, 'r', unit);

		/* write rectangle */
		unit = fmt->outunitp;
		wrangle(azmin, fmt->outunitp, fmt->outprecision, azmin_str, AZEL_STR_LEN);
		wrangle(azmax, fmt->outunitp, fmt->outprecision, azmax_str, AZEL_STR_LEN);
		wrangle(elmin, fmt->outunitp, fmt->outprecision, elmin_str, AZEL_STR_LEN);
		wrangle(elmax, fmt->outunitp, fmt->outprecision, elmax_str, AZEL_STR_LEN);
		fprintf(file, " %s %s %s %s\n", azmin_str, azmax_str, elmin_str, elmax_str);

		/* increment rectangle count */
		nrect++;
	}

	/* warn about rectangles with incorrect area */
	if (nbadarea > 0) {
		msg("%d rectangles have incorrect area, but kept\n", nbadarea);
	}

	/* advise */
	msg("%d rectangles written to %s\n",
			nrect, (file == stdout)? "output": filename);

	/* close file */
	if (file != stdout) fclose(file);

	return(nrect);
}

/*------------------------------------------------------------------------------
  Write mask data in polygon format.

	Input: filename = name of file to write to;
	"" or "-" means write to standard output.
	polys = polygons to write.
	npolys = number of polygons.
	npolyw = number of polygons to write.
  Return value: number of polygons written,
	or -1 if error occurred.
*/ 
int wr_poly(filename, polys, npolys, npolyw, selfdestruct)
		 char *filename;
		 int npolys, npolyw, selfdestruct;
#ifdef GCC
		 polygon *polys[npolys];
#else
		 polygon *polys[];
#endif
{
	int ier, ip, ipoly, nbadarea, npoly;
	double area, tol;
	FILE *file;
	char *polygon_fmt = "polygon %d ( %d caps, %.15lg weight, %.15lf str):\n";

	/* open filename for writing */
	if (!filename || strcmp(filename, "-") == 0) {
		file = stdout;
	} else {
		file = fopen(filename, "w");
		if (!file) {
	    fprintf(stderr, "wr_poly: cannot open %s for writing\n", filename);
	    return(-1);
		}
	}

	/* write number of polygons */
	fprintf(file, "%d polygons\n", npolyw);

	npoly = 0;
	nbadarea = 0;
	for (ipoly = 0; ipoly < npolys; ipoly++) {
		/* discard null polygons */
		if (!polys[ipoly]) continue;

		/* area of polygon */
		tol = mtol;
		ier = garea(polys[ipoly], &tol, &verb, &area);
		if (ier == -1) return(-1);
		if (ier) {
	    fprintf(stderr, "wr_poly: area of polygon %d is incorrect\n", polys[ipoly]->id);
	    nbadarea++;
		}

		/* number of caps, weight, and area of polygon */
		fprintf(file, polygon_fmt,
						polys[ipoly]->id, polys[ipoly]->np, polys[ipoly]->weight, area);

		/* write boundaries of polygon */
		for (ip = 0; ip < polys[ipoly]->np; ip++) {
	    fprintf(file, " %19.16f %19.16f %19.16f %.16g\n",
							polys[ipoly]->rp_(0, ip), polys[ipoly]->rp_(1, ip), polys[ipoly]->rp_(2, ip), polys[ipoly]->cm[ip]);
		}

    if(selfdestruct==1) free_poly(polys[ipoly]);

		/* increment polygon count */
		npoly++;
	}

	/* warn about polygons with incorrect area */
	if (nbadarea > 0) {
		msg("%d polygons have incorrect area, but kept\n", nbadarea);
	}

	/* advise */
	msg("%d polygons written to %s\n",
			npoly, (file == stdout)? "output": filename);

	/* close file */
	if (file != stdout) fclose(file);

	return(npoly);
}

/*------------------------------------------------------------------------------
  Write mask data in Max Tegmark's Region format.

	Input: filename = name of file to write to;
	"" or "-" means write to standard output.
	polys = polygons to write.
	npolys = number of polygons.
	npolyw = number of polygons to write.
  Return value: number of polygons written,
	or -1 if error occurred.
*/ 
int wr_Reg(filename, polys, npolys, npolyw)
		 char *filename;
		 int npolys;
#ifdef GCC
		 polygon *polys[npolys];
#else
		 polygon *polys[];
#endif
{
	int ip, ipoly, npoly;
	FILE *file;
	/* there are no holes in our treatment */
	int nholes = 0;
	char *Region_fmt = " Region %d ( %d caps, %d holes):\n";

	/* open filename for writing */
	if (!filename || strcmp(filename, "-") == 0) {
		file = stdout;
	} else {
		file = fopen(filename, "w");
		if (!file) {
	    fprintf(stderr, "wr_Reg: cannot open %s for writing\n", filename);
	    return(-1);
		}
	}

	/* write number of polygons */
	fprintf(file, " %d\n", npolyw);

	/* write number of caps */
	for (ipoly = 0; ipoly < npolys; ipoly++) {
		if (!polys[ipoly]) continue;
		fprintf(file, " %d", polys[ipoly]->np);
		if ((ipoly + 1) % 40 == 0 || ipoly == npolys - 1) fprintf(file, "\n");
	}

	/* write number of holes */
	for (ipoly = 0; ipoly < npolys; ipoly++) {
		if (!polys[ipoly]) continue;
		fprintf(file, " %d", nholes);
		if ((ipoly + 1) % 40 == 0 || ipoly == npolys - 1) fprintf(file, "\n");
	}

	/* write polygons */
	npoly = 0;
	for (ipoly = 0; ipoly < npolys; ipoly++) {
		/* discard null polygons */
		if (!polys[ipoly]) continue;

		/* blank line */
		fprintf(file, "\n");

		/* id line */
		fprintf(file, Region_fmt,
						polys[ipoly]->id, polys[ipoly]->np, nholes);

		/* write boundaries of polygon */
		for (ip = 0; ip < polys[ipoly]->np; ip++) {
	    fprintf(file, " %19.16f %19.16f %19.16f %19.16f\n",
							polys[ipoly]->rp_(0, ip), polys[ipoly]->rp_(1, ip), polys[ipoly]->rp_(2, ip), polys[ipoly]->cm[ip]);
		}

		/* increment polygon count */
		npoly++;
	}

	/* advise */
	msg("%d polygons written to %s\n",
			npoly, (file == stdout)? "output": filename);

	/* close file */
	if (file != stdout) fclose(file);

	return(npoly);
}

/*------------------------------------------------------------------------------
  Discard polygons with weight or area outside specified limits.

	Input: polys = polygons.
	npolys = number of polygons.
  Return value: number of polygons retained,
	or -1 if error occurred.
*/ 
int discard_poly(polys, npolys)
		 int npolys;
#ifdef GCC
		 polygon *polys[npolys];
#else
		 polygon *polys[];
#endif
{
	int discard, ier, ipoly, nbadarea, noutarea, noutweight, npoly;
	double area, tol;

	if (is_weight_min || is_weight_max || is_area_min || is_area_max) {
		noutweight = 0;
		nbadarea = 0;
		noutarea = 0;
		for (ipoly = 0; ipoly < npolys; ipoly++) {
	    discard = 0;

	    /* discard polygons with weights outside interval */
	    if (is_weight_min && is_weight_max) {
				/* min <= max */
				if (weight_min <= weight_max) {
					if (polys[ipoly]->weight < weight_min
							|| polys[ipoly]->weight > weight_max) {
						discard = 1;
					}
					/* min > max */
				} else {
					if (polys[ipoly]->weight < weight_min
							&& polys[ipoly]->weight > weight_max) {
						discard = 1;
					}
				}
	    } else if (is_weight_min) {
				if (polys[ipoly]->weight < weight_min) {
					discard = 1;
				}
	    } else if (is_weight_max) {
				if (polys[ipoly]->weight > weight_max) {
					discard = 1;
				}
	    }
	    if (discard) {
				noutweight++;
				free_poly(polys[ipoly]);
				polys[ipoly] = 0x0;
				continue;
	    }

	    /* area of polygon */
	    tol = mtol;
	    ier = garea(polys[ipoly], &tol, &verb, &area);
	    if (ier == -1) return(-1);
	    if (ier) {
				nbadarea++;

				/* discard polygons with areas outside interval */
	    } else if (is_area_min && is_area_max) {
				/* min <= max */
				if (area_min <= area_max) {
					if (area < area_min
							|| area > area_max) {
						discard = 1;
					}
					/* min > max */
				} else {
					if (area < area_min
							&& area > area_max) {
						discard = 1;
					}
				}
	    } else if (is_area_min) {
				if (area < area_min) {
					discard = 1;
				}
	    } else if (is_area_max) {
				if (area > area_max) {
					discard = 1;
				}
	    }
	    if (discard) {
				noutarea++;
				free_poly(polys[ipoly]);
				polys[ipoly] = 0x0;
				continue;
	    }

		}

		/* warn about discarded polygons */
		if (noutweight > 0) {
	    if (is_weight_min && is_weight_max) {
				if (weight_min < weight_max) {
					msg("%d polygons with weights outside [%g, %g] discarded\n",
							noutweight, weight_min, weight_max);
				} else {
					msg("%d polygons with weights inside (%g, %g) discarded\n",
							noutweight, weight_max, weight_min);
				}
	    } else if (is_weight_min) {
				msg("%d polygons with weights < %g discarded\n",
						noutweight, weight_min);
	    } else if (is_weight_max) {
				msg("%d polygons with weights > %g discarded\n",
						noutweight, weight_max);
	    }
		}
		if (noutarea > 0) {
	    if (is_area_min && is_area_max) {
				if (area_min < area_max) {
					msg("%d polygons with areas outside [%g, %g] discarded\n",
							noutarea, area_min, area_max);
				} else {
					msg("%d polygons with areas inside (%g, %g) discarded\n",
							noutarea, area_max, area_min);
				}
	    } else if (is_area_min) {
				msg("%d polygons with areas < %g discarded\n",
						noutarea, area_min);
	    } else if (is_area_max) {
				msg("%d polygons with areas > %g discarded\n",
						noutarea, area_max);
	    }
		}

	}

	/* count non-null polygons */
	npoly = 0;
	for (ipoly = 0; ipoly < npolys; ipoly++) {
		if (polys[ipoly]) npoly++;
	}

	return(npoly);
}
