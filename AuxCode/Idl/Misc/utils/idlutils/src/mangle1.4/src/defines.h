/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#ifndef DEFINES_H
#define DEFINES_H

#define	MAXINT		(((unsigned int)-1) / 2)

/* maximum number of polygons */
/*
  This is the only hard limit built into mangle.
  It's supposed to be a feature, not a bug.
  If you are making zillions of polygons, chances are it's a silly mistake.
*/
#define NPOLYSMAX	100000000

/* number of extra caps to allocate to polygon, to allow for later splitting */
#define DNP		4

#include "pi.h"
#define TWOPI		(2. * PI)
#define PIBYTWO		(PI / 2.)

/* list of possible units */
#define UNITS		"rdmsh"
/* #define UNITS	"rd°m'´s\"¨h" */

/* angular units in arcseconds */
#define RADIAN		(648000. / PI)
#define HOUR		54000.
#define DEGREE		3600.
#define MINUTE		60.
#define SECOND		1.

/* default maximum harmonic */
#define LMAX		0
/* default smoothing harmonic (0. = no smooth) */
#define LSMOOTH		0.
/* default smoothing exponent (2. = gaussian) */
#define ESMOOTH		2.
/* default snap angles for axis, latitude, and edge */
#define AXTOL		2.
#define BTOL		2.
#define THTOL		2.
/* default value of ytol */
#define YTOL		.01
/* default snap angle for multiple intersections */
#define MTOL		.00001
/* default input units of snap angles */
#define AXUNIT		's'
#define BUNIT		's'
#define THUNIT		's'
#define MUNIT		's'
/* default seed for random number generator */
#define SEED		1
/* default number of random points to generate */
#define NRANDOM		1
/* default number of points per edge */
#define	NVE		2
/* default input unit of polygon data is degrees */
#define INUNITP		'd'
/* default output unit of polygon data is degrees */
#define OUTUNITP	'd'
/* default input unit of az, el data is degrees */
#define INUNIT		'd'
/* default output unit of az, el data is degrees */
#define OUTUNIT		'd'
/* default output phase: '+' or '-' to make output azimuth in interval (-pi, pi] or [0, 2 pi) */
#define	OUTPHASE	'+'
/* identity transformation between angular frames */
#define AZN		0.
#define	ELN		90.
#define	AZP		180.
/* unit of transformation between angular frames
   must be 'd', since degrees is hard-wired into transformation routines */
#define TRUNIT		'd'
/* default number of characters of line of input data to skip */
#define SKIP		0
/* default last character to read from line of input data */
#define END		0

/* list of possible input formats */
#define RFMTS		"ceprRsv"
/* list of possible output formats */
#define WFMTS		"acegimprRsvw"

/* possible polygon file formats */
#define AREA		0
#define CIRCLE		1
#define EDGES		2
#define	GRAPHICS	3
#define	ID		4
#define	MIDPOINT	5
#define POLYGON		6
#define RECTANGLE	7
#define REGION		8
#define SPOLYGON	9
#define	VERTICES	10
#define	WEIGHT		11
#define BINARY		15

#endif	/* DEFINES_H */
