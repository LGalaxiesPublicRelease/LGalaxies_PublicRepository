#ifndef DEFINES_H
#define DEFINES_H

/* maximum number of polygons */
/*
  This is the only hard limit built into mangle.
  It's supposed to be a feature, not a bug.
  If you are making zillions of polygons, chances are it's a silly mistake.
*/
/* #define NPOLYSMAX	65536 */
/* a suitably huge alternative value */
/* #define NPOLYSMAX	524288  */
#define NPOLYSMAX	10000000

/* number of extra caps to allocate to polygon, to allow for later splitting */
#define DNP		4

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
#define MTOL		0.
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
#define RFMTS		"ceprRv"
/* list of possible output formats */
#define WFMTS		"ceprRv"

#endif	/* DEFINES_H */
