#include "defines.h"

/* maximum harmonic */
static int lmax = LMAX;
/* smoothing parameters */
static double lsmooth = LSMOOTH, esmooth = ESMOOTH;

/* name of file containing harmonics */
static char *Wlm_filename = 0x0;

/* name of survey */
static char *survey = 0x0;

/* option in -f<fopt> command line switch */
static char *fopt = 0x0;

/* seed for random number generator */
static unsigned int seed = SEED;
/* whether seed was read */
static int seed_read = 0;

/* number of random points to generate */
static int nrandom = NRANDOM;

/* initial data format */
static format fmt = {
	0x0,		/* keyword defining the input data format */
	0x0,		/* default format of the output data */
	SKIP,		/* skip first skip characters of each line of data */
	END,		/* last character to read from line of data */
	0,		/* keyword does not define precisely one polygon */
	0,		/* the number of thingys defined by keyword */
	0,		/* the number of numbers per thingy */
	NVE,		/* the input number of points per edge */
	NVE,		/* the output number of points per edge */
	0,		/* id number of current polygon */
	'o',		/* whether to use old or new id number */
	1.,		/* weight of current polygon */
	INUNITP,	/* default unit of angles in input polygon data */
	OUTUNITP,	/* default unit of angles in output polygon data */
	0,		/* angular frame of input az, el data */
	0,		/* angular frame of output az, el data */
	INUNIT,		/* default unit of input az, el data */
	OUTUNIT,	/* default unit of output az, el data */
	-1,		/* digits after decimal point in output angles */
	AZN,		/* default			       */
	ELN,		/*	transformation		       */
	AZP,		/* 		between angular frames */
	TRUNIT		/* unit of transformation angles */
};

/* global variables */

/* default is to be verbose */
#ifdef DEBUG
int verbose = 2;
#else
int verbose = 1;
#endif

/* snap angles for axis, latitude, and edge */
double axtol = AXTOL, btol = BTOL, thtol = THTOL;
/* units of snap angles for axis, latitude, and edge */
char axunit = AXUNIT, bunit = BUNIT, thunit = THUNIT;
/* snap edge to boundary only if closer than ytol x edge length */
double ytol = YTOL;
/* snap angle for multiple intersections */
double mtol = MTOL;
/* unit of snap angle for multiple intersections */
char munit = MUNIT;

/* whether min, max weight are turned on */
int is_weight_min = 0, is_weight_max = 0;
/* min, max weight to keep */
double weight_min, weight_max;

/* whether min, max area are turned on */
int is_area_min = 0, is_area_max = 0;
/* min, max area to keep */
double area_min, area_max;

/* whether to take intersection of polygons in input files */
int intersect = 0;

/* dictionary of keywords */
char *keywords[] = {
    "Region",
    "polygon",
    "circle",
    "edges",
    "rectangle",
    "vertices",
    "skip",
    "end",
    "unit",
    '\0'
};

/* dictionary of frames */
/* THE ORDER OF FRAMES MUST AGREE WITH THAT IN frames.par ! */
char *frames[] = {
    "unknown",
    "eqB1950",
    "eqJ2000",
    "galactic",
    "ecliptic",
    "ecliptic2k",
    "sdss",
    '\0'
};
