/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include "defines.h"

/* LOCAL VARIABLES */

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

/* write only summary to output */
static int summary = 0;

/* data format */
static format fmt = {
	0x0,		/* keyword defining the input data format */
	0x0,		/* default format of the output data */
	SKIP,		/* skip first skip characters of each line of data */
	END,		/* last character to read from line of data */
	0,		/* keyword does not define precisely one polygon */
	0,		/* the number of thingys defined by keyword */
	0,		/* the number of numbers per thingy */
	NVE,		/* the input number of points per edge */
	0,		/* controls interpetation of nve */
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
	-1,		/* digits after decimal point in output angles (-1 = automatic) */
	OUTPHASE,	/* '+' or '-' to make output azimuth in interval (-pi, pi] or [0, 2 pi) */
	AZN,		/* default			       */
	ELN,		/*	transformation		       */
	AZP,		/* 		between angular frames */
	TRUNIT		/* unit of transformation angles */
};

/* GLOBAL VARIABLES */

/* default is to be verbose */
#ifdef DEBUG
int verbose = 2;
#else
int verbose = 1;
#endif

/* tolerances */
double axtol = AXTOL;		/* snap angle for axis */
char axunit = AXUNIT;		/* unit of snap angle for axis */
double btol = BTOL;		/* snap angle for latitude */
char bunit = BUNIT;		/* unit of snap angle for latitude */
double thtol = THTOL;		/* snap angle for edge */
char thunit = THUNIT;		/* unit of snap angle for edge */
double ytol = YTOL;		/* edge to length tolerance */
double mtol = MTOL;		/* tolerance angle for multiple intersections */
char munit = MUNIT;		/* unit of tolerance angle for multiple intersections */

/* whether min, max weight are turned on */
int is_weight_min = 0;
int is_weight_max = 0;
/* min, max weight to keep */
double weight_min;
double weight_max;

/* whether min, max area are turned on */
int is_area_min = 0;
int is_area_max = 0;
/* min, max area to keep */
double area_min;
double area_max;

/* whether to take intersection of polygons in input files */
int intersect = 0;

/* dictionary of keywords */
/* THE NAMES OF FORMATS MUST AGREE WITH THE INDICES IN defines.h */
char *keywords[] = {
    "area",
    "circle",
    "edges",
    "graphics",
    "id",
    "midpoint",
    "polygon",
    "rectangle",
    "Region",
    "spolygon",
    "vertices",
    "weight",
    "skip",
    "end",
    "unit",
    "binary",
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
