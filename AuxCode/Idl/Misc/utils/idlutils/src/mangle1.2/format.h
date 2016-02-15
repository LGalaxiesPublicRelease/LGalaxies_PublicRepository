#ifndef FORMAT_H
#define FORMAT_H

/*
  Structure defining format of data.
  If your format needs more descriptives, shove 'em in here.
  DON'T FORGET TO SET THE DEFAULT VALUES IN defaults.h,
  AND TO UPDATE copy_format().
*/ 
typedef struct {
    char *in;		/* keyword defining the input data format */
    char *out;		/* keyword defining the output data format */
    size_t skip;	/* skip first skip characters of each line */
    size_t end;		/* read only up to end'th character of each line */
    char single;	/* keyword defines precisely one polygon */
    int n;		/* the number of thingys defined by keyword */
    int nn;		/* the number of thingys per thingy */
    int innve;		/* the input number of points per edge */
    int outnve;		/* the output number of points per edge */
    int id;		/* id number of current polygon */
    char newid;		/* whether to use old or new id number */
    double weight;	/* weight of current polygon */
    char inunitp;	/* angular units of input polygon data */
    char outunitp;	/* angular units of output polygon data */
    int inframe;	/* angular frame of input az, el data */
    int outframe;	/* angular frame of output az, el data */
    char inunit;	/* angular units of input az, el data */
    char outunit;	/* angular units of output az, el data */
    int outprecision;	/* digits after decimal point in output angles */
    double azn;		/* azimuth of new pole wrt original frame */
    double eln;		/* elevation of new pole wrt original frame
			 = elevation of original pole wrt new frame */
    double azp;		/* azimuth of original pole wrt new frame */
	char trunit;	/* angular units of transformation angles */
	char *linklist;	/* input file with links to check */
	char *parents;	/* output file with parents of balkans */
} format;

#define REGION		0
#define POLYGON		1
#define CIRCLE		2
#define EDGES		3
#define RECTANGLE	4
#define	VERTICES	5

#endif	/* FORMAT_H */
