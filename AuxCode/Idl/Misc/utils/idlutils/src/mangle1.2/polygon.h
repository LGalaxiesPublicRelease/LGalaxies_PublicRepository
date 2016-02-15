#ifndef POLYGON_H
#define POLYGON_H

/* How is one supposed to allocate memory dynamically for a multi-dimensional
   array inside a structure?  The following hack works around the problem,
   and is consistent with what fortran wants.
   The underline is because Sun goes for a loop without it.  */
#define rp_(i, j)		rp[(i) + 3 * (j)]

typedef struct {		/* polygon structure */
	int np;			/* number of caps of polygon */
	int npmax;			/* dimension of allocated rp and cm arrays */
	double *rp;			/* pointer to array rp(3, np) of axis coords */
	double *cm;			/* pointer to array cm[np] of 1 - cos(theta) */
	int id;			/* id number of polygon */
	double weight;		/* weight of polygon */
	int *parent_polys;  /* for mrb_balkanize, track the indices of the 
												 original set of polygons */
	int nparents;
} polygon;

#endif	/* POLYGON_H */
