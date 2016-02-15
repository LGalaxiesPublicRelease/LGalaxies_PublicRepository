#include <stdio.h>
#include <stdlib.h>
#include "polygon.h"

#define MEG(i)		(double)((i+999)/1000)/1000.

static long memory = 0, femory = 0;
static int mpoly = 0, fpoly = 0;

/* local functions */
int room_poly();
polygon *new_poly();
void free_poly();

/* external functions */
extern void copy_poly();
extern void msg(char *, ...);

/*------------------------------------------------------------------------------
  Allocate memory for polygon of np caps.

	Input: npmax = desired number of caps of polygon.
  Return value: pointer to a new polygon,
	or null if failed to allocate memory.
*/ 
polygon *new_poly(npmax)
		 int npmax;
{
	polygon *poly;

	/* allocate memory for new polygon */
	poly = (polygon *) malloc(sizeof(polygon));
	if (!poly) return(0x0);
	mpoly++;
	memory += sizeof(polygon);

	/* allocate new rp array */
	poly->rp = (double *) malloc(sizeof(double) * npmax * 3);
	if (!poly->rp) return(0x0);
	memory += sizeof(double) * npmax * 3;

	/* allocate new cm array */
	poly->cm = (double *) malloc(sizeof(double) * npmax);
	if (!poly->cm) return(0x0);
	memory += sizeof(double) * npmax;

	/* allocated number of caps of polygon */
	poly->npmax = npmax;

	/* parents */
	poly->parent_polys=0x0;
	poly->nparents=0;

	return(poly);
}

/*------------------------------------------------------------------------------
  Free polygon memory.
*/ 
void free_poly(poly)
		 polygon *poly;
{
	if (poly) {
		fpoly++;
		femory += sizeof(polygon);
		if (poly->rp) {
	    free(poly->rp);
	    poly->rp = 0x0;
	    femory += sizeof(double) * poly->npmax * 3;
		}
		if (poly->cm) {
	    free(poly->cm);
	    poly->cm = 0x0;
	    femory += sizeof(double) * poly->npmax;
		}
		if(poly->parent_polys) 
			free(poly->parent_polys);
		poly->parent_polys=0x0;  
		poly->nparents=0;
		free(poly);
	}
}

/*------------------------------------------------------------------------------
  Test whether a polygon contains enough space;
  if not, free the polygon, and allocate a new polygon with enough space,
  plus a bit extra to allow for subsequent expansion.

	Input: *poly = pointer to polygon.
	np = desired number of caps of polygon.
	dnp = number of extra caps to allocate.
	save = 0 to discard contents of poly when making new poly;
	1 to copy contents of poly to new poly.
	If poly already contains enough space, it is left unchanged.
  Return value: -1 if could not allocate new memory;
	here *poly is left unchanged, so use the return value,
	not whether *poly is null, to test for failure;
	0 if polygon contained enough space;
	1 if new polygon was allocated.
*/
int room_poly(poly, np, dnp, save)
		 polygon **poly;
		 int np, dnp, save;
{
	polygon *newpoly;

	/* polygon contains enough space */
	if (*poly && (*poly)->npmax >= np) return(0);

	/* allocate new polygon with np + dnp caps */
	newpoly = new_poly(np + dnp);
	if (!newpoly) return(-1);

	/* copy poly to new polygon */
	if (*poly && save) copy_poly(*poly, newpoly);

	/* free polygon */
	if (*poly) free_poly(*poly);
	(*poly)=0x0;

	/* point poly to new polygon */
	*poly = newpoly;

	return(1);
}

/*------------------------------------------------------------------------------
  Advise memory used.
*/ 
void memmsg()
{
	if (mpoly > 0)
		msg("%d polygons (%.3fMb) allocated, %d (%.3fMb) freed\n",
				mpoly, MEG(memory), fpoly, MEG(femory));
}
