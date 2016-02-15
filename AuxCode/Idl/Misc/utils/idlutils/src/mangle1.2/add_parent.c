#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "polygon.h"

/* add a parent to the polygon parent list */
void add_parent(polygon *poly, int parent) 
{
	int i;
	for(i=0;i<poly->nparents;i++) 
		if(poly->parent_polys[i]==parent)
			return;
	if(poly->parent_polys==0x0) 
		poly->parent_polys= (int *) malloc(1*sizeof(int));
	else 
		poly->parent_polys= (int *)
			realloc((void *) poly->parent_polys, (poly->nparents+1)*sizeof(int));
	poly->parent_polys[poly->nparents]=parent;
	poly->nparents++;
} /* end add_parent */

/* take a parent out of the polygon parent list */
void trim_parent(polygon *poly, int parent) 
{
	int i,*tmp_parent_polys, tmp_nparents;
	tmp_parent_polys=poly->parent_polys;
	tmp_nparents=poly->nparents;
	poly->parent_polys=0x0;
	poly->nparents=0;
	for(i=0;i<tmp_nparents;i++)
		if(tmp_parent_polys[i]!=parent)
			add_parent(poly,tmp_parent_polys[i]);
	free(tmp_parent_polys);
	tmp_parent_polys=0x0;
	tmp_nparents=0;
} /* end add_parent */

