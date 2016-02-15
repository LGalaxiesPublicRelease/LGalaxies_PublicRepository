#include <stdio.h>
#include <stdlib.h>
#include "polygon.h"

/* outside functions */
extern void add_parent();

/* local functions */
void copy_poly(), poly_poly(), poly_polyn();

/*------------------------------------------------------------------------------
  Copy poly1.

	Input: poly1 is a polygon.
  Output: poly is a copy of poly1.
*/
void copy_poly(poly1, poly)
		 polygon *poly1, *poly;
{
	int i, ip;
	
	for (ip = 0; ip < poly1->np; ip++) {
		for (i = 0; i < 3; i++) {
			poly->rp_(i, ip) = poly1->rp_(i, ip);
		}
		poly->cm[ip] = poly1->cm[ip];
	}
	poly->np = poly1->np;
	poly->id = poly1->id;
	poly->weight = poly1->weight;
	poly->nparents=0;
	for(i=0;i<poly1->nparents;i++)
		add_parent(poly,poly1->parent_polys[i]);
		
}

/*------------------------------------------------------------------------------
  Intersection of two polygons.

	Input: poly1, poly2 are two polygons.
  Ouptput: poly is the intersection of poly1 and poly2,
	The id number and weight of poly are taken equal to those of poly1.
*/
void poly_poly(poly1, poly2, poly)
		 polygon *poly1, *poly2, *poly;
{
	int i, ip;

	for (ip = 0; ip < poly1->np; ip++) {
		for (i = 0; i < 3; i++) {
	    poly->rp_(i, ip) = poly1->rp_(i, ip);
		}
		poly->cm[ip] = poly1->cm[ip];
	}
	for (ip = 0; ip < poly2->np; ip++) {
		for (i = 0; i < 3; i++) {
	    poly->rp_(i, ip + poly1->np) = poly2->rp_(i, ip);
		}
		poly->cm[ip + poly1->np] = poly2->cm[ip];
	}

	poly->np = poly1->np + poly2->np;
	poly->id = poly1->id;
	poly->weight = poly1->weight;
	poly->nparents=0;
	for(i=0;i<poly1->nparents;i++)
		add_parent(poly,poly1->parent_polys[i]);
	for(i=0;i<poly2->nparents;i++)
		add_parent(poly,poly2->parent_polys[i]);
}

/*------------------------------------------------------------------------------
  Intersection of a polygon with the n'th cap of another polygon.

	Input: poly1, poly2 are polygons.
  Output: poly is the intersection of poly1 and n'th cap of poly2.
	The id number and weight of poly are taken equal to those of poly1.
*/
void poly_polyn(poly1, poly2, n, poly)
		 polygon *poly1, *poly2, *poly;
		 int n;
{
	int i, ip;

	for (ip = 0; ip < poly1->np; ip++) {
		for (i = 0; i < 3; i++) {
	    poly->rp_(i, ip) = poly1->rp_(i, ip);
		}
		poly->cm[ip] = poly1->cm[ip];
	}
	ip = poly1->np;
	for (i = 0; i < 3; i++) {
		poly->rp_(i, ip) = poly2->rp_(i, n);
	}
	poly->cm[ip] = poly2->cm[n];

	poly->np = poly1->np + 1;
	poly->id = poly1->id;
	poly->weight = poly1->weight;
	poly->nparents=0;
	for(i=0;i<poly1->nparents;i++)
		add_parent(poly,poly1->parent_polys[i]);
	for(i=0;i<poly2->nparents;i++)
		add_parent(poly,poly2->parent_polys[i]);
}
