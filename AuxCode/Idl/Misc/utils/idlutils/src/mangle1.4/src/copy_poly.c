/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdlib.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Copy poly1.
  poly may be the same polygon as poly1.

   Input: poly1 is a polygon.
  Output: poly is a copy of poly1.
*/
void copy_poly(polygon *poly1, polygon *poly)
{
    int i, ip;

    for (ip = 0; ip < poly1->np; ip++) {
      for (i = 0; i < 3; i++) {
        poly->rp[ip][i] = poly1->rp[ip][i];
      }
      poly->cm[ip] = poly1->cm[ip];
    }
    poly->np = poly1->np;
    poly->id = poly1->id;
    poly->weight = poly1->weight;
    poly->nparents=0;
    poly->maxparents=0;
    if(poly->parent_polys!=0x0) {
      free(poly->parent_polys);
      poly->parent_polys=0x0;
    }
    for(i=0;i<poly1->nparents;i++)
      add_parent(poly,poly1->parent_polys[i]);
}

/*------------------------------------------------------------------------------
  Copy first np caps of poly1 into poly.
  poly may be the same polygon as poly1.

   Input: poly1 is a polygon.
	  np = number of caps of poly1 to copy into poly.
  Output: poly is a copy of the first np caps of poly1.
*/
void copy_polyn(polygon *poly1, int np, polygon *poly)
{
    int i, ip;

    for (ip = 0; ip < np; ip++) {
      for (i = 0; i < 3; i++) {
        poly->rp[ip][i] = poly1->rp[ip][i];
      }
      poly->cm[ip] = poly1->cm[ip];
    }
    poly->np = np;
    poly->id = poly1->id;
    poly->weight = poly1->weight;
    poly->nparents=0;
    poly->maxparents=0;
    if(poly->parent_polys!=0x0) {
      free(poly->parent_polys);
      poly->parent_polys=0x0;
    }
    for(i=0;i<poly1->nparents;i++)
      add_parent(poly,poly1->parent_polys[i]);
}

/*------------------------------------------------------------------------------
  Intersection of two polygons.
  poly may be the same polygon as poly1, but not the same polygon as poly2.

   Input: poly1, poly2 are two polygons.
  Ouptput: poly is the intersection of poly1 and poly2,
	   The id number and weight of poly are taken equal to those of poly1.
*/
void poly_poly(polygon *poly1, polygon *poly2, polygon *poly)
{
    int i, ip;

    for (ip = 0; ip < poly1->np; ip++) {
	for (i = 0; i < 3; i++) {
	    poly->rp[ip][i] = poly1->rp[ip][i];
	}
	poly->cm[ip] = poly1->cm[ip];
    }
    for (ip = 0; ip < poly2->np; ip++) {
	for (i = 0; i < 3; i++) {
	    poly->rp[ip + poly1->np][i] = poly2->rp[ip][i];
	}
	poly->cm[ip + poly1->np] = poly2->cm[ip];
    }

    poly->np = poly1->np + poly2->np;
    poly->id = poly1->id;
    poly->weight = poly1->weight;
    poly->nparents=0;
    poly->maxparents=0;
    if(poly->parent_polys!=0x0) {
      free(poly->parent_polys);
      poly->parent_polys=0x0;
    }
    for(i=0;i<poly1->nparents;i++)
      add_parent(poly,poly1->parent_polys[i]);
    for(i=0;i<poly2->nparents;i++)
      add_parent(poly,poly2->parent_polys[i]);
}

/*------------------------------------------------------------------------------
  Intersection of a polygon
  with the n'th cap, or the complement of the n'th cap, of another polygon.
  poly may be the same polygon as poly1, but not the same polygon as poly2.

   Input: poly1, poly2 are polygons.
	  n = intersect poly1 with n'th cap of poly2.
	  scm =  1 to intersect with n'th cap of poly2;
		-1 to intersect with complement of n'th cap of poly2.
  Output: poly is the intersection of poly1 and n'th cap of poly2.
	  The id number and weight of poly are taken equal to those of poly1.
*/
void poly_polyn(polygon *poly1, polygon *poly2, int n, int scm, polygon *poly)
{
    int i, ip;

    for (ip = 0; ip < poly1->np; ip++) {
	for (i = 0; i < 3; i++) {
	    poly->rp[ip][i] = poly1->rp[ip][i];
	}
	poly->cm[ip] = poly1->cm[ip];
    }
    ip = poly1->np;
    for (i = 0; i < 3; i++) {
	poly->rp[ip][i] = poly2->rp[n][i];
    }
    poly->cm[ip] = scm * poly2->cm[n];

    poly->np = poly1->np + 1;
    poly->id = poly1->id;
    poly->weight = poly1->weight;

    poly->nparents=0;
    poly->maxparents=0;
    if(poly->parent_polys!=0x0) {
      free(poly->parent_polys);
      poly->parent_polys=0x0;
    }
    for(i=0;i<poly1->nparents;i++)
      add_parent(poly,poly1->parent_polys[i]);
    for(i=0;i<poly2->nparents;i++)
      add_parent(poly,poly2->parent_polys[i]);
}

/*------------------------------------------------------------------------------
  Make group polygon gpoly.
  gpoly may be the same polygon as poly.

   Input: poly = pointer to polygon.
	  gp[ip], ip=1,np = group number of cap ip.
	  gpg = group to put in gpoly.
	  gpoly = pointer to group polgyon.
*/
void group_poly(polygon *poly, int gp[/*poly->np*/], int gpg, polygon *gpoly)
{
    int i, ip, np;

    /* number of caps of group polygon */
    np = 0;
    for (ip = 0; ip < poly->np; ip++) {
      if (gp[ip] == gpg) {
        np++;
      }
    }
    
    /* make group polygon */
    np = 0;
    for (ip = 0; ip < poly->np; ip++) {
      if (gp[ip] == gpg) {
        for (i = 0; i < 3; i++) {
          gpoly->rp[np][i] = poly->rp[ip][i];
        }
        gpoly->cm[np] = poly->cm[ip];
        np++;
      }
    }
    gpoly->np = np;
    gpoly->id = poly->id;
    gpoly->weight = poly->weight;
    gpoly->nparents=0;
    gpoly->maxparents=0;
    if(gpoly->parent_polys!=0x0) {
      free(gpoly->parent_polys);
      gpoly->parent_polys=0x0;
    }
    for(i=0;i<poly->nparents;i++)
      add_parent(gpoly,poly->parent_polys[i]);
    for(i=0;i<poly->nparents;i++)
      add_parent(gpoly,poly->parent_polys[i]);
}
