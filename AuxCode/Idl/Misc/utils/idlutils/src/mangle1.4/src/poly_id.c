/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Id numbers of polygons containing position az, el.

   Input: poly = array of pointers to npoly polygons.
	  npoly = number of polygons in poly array.
	  az, el = angular position in radians.
  Output: id_p = pointer to array containing id numbers of polygons;
		 the required memory is allocated.
  Return value: number of polygons that contain az, el position.
*/
int poly_id(int npoly, polygon *poly[/*npoly*/], double az, double el, int **id_p)
{
/* number of extra polygon id numbers to allocate, to allow for expansion */
#define DNID		16
    static int nidmax = 0;
    static int *id = 0x0;

    int ipoly, nid;
    double rp[3];

    /* unit vector corresponding to angular position az, el */
    rp[0] = cos(el) * cos(az);
    rp[1] = cos(el) * sin(az);
    rp[2] = sin(el);

    nid = 0;

    /* keep trying till the id array is big enough */
    do {
        /* make sure that allocated id array contain enough space */
        if (!id || nid > nidmax) {
            if (id) free(id);
            id = (int *) malloc(sizeof(int) * (nid + DNID));
            if (!id) {
                fprintf(stderr, "poly_id: failed to allocate memory for %d ints\n", nid + DNID);
                return(-1);
            }
	    nidmax = nid + DNID;
	}

	nid = 0;
	/* do each polygon in turn */
	for (ipoly = 0; ipoly < npoly; ipoly++) {
	    /* id number of each polygon that contains az, el position */
	    if (gptin(poly[ipoly], rp)) {
		if (nid < nidmax) id[nid] = poly[ipoly]->id;
		nid++;
	    }
	}

    } while (nid > nidmax);

    /* point id_p at id array */
    *id_p = id;

    /* number of polygons containing az, el position */
    return(nid);
}
