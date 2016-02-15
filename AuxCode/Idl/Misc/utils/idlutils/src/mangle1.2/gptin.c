#include "logical.h"
#include "polygon.h"

/* local functions */
int gptin();

/* external functions */
extern logical gptin_();

/*------------------------------------------------------------------------------
  Determine whether unit vector lies inside polygon.

  This is a c interface to fortran logical function gptin.

   Input: poly is a polygon.
	  rp = unit vector.
  Return value: 1 if in;
		0 if not in.
*/
int gptin(poly, rp)
double rp[3];
polygon *poly;
{
    /* fortran routine */
    if (gptin_(poly->rp, poly->cm, &poly->np, rp)) return(1);
    return(0);
}
