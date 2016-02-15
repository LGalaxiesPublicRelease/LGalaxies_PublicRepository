#include <stdio.h>
#include <string.h>
#include <math.h>
#include "export.h" 

/*
 * Checks separation distance
 */

#define RAD2DEG 57.29577951

double separation(double x1,
									double y1,
									double z1,
									double x2,
									double y2,
									double z2)
{
	double dx,dy,dz;
	
	dx=x1-x2;
	dy=y1-y2;
	dz=z1-z2;
	return(RAD2DEG*2.*asin(0.5*sqrt(dx*dx+dy*dy+dz*dz)));
} /* end separation */
