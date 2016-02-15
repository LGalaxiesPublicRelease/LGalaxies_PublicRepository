#ifndef SCALE_H
#define SCALE_H

#include "pi.h"
#define TWOPI		(2. * PI)

/* list of possible units */
#define UNITS		"rdmsh"
/* #define UNITS	"rd°m'´s\"¨h" */

/* angular units in arcseconds */
#define RADIAN		(648000. / PI)
#define HOUR		54000.
#define DEGREE		3600.
#define MINUTE		60.
#define SECOND		1.

#endif	/* SCALE_H */
