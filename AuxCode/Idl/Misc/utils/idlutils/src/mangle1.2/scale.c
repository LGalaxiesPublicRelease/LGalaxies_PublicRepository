#include "scale.h"
#include "vertices.h"

/* local functions */ 
void scale(), scale_azel(), scale_vert();

/*------------------------------------------------------------------------------
  Convert angle from specified unit to specified unit.
*/
void scale(angle, from, to)
double *angle;
char from, to;
{
    /* scale from specified unit to arcseconds */
    switch (from) {
    case 'r':
	*angle = *angle * RADIAN;
	break;
    case 'h':
	*angle = *angle * HOUR;
	break;
    case 'd':
    case '°':
	*angle = *angle * DEGREE;
	break;
    case 'm':
    case '\'':
    case '´':
	*angle = *angle * MINUTE;
	break;
    case 's':
    case '"':
    case '¨':
    default:
	*angle = *angle * SECOND;
	break;
    }

    /* scale from arcseconds to specified unit */
    switch (to) {
    case 'r':
	*angle = *angle / RADIAN;
	break;
    case 'h':
	*angle = *angle / HOUR;
	break;
    case 'd':
    case '°':
	*angle = *angle / DEGREE;
	break;
    case 'm':
    case '\'':
    case '´':
	*angle = *angle / MINUTE;
	break;
    case 's':
    case '"':
    case '¨':
    default:
	*angle = *angle / SECOND;
	break;
    }
}

/*------------------------------------------------------------------------------
  Scale azel structure to desired units.
*/
void scale_azel(v, from, to)
azel *v;
char from, to;
{
    scale(&v->az, from, to);
    scale(&v->el, (from == 'h')? 'd' : from, (to == 'h')? 'd' : to);
}

/*------------------------------------------------------------------------------
  Scale vertices structure to desired units.
*/
void scale_vert(vert, from, to)
vertices *vert;
char from, to;
{
    int iv;

    for (iv = 0; iv < vert->nv; iv++) {
	scale_azel(&vert->v[iv], from, to);
    }
}
