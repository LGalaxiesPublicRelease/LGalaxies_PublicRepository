#include <stdio.h>
#include <string.h>

/* local functions */ 
int rdangle();

/*------------------------------------------------------------------------------
  Read angle from word.

   Input: word = pointer to string.
	  unit = units of angle;
		 this is needed to detect 'h' unit [hms(RA) & dms(Dec)],
		 not to change the units of the angle, which is unchanged.
  Output: *next = pointer to character immediately after word;
		  *next may be the same pointer as word.
	  angle = angle.
  Return value: 1 = ok,
		otherwise error.
*/
int rdangle(word, next, unit, angle)
char *word, **next;
char unit;
double *angle;
{
    static char *blank = " \t";
    static char *number = "+-.0123456789eE";

    int i, ird, sgn;
    unsigned int deg, min;
    double sec;
    char *ch;

    /* find first non-blank character */
    ch = word;
    while (*ch && strchr(blank, *ch)) ch++;
    if (!*ch) return(0);

    /* units are hours, minutes, seconds (Right Ascension)
       and degrees, minutes, seconds (Declination) */
    if (unit == 'h') {
	sgn = 1;
	if (*ch == '-') {
	    sgn = -1;
	    ch++;
	}
	ird = sscanf(ch, "%d %d %lf", &deg, &min, &sec);
	*angle = sgn * ((sec/60. + min)/60. + deg);
	for (i = 0; i < ird; i++) {
	    while (*ch && strchr(blank, *ch)) ch++;
	    while (*ch && strchr(number, *ch)) ch++;
	}
	if (ird == 3) {
	   ird = 1;
	} else {
	   ird = 0;
	}
    } else {
	ird = sscanf(ch, "%lf", angle);
	while (*ch && strchr(number, *ch)) ch++;
    }

    /* point next at character after word */
    *next = ch;

    return(ird);
}
