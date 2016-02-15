/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include "manglefn.h"
#include <stdarg.h>

extern int verbose;

/*------------------------------------------------------------------------------
  Print messages.
*/
void msg(char *fmt, ...)
{
    va_list args;

    if (verbose) {
	va_start(args, fmt);
	vfprintf(stderr,fmt, args);
	va_end(args);
	fflush(stderr);
    }
}
