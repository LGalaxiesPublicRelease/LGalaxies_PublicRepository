#include <math.h>

/*------------------------------------------------------------------------------
  Round x to n decimal places.
*/
double places(x, n)
int n;
double x;
{
    int i;

    for (i = 0; i < n; i++) x = x * 10.;
    x = rint(x);
    for (i = 0; i < n; i++) x = x / 10.;

    return(x);
}
