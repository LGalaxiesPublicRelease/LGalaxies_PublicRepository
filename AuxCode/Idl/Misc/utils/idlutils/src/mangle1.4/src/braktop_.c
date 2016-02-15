/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include "manglefn.h"

/*------------------------------------------------------------------------------
  c interface to fortran subroutines in braktop.s.f
*/
void braktop(double aa, int *ia, double a[], int n, int l)
{
    braktop_(&aa, ia, a, &n, &l);
}

void brakbot(double aa, int *ia, double a[], int n, int l)
{
    brakbot_(&aa, ia, a, &n, &l);
}

void braktpa(double aa, int *ia, double a[], int n, int l)
{
    braktpa_(&aa, ia, a, &n, &l);
}

void brakbta(double aa, int *ia, double a[], int n, int l)
{
    brakbta_(&aa, ia, a, &n, &l);
}
