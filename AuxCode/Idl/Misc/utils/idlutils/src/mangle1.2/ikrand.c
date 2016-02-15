/* called as fortran integer function
      ikrand(ik)

   Pseudo-random number associated with integer ik.
*/
#include <stdlib.h>

int ikrand_(int *ik)
{
    int i;

    srandom((unsigned int) *ik);
    i = random() + 1;

    return(i);
}
