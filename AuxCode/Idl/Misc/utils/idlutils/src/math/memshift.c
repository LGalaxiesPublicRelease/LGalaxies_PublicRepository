#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "nr.h"

/******************************************************************************/
IDL_LONG memshift
  (int      argc,
   void *   argv[])
{
   unsigned char *  array;
   IDL_LONG    isrc;
   IDL_LONG    idest;
   IDL_LONG    nbyte;
   IDL_LONG    retval = 1;

   /* Allocate pointers from IDL */
   array = (unsigned char *)argv[0];
   isrc = *((IDL_LONG *)argv[1]);
   idest = *((IDL_LONG *)argv[2]);
   nbyte = *((IDL_LONG *)argv[3]);

   memmove(&array[idest], &array[isrc], nbyte);

   return retval;
}

/******************************************************************************/
