/*
 * metals.c
 *
 * Routines to handle the creation and addition of containers for metals.
 *
 * The routines are needed even if the metallicity is a single float,
 * because the lines for creation and addition of floats in the code
 * have been replaced by function calls (this to avoid a multitude of
 * #ifdef DETAILED_METALS_AND_MASS_RETURN statements throughout the code).
 *
 * In the description below metal is either:
 *   float
 *   struct metals
 * It would be relatively straight forward to add a metallicity array option also.
 *
 * metal metals_add(metal m1, metal m2, float fraction);
 *  Function returns m1+fraction*m2 for each metal component.
 *
 * metal metals_init();
 *  Function returns 0.0 for each metal component.
 *
 * void metals_print(char [s] ,metal m);
 *  Prints out the value of each metal component, preceded by string s.
 *
 * float metals_total(metal m);
 *  Function returns the total of all components of metal.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/* Code for metals array.
 * For simplicity, did away with non-array version. */

void metals_add(double *m1,
		double *m2,
		double fraction)
{
  int ii;
  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
      m1[ii]=m1[ii]+fraction*m2[ii];
  return;
}

void metals_init(double *m)
{
  int ii;
  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
      m[ii]=0.;
  return;
}

void metals_print(char s[],double *m)
{
  printf("%s.type1a [Msun] = %.2f\n",s,m[0]*1.0e10/Hubble_h);
  printf("%s.type2 [Msun]  = %.2f\n",s,m[1]*1.0e10/Hubble_h);
  printf("%s.agb  [Msun]   = %.2f\n",s,m[2]*1.0e10/Hubble_h);
  return;
}

double metals_total(double *m)
{
  int ii;
  double mTotal=0.;
  for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
      mTotal+=m[ii];
  return(mTotal);
}
