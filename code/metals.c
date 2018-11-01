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

#ifdef DETAILED_METALS_AND_MASS_RETURN

/* Defined in h_metals.h
 * struct metals
 * {
 *  float type2;
 *  float type1a;
 *  float agb;
 * };
 */

union metals_arr metals_add(union metals_arr m1, union metals_arr m2, double fraction)
{
  union metals_arr m;
  int ii;
  for(ii==0;ii<3;ii++)
    m.arr[ii]=m1.arr[ii]+fraction*m2.arr[ii];
  return(m1);
}

union metals_arr metals_init()
{
  union metals_arr m;
  int ii;
  for(ii=0;ii<3;ii++)
     m.arr[ii]=0.;
  return(m);
}

void metals_print(char s[],union metals_arr m)
{
  printf("%s.type1a [Msun] = %.2f\n",s,m.str.type1a*1.0e10/Hubble_h);
  printf("%s.type2 [Msun]  = %.2f\n",s,m.str.type2*1.0e10/Hubble_h);
  printf("%s.agb  [Msun]   = %.2f\n",s,m.str.agb*1.0e10/Hubble_h);
  return;
}

double metals_total(union metals_arr m)
{
  return(m.str.type1a+m.str.type2+m.str.agb);
}

#else

// The following mimics the original code with a single metallicity

double metals_add(double m1,
		 double m2,
		 double fraction)
{
  return(m1+fraction*m2);
}

double metals_init()
{
  return(0.);
}

void metals_print(char s[],double m)
{
  printf("%s=%f\n",s,m);
  return;
}

double metals_total(double m)
{
  return(m);
}

#endif
