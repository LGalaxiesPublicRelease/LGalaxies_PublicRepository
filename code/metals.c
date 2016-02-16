/*  Copyright (C) <2016>  <L-Galaxies>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/> */

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

/* Defined in allvars.h
 * struct metals
 * {
 *  float type1a;
 *  float type2;
 *  float agb;
 * };
 */

struct metals metals_add(struct metals m1,
			 struct metals m2,
			 float fraction)
{
  struct metals m;
  m.type1a=m1.type1a+fraction*m2.type1a;
  m.type2=m1.type2+fraction*m2.type2;
  m.agb=m1.agb+fraction*m2.agb;

  return(m);
}

struct metals metals_init()
{
  struct metals m;
  m.type1a=0.;
  m.type2=0.;
  m.agb=0.;
  return(m);
}

void metals_print(char s[],struct metals m)
{
  printf("%s.type1a [Msun] = %.2f\n",s,m.type1a*1.0e10/Hubble_h);
  printf("%s.type2 [Msun]  = %.2f\n",s,m.type2*1.0e10/Hubble_h);
  printf("%s.agb  [Msun]   = %.2f\n",s,m.agb*1.0e10/Hubble_h);
  return;
}

float metals_total(struct metals m)
{
  return(m.type1a+m.type2+m.agb);
}

#else

// The following mimics the original code with a single metallicity

float metals_add(float m1,
		 float m2,
		 float fraction)
{
  return(m1+fraction*m2);
}

float metals_init()
{
  return(0.);
}

void metals_print(char s[],float m)
{
  printf("%s=%f\n",s,m);
  return;
}

float metals_total(float m)
{
  return(m);
}

#endif


