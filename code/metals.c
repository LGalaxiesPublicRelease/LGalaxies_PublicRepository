/*
 * metals.c
 *
 * Routines to handle the creation and addition of containers for metals.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

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

