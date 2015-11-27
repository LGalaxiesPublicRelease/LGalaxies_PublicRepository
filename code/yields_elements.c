/*
 * elements.c
 *
 *  Created on: 20.01.2012
 *      Author: robyates
 *
 *  Where individual chemical element history arrays are created and dealt with (like metals.c for metal history arrays)
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"


#ifdef INDIVIDUAL_ELEMENTS

/* Defined in allvars.h
struct elements
{
  float H;
  float He;
#ifndef MAINELEMENTS
  float Cb; //NOTE: Carbon (C) is stored as Cb here
  float N;
#endif
  float O;
#ifndef MAINELEMENTS
  float Ne;
#endif
  float Mg;
#ifndef MAINELEMENTS
  float Si;
  float S;
  float Ca;
#endif
  float Fe;
};
*/

struct elements elements_add(struct elements ele1, struct elements ele2, float fraction)
{
  struct elements ele;
  ele.H=ele1.H+fraction*ele2.H;
  ele.He=ele1.He+fraction*ele2.He;
#ifndef MAINELEMENTS
  ele.Cb=ele1.Cb+fraction*ele2.Cb;
  ele.N=ele1.N+fraction*ele2.N;
#endif
  ele.O=ele1.O+fraction*ele2.O;
#ifndef MAINELEMENTS
  ele.Ne=ele1.Ne+fraction*ele2.Ne;
#endif
  ele.Mg=ele1.Mg+fraction*ele2.Mg;
#ifndef MAINELEMENTS
  ele.Si=ele1.Si+fraction*ele2.Si;
  ele.S=ele1.S+fraction*ele2.S;
  ele.Ca=ele1.Ca+fraction*ele2.Ca;
#endif
  ele.Fe=ele1.Fe+fraction*ele2.Fe;

  return(ele);
}

struct elements elements_init()
{
  struct elements ele;
  ele.H=0.0;
  ele.He=0.0;
#ifndef MAINELEMENTS
  ele.Cb=0.0;
  ele.N=0.0;
#endif
  ele.O=0.0;
#ifndef MAINELEMENTS
  ele.Ne=0.0;
#endif
  ele.Mg=0.0;
#ifndef MAINELEMENTS
  ele.Si=0.0;
  ele.S=0.0;
  ele.Ca=0.0;
#endif
  ele.Fe=0.0;
  return(ele);
}


void elements_print(char s[],struct elements ele)
{
  printf("%s.H [Msun]  = %.2f\n",s,ele.H);
  printf("%s.He [Msun] = %.2f\n",s,ele.He);
#ifndef MAINELEMENTS
  printf("%s.Cb [Msun] = %.2f\n",s,ele.Cb);
  printf("%s.N [Msun]  = %.2f\n",s,ele.N);
#endif
  printf("%s.O [Msun]  = %.2f\n",s,ele.O);
#ifndef MAINELEMENTS
  printf("%s.Ne [Msun] = %.2f\n",s,ele.Ne);
#endif
  printf("%s.Mg [Msun] = %.2f\n",s,ele.Mg);
#ifndef MAINELEMENTS
  printf("%s.Si [Msun] = %.2f\n",s,ele.Si);
  printf("%s.S [Msun]  = %.2f\n",s,ele.S);
  printf("%s.Ca [Msun] = %.2f\n",s,ele.Ca);
#endif
  printf("%s.Fe [Msun] = %.2f\n",s,ele.Fe);
  return;
}

double elements_total(struct elements ele)
{
#ifndef MAINELEMENTS
  return(ele.H+ele.He+ele.Cb+ele.N+ele.O+ele.Ne+ele.Mg+ele.Si+ele.S+ele.Ca+ele.Fe);
#else
  return(ele.H+ele.He+ele.O+ele.Mg+ele.Fe);
#endif
}

double metal_elements_total(struct elements ele)
{
#ifndef MAINELEMENTS
  return(ele.Cb+ele.N+ele.O+ele.Ne+ele.Mg+ele.Si+ele.S+ele.Ca+ele.Fe);
#else
  return(ele.O+ele.Mg+ele.Fe);
#endif
}

#endif //INDIVIDUAL_ELEMENTS
