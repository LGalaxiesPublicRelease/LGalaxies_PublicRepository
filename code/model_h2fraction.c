#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

void update_h2fraction(int p)
{
  int j, ii;
  //the central stellar surface density converted from (10^10M_sun/h)/(Mpc/h)^2 to (M_sun/pc^2)
  double SigmaHRings, inverse_ColdGas;

  Gal[p].H2fraction=0.;
  inverse_ColdGas = 1./Gal[p].ColdGas;

  for(j=0;j<RNUM;j++)
    {
      //KMT09
      if(H2FractionRecipe==0)
  	{
	  double metallicityr=0.;;

	  if(Gal[p].ColdGasRings[j]>1.0e-8)
	    {
	      for(ii=0;ii<NUM_METAL_CHANNELS;ii++)
		metallicityr+= Gal[p].MetalsColdGasRings[j][ii];
	      metallicityr /= (Gal[p].ColdGasRings[j]*0.0134);
	    }

	  if(metallicityr<0.01)
	    metallicityr=0.01;

	  //SigmaHRings = Gal[p].ColdGasRings[j] / (RingArea[j]*WARM_PHASE_FACTOR);
	  /* now convert from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2) -> *(1.e10/h)/(1.e6/h*1.e6/h)*/
	  //SigmaHRings = SigmaHRings * 0.01 * Hubble_h;
	  SigmaHRings = Gal[p].ColdGasRings[j] * InverseRingArea[j];

	  //update to clumping factor, as in Fu2013, to solve problems with low-z galaxies
	  if(metallicityr<1.0)
	    SigmaHRings=SigmaHRings*Clumpingfactor*pow(metallicityr,-0.7);

	  Gal[p].H2fractionRings[j] = update_H2fraction_KMT09(log10(SigmaHRings), metallicityr);
  	}
      //Blitz & Rosolowsky 2006, pressure recipe
      else if(H2FractionRecipe==1)
	{
	  double SigmaStarRings, alpha_p=0.92;
	  double SigmaStar0 = (Gal[p].DiskMassRings[0]/RingArea[0])*0.01*Hubble_h;

	  SigmaHRings = Gal[p].ColdGasRings[j] * InverseRingArea[j];
	  //from 10^10 M_sun/h / (Mpc/h)^2 to (M_sun/pc^2) */
	  SigmaStarRings= (Gal[p].DiskMassRings[j] / RingArea[j]) * 0.01*Hubble_h;

	  Gal[p].H2fractionRings[j]=1.38e-3*pow(SigmaHRings*(SigmaHRings+0.1*sqrt(SigmaStar0*SigmaStarRings)),alpha_p);
	  //Gal[p].H2fractionRings[j]=6.81e-3*pow(SigmaHRings*(SigmaHRings+0.1*sqrt(SigmaStar0*SigmaStarRings)),0.80);
	  if(Gal[p].H2fractionRings[j]<1.0e-8)
	    Gal[p].H2fractionRings[j]=0.0;
	  else
	    Gal[p].H2fractionRings[j]=1/(1+1/Gal[p].H2fractionRings[j]);
	}
      else Gal[p].H2fractionRings[j]=0;

      Gal[p].H2fraction += Gal[p].H2fractionRings[j] * Gal[p].ColdGasRings[j] * inverse_ColdGas;

    }

  if(Gal[p].ColdGas<1.0e-7) Gal[p].H2fraction=0.0;




}




void init_H2fraction_KMT09()
{
  int ii, jj;
  double khi, tau, s;

  H2Fraction_Zgrid[0]=0.01;
  for (ii=1; ii<LENZ; ii++)
    H2Fraction_Zgrid[ii] = H2Fraction_Zgrid[ii-1]*3.;

  for (jj=0; jj<LENSIGMAH; jj++)
    H2Fraction_SigmaHgrid[jj] = jj*0.1;

  for (ii=0; ii<LENZ; ii++)
    for (jj=0; jj<LENSIGMAH; jj++)
      {
	khi = 3.1/4.1 * ( 1.+3.1 * pow(H2Fraction_Zgrid[ii],0.365) );
	//khi = 3.1/4.1 * ( 1.+0.001 * pow(H2Fraction_Zgrid[ii],0.365) );
	//tau = 0.066 * pow(10.,H2Fraction_SigmaHgrid[jj]) * H2Fraction_Zgrid[ii];
	tau = 0.066 * pow(10.,H2Fraction_SigmaHgrid[jj]) * H2Fraction_Zgrid[ii];
	s = log(1.+0.6*khi+0.01*khi*khi)/(0.6*tau);
	H2Fraction[ii][jj] = max(4-2.*s,0.)/(4.+s);
      }
}


#define NJUMPTAB_H2 (LENSIGMAH-1)
int jumptab_H2[NJUMPTAB_H2];
double jumpfac_H2;

void init_jump_index_H2Fraction(void)
{
  double sigmaH;
  int i, idx;

  jumpfac_H2 = NJUMPTAB_H2 / (H2Fraction_SigmaHgrid[LENSIGMAH - 1] - H2Fraction_SigmaHgrid[1]);

  for(i = 0; i < NJUMPTAB_H2; i++)
    {
      sigmaH = H2Fraction_SigmaHgrid[1] + i / jumpfac_H2;
      idx = 1;
      while(H2Fraction_SigmaHgrid[idx + 1] < sigmaH)
	idx++;
      jumptab_H2[i] = idx;
    }
}


int get_jump_index_H2Fraction(double sigmaH)
{
	int idx;

	idx = (int) ((sigmaH - H2Fraction_SigmaHgrid[1]) * jumpfac_H2);
	if(idx<0)
		idx=0;
	return jumptab_H2[idx];
	  //return jumptab_H2[(int) ((sigmaH - H2Fraction_SigmaHgrid[1]) * jumpfac_H2)];
}


double update_H2fraction_KMT09(double logsigmah, double metallicity )
{
  int ii, jj;
  double mf1, mf2, mf;

  metallicity = max(metallicity, H2Fraction_Zgrid[0]);
  metallicity = min(metallicity, H2Fraction_Zgrid[LENZ-1]);
  //find interpolation point for metallicity
  for (ii=0; metallicity > H2Fraction_Zgrid[ii+1]; ii++ );

  if(logsigmah<=H2Fraction_SigmaHgrid[0])
    {
	  logsigmah=H2Fraction_SigmaHgrid[0];
	  jj=0;
    }
  else if(logsigmah>=H2Fraction_SigmaHgrid[LENSIGMAH-1])
    {
	  logsigmah= H2Fraction_SigmaHgrid[LENSIGMAH-1];
	  jj = LENSIGMAH-1;
    }
  else
    {
	  //logsigmah = max(logsigmah, H2Fraction_SigmaHgrid[0]);
	  //logsigmah = min(logsigmah, H2Fraction_SigmaHgrid[LENSIGMAH-1]);

	  //find interpolation point for sigmaH
	  jj = get_jump_index_H2Fraction(logsigmah);
	  while(H2Fraction_SigmaHgrid[jj+1] < logsigmah)
		  jj++;
    }

  mf1=H2Fraction[ii][jj]+ ( H2Fraction[ii+1][jj]-H2Fraction[ii][jj] ) * ( metallicity-H2Fraction_Zgrid[ii] ) / ( H2Fraction_Zgrid[ii+1]-H2Fraction_Zgrid[ii] );
  mf2=H2Fraction[ii][jj+1]+ ( H2Fraction[ii+1][jj+1]-H2Fraction[ii][jj+1] ) * ( metallicity-H2Fraction_Zgrid[ii] ) / ( H2Fraction_Zgrid[ii+1]-H2Fraction_Zgrid[ii] );
  mf=mf1+ ( mf2-mf1 ) * ( logsigmah-H2Fraction_SigmaHgrid[jj] ) / ( H2Fraction_SigmaHgrid[jj+1]-H2Fraction_SigmaHgrid[jj] );

  //mf+=0.1;

  return (mf);
}



