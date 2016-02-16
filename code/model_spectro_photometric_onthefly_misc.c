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
 *  Created in: 2009
 *      Author: Chiara Tonini & Bruno Henriques
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <stddef.h>

#include "allvars.h"
#include "proto.h"

/* double get_AbsAB_magnitude(double *FluxInputSSPInt, double *FluxFilterInt, int BandNum, double redshift)
 * double get_area (double redshift)
 * double lum_distance(double redshift)
 * double error_on_mag(double AbsMag)
 *
 * void create_grid (double WaveMin, double WaveMax,double *LambdaInputSSP_SingleAge, double *lgrid,
 * int *Min_Wave_Grid, int *Max_Wave_Grid, int *Grid_Length)
 *
 * void locate(double *xx, int n, double x, int *j)
 * void interpolate(double *lgrid, int Grid_Length, double *lambda, int nlambda, double *flux, double *FluxOnGrid)
 * double integrate(double *flux, int Grid_Length)
 * void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
 */
double get_AbsAB_magnitude(double FluxInputSSPInt, double FluxFilterInt, double redshift)
{
  double zeropoint, distance_cm;
  double AbsAB,area;

//it needs to be converted to cm since the units of
//the original InputSSP spectra are (erg.s^-1.AA^-1) which
//was converted to (erg.s^-1.Hz^-1) by doing Flux*Lambda^2/Clight*1e8
//we need 3631Jy or 3631*10^-23 erg.s^-1.Hz^-1 cm-2

  area=get_area(redshift);

  distance_cm=10.0*3.08568025e18;
  //4*pi*10pc^2*3631jy*erg.s^-1.Hz^-1

#ifdef APP
  zeropoint=+48.6-2.5*area;
#else
  zeropoint=-2.5*log10(4.0*M_PI*distance_cm*distance_cm*3631.0*1.0e-23);
#endif

  AbsAB=-2.5*(log10(FluxInputSSPInt)
      -log10(FluxFilterInt))-zeropoint;

  return AbsAB;
}


//use for apparent magnitudes
double get_area (double redshift)
{
  double dist, area;

  dist=lum_distance(redshift);

  //when calculating apparent magnitudes, minimum dist set to 10pc
  if (redshift<0.00000001) dist=1e-5;
  if(dist != 0.0)
    area=log10(4.*M_PI)+2.*log10(dist*3.08568025e24);  //in cm (dl in Mpc)
  else area=0.0;

  return area;
}


//LUMINOSITY DISTANCE

double lum_distance(double redshift)
{
  int i, k, Npoints=1000;
  double x[1000];
  double sum[2], I[3], f[4];
  double h, integral, dl;

  for(i=0;i<2;i++)sum[i]=0.0;
  for(i=0;i<3;i++)I[i]=0.0;
  for(i=0;i<4;i++)f[i]=0.0;

  h=redshift/(Npoints-1);
  for(i=0;i<Npoints;i++) x[i]=h*(i-1);


  for (i=0;i<Npoints/2;i++)
    {
      k=2*i-1;
      f[2]=1./sqrt((1.+x[k])*(1.+x[k])*(1.+Omega*x[k])-x[k]*OmegaLambda*(2.+x[k]));
      sum[0]=sum[0]+f[2];
    }
  I[1]=sum[0]*4./3.;

  for (i=0;i<Npoints/2-1;i++)
    {
      k=2*i;
      f[3]=1./sqrt((1.+x[k])*(1.+x[k])*(1.+Omega*x[k])-x[k]*OmegaLambda*(2.+x[k]));
      sum[1]=sum[1]+f[3];
    }
  I[2]=sum[1]*2./3.;

  f[1]=1./sqrt((1.+x[0])*(1.+x[0])*(1.+Omega*x[0])-x[0]*OmegaLambda*(2.+x[0]));
  f[2]=1./sqrt((1.+x[Npoints-1])*(1.+x[Npoints-1])*(1.+Omega*x[Npoints-1])
	       -x[Npoints-1]*OmegaLambda*(2.+x[Npoints-1]));

  I[0]=(f[0]+f[1])/3.;

  integral=h*(I[0]+I[1]+I[2]);

  dl=integral/1000.0;    //!Mpc

  dl*=(1.+redshift)*(C/100.)/(Hubble_h*100.);   //in Mpc

  return dl;

  }


/* creates a grid of points based on the input SSP. The first point on the grid is the
 * first point for which  LambdaInputSSP>FilterWaveMin. The last point is the largest
 * wavelength for which LambdaInputSSP<FilterWaveMax*/

double* create_grid (double WaveMin, double WaveMax,int AgeLoop, double redshift, double LambdaInputSSP[SSP_NAGES][SSP_NLambda],
		                      int *Min_Wave_Grid, int *Max_Wave_Grid, int *Grid_Length)
{
  double x0, x1, h;
  int i, min, max;
  double *grid;

  min=min(WaveMin, WaveMax);
  max=max(WaveMin, WaveMax);

  *Grid_Length=0;

  //get minimum of grid
  for(i=0;i<SSP_NLambda;i++)
    if((1+redshift)*LambdaInputSSP[AgeLoop][i]>=min)
      {
	*Min_Wave_Grid=i;
	break;
      }

  for(i=0;i<SSP_NLambda;i++)
    if((1+redshift)*LambdaInputSSP[AgeLoop][i]>=min)
      {
	*Grid_Length+=1;

	//point at maximum range or out of it, set maximum
	if((1+redshift)*LambdaInputSSP[AgeLoop][i]>=max)
	  {
	    if((1+redshift)*LambdaInputSSP[AgeLoop][i]==max)
	      *Max_Wave_Grid=i;
	    else //if point out of range, set max to previous
	      {
		*Max_Wave_Grid=i-1;
		*Grid_Length-=1;
	      }
	    break;
	  }
      }

  grid = malloc(sizeof(double) * *Grid_Length);
  for(i=0;i<*Grid_Length;i++)
    grid[i]=(1+redshift)*LambdaInputSSP[AgeLoop][*Min_Wave_Grid+i];

  return grid;
}


//*****************************************************
//interpolate filters on integral grid
//*****************************************************

void interpolate(double *lgrid, int Grid_Length, double *lambda, int nlambda, double *flux, double *FluxOnGrid)
{
  int kk=0, nn=0, m=2, i;

  for(i=0;i<Grid_Length;i++)
    {
      if (lgrid[i] < lambda[0] || lgrid[i] > lambda[nlambda-1]) FluxOnGrid[i]=0;
      //outside filter range, transmission is 0
      else
	{
	  //finds where wavelenght is in respect to the grid
	  locate(lambda,nlambda-1,lgrid[i],&nn);
	  kk=min(max(nn-(m-1)/2,1),nlambda+1-m);
	  FluxOnGrid[i]=flux[kk];
	}
    }
}

