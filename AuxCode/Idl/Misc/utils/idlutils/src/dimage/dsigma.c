#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dimage.h"
#include "export.h"

/*
 * dsigma.c
 *
 * Simple guess at the sky sigma
 *
 * Mike Blanton
 * 1/2006 */

#define PI 3.14159265358979

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

static float *diff=NULL;

float dselip(unsigned long k, unsigned long n, float *arr);

int dsigma(float *image, 
					 int nx, 
					 int ny,
					 int sp,
					 float *sigma)
{
	float tot;
  int i,j,dx,dy, ndiff;

	if(nx==1 && ny==1) {
		(*sigma)=0.;
		return(0);
	}

	dx=50;
	if(dx>nx/4) dx=nx/4;
	if(dx<=0) dx=1;

	dy=50;
	if(dy>ny/4) dy=ny/4;
	if(dy<=0) dy=1;
	
	diff=(float *) malloc(2*nx*ny*sizeof(float));
	ndiff=0;
	for(j=0;j<ny;j+=dy) {
		for(i=0;i<nx;i+=dx) {
			if(i<nx-sp) {
				diff[ndiff]=fabs(image[i+j*nx]-image[i+sp+j*nx]);
				ndiff++;
			}
			if(j<ny-sp) {
				diff[ndiff]=fabs(image[i+j*nx]-image[i+(j+sp)*nx]);
				ndiff++;
			}
		}
	}

	if(ndiff<=1) {
		(*sigma)=0.;
		return(0);
	}

	if(ndiff<=10) {
		tot=0.;
		for(i=0;i<ndiff;i++)
			tot+=diff[i]*diff[i];
		(*sigma)=sqrt(tot/(float) ndiff);
		return(0);
	}

	(*sigma)=(dselip((int) floor(ndiff*0.68),ndiff,diff))/sqrt(2.);
	
	FREEVEC(diff);

	return(1);
} /* end dsigma */
