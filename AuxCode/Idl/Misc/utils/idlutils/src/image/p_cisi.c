#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "ph.h"
#define EPS 6.0e-8
#define EULER 0.57721566
#define MAXIT 100
#define PIBY2 1.5707963
#define FPMIN 1.0e-30
#define TMIN 2.0
#define ONE Complex(1.0,0.0)

/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file complex.c.  Do not confuse this file with the same-named
   file complex.c that is supplied in the same subdirectory or archive
   as the header file complex.h.  *That* file contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */

#include <math.h>

typedef struct FCOMPLEX {float r,i;} fcomplex;

fcomplex Cadd(fcomplex a, fcomplex b)
{
	fcomplex c; 
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}

fcomplex Csub(fcomplex a, fcomplex b)
{
	fcomplex c;
	c.r=a.r-b.r;
	c.i=a.i-b.i;
	return c;
}


fcomplex Cmul(fcomplex a, fcomplex b)
{
	fcomplex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}

fcomplex Complex(float re, float im)
{
	fcomplex c;
	c.r=re;
	c.i=im;
	return c;
}

fcomplex Conjg(fcomplex z)
{
	fcomplex c;
	c.r=z.r;
	c.i = -z.i;
	return c;
}

fcomplex Cdiv(fcomplex a, fcomplex b)
{
	fcomplex c;
	float r,den;
	if (fabs(b.r) >= fabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}

float Cabs(fcomplex z)
{
	float x,y,ans,temp;
	x=fabs(z.r);
	y=fabs(z.i);
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}

fcomplex Csqrt(fcomplex z)
{
	fcomplex c;
	float x,y,w,r;
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r=0.0;
		c.i=0.0;
		return c;
	} else {
		x=fabs(z.r);
		y=fabs(z.i);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r >= 0.0) {
			c.r=w;
			c.i=z.i/(2.0*w);
		} else {
			c.i=(z.i >= 0) ? w : -w;
			c.r=z.i/(2.0*c.i);
		}
		return c;
	}
}

fcomplex RCmul(float x, fcomplex a)
{
	fcomplex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}
/* (C) Copr. 1986-92 Numerical Recipes Software "!15L1. */

void p_cisi(float x, float *ci, float *si)
{
	int i,k,odd;
	float a,err,fact,sign,sum,sumc,sums,t,term;
	fcomplex h,b,c,d,del;

	t=fabs(x);
	if (t == 0.0) {
		*si=0.0;
		*ci = -1.0/FPMIN;
		return;
	}
	if (t > TMIN) {
		b=Complex(1.0,t);
		c=Complex(1.0/FPMIN,0.0);
		d=h=Cdiv(ONE,b);
		for (i=2;i<=MAXIT;i++) {
			a = -(i-1)*(i-1);
			b=Cadd(b,Complex(2.0,0.0));
			d=Cdiv(ONE,Cadd(RCmul(a,d),b));
			c=Cadd(b,Cdiv(Complex(a,0.0),c));
			del=Cmul(c,d);
			h=Cmul(h,del);
			if (fabs(del.r-1.0)+fabs(del.i) < EPS) break;
		}
		if (i > MAXIT) {
      fprintf(stderr,"cf failed in cisi\n");
      exit(1);
    }
		h=Cmul(Complex(cos(t),-sin(t)),h);
		*ci = -h.r;
		*si=PIBY2+h.i;
	} else {
		if (t < sqrt(FPMIN)) {
			sumc=0.0;
			sums=t;
		} else {
			sum=sums=sumc=0.0;
			sign=fact=1.0;
			odd=TRUE;
			for (k=1;k<=MAXIT;k++) {
				fact *= t/k;
				term=fact/k;
				sum += sign*term;
				err=term/fabs(sum);
				if (odd) {
					sign = -sign;
					sums=sum;
					sum=sumc;
				} else {
					sumc=sum;
					sum=sums;
				}
				if (err < EPS) break;
				odd=!odd;
			}
			if (k > MAXIT) {
        fprintf(stderr,"maxits exceeded in cisi %e\n",x);
        exit(1);
      }
		}
		*si=sums;
		*ci=sumc+log(t)+EULER;
	}
	if (x < 0.0) *si = -(*si);
}
#undef EPS
#undef EULER
#undef MAXIT
#undef PIBY2
#undef FPMIN
#undef TMIN
#undef TRUE
#undef ONE
/* (C) Copr. 1986-92 Numerical Recipes Software "!15L1. */
