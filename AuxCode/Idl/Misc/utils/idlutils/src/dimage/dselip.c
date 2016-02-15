#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define M 64
#define BIG 1.0e30

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
#define FREEALL FREEVEC(sel);FREEVEC(isel)

float dselip(unsigned long k, unsigned long n, float *arr)
{
	void dshell(unsigned long n, float *a);
	unsigned long i,j,jl,jm,ju,kk,mm,nlo,nxtmm,*isel;
	float ahi,alo,sum,*sel;

	if (k < 1 || k > n || n <= 0) {
    printf("bad input to selip");
    exit(1);
  }
	isel=(unsigned long *) malloc(sizeof(unsigned long)*(M+2));
	sel=(float *) malloc(sizeof(float)*(M+2));
	kk=k+1;
	ahi=BIG;
	alo = -BIG;
	for (;;) {
		mm=nlo=0;
		sum=0.0;
		nxtmm=M+1;
		for (i=1;i<=n;i++) {
			if (arr[i-1] >= alo && arr[i-1] <= ahi) {
				mm++;
				if (arr[i-1] == alo) nlo++;
				if (mm <= M) sel[mm-1]=arr[i-1];
				else if (mm == nxtmm) {
					nxtmm=mm+mm/M;
					sel[1 + ((i+mm+kk) % M) - 1]=arr[i-1];
				}
				sum += arr[i-1];
			}
		}
		if (kk <= nlo) {
			FREEALL
			return alo;
		}
		else if (mm <= M) {
			dshell(mm,sel);
			ahi = sel[kk-1];
			FREEALL
			return ahi;
		}
		sel[M+1-1]=sum/mm;
		dshell(M+1,sel);
		sel[M+2-1]=ahi;
		for (j=1;j<=M+2;j++) isel[j-1]=0;
		for (i=1;i<=n;i++) {
			if (arr[i-1] >= alo && arr[i-1] <= ahi) {
				jl=0;
				ju=M+2;
				while (ju-jl > 1) {
					jm=(ju+jl)/2;
					if (arr[i-1] >= sel[jm-1]) jl=jm;
					else ju=jm;
				}
				isel[ju-1]++;
			}
		}
		j=1;
		while (kk > isel[j-1]) {
			alo=sel[j-1];
			kk -= isel[j-1];
      j++;
		}
		ahi=sel[j-1];
	}
}
#undef M
#undef BIG
#undef FREEALL
void dshell(unsigned long n, float *a)
{
	unsigned long i,j,inc;
	float v;
	inc=1;
	do {
		inc *= 3;
		inc++;
	} while (inc <= n);
	do {
		inc /= 3;
		for (i=inc+1;i<=n;i++) {
			v=a[i-1];
			j=i;
			while (a[j-inc-1] > v) {
				a[j-1]=a[j-inc-1];
				j -= inc;
				if (j <= inc) break;
			}
			a[j-1]=v;
		}
	} while (inc > 1);
}
