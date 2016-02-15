#define NRANSI
#include "nrutil.h"
#define NR_M 64
#define BIG 1.0e30
#define FREEALL nr_free_vector(sel,1,NR_M+2);nr_free_lvector(isel,1,NR_M+2);

float selip(unsigned long k, unsigned long n, float arr[])
{
	void shell(unsigned long n, float a[]);
	unsigned long i,j,jl,jm,ju,kk,mm,nlo,nxtmm,*isel;
	float ahi,alo,sum,*sel;

	if (k < 1 || k > n || n <= 0) nrerror("bad input to selip");
	isel=nr_lvector(1,NR_M+2);
	sel=nr_vector(1,NR_M+2);
	kk=k;
	ahi=BIG;
	alo = -BIG;
	for (;;) {
		mm=nlo=0;
		sum=0.0;
		nxtmm=NR_M+1;
		for (i=1;i<=n;i++) {
			if (arr[i] >= alo && arr[i] <= ahi) {
				mm++;
				if (arr[i] == alo) nlo++;
				if (mm <= NR_M) sel[mm]=arr[i];
				else if (mm == nxtmm) {
					nxtmm=mm+mm/NR_M;
					sel[1 + ((i+mm+kk) % NR_M)]=arr[i];
				}
				sum += arr[i];
			}
		}
		if (kk <= nlo) {
			FREEALL
			return alo;
		}
		else if (mm <= NR_M) {
			shell(mm,sel);
			ahi = sel[kk];
			FREEALL
			return ahi;
		}
		sel[NR_M+1]=sum/mm;
		shell(NR_M+1,sel);
		sel[NR_M+2]=ahi;
		for (j=1;j<=NR_M+2;j++) isel[j]=0;
		for (i=1;i<=n;i++) {
			if (arr[i] >= alo && arr[i] <= ahi) {
				jl=0;
				ju=NR_M+2;
				while (ju-jl > 1) {
					jm=(ju+jl)/2;
					if (arr[i] >= sel[jm]) jl=jm;
					else ju=jm;
				}
				isel[ju]++;
			}
		}
		j=1;
		while (kk > isel[j]) {
			alo=sel[j];
			kk -= isel[j++];
		}
		ahi=sel[j];
	}
}
#undef NR_M
#undef BIG
#undef FREEALL
#undef NRANSI
void shell(unsigned long n, float a[])
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
			v=a[i];
			j=i;
			while (a[j-inc] > v) {
				a[j]=a[j-inc];
				j -= inc;
				if (j <= inc) break;
			}
			a[j]=v;
		}
	} while (inc > 1);
}
