#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"

void idl_nnls
  (int         argc,
   void    *   argv[])
{
  extern void nnls_();/* Fortran routine */
  int *mda,*m,*n,*indx,*mde;
  float *a,*b,*x,*resnorm,*w,*zz;
  a = (float *) argv[0];
  mda = * (int *) argv[1];
  m = * (int *) argv[2];
  n = * (int *) argv[3];
  b = (float *) argv[4];
  x = (float *) argv[5];
  resnorm = * (float *) argv[6];
  w = (float *) argv[7];
  zz = (float *) argv[8] ;
  indx = (int *) argv[9] ;
  mde = (int *) argv[10] ;
  for(i=0;i<mda[0];i++){
    for(j=0;j<n[0];j++){
      printf("%d %d %lf\n",i,j,a[i][j]) ;
    }
  }
  nnls_(a,mda,m,n,b,x,resnorm,w,zz,indx,mde);
}
