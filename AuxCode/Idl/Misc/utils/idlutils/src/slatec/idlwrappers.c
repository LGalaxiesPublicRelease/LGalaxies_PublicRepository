#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "export.h"

void bvalu_(float *t, float *a, IDL_LONG *n, IDL_LONG *k, 
         IDL_LONG *ideriv, float *x, IDL_LONG *inbv, float *w, float *retval);

void efc_(IDL_LONG *ndata, float *xdata, float *ydata, float *sddata, 
         IDL_LONG *nord, IDL_LONG *nbkpt, float *bkpt, IDL_LONG *mdein, 
         IDL_LONG *mdeout, float *coeff, IDL_LONG *lw, float *w);

void efcmn_(IDL_LONG *, float *, float *, float *, IDL_LONG *, IDL_LONG *, 
         float *, IDL_LONG *, IDL_LONG *, float *, float *, float *, float *,
         float *, float *, IDL_LONG *, float *, IDL_LONG *, IDL_LONG *);

void bndacc_(float *g, IDL_LONG *mdg, IDL_LONG *nb, 
       IDL_LONG *ip, IDL_LONG *ir, IDL_LONG *mt, IDL_LONG *jt); 

void bndsol_(IDL_LONG *mode, float *g, IDL_LONG *mdg, IDL_LONG *nb, 
       IDL_LONG *ip, IDL_LONG *ir, float *x, IDL_LONG *n, float *rnorm);

IDL_LONG bndacc_idl
 (int      argc,
   void *   argv[])
{
  float *g;
  IDL_LONG *mdg;
  IDL_LONG *nb;
  IDL_LONG *ip;
  IDL_LONG *ir;
  IDL_LONG *mt;
  IDL_LONG *jt;

  int argct = 0;
  g      = (float *)argv[argct++];
  mdg    = (IDL_LONG *)argv[argct++];
  nb     = (IDL_LONG *)argv[argct++];
  ip     = (IDL_LONG *)argv[argct++];
  ir     = (IDL_LONG *)argv[argct++];
  mt     = (IDL_LONG *)argv[argct++];
  jt     = (IDL_LONG *)argv[argct++];

  bndacc_(g, mdg, nb, ip, ir, mt, jt);

  return 0;
}

IDL_LONG bndsol_idl
 (int      argc,
   void *   argv[])
{
  IDL_LONG *mode;
  float *g;
  IDL_LONG *mdg;
  IDL_LONG *nb;
  IDL_LONG *ip;
  IDL_LONG *ir;
  float *x;
  IDL_LONG *n;
  float *rnorm;

  int argct = 0;
  mode   = (IDL_LONG *)argv[argct++];
  g      = (float *)argv[argct++];
  mdg    = (IDL_LONG *)argv[argct++];
  nb     = (IDL_LONG *)argv[argct++];
  ip     = (IDL_LONG *)argv[argct++];
  ir     = (IDL_LONG *)argv[argct++];
  x      = (float *)argv[argct++];
  n      = (IDL_LONG*)argv[argct++];
  rnorm  = (float *)argv[argct++];

  bndsol_(mode, g, mdg, nb, ip, ir, x, n, rnorm);

  return 0;
}


IDL_LONG efcmn_idl
 (int      argc,
   void *   argv[])
{
  IDL_LONG *ndata;
  float *xdata;
  float *ydata; 
  float *sddata;
  IDL_LONG *nord;
  IDL_LONG *nbkpt;
  float  *bkptin;
  IDL_LONG *mdein;
  IDL_LONG *mdeout;
  float *coeff;
  float *bf;
  float *xtemp;
  float *ptemp;
  float  *bkpt;
  float *g;
  IDL_LONG *mdg;
  float  *w;
  IDL_LONG *mdw;
  IDL_LONG *lw;

  int argct = 0;
  ndata  = (IDL_LONG *)argv[argct++];
  xdata  = (float *)argv[argct++];
  ydata  = (float *)argv[argct++];
  sddata = (float *)argv[argct++];
  nord   = (IDL_LONG*)argv[argct++];
  nbkpt  = (IDL_LONG*)argv[argct++];
  bkptin = (float *)argv[argct++];
  mdein  = (IDL_LONG*)argv[argct++];
  mdeout = (IDL_LONG*)argv[argct++];
  coeff  = (float *)argv[argct++];
  bf     = (float *)argv[argct++];
  xtemp  = (float *)argv[argct++];
  ptemp  = (float *)argv[argct++];
  bkpt   = (float *)argv[argct++];
  g      = (float *)argv[argct++];
  mdg    = (IDL_LONG*)argv[argct++];
  w      = (float *)argv[argct++];
  mdw    = (IDL_LONG*)argv[argct++];
  lw     = (IDL_LONG*)argv[argct++];

  efcmn_(ndata, xdata, ydata, sddata, nord, nbkpt, bkptin, mdein, mdeout, 
       coeff, bf, xtemp, ptemp, bkpt, g, mdg, w, mdw, lw);

  return *mdeout;

} 
  
IDL_LONG efc_idl
 (int      argc,
   void *   argv[])
{
  IDL_LONG *ndata;
  float *xdata;
  float *ydata; 
  float *sddata;
  IDL_LONG *nord;
  IDL_LONG *nbkpt;
  float  *bkpt;
  IDL_LONG *mdein;
  IDL_LONG *mdeout;
  float *coeff;
  IDL_LONG *lw;
  float  *w;

  int argct = 0;
  ndata  = (IDL_LONG *)argv[argct++];
  xdata  = (float *)argv[argct++];
  ydata  = (float *)argv[argct++];
  sddata = (float *)argv[argct++];
  nord   = (IDL_LONG*)argv[argct++];
  nbkpt  = (IDL_LONG*)argv[argct++];
  bkpt   = (float *)argv[argct++];
  mdein  = (IDL_LONG*)argv[argct++];
  mdeout = (IDL_LONG*)argv[argct++];
  coeff  = (float *)argv[argct++];
  lw     = (IDL_LONG*)argv[argct++];
  w      = (float *)argv[argct++];
  efc_(ndata, xdata, ydata, sddata, nord, nbkpt, bkpt, mdein, mdeout, 
       coeff, lw, w);

  return *mdeout;

} 
  
IDL_LONG bvalu_idl
 (int      argc,
   void *   argv[])
{
  float *t;
  float *a; 
  IDL_LONG *n;
  IDL_LONG *k;
  IDL_LONG *ideriv;
  float *x; 
  IDL_LONG *ndata;
  IDL_LONG *inbv;
  float  *w;
  float *value;
  int i;

  int argct = 0;
  t = (float *)argv[argct++];
  a = (float *)argv[argct++];
  n = (IDL_LONG *)argv[argct++];
  k = (IDL_LONG *)argv[argct++];
  ideriv = (IDL_LONG *)argv[argct++];
  x      = (float *)argv[argct++];
  ndata  = (IDL_LONG *)argv[argct++];
  inbv   = (IDL_LONG *)argv[argct++];
  w      = (float *)argv[argct++];
  value   = (float *)argv[argct++];

  for (i=0;i<*ndata;i++)    
    bvalu_(t, a, n, k, ideriv, &x[i], inbv, w, &value[i]);

  return 0;
} 
  



