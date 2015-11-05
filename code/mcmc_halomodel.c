
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>

#include "allvars.h"
#include "mcmc_vars.h"
#include "mcmc_halomodel.h"
#include "mcmc_mpfit.h"

#define WORKSIZE 100000
#define PI 3.14159
#ifdef MCRIT
#define xrvir 6
#else
#define xrvir 3.16
#endif

void halomodel(double* r_arr,double* proj_arr,float masslimit_low, float masslimit_high,int snap) {
  int i;
  int NK=60;
  int extralin=16;
  double k,p1h,p2h,m10,m11;
  double kbin=6./79.;
  double r,corrtmp;
  double *rCorrTable,*CorrTable;
  gsl_set_error_handler_off();
  cutoff_low=malloc(6*sizeof(double));
  cutoff_high=malloc(6*sizeof(double));
  init_numgal(masslimit_low,masslimit_high,snap);
  ngal_mean=ngal_mean_calc(0);
  PowerTable=malloc((NK+extralin)*sizeof(double));
  kPowerTable=malloc((NK+extralin)*sizeof(double));
  m10=Mcensat(1,1,0);
  for (i=0; i<NK+extralin; i++) {
    if (i>2*extralin) k=pow(10.,(i-extralin)*kbin-2.);
    else k=pow(10.,i*0.5*kbin-2.);
    if (i>0 && kPowerTable[i-1]>2) k=pow(10.,kPowerTable[i-1]+2.044/6.);
    p1h=2*Mcensat(k,0,1)+Mcensat(k,0,2);
    m11=Mcensat(k,1,1);
    p2h=PowerSpec(k)*(m10*m10+2*m10*m11+m11*m11);
    kPowerTable[i]=log10(k);
    PowerTable[i]=log10((p1h+p2h)/(gsl_spline_eval(ellipSpline,kPowerTable[i],ellipAcc)+1));
  } //for
  gsl_spline_free(pcSpline);
  gsl_interp_accel_free(pcAcc);
  gsl_spline_free(pbSpline);
  gsl_interp_accel_free(pbAcc);
  gsl_spline_free(paSpline);
  gsl_interp_accel_free(paAcc);
  for (i=5; i>=0; --i) {
    gsl_spline_free(NgalSpline[i]);
    gsl_interp_accel_free(NgalAcc[i]);
  } //for
  free(NgalSpline);
  free(NgalAcc);
  NewpowAcc=gsl_interp_accel_alloc();
  NewpowSpline=gsl_spline_alloc(gsl_interp_cspline,(NK+extralin));
  gsl_spline_init(NewpowSpline,kPowerTable,PowerTable,(NK+extralin));
#if defined(OUTPUTCORR) || defined(OUTPUTPOW)
  FILE *fd;
  char buf[500];
  float mingalmass,maxgalmass;
  mingalmass=8.77+(ThisTask%6)*0.5;
  maxgalmass=8.77+(ThisTask%6+1)*0.5;
#endif
#ifdef OUTPUTPOW
  sprintf(buf,"pow_%.2f-%.2f_%d.dat",mingalmass,maxgalmass,snap);
  fd=fopen(buf,"w");
  for (i=0; i<1000; ++i) {
    fprintf(fd,"%g   %g\n",pow(10.,0.001*i*5-2),pow(10.,gsl_spline_eval(NewpowSpline,0.001*i*5-2,NewpowAcc)));
  } //for
  fclose(fd);
#endif
  CorrTable=malloc(NR*10*sizeof(double));
  rCorrTable=malloc(NR*10*sizeof(double));
  for (i=0; i<NR*10; ++i) {
    rCorrTable[i]=pow(10.,(i+0.5)*(log10(610.)+3.1)/(float)(NR*10)-3.1);
    CorrTable[i]=corr_qawo(rCorrTable[i],0.01,0.99)+corr_qawo(rCorrTable[i],1.,99.)+corr_qawo(rCorrTable[i],100.,900.)+corr_qawo(rCorrTable[i],1e3,5e3)+corr_qawo(rCorrTable[i],6e3,5e3);
  } //for
  CorrAcc=gsl_interp_accel_alloc();
  CorrSpline=gsl_spline_alloc(gsl_interp_cspline,NR*10);
  gsl_spline_init(CorrSpline,rCorrTable,CorrTable,NR*10);
#ifdef OUTPUTCORR
  sprintf(buf,"corr_%.2f-%.2f_%d.dat",mingalmass,maxgalmass,snap);
  fd=fopen(buf,"w");
  for (i=0; i<1000; ++i) {
    fprintf(fd,"%g   %g\n",pow(10.,0.001*i*5.785-3),gsl_spline_eval(CorrSpline,pow(10.,0.001*i*5.785-3),CorrAcc));
  } //for
  fclose(fd);
#endif
  for (i=0; i<NR; ++i) {
    r=pow(10.,(i+0.5)*(log10(80.)+2.1)/(float)(NR-1)-2.3);
    corrtmp=proj_corr(r);
    r_arr[i]=r/Hubble_h;
    proj_arr[i]=corrtmp/Hubble_h;
  } //for
#ifdef OUTPUTPROJ
  sprintf(buf,"proj_%.2f-%.2f_%d.dat",mingalmass,maxgalmass,snap);
  fd=fopen(buf,"w");
  for (i=0; i<NR; ++i) {
    fprintf(fd,"%g   %g\n",r_arr[i],proj_arr[i]);
  } //for
  fclose(fd);
#endif
  gsl_spline_free(CorrSpline);
  gsl_interp_accel_free(CorrAcc);
  gsl_spline_free(NewpowSpline);
  gsl_interp_accel_free(NewpowAcc);
  free(kPowerTable);
  free(PowerTable);
  free(cutoff_high);
  free(cutoff_low);
} //halomodel

double NewPowerSpec(double k) {
  return pow(10.,gsl_spline_eval(NewpowSpline,log10(k),NewpowAcc));
} //NewPowerSpec

double corr_qawo_func(double k,void *params) {
  return k*NewPowerSpec(k)*exp(-k*k/1e7); //exponential needed for convergence
} //corr_qawo_func

double corr_qawo(double r,double a,double L) {
  double result=0,abserr;
  gsl_function F;
  int status;
  F.function=&corr_qawo_func;
  F.params=0;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  gsl_integration_qawo_table *t=gsl_integration_qawo_table_alloc(r,L,GSL_INTEG_SINE,40);
  status=gsl_integration_qawo(&F,a,0,1.0e-3,WORKSIZE,w,t,&result,&abserr);
  gsl_integration_qawo_table_free(t);
  gsl_integration_workspace_free(w);
  return result/(2.*PI*PI*r);
} //corr_qawo

double proj_corr_func(double r,void *params) {
  double sigma=*(double *) params;
  return 2*r*gsl_spline_eval(CorrSpline,r,CorrAcc)/sqrt(pow(r,2)-pow(sigma,2));
} //proj_corr_func

double proj_corr(double sigma) {
  double result=0,abserr;
  gsl_function F;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  F.function=&proj_corr_func;
  F.params=&sigma;
  gsl_integration_qag(&F,sigma,600.,0,1.0e-3,WORKSIZE,GSL_INTEG_GAUSS51,w,&result,&abserr);
  gsl_integration_workspace_free(w);
  return result;
} //proj_corr

double Radius(double m) {
  return pow(m/rho_mean*3./4./PI,1./3.);
} //Radius

double Sigma2(double m) {
  return pow(10.,gsl_spline_eval(SigmaSpline,log10(m),SigmaAcc));
} //Sigma2

double PowerSpec(double k) {
  double lk,lp;
  lk=log10(k);
  lp=gsl_spline_eval(PowSpline,lk,PowAcc);
  return Norm*pow(10.,lp-3*lk)*2.*PI*PI;
} //PowerSpec

double TopHatWindow(double kr) {
  double kr2,kr3;
  kr2=kr*kr;
  kr3=kr2*kr;
  if (kr<1e-8) return 0;
  return 3.*(sin(kr)/kr3-cos(kr)/kr2);
} //TopHatWindow

double sigma2_func(double k,void *params) {
  double r_tophat=*(double *) params;
  double W=TopHatWindow(r_tophat*k);
  return 1./(2.*PI*PI)*k*k*W*W*PowerSpec(k);
} //sigma2_func

double TopHatSigma2(double R) {
  double result=0,abserr;
  gsl_function F;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  F.function=&sigma2_func;
  F.params=&R;
  gsl_integration_qag(&F,2*PI/500.,500.,0,1.0e-5,WORKSIZE,GSL_INTEG_GAUSS41,w,&result,&abserr);
  gsl_integration_workspace_free(w);
  return result;
} //TopHatSigma2

double Mass(double R) {
  return 4.0*PI*R*R*R/3.0*Omega*rho_c;
} //Mass

double lnSigma2(double lR,void *params) {
  double M=Mass(exp(lR));
  return log(Sigma2(M));
} //lnSigma2

double lnSigma2lnM(double lm,void *params) {
  return log(Sigma2(exp(lm)));
} //lnSigma2lnM

double gammaM(double m) {
  gsl_function F;
  double result,abserr,lm;
  lm=log(m);
  F.function=&lnSigma2lnM;
  F.params=0;
  gsl_deriv_central(&F,lm,1e-8,&result,&abserr);
  return -result;
} //gammaM

double nbargal(double m) {
  double lm=log10(m);
  double value;
  if (lm<cutoff_fof_low || lm>cutoff_fof_high) return 0.;
  else {
    value=gsl_spline_eval(FofSpline,lm,FofAcc);
    if (value>0) return value;
    else return 0.;
  } //else
} //nbargal

double b(double m,int i) {
  double nu2,delta_c2,b1;
  delta_c2=delta_c*delta_c;
  nu2=delta_c2/Sigma2(m);
  //Tinker et al. (2010)
  double AT,aT,BT,bT,CT,cT;
  double y=log10(Delta);
  AT=1.0+0.24*y*exp(-pow(4./y,4));
  aT=0.44*y-0.88;
  BT=0.183;
  bT=1.5;
  CT=0.019+0.107*y+0.19*exp(-pow(4./y,4));
  cT=2.4;
  b1=1-AT*pow(nu2,0.5*aT)/(pow(nu2,0.5*aT)+pow(delta_c2,0.5*aT))+BT*pow(nu2,0.5*bT)+CT*pow(nu2,0.5*cT);
  switch(i) {
  case 0:
    return 1.0;
    break;
  case 1:
    return b1;
    break;
  } //switch
  return 1.0;
} //b

double mugal_qawo_func(double r,void *p) {
  struct mugal_qawo_params *params=(struct mugal_qawo_params *)p;
  double pa=(params->pa);
  double pb=(params->pb);
  double pc=(params->pc);
  double nm=pow(r/pb,pa*pc)/pow(r,3)*exp(-pow(r/pb,pc));
  return r*nm;
} //mugal_qawo_func

double mugal_qawo(double k,double m) {
  double rvir=Radius(m)*Delta_invth;
  double rvir3=rvir*rvir*rvir;
  double pa=pow(10.,pa_eval(log10(m)));
  double pb=pow(10.,pb_eval(log10(m)));
  double pc=pow(10.,pc_eval(log10(m)));
  double norm=pc/(rvir3*4*PI*exp(gsl_sf_lngamma(pa)+log(gsl_sf_gamma_inc_P(pa,pow(xrvir/pb,pc)))));
  struct mugal_qawo_params params={ pa,pb,pc };
  double result=0,abserr;
  gsl_function F;
  int status;
  F.function=&mugal_qawo_func;
  F.params=&params;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  gsl_integration_qawo_table *t=gsl_integration_qawo_table_alloc(k*rvir,xrvir,GSL_INTEG_SINE,25);
  status=gsl_integration_qawo(&F,0.0,0,1.0e-3,WORKSIZE,w,t,&result,&abserr);
  gsl_integration_qawo_table_free(t);
  gsl_integration_workspace_free(w);
  return norm*4*PI*rvir*rvir*result/k;
} //mugal_qawo

double NgalF(double m,int j) { //0: ngal, 1: ngal*(ngal-1), 2: ncen, 3: nsat, 4:ncen*nsat, 5: nsat*(nsat-1)
  double lm=log10(m);
  double value;
  if (lm<cutoff_low[j] || lm>cutoff_high[j]) return 0.;
  else {
    value=gsl_spline_eval(NgalSpline[j],lm,NgalAcc[j]);
    if (value>0) return value;
    else return 0.;
  } //else
} //NgalF

double pa_eval(double m) {
  double value;
  if (m<parscutoff_low) return gsl_spline_eval(paSpline,parscutoff_low,paAcc);
  else if (m>parscutoff_high) return gsl_spline_eval(paSpline,parscutoff_high,paAcc);
  else {
    value=gsl_spline_eval(paSpline,m,paAcc);
    if (value>=-2. && value<=2.) return value;
    else if (value<-2.) return -2.;
    else return 2.;
  } //else
} //pa_eval

double pb_eval(double m) {
  double value;
  if (m<parscutoff_low) return gsl_spline_eval(pbSpline,parscutoff_low,pbAcc);
  else if (m>parscutoff_high) return gsl_spline_eval(pbSpline,parscutoff_high,pbAcc);
  else {
    value=gsl_spline_eval(pbSpline,m,pbAcc);
    if (value>=-1. && value<=0.5) return value;
    else if (value<-1.) return -1.;
    else return 0.5;
  } //else
} //pb_eval

double pc_eval(double m) {
  double value;
  if (m<parscutoff_low) return gsl_spline_eval(pcSpline,parscutoff_low,pcAcc);
  else if (m>parscutoff_high) return gsl_spline_eval(pcSpline,parscutoff_high,pcAcc);
  else {
    value=gsl_spline_eval(pcSpline,m,pcAcc);
    if (value>=-2. && value<=2.) return value;
    else if (value<-2.) return -2.;
    else return 2.;
  } //else
} //pc_eval

double Mcensat_func(double lm,void *p) {
  struct M_params *params=(struct M_params *) p;
  double k=(params->k);
  int i=(params->i);
  int j=(params->j);
  double m=exp(lm);
  if (i==0) { //counterterm due to Valageas & Nishimichi (2011)
    if (j==0) return 0.;
#ifdef MCRIT
    else if (j==1) return nbargal(m)*NgalF(m,4)/pow(ngal_mean,2)*mugal_qawo(k,m)*m;
    else return nbargal(m)*NgalF(m,5)/pow(ngal_mean,2)*pow(mugal_qawo(k,m),2)*m;
#else
    else if (j==1) return nbargal(m)*NgalF(m,4)/pow(ngal_mean,2)*(mugal_qawo(k,m)-TopHatWindow(k*Radius(m)))*m;
    else return nbargal(m)*NgalF(m,5)/pow(ngal_mean,2)*(pow(mugal_qawo(k,m),2)-pow(TopHatWindow(k*Radius(m)),2))*m;
#endif
  } //if
  else {
    if (j==0) return nbargal(m)*b(m,1)*NgalF(m,2)/ngal_mean*m;
    else return nbargal(m)*b(m,1)*NgalF(m,3)/ngal_mean*mugal_qawo(k,m)*m;
  } //else
} //Mcensat_func

double Mcensat(double k,int i,int j) { //k=wavenumber, i=which haloterm [0/1], j=central/satellite [0/1/2]
  double result=0,abserr;
  struct M_params params={ k,i,j };
  gsl_function F;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  F.function=&Mcensat_func;
  F.params=&params;
  gsl_integration_qag(&F,log(pow(10.,cutoff_low[3-i+j])),log(pow(10.,cutoff_high[3-i+j])),0,1.0e-3,WORKSIZE,GSL_INTEG_GAUSS41,w,&result,&abserr);
  gsl_integration_workspace_free(w);
  return result;
} //Mcensat

double ngal_mean_func(double lm,void *p) {
  struct N_params *params=(struct N_params *) p;
  int j=(params->j);
  double m=exp(lm);
  return nbargal(m)*NgalF(m,j)*m; //extra m for logarithmic integration
} //ngal_mean_func

double ngal_mean_calc(int j) {
  double result=0,abserr;
  struct N_params params={ j };
  gsl_function F;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  F.function=&ngal_mean_func;
  F.params=&params;
  gsl_integration_qag(&F,log(pow(10.,cutoff_low[j])),log(pow(10.,cutoff_high[j])),0,1.0e-3,WORKSIZE,GSL_INTEG_GAUSS41,w,&result,&abserr);
  gsl_integration_workspace_free(w);
  return result;
} //ngal_mean_calc

void init_power() {
  FILE *fd;
  char buf[500];
  double k,p;
  int NPowerTable=0;
  sprintf(buf,"%s/powrealized_rebin_corrected.dat",MCMCHaloModelDir);
  if (!(fd=fopen(buf,"r"))) {
    char sbuf[1000];
    sprintf(sbuf,"Can't read input spectrum in file '%s'.\n",buf);
    terminate(sbuf);
  } //if
  do {
    if (fscanf(fd," %lg %lg ",&k,&p)==2) NPowerTable++;
    else break;
  } //do
  while(1);
  fclose(fd);
  PowerTable=malloc(NPowerTable*sizeof(double));
  kPowerTable=malloc(NPowerTable*sizeof(double));
  fd=fopen(buf,"r");
  NPowerTable=0;
  do { //k and Delta
    if (fscanf(fd," %lg %lg ",&k,&p)==2) {
      kPowerTable[NPowerTable]=k-log10(ScalePos);
      PowerTable[NPowerTable]=p;
      NPowerTable++;
    } //if
    else break;
  } //do
  while(1);
  fclose(fd);
  PowAcc=gsl_interp_accel_alloc();
  PowSpline=gsl_spline_alloc(gsl_interp_cspline,NPowerTable);
  gsl_spline_init(PowSpline,kPowerTable,PowerTable,NPowerTable);
  free(kPowerTable);
  free(PowerTable);
  NPowerTable=0;
  sprintf(buf,"%s/ellip_corr.dat",MCMCHaloModelDir);
  if (!(fd=fopen(buf,"r"))) {
    char sbuf[1000];
    sprintf(sbuf,"Can't read correction spectrum in file '%s'.\n",buf);
    terminate(sbuf);
  } //if
  do {
    if (fscanf(fd," %lg %lg ",&k,&p)==2) NPowerTable++;
    else break;
  } //do
  while(1);
  fclose(fd);
  PowerTable=malloc(NPowerTable*sizeof(double));
  kPowerTable=malloc(NPowerTable*sizeof(double));
  fd=fopen(buf,"r");
  NPowerTable=0;
  do { //k and Delta
    if (fscanf(fd," %lg %lg ",&k,&p)==2) {
      kPowerTable[NPowerTable]=k-log10(ScalePos);
      PowerTable[NPowerTable]=p;
      NPowerTable++;
    } //if
    else break;
  } //do
  while(1);
  fclose(fd);
  ellipAcc=gsl_interp_accel_alloc();
  ellipSpline=gsl_spline_alloc(gsl_interp_cspline,NPowerTable);
  gsl_spline_init(ellipSpline,kPowerTable,PowerTable,NPowerTable);
  free(kPowerTable);
  free(PowerTable);
} //init_power

void init_sigma() {
  int i,NSigmaTable=500;
  double MSigmaTable[500],SigmaTable[500];
  double R;
  for (i=0; i<NSigmaTable; i++) {
    MSigmaTable[i]=i*15./(float)NSigmaTable+5.;
    R=Radius(pow(10.,MSigmaTable[i]));
    SigmaTable[i]=log10(TopHatSigma2(R));
  } //for
  SigmaAcc=gsl_interp_accel_alloc();
  SigmaSpline=gsl_spline_alloc(gsl_interp_cspline,NSigmaTable);
  gsl_spline_init(SigmaSpline,MSigmaTable,SigmaTable,NSigmaTable);
} //init_sigma

void init_numgal(float masslimit_low, float masslimit_high, int snap) {
  int i,j,jj;
  double *MassTable;
  int mbin,mbin2,k,found,ncen,nsat;
  int nsat_tot=0;
  double *masstmp,*patmp,*pbtmp,*pctmp;
  double **NgalTable;
  int NMassTable2=5; //5
  double borders[6]; //NMassTable2+1
  double p[3];
  double perror[3];
  double *pberr;
  double *paerr;
  double *pcerr;
  double *x;
  double *y;
  struct fitvars v;
  int status;
  mp_par pars[3];
  mp_result result;
  mp_config config;
  int NRadiusTable=25+(int)(log10(xrvir)/0.1+0.5);
  double minlogr=-2.5;
  int rbin,NgalTotalTotal,mstart,mend,NValid;
  int *NgalTotal,*NfofTotal,**NProfileTable,*usedbymass;
  double rbinsize=0.1;
  double boxsize=BoxSize;
  double massoffset=log10(2*(G/1e10)/(Delta*Omega*1e4));
  double logrvir,logr,relx,rely,relz;
  double *MassTable2,*MassVarTable2,**RadiusTable,**ProfileTable,**ProfileErrTable,*pa_m,*pb_m,*pc_m,*weights;
  for (i=0; i<6; ++i) {
    cutoff_low[i]=0;
    cutoff_high[i]=0;
  } //for
  NgalAcc=malloc(6*sizeof(gsl_interp_accel*));
  NgalSpline=malloc(6*sizeof(gsl_spline*));
  usedbymass=malloc(massbins*sizeof(int));
  for (i=0; i<massbins; ++i) usedbymass[i]=0;
  for (i=0; i<6; ++i) {
    NgalAcc[i]=gsl_interp_accel_alloc();
    NgalSpline[i]=gsl_spline_alloc(gsl_interp_cspline,massbins);
  } //for
  MassTable=malloc(massbins*sizeof(double));
  NgalTable=malloc(6*sizeof(double*));
  for (i=0; i<6; ++i) NgalTable[i]=malloc(massbins*sizeof(double));
  for (j=0; j<massbins; ++j) {
    MassTable[j]=0;
    for (i=0; i<6; ++i) NgalTable[i][j]=0;
  } //for
  for (mbin=0; mbin<massbins; ++mbin) {
    for (j=0; j<UsedFofsInSample[snap]; ++j) {
      ncen=0;
      nsat=0;
      //if (MCMC_FOF2[j].M_Crit200[snap]>minfofmass+mbin*(maxfofmass-minfofmass)/(double)massbins && MCMC_FOF2[j].M_Crit200[snap]<minfofmass+(mbin+1)*(maxfofmass-minfofmass)/(double)massbins) {
      if (MCMC_FOF2[j].M_Mean200[snap]>minfofmass+mbin*(maxfofmass-minfofmass)/(double)massbins &&
      		MCMC_FOF2[j].M_Mean200[snap]<minfofmass+(mbin+1)*(maxfofmass-minfofmass)/(double)massbins)
      {
        i=MCMC_FOF2[j].IndexOfCentralGal[snap];
        if (MCMC_GAL[HashTable[i]].StellarMass[snap]>masslimit_low && MCMC_GAL[HashTable[i]].StellarMass[snap]<=masslimit_high) ncen=1;
        for (k=1; k<MCMC_GAL[HashTable[i]].ngal[snap]; ++k) {
          if (MCMC_GAL[HashTable[i+k]].StellarMass[snap]>masslimit_low && MCMC_GAL[HashTable[i+k]].StellarMass[snap]<=masslimit_high) nsat++;
        } //for
        usedbymass[mbin]++;
      } //if
      NgalTable[0][mbin]+=ncen+nsat;
      NgalTable[1][mbin]+=(ncen+nsat)*max((ncen+nsat)-1,0);
      NgalTable[2][mbin]+=ncen;
      NgalTable[3][mbin]+=nsat;
      NgalTable[4][mbin]+=ncen*nsat;
      NgalTable[5][mbin]+=nsat*max(nsat-1,0);
      nsat_tot+=nsat;
    } //for
  } //for
  for (j=0; j<massbins; ++j) MassTable[j]=minfofmass+(j+0.5)*(maxfofmass-minfofmass)/(double)massbins;
  nsat=NgalTable[3][0];
  //borders[0]=10.;
  if (NgalTable[0][0]>0) {
    char sbuf[1000];
    sprintf(sbuf,"First Ngal mass bin not equal to zero, divergence will ensue.\n");
    terminate(sbuf);
  } //if

  borders[0]=minfofmass;
  borders[NMassTable2]=maxfofmass;

  j=1;
  for (mbin=1; mbin<massbins; ++mbin) {
    if (nsat<=j*nsat_tot/(double)NMassTable2 && nsat+NgalTable[3][mbin]>j*nsat_tot/(double)NMassTable2) {
      borders[j]=MassTable[mbin]+(MassTable[mbin]-MassTable[mbin-1])/NgalTable[3][mbin]*(j*nsat_tot/(double)NMassTable2-nsat);
      j++;
      if (j==NMassTable2) break;
    } //if
    nsat+=NgalTable[3][mbin];
  } //for
  for (j=0; j<massbins; ++j) {
    if (usedbymass[j]>0) {
      for (i=0; i<6; ++i) NgalTable[i][j]/=usedbymass[j];
    } //if
  } //for
  for (i=0; i<6; ++i) {
    for (mbin=0; mbin<massbins-1 && cutoff_low[i]==0; ++mbin) {
      if (NgalTable[i][mbin]==0 && NgalTable[i][mbin+1]>0) {
	    cutoff_low[i]=MassTable[mbin];
      } //if
    } //for
    for (mbin=massbins-1; mbin>0 && cutoff_high[i]==0; --mbin) {
      if (NgalTable[i][mbin-1]>0 && NgalTable[i][mbin]==0) {
	    cutoff_high[i]=MassTable[mbin];
      } //if
    } //for
  } //for
  for (i=0; i<6; ++i) gsl_spline_init(NgalSpline[i],MassTable,NgalTable[i],massbins);
  for (i=5; i>=0; --i) free(NgalTable[i]);
  free(NgalTable);
  free(MassTable);
  MassTable2=malloc(NMassTable2*sizeof(double));
  MassVarTable2=malloc(NMassTable2*sizeof(double));
  RadiusTable=malloc(NMassTable2*sizeof(double*));
  NProfileTable=malloc(NMassTable2*sizeof(int*));
  ProfileTable=malloc(NMassTable2*sizeof(double*));
  ProfileErrTable=malloc(NMassTable2*sizeof(double*));
  for (i=0; i<NMassTable2; ++i) {
    RadiusTable[i]=malloc(NRadiusTable*sizeof(double));
    NProfileTable[i]=malloc(NRadiusTable*sizeof(int));
    ProfileTable[i]=malloc(NRadiusTable*sizeof(double));
    ProfileErrTable[i]=malloc(NRadiusTable*sizeof(double));
    for (j=0; j<NRadiusTable; ++j) {
      RadiusTable[i][j]=0;
      NProfileTable[i][j]=0;
      ProfileTable[i][j]=0;
      ProfileErrTable[i][j]=0;
    } //for
  } //for
  NgalTotal=malloc(NMassTable2*sizeof(int));
  NfofTotal=malloc(NMassTable2*sizeof(int));
  for (i=0; i<NMassTable2; ++i) {
    NgalTotal[i]=0;
    NfofTotal[i]=0;
    MassTable2[i]=0;
    MassVarTable2[i]=0;
  } //for
  for (j=0; j<UsedFofsInSample[snap]; ++j) {
    i=MCMC_FOF2[j].IndexOfCentralGal[snap];
    for (jj=1; jj<=NMassTable2; ++jj) {
      //if (MCMC_GAL[HashTable[i]].M_Crit200[snap]<borders[jj]) {
    	if (MCMC_GAL[HashTable[i]].M_Mean200[snap]<borders[jj]) {
        mbin2=jj-1;
        break;
      } //if
    } //for
    if (mbin2>=0 && mbin2<NMassTable2 && MCMC_GAL[HashTable[i]].ngal[snap]>1) {
      found=0;
      //logrvir=1./3.*(MCMC_GAL[HashTable[i]].M_Crit200[snap]+massoffset);
      logrvir=1./3.*(MCMC_GAL[HashTable[i]].M_Mean200[snap]+massoffset);
      for (jj=1; jj<MCMC_GAL[HashTable[i]].ngal[snap]; ++jj) {
        if (MCMC_GAL[HashTable[i+jj]].StellarMass[snap]>masslimit_low && MCMC_GAL[HashTable[i+jj]].StellarMass[snap]<=masslimit_high) {
          found++;
          relx=min(fabs(MCMC_GAL[HashTable[i+jj]].x[snap]-MCMC_GAL[HashTable[i]].x[snap]),fabs(boxsize-fabs(MCMC_GAL[HashTable[i+jj]].x[snap]-MCMC_GAL[HashTable[i]].x[snap])));
          rely=min(fabs(MCMC_GAL[HashTable[i+jj]].y[snap]-MCMC_GAL[HashTable[i]].y[snap]),fabs(boxsize-fabs(MCMC_GAL[HashTable[i+jj]].y[snap]-MCMC_GAL[HashTable[i]].y[snap])));
          relz=min(fabs(MCMC_GAL[HashTable[i+jj]].z[snap]-MCMC_GAL[HashTable[i]].z[snap]),fabs(boxsize-fabs(MCMC_GAL[HashTable[i+jj]].z[snap]-MCMC_GAL[HashTable[i]].z[snap])));
          logr=0.5*log10(relx*relx+rely*rely+relz*relz);
          rbin=floor((logr-logrvir-minlogr)/rbinsize);
          if (rbin>=0 && rbin<NRadiusTable) {
            RadiusTable[mbin2][rbin]+=logr-logrvir;
            NProfileTable[mbin2][rbin]++;
          } //if
        } //if
      } //for
      if (found>0) {
        NgalTotal[mbin2]+=found;
        NfofTotal[mbin2]++;
        //MassTable2[mbin2]+=found*MCMC_GAL[HashTable[i]].M_Crit200[snap];
        //MassVarTable2[mbin2]+=found*SQR(MCMC_GAL[HashTable[i]].M_Crit200[snap]);
        MassTable2[mbin2]+=found*MCMC_GAL[HashTable[i]].M_Mean200[snap];
        MassVarTable2[mbin2]+=found*SQR(MCMC_GAL[HashTable[i]].M_Mean200[snap]);
      } //if
    } //if
  } //for
  NgalTotalTotal=0;
  for (mbin2=0; mbin2<NMassTable2; ++mbin2) {
    if (NgalTotal[mbin2]>0) MassTable2[mbin2]/=NgalTotal[mbin2];
    else MassTable2[mbin2]=0.5*(borders[mbin2]+borders[mbin2+1]);
    MassVarTable2[mbin2]-=NgalTotal[mbin2]*SQR(MassTable2[mbin2]);
    MassVarTable2[mbin2]/=(NgalTotal[mbin2]-1);
    NgalTotalTotal+=NgalTotal[mbin2];
    for (rbin=0; rbin<NRadiusTable; ++rbin) {
      if (NProfileTable[mbin2][rbin]>0) {
	    RadiusTable[mbin2][rbin]/=(double)NProfileTable[mbin2][rbin];
	    ProfileTable[mbin2][rbin]=(double)NProfileTable[mbin2][rbin]/(NgalTotal[mbin2]*4./3.*PI*(pow(10.,3*((rbin+1)*rbinsize+minlogr))-pow(10.,3*(rbin*rbinsize+minlogr)))*pow(10.,MassTable2[mbin2]+massoffset));
	    ProfileErrTable[mbin2][rbin]=sqrt(NProfileTable[mbin2][rbin])/(NgalTotal[mbin2]*4./3.*PI*(pow(10.,3*((rbin+1)*rbinsize+minlogr))-pow(10.,3*(rbin*rbinsize+minlogr)))*pow(10.,MassTable2[mbin2]+massoffset));
      } //if
      else RadiusTable[mbin2][rbin]=(rbin+0.5)*rbinsize+minlogr;
    } //for
  } //for
  pa_m=malloc(NMassTable2*sizeof(double));
  pb_m=malloc(NMassTable2*sizeof(double));
  pc_m=malloc(NMassTable2*sizeof(double));
  paerr=malloc(NMassTable2*sizeof(double));
  pberr=malloc(NMassTable2*sizeof(double));
  pcerr=malloc(NMassTable2*sizeof(double));
  memset(&result,0,sizeof(result));
  memset(&config,0,sizeof(config));
  memset(&pars[0],0,sizeof(pars));
  pars[0].limited[0]=1;
  pars[0].limits[0]=-2.;
  pars[0].limited[1]=1;
  pars[0].limits[1]=2.;
  pars[1].limited[0]=1;
  pars[1].limits[0]=-1.;
  pars[1].limited[1]=1;
  pars[1].limits[1]=0.5;
  pars[2].limited[0]=1;
  pars[2].limits[0]=-2.;
  pars[2].limited[1]=1;
  pars[2].limits[1]=2.;
  for (mbin2=0; mbin2<NMassTable2; ++mbin2) {
    NValid=0;
    for (rbin=0; rbin<NRadiusTable; ++rbin) if (NProfileTable[mbin2][rbin]>0) NValid++;
    if (NValid>2) {
      x=malloc(NValid*sizeof(double));
      y=malloc(NValid*sizeof(double));
      weights=malloc(NValid*sizeof(double));
      i=0;
      for (rbin=0; rbin<NRadiusTable; ++rbin) {
        if (NProfileTable[mbin2][rbin]>0) {
	  x[i]=pow(10.,RadiusTable[mbin2][rbin]);
	  y[i]=log10(ProfileTable[mbin2][rbin]);
	  weights[i]=max(1./sqrt(MassVarTable2[mbin2]+SQR(0.5*(log10(ProfileTable[mbin2][rbin]+ProfileErrTable[mbin2][rbin])-log10(ProfileTable[mbin2][rbin]-ProfileErrTable[mbin2][rbin])))),1./sqrt(MassVarTable2[mbin2]+SQR(3*(log10(ProfileTable[mbin2][rbin]+ProfileErrTable[mbin2][rbin])-log10(ProfileTable[mbin2][rbin])))));
	  if (isnan(weights[i])) {
	    printf("WARNING: Task %d found a NaN for mbin2=%d and i=%d, setting weight to 0...\n",ThisTask,mbin2,i);
	    weights[i]=0.;
	  } //if
	  i++;
        } //if
      } //for
      result.xerror=perror;
      v.x=x;
      v.y=y;
      v.w=weights;
      fitconst=MassTable2[mbin2]+massoffset;
      p[0]=0.25;
      p[1]=-0.64;
      p[2]=-0.05;
      config.ftol=1e-12;
      config.xtol=1e-12;
      config.maxiter=2000;
      status=mpfit(gammafit,NValid,3,p,pars,&config,(void *) &v,&result);
      pa_m[mbin2]=p[0];
      pb_m[mbin2]=p[1];
      pc_m[mbin2]=p[2];
      paerr[mbin2]=perror[0];
      pberr[mbin2]=perror[1];
      pcerr[mbin2]=perror[2];
      if (paerr[mbin2]==0) paerr[mbin2]=0.5;
      if (pberr[mbin2]==0) pberr[mbin2]=0.5;
      if (pcerr[mbin2]==0) pcerr[mbin2]=0.5;
      if (paerr[mbin2]<0.1) paerr[mbin2]=0.1;
      if (pberr[mbin2]<0.1) pberr[mbin2]=0.1;
      if (pcerr[mbin2]<0.1) pcerr[mbin2]=0.1;
      free(weights);
      free(y);
      free(x);
    } //if
  } //for
  mstart=0;
  mend=NMassTable2;
  NValid=0;
  for (mbin2=0; NValid<3 && mbin2<NMassTable2; ++mbin2) {
    NValid=0;
    for (rbin=0; rbin<NRadiusTable; ++rbin) if (NProfileTable[mbin2][rbin]>0) NValid++;
    if (NValid<3) mstart++;
  } //for
  NValid=0;
  for (mbin2=NMassTable2-1; NValid<3 && mbin2>=0; --mbin2) {
    NValid=0;
    for (rbin=0; rbin<NRadiusTable; ++rbin) if (NProfileTable[mbin2][rbin]>0) NValid++;
    if (NValid<3) mend--;
  } //for
  masstmp=malloc((NMassTable2+2)*sizeof(double));
  patmp=malloc((NMassTable2+2)*sizeof(double));
  pbtmp=malloc((NMassTable2+2)*sizeof(double));
  pctmp=malloc((NMassTable2+2)*sizeof(double));
  masstmp[0]=minfofmass;
  masstmp[NMassTable2+1]=maxfofmass;
  patmp[0]=pa_m[0];
  pbtmp[0]=pb_m[0];
  pctmp[0]=pc_m[0];
  patmp[NMassTable2+1]=pa_m[NMassTable2-1];
  pbtmp[NMassTable2+1]=pb_m[NMassTable2-1];
  pctmp[NMassTable2+1]=pc_m[NMassTable2-1];
  for (i=1; i<=NMassTable2; ++i) {
    masstmp[i]=MassTable2[i-1];
    patmp[i]=pa_m[i-1];
    pbtmp[i]=pb_m[i-1];
    pctmp[i]=pc_m[i-1];
  } //for
  paAcc=gsl_interp_accel_alloc();
  paSpline=gsl_spline_alloc(gsl_interp_akima,NMassTable2+2);
  pbAcc=gsl_interp_accel_alloc();
  pbSpline=gsl_spline_alloc(gsl_interp_akima,NMassTable2+2);
  pcAcc=gsl_interp_accel_alloc();
  pcSpline=gsl_spline_alloc(gsl_interp_akima,NMassTable2+2);
  gsl_spline_init(paSpline,masstmp,patmp,NMassTable2+2);
  gsl_spline_init(pbSpline,masstmp,pbtmp,NMassTable2+2);
  gsl_spline_init(pcSpline,masstmp,pctmp,NMassTable2+2);
  parscutoff_low=10.;
  parscutoff_high=16.;
  free(pctmp);
  free(pbtmp);
  free(patmp);
  free(masstmp);
  free(pcerr);
  free(pberr);
  free(paerr);
  free(pc_m);
  free(pb_m);
  free(pa_m);
  free(NfofTotal);
  free(NgalTotal);
  for (i=NMassTable2-1; i>=0; --i) {
    free(ProfileErrTable[i]);
    free(ProfileTable[i]);
    free(NProfileTable[i]);
    free(RadiusTable[i]);
  } //for
  free(ProfileErrTable);
  free(ProfileTable);
  free(NProfileTable);
  free(RadiusTable);
  free(MassTable2);
} //init_numgal

void initialize_halomodel() {
  int i,NnuTable=100;
  double res;
  double MnuTable[100],nuTable[100];
  FILE *fd;
  char buf[500];
  double m,h;
  double FofmassTable[massbins],FofnumTable[massbins];
  rho_c=RhoCrit*1e10;
  delta_c=(3./5.*pow(0.5*3*PI,2./3.))*(1-0.0123*log10(1+(1./Omega-1)));
  Delta_invth=1./pow(Delta,1./3.);
  rho_mean=Omega*rho_c;
  Norm=1.0;
  init_power();
  res=TopHatSigma2(8.);
  Norm=Sigma8*Sigma8/res;
  res=TopHatSigma2(8.);
  init_sigma();
  for (i=0; i<NnuTable; i++) {
    MnuTable[i]=log(pow(10,i*3./(float)NnuTable+11.));
    nuTable[i]=delta_c/sqrt(Sigma2(exp(MnuTable[i])));
  } //for
  nuAcc=gsl_interp_accel_alloc();
  nuSpline=gsl_spline_alloc(gsl_interp_cspline,NnuTable);
  gsl_spline_init(nuSpline,nuTable,MnuTable,NnuTable);
  Mstar=exp(gsl_spline_eval(nuSpline,1.0,nuAcc));
  gsl_spline_free(nuSpline);
  gsl_interp_accel_free(nuAcc);
  Rstar=Radius(Mstar);
  //sprintf(buf,"%s/fofnum.dat",MCMCHaloModelDir);
#ifdef MCRIT
  sprintf(buf,"%s/fofnum_m200_z0.02.dat",MCMCHaloModelDir);
#else
  sprintf(buf,"%s/fofnum_m200mean.dat",MCMCHaloModelDir);
#endif
  if (!(fd=fopen(buf,"r"))) {
    char sbuf[1000];
    sprintf(sbuf,"Can't read input FoF mass function in file '%s'.\n",buf);
    terminate(sbuf);
  } //if
  i=0;
  do {
    if (fscanf(fd," %lg %lg ",&m,&h)==2) i++;
    else break;
    FofmassTable[i-1]=m+log10(ScaleMass);
    FofnumTable[i-1]=h-3*log10(ScalePos);
  } //do
  while(1);
  fclose(fd);
  if (FofnumTable[0]>0) {
    char sbuf[1000];
    sprintf(sbuf,"First FoF mass bin not equal to zero, divergence will ensue.\n");
    terminate(sbuf);
  } //if
  cutoff_fof_low=0.;
  cutoff_fof_high=0.;
  for (i=0; i<massbins && cutoff_fof_low==0; ++i) {
    if (FofnumTable[i]==0 && FofnumTable[i+1]>0) {
      cutoff_fof_low=FofmassTable[i];
    } //if
  } //for
  for (i=massbins-1; i>0 && cutoff_fof_high==0; --i) {
    if (FofnumTable[i-1]>0 && FofnumTable[i]==0) {
      cutoff_fof_high=FofmassTable[i];
    } //if
  } //for
  for (i=0; i<massbins; ++i) FofnumTable[i]/=pow(10.,FofmassTable[i])*log(10.);
  FofAcc=gsl_interp_accel_alloc();
  FofSpline=gsl_spline_alloc(gsl_interp_cspline,massbins);
  gsl_spline_init(FofSpline,FofmassTable,FofnumTable,massbins);
} //initialize_halomodel

int gammafit(int m,int n,double *p,double *dy,double **dvec,void *vars) {
  int i;
  double lognorm=p[2]-log10(4*PI)-gsl_sf_lngamma(pow(10.,p[0]))/log(10.)-log10(gsl_sf_gamma_inc_P(pow(10.,p[0]),pow(xrvir/pow(10.,p[1]),pow(10.,p[2]))))-fitconst;
  struct fitvars *v=(struct fitvars *) vars;
  double *x,*y,*w,f;
  x=v->x; //should be r/rvir (NOT log)
  y=v->y; //should be log10(profile)
  w=v->w; //should be a standard deviation in log space
  for (i=0; i<m; i++) {
    f=lognorm+pow(10.,p[0]+p[2])*(log10(x[i])-p[1])-3*log10(x[i])-pow(x[i]/pow(10.,p[1]),pow(10.,p[2]))/log(10.);
    dy[i]=(y[i]-f)*w[i];
  } //for
  return 0;
} //gammafit
