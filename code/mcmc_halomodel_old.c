
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

void halomodel(double* r_arr,double* proj_arr,float masslimit_low, float masslimit_high, int snap)
{
  FILE *fd;
  int i,j;
  int NK=60;
  int extralin=16;
  double k,p1h,p2h,m10,m11;
  double lk;
  double kbin=6./79.;
  double hubble=0.73;
  double r,corrtmp;
  double *rCorrTable,*CorrTable;

  gsl_set_error_handler_off();
  cutoff_low=malloc(6*sizeof(double));
  cutoff_high=malloc(6*sizeof(double));

  init_numgal(masslimit_low, masslimit_high, snap);

  ngal_mean=ngal_mean_calc(0);
  ncen_mean=ngal_mean_calc(2);
  PowerTable=malloc((NK+extralin)*sizeof(double));
  kPowerTable=malloc((NK+extralin)*sizeof(double));

  //printf("\n\npower\n");

	/* char buf[1000];
	 FILE *fa;
	 //sprintf(buf, "/galformod/scratch/bmh20/plots/data/correlation_3.txt");
	 sprintf(buf, "/galformod/scratch/bmh20/plots/data/power_100_%0.2f.txt",masslimit_low);
	 if(!(fa = fopen(buf, "w")))
	 {
		 char sbuf[1000];
		 sprintf(sbuf, "can't open file `%s'\n", buf);
		 terminate(sbuf);
	 }*/



  for (i=0; i<NK+extralin; i++)
  {
    if (i>2*extralin)
    	k=pow(10.,(i-extralin)*kbin-2.);
    else
    	k=pow(10.,i*0.5*kbin-2.);

    if (i>0 && kPowerTable[i-1]>2)
    	k=pow(10.,kPowerTable[i-1]+2.044/6.);

    p1h=2*Mcensat(k,0,1)+Mcensat(k,0,2);
    m10=Mcensat(k,1,0);
    m11=Mcensat(k,1,1);
    p2h=PowerSpec(k)*(m10*m10+2*m10*m11+m11*m11);
    kPowerTable[i]=log10(k);
    PowerTable[i]=log10((p1h+p2h)/(gsl_spline_eval(ellipSpline,kPowerTable[i],ellipAcc)+1));
   // fprintf(fa,"%f %f\n",kPowerTable[i],PowerTable[i]);
  } //for

	// fclose(fa);

  gsl_spline_free(rsSpline);
  gsl_interp_accel_free(rsAcc);
  gsl_spline_free(alphaSpline);
  gsl_interp_accel_free(alphaAcc);
  for (i=5; i>=0; --i) {
    gsl_spline_free(NgalSpline[i]);
    gsl_interp_accel_free(NgalAcc[i]);
  } //for
  free(NgalSpline);
  free(NgalAcc);

  NewpowAcc=gsl_interp_accel_alloc();
  NewpowSpline=gsl_spline_alloc(gsl_interp_cspline,(NK+extralin));
  gsl_spline_init(NewpowSpline,kPowerTable,PowerTable,(NK+extralin));

  CorrTable=malloc((NR+8)*8*sizeof(double));
  rCorrTable=malloc((NR+8)*8*sizeof(double));

  for (i=0; i<(NR+8)*8; ++i)
  {
    rCorrTable[i]=pow(10.,(i+0.5)*(log10(150.)+2.3)/(float)(NR*8)-2.3);
    CorrTable[i]=corr_qawo(rCorrTable[i],0.01,0.99) + corr_qawo(rCorrTable[i],1.,99.) + corr_qawo(rCorrTable[i],100.,900.) +
    			       corr_qawo(rCorrTable[i],1e3,5e3) + corr_qawo(rCorrTable[i],6e3,5e3);
  } //for

  CorrAcc=gsl_interp_accel_alloc();
  CorrSpline=gsl_spline_alloc(gsl_interp_cspline,(NR+8)*8);
  gsl_spline_init(CorrSpline,rCorrTable,CorrTable,(NR+8)*8);

#ifdef OUTPUTCORR
  char buf[500];
  float mingalmass,maxgalmass;
  mingalmass=8.77+(ThisTask%6)*0.5;
  maxgalmass=8.77+(ThisTask%6+1)*0.5;
  sprintf(buf,"corr_%.2f-%.2f.dat",mingalmass,maxgalmass);
  fd=fopen(buf,"w");
  for (i=0; i<1000; ++i) {
    fprintf(fd,"%g   %g\n",pow(10.,0.001*i*5.7-3),gsl_spline_eval(CorrSpline,pow(10.,0.001*i*5.7-3),CorrAcc));
  } //for
  fclose(fd);
#endif


  for (i=0; i<NR; ++i)
  {
    r=pow(10.,(i+0.5)*(log10(80.)+2.1)/(float)(NR-1)-2.3);
    corrtmp=proj_corr(r);
    r_arr[i]=r/hubble;
    proj_arr[i]=corrtmp/hubble;
  } //for

  gsl_spline_free(CorrSpline);
  gsl_interp_accel_free(CorrAcc);
  gsl_spline_free(NewpowSpline);
  gsl_interp_accel_free(NewpowAcc);
  free(kPowerTable);
  free(PowerTable);
  free(cutoff_high);
  free(cutoff_low);
} //halomodel





void init_numgal(float masslimit_low, float masslimit_high, int snap) {
  FILE *fd;
  int i,j,jj;
  double *MassTable;
  int mbin,mbin2,k,found,ncen,nsat;
  double **NgalTable;
  double loghubble=log10(Hubble_h);
  char buf[500];
  int NMassTable2=8;
  double p[2];
  double perror[2];
  double *alphaerr;
  double *rserr;
  double *x;
  double *y;
  struct fitvars v;
  int status;
  mp_par pars[2];
  mp_result result;
  mp_config config;
  int NRadiusTable=25;
  int NParams=400;
  int NParams2=800;
  int rbin,besta,bestb,bestr,bestc,NgalTotalTotal,mstart,mend,NValid;
  int *NgalTotal,*NfofTotal,*usedbymass, *besti,*bestj,**NProfileTable;
  double rbinsize=0.1;
  double boxsize=BoxSize;
  double G=4.302e-9;
  double massoffset=log10(2*G/(Delta*Omega*1e4));
  double logrvir,logr,relx,rely,relz,sigma1,minsig1,sigma2,minsig2,lognorm,mbinsize,minfofmass2,maxfofmass2;
  double *MassTable2,*MassVarTable2,**RadiusTable,**ProfileTable,**ProfileErrTable,*alpha_m,*rs_m,*beta,*gamma,*weights;


  for (i=0; i<6; ++i)
  {
    cutoff_low[i]=0;
    cutoff_high[i]=0;
  } //for


  NgalAcc=malloc(6*sizeof(gsl_interp_accel*));
  NgalSpline=malloc(6*sizeof(gsl_spline*));

  usedbymass=malloc(massbins*sizeof(int));
  for(i=0;i<massbins;i++)
  	usedbymass[i]=0;

  for (i=0; i<6; ++i)
  {
    NgalAcc[i]=gsl_interp_accel_alloc();
    NgalSpline[i]=gsl_spline_alloc(gsl_interp_cspline,massbins);
  } //for


  MassTable=malloc(massbins*sizeof(double));
  NgalTable=malloc(6*sizeof(double*));

  for (i=0; i<6; ++i)
  	NgalTable[i]=malloc(massbins*sizeof(double));
  for (j=0; j<massbins; ++j)
  {
    MassTable[j]=0;
    for (i=0; i<6; ++i)
    	NgalTable[i][j]=0;
  } //for

  minfofmass2=minfofmass;
  maxfofmass2=maxfofmass;

  for (mbin=0; mbin<massbins; ++mbin)
  {
  	for (j=0; j<NFofsInSample[snap]; ++j)
  	{
  		ncen=0;
  		nsat=0;

  		//printf("mass bins=%f %f\n",minfofmass2+mbin*(maxfofmass2-minfofmass2)/(double)massbins,
  		//		minfofmass2+(mbin+1)*(maxfofmass2-minfofmass2)/(double)massbins);
  		if(MCMC_FOF2[j].M_Crit200[snap]>minfofmass+mbin*(maxfofmass-minfofmass)/(double)massbins &&
  				MCMC_FOF2[j].M_Crit200[snap]<minfofmass+(mbin+1)*(maxfofmass-minfofmass)/(double)massbins)
  		{
  			i=MCMC_FOF2[j].IndexOfCentralGal[snap];

  			if (MCMC_GAL[HashTable[i]].StellarMass[snap]>masslimit_low && MCMC_GAL[HashTable[i]].StellarMass[snap]<=masslimit_high)
  				ncen=1;

  			for (k=1; k<MCMC_GAL[HashTable[i]].ngal[snap]; ++k)
  			{
  				if (MCMC_GAL[HashTable[i+k]].StellarMass[snap]>masslimit_low && MCMC_GAL[HashTable[i+k]].StellarMass[snap]<=masslimit_high)
  					nsat++;
  			} //for
  			usedbymass[mbin]++;
  		}
  		NgalTable[0][mbin]+=ncen+nsat;
  		NgalTable[1][mbin]+=(ncen+nsat)*max((ncen+nsat)-1,0);
  		NgalTable[2][mbin]+=ncen;
  		NgalTable[3][mbin]+=nsat;
  		NgalTable[4][mbin]+=ncen*nsat;
  		NgalTable[5][mbin]+=nsat*max(nsat-1,0);
  		if (MCMC_FOF2[j].M_Crit200[snap]<minfofmass2 && nsat>0)	minfofmass2=MCMC_FOF2[j].M_Crit200[snap];
  		if (MCMC_FOF2[j].M_Crit200[snap]>maxfofmass2 && nsat>0)	maxfofmass2=MCMC_FOF2[j].M_Crit200[snap];
    } //for (loop on FoFID)
  	//printf("%f %f %f", minfofmass+(mbin+0.5)*(maxfofmass-minfofmass)/(double)massbins, NgalTable[2][mbin],NgalTable[3][mbin]);
  } //for (loop on halo mass bins)


  //printf("\n\n\nMass function\n");
  for (j=0; j<massbins; ++j)
  {
  	if (usedbymass[j]>0)
    {
      //MassTable[j]/=usedbymass[j];
      for (i=0; i<6; ++i)
      	NgalTable[i][j]/=usedbymass[j];
    } //if
   // else
    	MassTable[j]=minfofmass+(j+0.5)*(maxfofmass-minfofmass)/(double)massbins;
    	//printf("%f %f %f %d\n", MassTable[j], NgalTable[2][j],NgalTable[3][j], usedbymass[j]);
  } //for

  //printf("\n\n\n");

  for (i=0; i<6; ++i) {
    for (mbin=0; mbin<massbins-1 && cutoff_low[i]==0; ++mbin) {
      if (NgalTable[i][mbin]==0 && NgalTable[i][mbin+1]>0)
      {
      	cutoff_low[i]=MassTable[mbin];
      } //if
    } //for
    for (mbin=massbins-1; mbin>0 && cutoff_high[i]==0; --mbin) {
      if (NgalTable[i][mbin-1]>0 && NgalTable[i][mbin]==0)
      {
      	cutoff_high[i]=MassTable[mbin];
      } //if
    } //for
  } //for

  for (i=0; i<6; ++i)
  	gsl_spline_init(NgalSpline[i],MassTable,NgalTable[i],massbins);
  for (i=5; i>=0; --i)
  	free(NgalTable[i]);

  free(NgalTable);
  free(MassTable);


  MassTable2=malloc(NMassTable2*sizeof(double));
  MassVarTable2=malloc(NMassTable2*sizeof(double));
  RadiusTable=malloc(NMassTable2*sizeof(double*));
  NProfileTable=malloc(NMassTable2*sizeof(int*));
  ProfileTable=malloc(NMassTable2*sizeof(double*));
  ProfileErrTable=malloc(NMassTable2*sizeof(double*));

  for (i=0; i<NMassTable2; ++i)
  {
    RadiusTable[i]=malloc(NRadiusTable*sizeof(double));
    NProfileTable[i]=malloc(NRadiusTable*sizeof(int));
    ProfileTable[i]=malloc(NRadiusTable*sizeof(double));
    ProfileErrTable[i]=malloc(NRadiusTable*sizeof(double));

    for (j=0; j<NRadiusTable; ++j)
    {
      RadiusTable[i][j]=0;
      NProfileTable[i][j]=0;
      ProfileTable[i][j]=0;
      ProfileErrTable[i][j]=0;
    } //for

  } //for

  mbinsize=(maxfofmass2-minfofmass2)/(double)NMassTable2;
  NgalTotal=malloc(NMassTable2*sizeof(int));
  NfofTotal=malloc(NMassTable2*sizeof(int));

  for (i=0; i<NMassTable2; ++i)
  {
    NgalTotal[i]=0;
    NfofTotal[i]=0;
    MassTable2[i]=0;
    MassVarTable2[i]=0;
  } //for

 // for (mbin=0; mbin<massbins; ++mbin)
 // {
    for (j=0; j<NFofsInSample[snap]; ++j)
    {

    //	if(MCMC_FOF2[j].M_Crit200[snap]>minfofmass2+mbin*(maxfofmass2-minfofmass2)/(double)NMassTable2 &&
    //  				MCMC_FOF2[j].M_Crit200[snap]<minfofmass2+(mbin+1)*(maxfofmass2-minfofmass2)/(double)NMassTable2)
   //  {

    	i=MCMC_FOF2[j].IndexOfCentralGal[snap];

    	mbin2=floor((MCMC_GAL[HashTable[i]].fofmass[snap]-minfofmass2)/mbinsize);
    	//printf("mbin2=%d fofmass=%f FofId=%d i=%d\n",mbin2,MCMC_GAL[HashTable[i]].fofmass[snap],j,i);

      	if (mbin2==NMassTable2)
      		mbin2--; //because of maxfofmass2/roundoff

      	if (mbin2>=0 && mbin2<NMassTable2 && MCMC_GAL[HashTable[i]].ngal[snap]>1)
      	{
      		found=0;
      		logrvir=1./3.*(MCMC_GAL[HashTable[i]].fofmass[snap]+massoffset);

      		for (jj=1; jj<MCMC_GAL[HashTable[i]].ngal[snap]; ++jj)
      		{
      			if (MCMC_GAL[HashTable[i+jj]].StellarMass[snap]>masslimit_low && MCMC_GAL[HashTable[i+jj]].StellarMass[snap]<=masslimit_high)
      			{
      				found++;
      				//printf("Found satellite FOFID=%d SatType=%d SatMass=%f\n",
      				//		MCMC_GAL[HashTable[i]].fofid[snap], MCMC_GAL[HashTable[i+jj]].Type[snap], MCMC_GAL[HashTable[i+jj]].StellarMass[snap]);

      				relx=min(fabs(MCMC_GAL[HashTable[i+jj]].x[snap]-MCMC_GAL[HashTable[i]].x[snap]),
      						fabs(boxsize-fabs(MCMC_GAL[HashTable[i+jj]].x[snap]-MCMC_GAL[HashTable[i]].x[snap])));
      				rely=min(fabs(MCMC_GAL[HashTable[i+jj]].y[snap]-MCMC_GAL[HashTable[i]].y[snap]),
      						fabs(boxsize-fabs(MCMC_GAL[HashTable[i+jj]].y[snap]-MCMC_GAL[HashTable[i]].y[snap])));
      				relz=min(fabs(MCMC_GAL[HashTable[i+jj]].z[snap]-MCMC_GAL[HashTable[i]].z[snap]),
      						fabs(boxsize-fabs(MCMC_GAL[HashTable[i+jj]].z[snap]-MCMC_GAL[HashTable[i]].z[snap])));

      				logr=0.5*log10(relx*relx+rely*rely+relz*relz);
      				rbin=floor((logr-logrvir+2.5)/rbinsize);

      				if (rbin>=0 && rbin<NRadiusTable)
      				{
      					RadiusTable[mbin2][rbin]+=logr-logrvir;
      					NProfileTable[mbin2][rbin]++;
      				} //if
      			} //if
      		} //for

      		if (found>0)
      		{
      			NgalTotal[mbin2]+=found;
      			NfofTotal[mbin2]++;
      			MassTable2[mbin2]+=found*MCMC_GAL[HashTable[i]].fofmass[snap];
      			MassVarTable2[mbin2]+=found*SQR(MCMC_GAL[HashTable[i]].fofmass[snap]);
      		} //if

      	} //if

   //   }//if (FOF is in current mass bin)
    } //for (loop on FOFIds)
 // } //for (loop on halo mass bins)



  NgalTotalTotal=0;

  for (mbin2=0; mbin2<NMassTable2; ++mbin2)
  {
    if (NgalTotal[mbin2]>0)
    	MassTable2[mbin2]/=NgalTotal[mbin2];
    else
    	MassTable2[mbin2]=minfofmass2+(mbin2+0.5)*mbinsize;

    MassVarTable2[mbin2]-=NgalTotal[mbin2]*SQR(MassTable2[mbin2]);
    MassVarTable2[mbin2]/=(NgalTotal[mbin2]-1);

    NgalTotalTotal+=NgalTotal[mbin2];

    for (rbin=0; rbin<NRadiusTable; ++rbin)
    {
      if (NProfileTable[mbin2][rbin]>0)
      {
      	RadiusTable[mbin2][rbin]/=(double)NProfileTable[mbin2][rbin];
      	ProfileTable[mbin2][rbin]=(double)NProfileTable[mbin2][rbin]/(NgalTotal[mbin2] * 4./3.*PI*(pow(10.,3*((rbin+1)*rbinsize-2.5))-pow(10.,3*(rbin*rbinsize-2.5)))*pow(10.,MassTable2[mbin2]+massoffset));
      	ProfileErrTable[mbin2][rbin]=sqrt(NProfileTable[mbin2][rbin])/(NgalTotal[mbin2]*4./3.*PI*(pow(10.,3*((rbin+1)*rbinsize-2.5))-pow(10.,3*(rbin*rbinsize-2.5)))*pow(10.,MassTable2[mbin2]+massoffset));
      } //if
      else
      	RadiusTable[mbin2][rbin]=(rbin+0.5)*rbinsize-2.5;

    } //for
  } //for

  alpha_m=malloc(NMassTable2*sizeof(double));
  rs_m=malloc(NMassTable2*sizeof(double));
  alphaerr=malloc(NMassTable2*sizeof(double));
  rserr=malloc(NMassTable2*sizeof(double));

  memset(&result,0,sizeof(result));
  memset(&config,0,sizeof(config));
  memset(&pars[0],0,sizeof(pars));

  pars[0].limited[0]=1;
  pars[0].limits[0]=-2.3;
  pars[0].limited[1]=1;
  pars[1].limited[0]=1;
  pars[1].limits[0]=-9;
  pars[1].limited[1]=1;

  for (mbin2=0; mbin2<NMassTable2; ++mbin2)
  {
    NValid=0;
    for (rbin=0; rbin<NRadiusTable; ++rbin)
    	if (NProfileTable[mbin2][rbin]>0)
    		NValid++;

    if (NValid>2)
    {
      x=malloc(NValid*sizeof(double));
      y=malloc(NValid*sizeof(double));
      weights=malloc(NValid*sizeof(double));
      i=0;
      for (rbin=0; rbin<NRadiusTable; ++rbin)
      {
        if (NProfileTable[mbin2][rbin]>0)
        {
        	x[i]=pow(10.,RadiusTable[mbin2][rbin]);
        	y[i]=log10(ProfileTable[mbin2][rbin]);
        	weights[i]=max(1./sqrt(MassVarTable2[mbin2]+SQR(0.5*(log10(ProfileTable[mbin2][rbin]+ProfileErrTable[mbin2][rbin])-log10(ProfileTable[mbin2][rbin]-ProfileErrTable[mbin2][rbin])))),1./sqrt(MassVarTable2[mbin2]+SQR(3*(log10(ProfileTable[mbin2][rbin]+ProfileErrTable[mbin2][rbin])-log10(ProfileTable[mbin2][rbin])))));

        	if (isnan(weights[i]))
        	{
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

      p[0]=-1.;
      p[1]=-2.;

      config.ftol=1e-12;
      config.xtol=1e-12;
      config.maxiter=1000;

      status=mpfit(einastofit,NValid,2,p,pars,&config,(void *) &v,&result);

      alpha_m[mbin2]=p[0];
      rs_m[mbin2]=p[1];
      alphaerr[mbin2]=perror[0];
      rserr[mbin2]=perror[1];

      if (alphaerr[mbin2]==0) alphaerr[mbin2]=0.5;
      if (rserr[mbin2]==0) rserr[mbin2]=2.;
      if (alphaerr[mbin2]<0.1) alphaerr[mbin2]=0.1;
      if (rserr[mbin2]<0.3) rserr[mbin2]=0.3;

      free(weights);
      free(y);
      free(x);

    } //if
  } //for

  mstart=0;
  mend=NMassTable2;
  NValid=0;

  for (mbin2=0; NValid<3 && mbin2<NMassTable2; ++mbin2)
  {
    NValid=0;
    for (rbin=0; rbin<NRadiusTable; ++rbin)
    	if (NProfileTable[mbin2][rbin]>0)
    		NValid++;
    if (NValid<3) mstart++;
  } //for

  NValid=0;
  for (mbin2=NMassTable2-1; NValid<3 && mbin2>=0; --mbin2)
  {
    NValid=0;
    for (rbin=0; rbin<NRadiusTable; ++rbin)
    	if (NProfileTable[mbin2][rbin]>0)
    		NValid++;
    if (NValid<3) mend--;
  } //for

  while (mend-mstart<5) {
    if (NMassTable2-mend>mstart) {
      mend++;
      MassTable2[mend]=MassTable2[mend-1]+0.1;
      alpha_m[mend]=alpha_m[mend-1];
      alphaerr[mend]=alphaerr[mend-1]*10;
      rs_m[mend]=rs_m[mend-1];
      rserr[mend]=rserr[mend-1]*10;
    } //if
    else {
      mstart--;
      MassTable2[mstart]=MassTable2[mstart+1]-0.1;
      alpha_m[mstart]=alpha_m[mstart+1];
      alphaerr[mstart]=alphaerr[mstart+1]*10;
      rs_m[mstart]=rs_m[mstart+1];
      rserr[mstart]=rserr[mstart+1]*10;
    } //else
  } //while



  alphaAcc=gsl_interp_accel_alloc();
  alphaSpline=gsl_spline_alloc(gsl_interp_akima,mend-mstart);
  rsAcc=gsl_interp_accel_alloc();
  rsSpline=gsl_spline_alloc(gsl_interp_akima,mend-mstart);

  gsl_spline_init(alphaSpline,&MassTable2[mstart],&alpha_m[mstart],mend-mstart);
  gsl_spline_init(rsSpline,&MassTable2[mstart],&rs_m[mstart],mend-mstart);

  alpharscutoff_low=MassTable2[mstart];
  alpharscutoff_high=MassTable2[mend-1];

  free(rserr);
  free(alphaerr);
  free(rs_m);
  free(alpha_m);
  free(NfofTotal);
  free(NgalTotal);

  for (i=NMassTable2-1; i>=0; --i)
  {
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


void init_power()
{
  FILE *fd;
  char buf[500];
  double k,p;
  int NPowerTable=0;

  sprintf(buf,"%s/powrealized_rebin_corrected.dat",MCMCHaloModelDir);
  if(!(fd=fopen(buf,"r"))) {
  	char sbuf[1000];
  	sprintf(sbuf, "can't open %s \n", buf);
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
    if (fscanf(fd," %lg %lg ",&k,&p)==2)
    {
      kPowerTable[NPowerTable]=k;
      PowerTable[NPowerTable]=p;
      NPowerTable++;
    } //if
    else
    	break;
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
  if(!(fd=fopen(buf,"r")))
  {
  	char sbuf[1000];
  	sprintf(sbuf, "can't open %s \n", buf);
   	terminate(sbuf);
   } //if

  do {
    if (fscanf(fd," %lg %lg ",&k,&p)==2)
    	NPowerTable++;
    else
    	break;
  } //do
  while(1);
  fclose(fd);

  PowerTable=malloc(NPowerTable*sizeof(double));
  kPowerTable=malloc(NPowerTable*sizeof(double));

  fd=fopen(buf,"r");
  NPowerTable=0;
  do { //k and Delta
    if (fscanf(fd," %lg %lg ",&k,&p)==2) {
      kPowerTable[NPowerTable]=k;
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
    MSigmaTable[i]=pow(10,i*15./(float)NSigmaTable+5.);
    R=Radius(MSigmaTable[i]);
    SigmaTable[i]=TopHatSigma2(R);
  } //for
  SigmaAcc=gsl_interp_accel_alloc();
  SigmaSpline=gsl_spline_alloc(gsl_interp_cspline,NSigmaTable);
  gsl_spline_init(SigmaSpline,MSigmaTable,SigmaTable,NSigmaTable);
} //init_sigma



void initialize_halomodel() {
  int i,NnuTable=100;
  double res;
  double MnuTable[100],nuTable[100];
  rho_c=2.775e11; //recalculated value (3*H^2/8*PI*G in units of (Msun/h)/(Mpc/h)^3)
  Delta=180.;
  Delta_invth=1./pow(Delta,1./3.);
  rho_mean=Omega*rho_c;
  delta_c=1.674; //Millennium cosmology
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
} //initialization

int einastofit(int m,int n,double *p,double *dy,double **dvec,void *vars) {
  int i;
  double lognorm=-((2.-3./pow(10.,p[0]))*log10(2.)+(2./pow(10.,p[0]))/log(10.)+log10(PI)+3*log10(pow(10.,p[1]))+(3./pow(10.,p[0])-1)*log10(pow(10.,p[0]))+gsl_sf_lngamma(3./pow(10.,p[0]))/log(10.)+log10(gsl_sf_gamma_inc_P(3./pow(10.,p[0]),2*pow(pow(10.,p[1]),-pow(10.,p[0]))/pow(10.,p[0])))+fitconst);
  struct fitvars *v=(struct fitvars *) vars;
  double *x,*y,*w,f;
  x=v->x; //should be r/rvir (NOT log)
  y=v->y; //should be log10(profile)
  w=v->w; //should be a standard deviation in log space
  for (i=0; i<m; i++) {
    f=lognorm-2./pow(10.,p[0])*(pow(x[i]/pow(10.,p[1]),pow(10.,p[0]))-1)/log(10.);
    dy[i]=(y[i]-f)*w[i];
  } //for
  return 0;
} //einastofit



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
  gsl_integration_qawo_table *t=gsl_integration_qawo_table_alloc(r,L,GSL_INTEG_SINE,30);
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
  gsl_integration_qag(&F,sigma,500.,0,1.0e-3,WORKSIZE,GSL_INTEG_GAUSS51,w,&result,&abserr);
  gsl_integration_workspace_free(w);
  return result;
} //proj_corr

double Radius(double m) {
  return pow(m/rho_mean*3./4./PI,1./3.);
} //Radius

double Sigma2(double m) {
  return gsl_spline_eval(SigmaSpline,m,SigmaAcc);
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
  return rho_mean/(m*m)*0.5*gammaM(m)*0.315*exp(-pow(fabs(0.61-0.5*log(Sigma2(m))),3.8));
} //nbargal

double b(double m,int i) {
  double nu2,delta_c2,b1;
  delta_c2=delta_c*delta_c;
  nu2=delta_c2/Sigma2(m);
  //Tinker et al. (2010)
  double AT,aT,BT,bT,CT,cT;
  double y=log10(180.);
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
  double rs=(params->rs);
  double alpha=(params->alpha);
  double nm=exp(-2./alpha*(pow(r/rs,alpha)-1));
  return r*nm;
} //mugal_qawo_func

double mugal_qawo(double k,double m) {
  double rvir=Radius(m)*Delta_invth;
  double rvir3=rvir*rvir*rvir;
  double alpha=pow(10.,alpha_eval(log10(m)));
  double rs=pow(10.,rs_eval(log10(m)));
  double norm=1./(rvir3*PI*exp((2.-3./alpha)*log(2.)+2./alpha+3*log(rs)+(3./alpha-1)*log(alpha)+gsl_sf_lngamma(3./alpha)+log(gsl_sf_gamma_inc_P(3./alpha,2*pow(rs,-alpha)/alpha))));
  struct mugal_qawo_params params={ rs,alpha };
  double result=0,abserr;
  gsl_function F;
  int status;
  F.function=&mugal_qawo_func;
  F.params=&params;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(WORKSIZE);
  gsl_integration_qawo_table *t=gsl_integration_qawo_table_alloc(k*rvir,1.0,GSL_INTEG_SINE,25);
  status=gsl_integration_qawo(&F,0.0,0,1.0e-3,WORKSIZE,w,t,&result,&abserr);
  gsl_integration_qawo_table_free(t);
  gsl_integration_workspace_free(w);
  return norm*4*PI*rvir*rvir*result/k;
} //mugal_qawo

double NgalF(double m,int j) {
  double lm=log10(m);
  double value;
  if (lm<cutoff_low[j] || lm>cutoff_high[j]) return 0.;
  else {
    value=gsl_spline_eval(NgalSpline[j],lm,NgalAcc[j]);
    if (value>0) return value;
    else return 0.;
  } //else
} //NgalF

double alpha_eval(double m) {
  double value;
  if (m<alpharscutoff_low) return gsl_spline_eval(alphaSpline,alpharscutoff_low,alphaAcc);
  else if (m>alpharscutoff_high) return gsl_spline_eval(alphaSpline,alpharscutoff_high,alphaAcc);
  else {
    value=gsl_spline_eval(alphaSpline,m,alphaAcc);
    if (value>-2.3) return value;
    else return -2.3;
  } //else
} //alpha_eval

double rs_eval(double m) {
  double value;
  if (m<alpharscutoff_low) return gsl_spline_eval(rsSpline,alpharscutoff_low,rsAcc);
  else if (m>alpharscutoff_high) return gsl_spline_eval(rsSpline,alpharscutoff_high,rsAcc);
  else {
    value=gsl_spline_eval(rsSpline,m,rsAcc);
    if (value>-9.) return value;
    else return -9.;
  } //else
} //rs_eval

double Mcensat_func(double lm,void *p) {
  struct M_params *params=(struct M_params *) p;
  double k=(params->k);
  int i=(params->i);
  int j=(params->j);
  double m=exp(lm);
  if (i==0) {
    if (j==0) return 0.;
    else if (j==1) return nbargal(m)*NgalF(m,4)/pow(ngal_mean,2)*mugal_qawo(k,m)*m;
    else return nbargal(m)*NgalF(m,5)/pow(ngal_mean,2)*pow(mugal_qawo(k,m),2)*m;
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

