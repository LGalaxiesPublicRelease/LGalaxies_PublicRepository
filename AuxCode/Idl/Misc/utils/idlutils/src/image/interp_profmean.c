#include <stdio.h>
#include <string.h>
#include <math.h>

/* Interpolate a profmean radial profile */

#define GAMMA 2.5
#define PI 3.14159265358979

/* interpolation info */
static int nprof;
static float gmma=GAMMA, *brek=NULL, *s=NULL, *coef=NULL;
static int ltaut, ktaut, iflag, jderiv;

float ppvalu_(float *brek, float *coef, int *l, int *k,
							 float *x, int *jderiv);
void tautsp_(float r[], float f[], int *n, float *gmma,
						 float *s, float *brek, float *coef, int *l,
						 int *k, int *iflag);

/* radii of annuli for radial profile */
static int nrad;
static float rscale,pmscale;
static float *radius=NULL;
static float *aaradius=NULL;

static float *profmean=NULL;
static float *profcum=NULL;
static float *aaprofcum=NULL;

/* sets up profile interpolation */
int init_interp_profmean(int in_nprof,
                          float *in_profmean,
                          float *profradius,
                          int in_nrad,
                          float *maxradius, 
                          float in_rscale, 
                          float in_pmscale)
{
	int i,j,initialized_profinterp;

	/* initialize some variables if necessary */
  rscale=in_rscale;
  pmscale=in_pmscale;
	nprof=in_nprof;
	nrad=in_nrad;
  initialized_profinterp=0;
	if(!initialized_profinterp) {
    radius=(float *) malloc((nrad+4)*sizeof(float));
    aaradius=(float *) malloc((nrad+4)*sizeof(float));
    profmean=(float *) malloc((nrad+4)*sizeof(float));
    profcum=(float *) malloc((nrad+4)*sizeof(float));
    aaprofcum=(float *) malloc((nrad+4)*sizeof(float));
		/* put radius in arcsec; asinh for interpolation later */
    radius[0]=aaradius[0]=0.;
		for(i=1;i<=nrad;i++) {
			radius[i]=profradius[i-1];
			aaradius[i]=asinh(radius[i]/rscale);
		} /* end for i */
		
		/* allocate some memory */
		s=(float *) malloc((nrad)*6*sizeof(float));
		brek=(float *) malloc(4*(nrad)*sizeof(float));
		coef=(float *) malloc((nrad)*4*(nrad)*sizeof(float));
		initialized_profinterp=1;
	} /* end if */
	
	if(nprof<=0) {  /* if no profile exists, flux is necessarily zero */
		(*maxradius)=0.;
	} else {  /* else we do the interpolation */
		
		/* copy profile */
		for(j=0;j<nprof;j++) 
			profmean[j]=in_profmean[j];
		for(j=nprof;j<nrad+4;j++) 
			profmean[j]=0.;
    
		/* return maximum radius to look at */
		(*maxradius)=radius[nprof]; 
    
		/* make sure extrapolation is zero */
		for(j=nprof;j<nprof+1;j++)
			profmean[j]=0.;
		nprof+=1;
		
		/* create cumulative profile */
		profcum[0]=0.;
		for(j=1;j<nprof;j++) 
			profcum[j]=profcum[j-1]+PI*profmean[j-1]
				*(radius[j]*radius[j]-radius[j-1]*radius[j-1]);
		for(j=0;j<nprof;j++) 
      aaprofcum[j]=asinh(profcum[j]/pmscale);
		
		/* create spline */
		ktaut=nprof;
		ltaut=4*nprof;
		tautsp_(aaradius,aaprofcum,&nprof,&gmma,s,
						brek,coef,&ltaut,&ktaut,&iflag);
	} /* end if */

  return(1);
} /* end initProfInterp */

/* returns cumulative maggies */
float interp_prof(float radius) 
{
	float val,ar;
	
	if(nprof==0) return(0.);
	jderiv=0;
	ar=asinh(radius/rscale);
	val=ppvalu_(brek,coef,&ltaut,&ktaut,&ar,&jderiv);
  val=pmscale*sinh(val);
	return(val);
} /* end flux */

