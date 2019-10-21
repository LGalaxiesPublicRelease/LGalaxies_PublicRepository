
#define SQR(a) ((a)*(a))
#define CUB(a) ((a)*(a)*(a))

void halomodel(double* r_arr,double* proj_arr,float masslimit_low, float masslimit_high, int snap);
double NewPowerSpec(double k);
double corr_qawo_func(double k,void *params);
double corr_qawo(double r,double a,double L);
double proj_corr_func(double r,void *p);
double proj_corr(double sigma);
double Radius(double m);
double Sigma2(double m);
double PowerSpec(double k);
double sigma2_func(double k,void *params);
double TopHatSigma2(double R);
double TopHatWindow(double kr);
double Mass(double R);
double lnSigma2(double lR,void *params);
double lnSigma2lnM(double lm,void *params);
double gammaM(double m);
double nbargal(double m);
double b(double m,int i);
double mugal_qawo_func(double r,void *p);
double mugal_qawo(double k,double m);
double NgalF(double m,int j);
double alpha_eval(double m);
double rs_eval(double m);
double Mcensat_func(double lm,void *p);
double Mcensat(double k,int i,int j);
double ngal_mean_func(double lm,void *p);
double ngal_mean_calc(int j);
void init_power();
void init_sigma();
void init_numgal(float masslimit_low, float masslimit_high, int snap);
void initialize_halomodel();
int einastofit(int m,int n,double *p,double *dy,double **dvec,void *vars);

struct fitvars {
  double *x;
  double *y;
  double *w;
};
