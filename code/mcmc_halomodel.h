#define SQR(a) ((a)*(a))
#define CUB(a) ((a)*(a)*(a))

void halomodel(double* r_arr,double* proj_arr,float masslimit_low,float
masslimit_high,int snap);
double TwoPowerSpec(double k,int censat);
double pconv_W_P_func(double theta,void *p);
double pconv_W_func(double lq,void *p);
double pconv_W(double k,double R,int censat);
double calc_mean_rhalo_simple_func(double lm,void *p);
double calc_mean_rhalo_simple(int censat,int norm);
double calc_mean_rsat_func(double lr,void *p);
double calc_mean_rsat(double lm,int extrar);
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
double nbargal(double m);
double b(double m,int i);
double mugal_qawo_func(double r,void *p);
double mugal_qawo(double k,double m);
double NgalF(double m,int j);
double pa_eval(double m);
double pb_eval(double m);
double pc_eval(double m);
double Mcensat_func(double lm,void *p);
double Mcensat(double k,int i,int j);
double ngal_mean_func(double lm,void *p);
double ngal_mean_calc(int j);
void init_power();
void init_sigma();
void init_numgal(float masslimit_low,float masslimit_high,int snap);
void initialize_halomodel();
double my_f(const gsl_vector *v,void *params);
void my_df(const gsl_vector *v,void *params,gsl_vector *df);
void my_fdf(const gsl_vector *x,void *params,double *f,gsl_vector *df);
void paramerror(double *x,double *p,double *perror);
int poissonfit(int m,int n,double *p,double *dy,double **dvec,void *vars);