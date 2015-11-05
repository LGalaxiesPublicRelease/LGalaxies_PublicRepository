
// Routines used in the MCMC sampling
void Senna (void);
void print_parameters (int AcceptanceLogic, FILE *fmcmc);
double SAM(int filenr);
void initialize_mcmc_par_and_lhood (FILE *fmcmc);
void read_mcmc_par (int snapnum);
void read_sample_info(void);
void read_observations(void);
void open_files_with_comparison_to_observations(void);
void close_files_with_comparison_to_observations(void);
double propose_new_parameters(void);
void save_galaxy_for_mcmc(int gal_index);
#ifdef MR_PLUS_MRII
void change_dark_matter_sim(char SimName[]);
#endif

double get_likelihood (void);
void bin_function(int ObsNr, double *binsamdata, double *samdata, int snap);
void bin_red_fraction(int ObsNr, double *binredfraction, int snap);
void bin_passive_fraction(int ObsNr, double *binpassivefraction, int snap);
void bin_bulge_fraction(int ObsNr, double *binbulgefraction, int snap);
void bin_ColdGasFractionvsStellarMass(int ObsNr, double *bingasfraction, int snap);
void bin_color_hist(int ObsNr, double *bincolorhist, int snap);
void bin_bhbm(double *binblackholeup, double *binblackholedown, int snap);
void correct_for_correlation(int snap);
void compute_correlation_func(int ObsNr, double *binsamdata, int snap, float mingalmass, float maxgalmass);
double chi_square_probability(int ObsNr, double *samdata, int snap);
double maximum_likelihood_probability(int ObsNr, double *samfract, int snap);
double binomial_probability(int ObsNr, double *samup, double *samdown, int snap);


double gammp (double a, double x);
double gammq (double a, double x);
double gser (double a, double x);
double gcf(double a, double x);
double gammpapprox(double a, double x, int psig);
double gammln(double xx);
double betai(double a, double b, double x);
double betacf(double a,double  b, double x);
double ran3(long *idum);
double gassdev(long *idum);

