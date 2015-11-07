
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
#include <gsl/gsl_deriv.h>



// Variables for the MCMC sampling

int MCMCNpar; //Number of parameters to sample
#define MCMCNConstraints 27  //Nr of Observational Constraints
#define MCMCMaxObsBins 100 //maximum number of bins per observation per redshift

double MCMCConstraintsZZ[NOUT];
//Nr of SAM galaxies used to compare with data
int TotMCMCGals[NOUT];
//To allocate structure with SAM Galaxies
#define MCMCAllocFactor 200000
//#define MCMCAllocFactor 1000000
long MCMCseed;
int Nbins[NOUT][MCMCNConstraints]; //bins on each obs test
double lhood1;

int NFofsInSample[NOUT];
int UsedFofsInSample[NOUT];

#ifdef MR_PLUS_MRII
int NTrees_Switch_MR_MRII;
int Switch_MR_MRII;
#endif
int CurrentMCMCStep;
int GlobalMCMCStep;

time_t global_starting_time;

FILE *FILE_MCMC_LIKELIHOOD;
FILE *FILE_MCMC_PredictionsPerStep[NOUT][MCMCNConstraints];

//READ FROM input.par
char MCMCParameterPriorsAndSwitches[512];
char MCMCObsConstraints[512];
char MCMCWeightsObsConstraints[512];
char ObsConstraintsDir[512];
char MCMCSampleDir[512];
char MCMCSampleFilePrefix[512];
int  MCMCSampleFile;
#ifdef MR_PLUS_MRII
char MCMCSampleFilePrefix_MR[512];
char MCMCSampleFilePrefix_MRII[512];
int  MCMCSampleFile_MR;
int  MCMCSampleFile_MRII;
#endif
int  MCMCTreeSampleFile;
int  ChainLength;
int  FirstChainNumber;
int  Time_Dependant_PhysPar;
int  MCMCMode;
double MCMC_LogStep_Size;
double MCMC_Initial_Par_Displacement;
double MCMC_Minimum_Obs_Error;
double MachineTimeOut;
char JobSubmitCommand[512];
char JobSubmitPipe[512];
char JobSubmitFile[512];

//Structures for the MCMC
struct MCMC_OBSCONSTRAINTS
{
	char Name[1000];
	char TestType[1000];
	//for chi-square like tests
	double Bin_low[NOUT][MCMCMaxObsBins];
	double Bin_high[NOUT][MCMCMaxObsBins];
	double Obs[NOUT][MCMCMaxObsBins];
	double Error[NOUT][MCMCMaxObsBins];
  //for binomial like tests
	double ObsUp[NOUT][MCMCMaxObsBins];
	double ObsDown[NOUT][MCMCMaxObsBins];
	int ObsTest_Switch_z[NOUT];
	double ObsTest_Weight_z[NOUT];
} *MCMC_Obs;

struct MCMC_GALAXY
{
  float StellarMass[NOUT];
  float ColdGas[NOUT];
  float BulgeMass[NOUT];
  float BlackHoleMass[NOUT];
  float Sfr[NOUT];
  float MagU[NOUT];
  float MagB[NOUT];
  float MagV[NOUT];
  float MagJ[NOUT];
  float MagK[NOUT];
  float Magu[NOUT];
  float Magg[NOUT];
  float Magr[NOUT];
  float Magi[NOUT];
  float Magz[NOUT];
  float Weight[NOUT];
} *MCMC_GAL;


struct MCMC_PAR
{
    char   Name[1000];
  	double Value[NOUT];
  	double PropValue[NOUT];
  	double PriorMin;
  	double PriorMax;
    int    Sampling_Switch;
} *MCMC_PAR;

struct MCMC_FOF_struct
{
	//values for the sample of FoF groups
	double Weight[NOUT];
	long long FoFID[NOUT];
} *MCMC_FOF;

