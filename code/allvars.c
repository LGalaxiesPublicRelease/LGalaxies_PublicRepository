
#include "allvars.h"


struct GALAXY			/* Galaxy data */
 *Gal, *HaloGal;

struct halo_data *Halo, *Halo_Data;

struct halo_aux_data		/* auxiliary halo data */
 *HaloAux;

struct halo_ids_data *HaloIDs, *HaloIDs_Data;


int FirstFile;			/* first and last file for processing */
int LastFile;


double AllocValue_MaxHaloGal;
double AllocValue_MaxGal;
double AllocValue_MaxGalTree;

int Ntrees;			/* number of trees in current file */

int MaxGal;
int NHaloGal, MaxHaloGal;
int NGalTree, MaxGalTree;
int *HaloGalHeap;
int IndexStored;

char SpecPhotDir[512];
char PhotPrefix[50];
char SpecPhotIMF[50];
char McFile[512];
char FileWithFilterNames[512];
char CoolFunctionsDir[512];
char OutputDir[512];
char FinalOutputDir[512];
char FileNameGalaxies[512];
char SimulationDir[512];
char FileWithOutputRedshifts[512];
char FileWithZList[512];
//variables used to scale to a different cosmology
char FileWithZList_OriginalCosm[512];
#ifdef MR_PLUS_MRII  //OPTION for MCMC
char FileWithZList_MR[512];
char FileWithZList_OriginalCosm_MR[512];
char FileWithZList_MRII[512];
char FileWithZList_OriginalCosm_MRII[512];
#endif

double ScalePos, ScaleMass;

#ifdef SPECIFYFILENR
char FileNrDir[512];
int ListInputFilrNr[111];
#endif


int TotHalos;
int TotGalaxies[NOUT];
int *TreeNgals[NOUT];

int LastSnapShotNr;
int LastDarkMatterSnapShot;
#ifdef MR_PLUS_MRII  //OPTION for MCMC
int LastDarkMatterSnapShot_MR;
int LastDarkMatterSnapShot_MRII;
#endif


int *FirstHaloInSnap;
int *TreeNHalos;
int *TreeFirstHalo;

double MaxMemSize;

size_t AllocatedBytes;
size_t HighMarkBytes;
size_t FreeBytes;

int ThisTask, NTask;


#ifdef GALAXYTREE
int GalCount;
int TotGalCount;
struct galaxy_tree_data *GalTree;
#endif

size_t HighMark;


/* cosmological parameters */
double BaryonFrac;
double Sigma8;
double Omega;
double OmegaLambda;
double Hubble_h;
double Omega_OriginalCosm;
double OmegaLambda_OriginalCosm;
double Hubble_h_OriginalCosm;
//SIMULATION RELATED
double PartMass;
double BoxSize;
double PartMass_OriginalCosm;
double BoxSize_OriginalCosm;

#ifdef MR_PLUS_MRII  //OPTION for MCMC
double PartMass_MR;
double BoxSize_MR;
double PartMass_OriginalCosm_MR;
double BoxSize_OriginalCosm_MR;
double PartMass_MRII;
double BoxSize_MRII;
double PartMass_OriginalCosm_MRII;
double BoxSize_OriginalCosm_MRII;
#endif


/* flags */
int ReionizationModel;
int DiskRadiusModel;
int StarFormationModel;
int FeedbackReheatingModel;
int FeedbackEjectionModel;
int FateOfSatellitesGas;
int ReIncorporationModel;
int AGNRadioModeModel;
int DiskInstabilityModel;
int BHGrowthInDiskInstabilityModel;
int HotGasStripingModel;
int DisruptionModel;
int StarBurstRecipe;
int BulgeFormationInMinorMergersOn;
int MetallicityOption;

/* parameters */
double Reionization_z0;
double Reionization_zr;
double RamPressureStrip_CutOffMass;
double SfrEfficiency;
double SfrColdCrit;
double SfrBurstEfficiency;
double SfrBurstSlope;
double Yield;
double RecycleFraction;
double ThreshMajorMerger;
double MergerTimeMultiplier;
double AgnEfficiency;
double BlackHoleGrowthRate;
double BlackHoleSeedMass;
double BlackHoleCutoffVelocity;
double FeedbackReheatingEpsilon;
double ReheatPreVelocity;
double ReheatSlope;
double FeedbackEjectionEfficiency;
double EjectPreVelocity;
double EjectSlope;
double ReIncorporationFactor;
double EnergySNcode, EnergySN;
double EtaSNcode, EtaSN;

double UnitTime_in_s,
  UnitPressure_in_cgs,
  UnitDensity_in_cgs,
  UnitCoolingRate_in_cgs,
  UnitEnergy_in_cgs,
  UnitTime_in_Megayears,
  UnitTime_in_years,
  G,
  Hubble,
  a0, ar;

int ListOutputSnaps[NOUT];
float ListOutputRedshifts[NOUT];


double ZZ[MAXSNAPS];
double AA[MAXSNAPS];
//variable used to scale to a different cosmology
double AA_OriginalCosm[MAXSNAPS];

double Age[MAXSNAPS];

int Zlistlen;

gsl_rng *random_generator;

int NumMergers;

/*  tabulated stuff */

/* fixed-metallicity spectrophotometric model */
/* tables hold magnitues of starburst population as a function of age */

#ifdef STAR_FORMATION_HISTORY
double SFH_t[MAXSNAPS][STEPS][SFH_NBIN];
double SFH_dt[MAXSNAPS][STEPS][SFH_NBIN];
int SFH_Nbins[MAXSNAPS][STEPS][SFH_NBIN];
int SFH_ibin[MAXSNAPS][STEPS];
#endif //STAR_FORMATION_HISTORY

#ifdef COMPUTE_SPECPHOT_PROPERTIES
//SSP PHOT TABLES
float SSP_logMetalTab[SSP_NMETALLICITES];
float SSP_logAgeTab[SSP_NAGES];
float RedshiftTab[MAXSNAPS];
float LumTables[NMAG][SSP_NMETALLICITES][MAXSNAPS][SSP_NAGES];
float FilterLambda[NMAG+1];	//wavelength of each filter + 1 for V-band
#ifdef SPEC_PHOTABLES_ON_THE_FLY
int NLambdaFilter[NMAG];
#endif


// dust
long mu_seed;
#endif

void *TreeAuxData;

#ifdef UPDATETYPETWO
int NtotHalos, TotIds, Nids, TotSnaps, OffsetIDs;
int *CountIDs_halo, *OffsetIDs_halo, *CountIDs_snaptree, *OffsetIDs_snaptree;
long long *IdList;
float *PosList, *VelList;
#endif


int Hashbits;
int NumWrittenInParallel;
double ScaleFactor;


#ifdef USE_MEMORY_TO_MINIMIZE_IO
char *ptr_auxdata, *ptr_treedata, *ptr_dbids, *ptr_galaxydata, *ptr_galsnapdata[NOUT];
size_t offset_auxdata, offset_treedata, offset_dbids;
size_t offset_galaxydata, maxstorage_galaxydata, filled_galaxydata;
size_t offset_galsnapdata[NOUT], maxstorage_galsnapdata[NOUT], filled_galsnapdata[NOUT];
#endif

/* reionization Okamoto et al. 2008*/
float Reion_z[46],Reion_Mc[46];

FILE *tree_file;
FILE *treeaux_file;
FILE *treedbids_file;
FILE *FdGalTree;
FILE *FdGalTreeSFH;
FILE *FdGalDumps[NOUT];




