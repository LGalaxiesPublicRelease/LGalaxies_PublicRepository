// Variables but some parameters hidden in here too.
// May be better to split by function

// I have no idea why the following variable is not declared extern,
// but compilation breaks if you do so.
char *inputFile; // Name of the input file
extern double MinGalOutputMass;

extern int FirstFile;		/* first and last file for processing */
extern int LastFile;

extern int Ntrees;		/* number of trees in current file */
extern double AllocValue_MaxHaloGal;
extern double AllocValue_MaxGal;
extern double AllocValue_MaxGalTree;

extern int MaxGal;		/* Maximum number of galaxies allowed for Gal[] array */
extern int NHaloGal, MaxHaloGal;
extern int NGalTree, MaxGalTree;
extern int *HaloGalHeap;
extern int IndexStored;

extern int LastSnapShotNr;

extern int LastDarkMatterSnapShot;
#ifdef MR_PLUS_MRII //OPTION for MCMC
extern int LastDarkMatterSnapShot_MR;
extern int LastDarkMatterSnapShot_MRII;
#endif


extern char SpecPhotDir[512];
extern char PhotPrefix[50];
extern char SpecPhotIMF[50];
extern char McFile[512];
extern char FileWithFilterNames[512];
extern char CoolFunctionsDir[512];
extern char CosmologyTablesDir[512];
extern char OutputDir[512];
/* in case a second parameter is given as argument to the code, this will be taken as a
 * temporary outputdir to allow fast I/O. OutputDir will be replaced by this directory
 * and in the end everything will be moved to the FinalOutputDir (original OutputDir
 * given in input.par )*/
extern char FinalOutputDir[512];
extern char FileNameGalaxies[512];
extern char SimulationDir[512];
extern char FileWithOutputRedshifts[512];

extern char FileWithZList[512];
//variable used to scale to a different cosmology
extern char FileWithZList_OriginalCosm[512];
#ifdef MR_PLUS_MRII  //OPTION for MCMC
extern char FileWithZList_MR[512];
extern char FileWithZList_OriginalCosm_MR[512];
extern char FileWithZList_MRII[512];
extern char FileWithZList_OriginalCosm_MRII[512];
#endif

extern double ScalePos;
extern double ScaleMass;

#ifdef SPECIFYFILENR
extern char   FileNrDir[512];
extern int    ListInputFilrNr[111];
#endif

extern int TotHalos;
extern int TotGalaxies[NOUT];
extern int *TreeNgals[NOUT];

extern int *FirstHaloInSnap;

extern int *TreeNHalos;
extern int *TreeFirstHalo;

extern void *TreeAuxData;


extern double MaxMemSize;

extern size_t AllocatedBytes;
extern size_t HighMarkBytes;
extern size_t FreeBytes;

extern int ThisTask, NTask;

#ifdef GALAXYTREE
extern int GalCount;
extern int TotGalCount;
#ifdef NORMALIZEDDB
extern int TotGalSFHBinCount;
#endif
#endif

/* Cosmological parameters */
extern double BaryonFrac;
extern double Sigma8;
extern double Omega;
extern double OmegaLambda;
extern double Hubble_h;
extern double Omega_OriginalCosm;
extern double OmegaLambda_OriginalCosm;
extern double Hubble_h_OriginalCosm;
//SIMULATION RELATED
extern double PartMass;
extern double BoxSize;
extern double PartMass_OriginalCosm;
extern double BoxSize_OriginalCosm;
#ifdef MR_PLUS_MRII  //OPTION for MCMC
extern double PartMass_MR;
extern double BoxSize_MR;
extern double PartMass_OriginalCosm_MR;
extern double BoxSize_OriginalCosm_MR;
extern double PartMass_MRII;
extern double BoxSize_MRII;
extern double PartMass_OriginalCosm_MRII;
extern double BoxSize_OriginalCosm_MRII;
#endif


/* flags */
extern int StarFormationModel;
#ifdef H2_AND_RINGS
extern int H2FractionRecipe;
extern int SFRtdyn;
#endif
extern int FeedbackReheatingModel;
extern int FeedbackReheatingDeansityScaling;
extern int FateOfSatellitesGas;
extern int ReIncorporationModel;
#ifdef EXCESS_MASS
extern int InfallModel;
#endif
extern int ReionizationModel;
extern int BlackHoleGrowth;
extern int AGNRadioModeModel;
extern int DiskRadiusModel;
extern int DiskInstabilityModel;
extern int BHGrowthInDiskInstabilityModel;
extern int HotGasStripingModel;
extern int HotGasOnType2Galaxies;
extern int StarBurstModel;
extern int BulgeFormationInMinorMergersOn;
extern int MetallicityOption;

/* parameters */
extern double Reionization_z0;
extern double Reionization_zr;
extern double Yield;
extern double RecycleFraction;
extern double ThreshMajorMerger;
extern double MergerTimeMultiplier;
extern double RamPressureStrip_CutOffMass;
extern double RamPressureRadiusThreshold;
extern double SfrEfficiency;
extern double SfrColdCrit;
extern double SfrBurstEfficiency;
extern double SfrBurstSlope;
extern double AgnEfficiency;
extern double BlackHoleGrowthRate;
extern double BlackHoleDisruptGrowthRate;
extern double BlackHoleSeedMass;
extern double BlackHoleAccretionRate;
extern double BlackHoleCutoffVelocity;
extern double FeedbackReheatingEpsilon;
extern double ReheatPreVelocity;
extern double ReheatSlope;
extern double FeedbackEjectionEfficiency;
extern double EjectPreVelocity;
extern double EjectSlope;
extern double ReIncorporationFactor;
extern double ReincZpower;
extern double ReincVelocitypower;
extern double FracZSNIItoHot;
extern double FracZSNIatoHot;
#ifdef H2_AND_RINGS
extern double Clumpingfactor;
extern double Warmphasefactor;
extern double GasInflowVel;
#endif
#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN
extern double EnergySNcode, EnergySN;
extern double EnergySNIIcode, EnergySNII;
extern double EnergySNIacode, EnergySNIa;
extern double EnergyAGBcode, EnergyAGB;
#else
extern double EnergySNcode, EnergySN;
#endif
extern double EtaSNcode, EtaSN;

#ifdef H2_AND_RINGS
extern double RingRadius[RNUM];
extern double RingArea[RNUM];
extern double InverseRingArea[RNUM];
#endif

extern double
	UnitLength_in_cm,
	UnitTime_in_s,
	UnitVelocity_in_cm_per_s,
	UnitMass_in_g,
	RhoCrit,
	UnitPressure_in_cgs,
	UnitDensity_in_cgs,
	UnitCoolingRate_in_cgs,
	UnitEnergy_in_cgs,
	UnitTime_in_Megayears, //Using time as stored in the code, this gives Myr/h
	UnitTime_in_years,
	G,
	Hubble,
	a0, ar;

extern int ListOutputSnaps[NOUT];
extern float ListOutputRedshifts[NOUT];

extern double ZZ[MAXSNAPS];
extern double AA[MAXSNAPS];
//variable used to scale to a different cosmology
extern double AA_OriginalCosm[MAXSNAPS];

extern double Age[MAXSNAPS];

extern int    Zlistlen;

extern gsl_rng *random_generator;


extern int    NumMergers;


/*  tabulated stuff */
#ifdef STAR_FORMATION_HISTORY
/* SFH_ is the reference structure for storing the star formation histories in
 * logarithmic bins. It is computed in init.c generating a binning structure for
 * each snapshot/time step. In the code galaxy structures are adjusted with respect
 * to this structure at each step. */
extern double SFH_t[MAXSNAPS][STEPS][SFH_NBIN]; //Time to present at the lower edge of the bin (code units)
extern double SFH_dt[MAXSNAPS][STEPS][SFH_NBIN]; //Time width of the bin (code units)
extern int SFH_Nbins[MAXSNAPS][STEPS][SFH_NBIN]; //Number of bins merged in each bin (only useful for the merging algorithm)
extern int SFH_ibin[MAXSNAPS][STEPS]; //Last active bin
#ifdef DETAILED_METALS_AND_MASS_RETURN
extern double tau_t[STEPS*MAXSNAPS]; //Time-to-z=0 of every timestep in the code. (Used for SNe rates in yield_integrals.c)
extern double tau_dt[STEPS*MAXSNAPS];//Width of every timestep in the code. (Used for SNe rates in yield_integrals.c)
#endif
#endif //STAR_FORMATION_HISTORY

#ifdef DETAILED_METALS_AND_MASS_RETURN

//Number of interpolated points within the mass ranges for the four types of yield table:
#define LIFETIME_MASS_NUM 150
#define LIFETIME_Z_NUM 6
#define AGB_MASS_NUM 59 //55 //ROB: 59, when going from 0.85 to 7 Msun
#define AGB_Z_NUM 3
#ifdef PORTINARI
#define SNII_MASS_NUM 85  //ROB: 85, from 6 <= M[Msun] <= 120. Change SNII_MIN_MASS and SNII_MAX_MASS for shorter ranges.
#define SNII_Z_NUM 5
#endif
#ifdef CHIEFFI
#define SNII_MASS_NUM 81 //ROB: 56 if 7 <= M[Msun] <= 50. 81 if 7 <= M[Msun] <= 120. (NB: You can set SNII_MASS_NUM 81, and SNII_MAX_MASS 50. But DON"T put SNII_MASS_NUM > 81 ever!)
#define SNII_Z_NUM 6
#endif
#define SNIA_MASS_NUM 83 //48 //Number increased after extending range to cover M2 masses (07-02-12)

//Mass ranges for the different modes of ejection:
#define AGB_MIN_MASS 0.85
#define AGB_MAX_MASS 7.0 //6.0
#define SNIA_MIN_MASS 3.0
#define SNIA_MAX_MASS 16.0
#define SNIA_MIN_TIME 35.0*1.0e6
#define SNIA_MAX_TIME 21.0*1.0e9
#ifdef PORTINARI
#define SNII_MIN_MASS 7.0 //6.0
#define SNII_MAX_MASS 120.0
#endif
#ifdef CHIEFFI
#define SNII_MIN_MASS 7.0
#define SNII_MAX_MASS 120.0 //50.0
#endif

int ELETOBIGCOUNTA;
int FRACCOUNTA;

//Arrays that yield tables are written to:
double lifetimeMasses[LIFETIME_MASS_NUM];
double lifetimeMetallicities[LIFETIME_Z_NUM];
double lifetimes[LIFETIME_Z_NUM][LIFETIME_MASS_NUM];
double AGBMasses[AGB_MASS_NUM]; //Initial star masses [Msun]
double AGBMetallicities[AGB_Z_NUM]; //Initial star metallicities [Msun]
double AGBEjectedMasses[AGB_Z_NUM][AGB_MASS_NUM]; //Total mass ejected [Msun]
double AGBTotalMetals[AGB_Z_NUM][AGB_MASS_NUM]; //Total metal YIELD ejected [Msun]
double AGBYields[AGB_Z_NUM][11][AGB_MASS_NUM]; //YIELD ejected, for each element [Msun]
double SNIIMasses[SNII_MASS_NUM];
double SNIIMetallicities[SNII_Z_NUM];
double SNIIEjectedMasses[SNII_Z_NUM][SNII_MASS_NUM];
double SNIITotalMetals[SNII_Z_NUM][SNII_MASS_NUM];
double SNIIYields[SNII_Z_NUM][11][SNII_MASS_NUM];
#ifndef DTD
double SNIaMasses[SNIA_MASS_NUM];
double SNIaEjectedMasses[SNIA_MASS_NUM];
double SNIaTotalMetals[SNIA_MASS_NUM];
double SNIaYields[42][SNIA_MASS_NUM];
#else
double SNIaYields[42];
#endif

//Integrated yields arrays:
double NormSNIIMassEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
double NormSNIIMetalEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
#ifdef INDIVIDUAL_ELEMENTS
double NormSNIIYieldRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM][NUM_ELEMENTS];
#endif
double NormAGBMassEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
double NormAGBMetalEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
#ifdef INDIVIDUAL_ELEMENTS
double NormAGBYieldRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM][NUM_ELEMENTS];
#endif
double NormSNIaMassEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
double NormSNIaMetalEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
#ifdef INDIVIDUAL_ELEMENTS
double NormSNIaYieldRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM][NUM_ELEMENTS];
#endif

//Arrays used to plot SNe rates from SFH bins (yield_integrals.c):
double TheSFH[SFH_NBIN];
double SNIIRate[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
double SNIaRate[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
double AGBRate[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
//Arrays used to plot SNe rates from SFH-timesteps (calc_SNe_rates.c):
double TheSFH2[STEPS*MAXSNAPS];
double SNIIRate2[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
double SNIaRate2[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
double AGBRate2[STEPS*MAXSNAPS][LIFETIME_Z_NUM];


//IMF parameters (for chemical enrichment):
#define IMF_SLOPE 2.3 //2.15 //High-mass slope of Chabrier IMF. (2.3 = normal Chabrier IMF. <2.3 = top-heavy, >2.3 = bottom-heavy.)

//SNIa parameters:
#define A_FACTOR 0.04 //0.028 //Fraction of mass from all objects between SNIa_MIN_MASS and SNIA_MAX_MASS that comes from SN-Ia. //0.028 preferred in Yates+13.
//#define FRAC2HOT 0.9 //Fraction of material released by disk stars that goes straight into the HotGas. Res goes in ColdGas.
#ifdef DTD
//#define KALPHA 1.4765 //1.59203 //Now set in yield_integrals.c
//#define	F316 0.0384 //Integral of the IMF (by number) from 3.0 - 16.0 Msun //Now set in yield_integrals.c
#define SNIAEJECMASS 1.2300971 //Total mass (and total metals) ejected by a SNIa explosion in Msun //Value form original yield table (42 elements): 1.3740855. //Value when only considering 11 elements: 1.2300971
#ifdef BIMODALDTD
	//#define DTD_NORM 0.903206 //(26Myrs - 21Gyrs)
	//#define DTD_NORM 0.896668 //For P98 Z=0.02 lifetimes (35Myrs - 17Gyrs)
	#define DTD_NORM 0.900348 //(35Myrs - 21Gyrs)
#endif
#ifdef CUSTOMDTD
	//#define DTD_NORM 0.524836 //(26Myrs - 21Gyrs)
	//#define DTD_NORM 0.606746 //For P98 Z=0.02 lifetimes (35Myrs - 17Gyrs, for ~32% in prompt component)
	#define DTD_NORM 0.610431 //(35Myrs - 21Gyrs, for ~32% in prompt component)
#endif
#ifdef GAUSSIANDTD
	#define DTD_NORM = 1.0 //For P98 Z=0.02 lifetimes (35Myrs - 17Gyrs) and (35Myrs - 21Gyrs)
	#define TAUCHARAC 1.0 //Characteristic delay time for SNe-Ia (i.e. peak of Gaussian distribution) in Gyrs //default: 2.0
	#define SIGMA_TD 0.2*TAUCHARAC //0.2 for narrow-DTD, 0.5 for wide_DTD
#endif
#ifdef POWERLAWDTD
	//#define DTD_NORM 7.21863 // (26Myrs - 21Gyrs)
	#define DTD_NORM 6.72574 //6.72544 // (35Myrs - 21Gyrs)
	//#define DTD_NORM 6.56087 //For P98 Z=0.02 lifetimes (35Myrs - 17Gyrs)
	//#define DTD_NORM 6.22432 // (35Myrs - 11Gyrs)
	//#define DTD_NORM 6.02197 // (50Myrs - 17Gyrs)
	//#define DTD_NORM 6.35503 // (40Myrs - 17Gyrs)
	#define DTD_SLOPE -1.12 //Slope of power law, according to Maoz et al. (2012)
#endif
#ifdef RUITERDTD
	//#define DTD_NORM 1.09545 //(26Myrs - 21Gyrs)
	#define DTD_NORM 1.09545 //(35Myrs - 21Gyrs)
	//#define DTD_NORM 1.08422 //For P98 Z=0.02 lifetimes (35Myrs - 17Gyrs)
	#define TAUCHARAC 0.5 //Peak of Gaussian (prompt) component [in Gyrs]
	#define SIGMA_TD 0.2*TAUCHARAC //Width of Gaussian (prompt) component
	#define DTD_SLOPE -2.0 //Slope of power law (delayed) component (see Ruiter et al. 2012)
#endif
#endif


#endif //DETAILED_METALS_AND_MASS_RETURN

#ifdef COMPUTE_SPECPHOT_PROPERTIES
// SSP PHOT_TABLES - magnitues of starburst population as a function of age

#ifdef M05
#define SSP_NAGES 220		// Age grid of the SSP tables
#define SSP_NMETALLICITES 4			// Number of Metalicities used
#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define SSP_NLambda 1221
#endif
#endif

#ifdef BC03
#define SSP_NAGES 221		// Age grid of the SSP tables
#define SSP_NMETALLICITES 6			// Number of Metalicities used
#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define SSP_NLambda 1221
#endif
#endif

#ifdef CB07
#define SSP_NAGES 221		// Age grid of the SSP tables
#define SSP_NMETALLICITES 6			// Number of Metalicities used
#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define SSP_NLambda 1238
#endif
#endif

//table containing the Metallicity grid of the SSP tables (converted to log10)
extern double SSP_logMetalTab[SSP_NMETALLICITES];
//table containing the Age grid of the SSP tables (originally in years, converted to log10(internal time units 1e12 Yrs/h))
extern double SSP_logAgeTab[SSP_NAGES];
//table containing redshift (different from the one in the code when scaling to future times)
extern double RedshiftTab[MAXSNAPS];
extern double LumTables[NMAG][SSP_NMETALLICITES][MAXSNAPS][SSP_NAGES];
extern double FilterLambda[NMAG+1];//wavelength of each filter + 1 for V-band

#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define MAX_NLambdaFilter 1000
extern int NLambdaFilter[NMAG];
//VEGA
#define NLambdaVega 3303
#endif

//DUST EXTINCTION
#define ExpTauBCBulge 0.5	// constant extinction for young stars in bulges.
#define MUWIDTH  0.2
#define MUCENTER 0.3
extern long mu_seed;

#endif //COMPUTE_SPECPHOT_PROPERTIES


extern size_t HighMark;

#ifdef UPDATETYPETWO
extern int NtotHalos, TotIds, Nids, TotSnaps, OffsetIDs;
extern int *CountIDs_halo, *OffsetIDs_halo, *CountIDs_snaptree, *OffsetIDs_snaptree;
extern long long *IdList;
extern float *PosList, *VelList;
#endif


extern int Hashbits;
extern int NumWrittenInParallel;
extern double ScaleFactor;	// factor by which to multiply a position to get its ph index (after floring)


#ifdef USE_MEMORY_TO_MINIMIZE_IO
extern char *ptr_auxdata, *ptr_treedata, *ptr_dbids, *ptr_galaxydata, *ptr_galsnapdata[NOUT];
extern size_t offset_auxdata, offset_treedata, offset_dbids;
extern size_t offset_galaxydata, maxstorage_galaxydata, filled_galaxydata;
extern size_t offset_galsnapdata[NOUT], maxstorage_galsnapdata[NOUT], filled_galsnapdata[NOUT];
#endif


extern float Reion_z[46],Reion_Mc[46];

extern FILE *tree_file;
extern FILE *treeaux_file;
extern FILE *treedbids_file;
extern FILE *FdGalTree;
extern FILE *FdGalTreeSFH;
extern FILE *FdGalDumps[NOUT];

/*H2 fraction table*/
#ifdef H2_AND_RINGS
//#define RHO_LEN 101
//#define Z_LEN 13
#define RHO_LEN 420
#define Z_LEN 6
extern double H2Fraction[LENZ][LENSIGMAH];
extern double H2Fraction_Zgrid[LENZ];
extern double H2Fraction_SigmaHgrid[LENSIGMAH];
#endif

#ifdef HDF5_OUTPUT
int b[NOUT];
struct GALAXY_OUTPUT galaxy_output_hdf5[NOUT][NRECORDS_APP];
#endif //HDF5_OUTPUT

#ifdef DEBUG
extern FILE *FdGalDebug;
#ifdef H2_AND_RINGS
#define NDebugProps 23
extern char DebugProperties[NDebugProps][100];
#else
#define NDebugProps 11
extern char DebugProperties[NDebugProps][100];
#endif
#endif
