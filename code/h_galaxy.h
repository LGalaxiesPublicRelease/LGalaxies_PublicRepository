/*Structure with all the data associated with galaxies (this is not the same as the output!)*/
struct GALAXY			/* Galaxy data */
{
  int HeapIndex;
  int GalTreeIndex;
  int NextGalaxy;
#ifdef GALAXYTREE
  int FirstProgGal;
#endif
  int Type;
  int HaloNr;
  long long MostBoundID;
  int SnapNum;
  int CentralGal;        /* Galaxy at the centre of this subhalo.
			  * Own ID for types 0 & 1. */
  int FOFCentralGal;     // Galaxy at centre of FOF group; ie the type 0.
  int MergerCentralGal;  /* Galaxy this galaxy will merge onto.  Own ID for
			  * types 0 and 1, unless 1's have little dark matter
			  * and already merging to type 0. */
  double CentralMvir;
  /* properties of subhalo at the last time this galaxy was a central galaxy */
  double Pos[3];
  double MergCentralPos[3];
  double Vel[3];
  double Pos_notupdated[3];
  double Vel_notupdated[3];
#ifdef HALOPROPERTIES
  double HaloM_Mean200;
  double HaloM_Crit200;
  double HaloM_TopHat;
  double HaloPos[3];
  double HaloVel[3];
  double HaloVelDisp;
  double HaloVmax;
#endif
  double HaloSpin[3];
  double ColdGasSpin[3];
  double DiskSpin[3];
  int   Len;   
  double Mvir;
  double Rvir;
  double Vvir;
  double Vmax;
  double InfallVmax;
  double InfallVmaxPeak; // km/s - Max previous Vmax at infall
  int InfallSnap;
  double InfallHotGas;
  double InfallHotGasRadius;
  double HotRadius;
  /* baryonic reservoirs */
  double ColdGas;
#ifdef H2_AND_RINGS
 //double the make higher precisions for multirings
  double ColdGasRings[RNUM];
  double H2fraction;
  double H2fractionRings[RNUM];
#endif
  double DiskMass;
  double BulgeMass; 
#ifdef H2_AND_RINGS
  double DiskMassRings[RNUM];
#ifdef RINGS_IN_BULGES
  double BulgeMassRings[RNUM];
#endif
#endif
  double HotGas;
  //double ReheatedGas;
  double EjectedMass;
#ifdef EXCESS_MASS
  double ExcessMass;
#endif
  double BlackHoleMass;
  double BlackHoleGas;
  double ICM;

  double MetalsColdGas[NUM_METAL_CHANNELS];
#ifdef H2_AND_RINGS
  double MetalsColdGasRings[RNUM][NUM_METAL_CHANNELS];
#endif
  double MetalsDiskMass[NUM_METAL_CHANNELS];
  double MetalsBulgeMass[NUM_METAL_CHANNELS];
#ifdef H2_AND_RINGS
  double MetalsDiskMassRings[RNUM][NUM_METAL_CHANNELS];
#ifdef RINGS_IN_BULGES
  double MetalsBulgeMassRings[RNUM][NUM_METAL_CHANNELS];
#endif
#endif
  double MetalsHotGas[NUM_METAL_CHANNELS];
  //double MetalsReheatedGas;
  double MetalsEjectedMass[NUM_METAL_CHANNELS];
#ifdef EXCESS_MASS
  double MetalsExcessMass[NUM_METAL_CHANNELS];
#endif
#ifdef METALS_SELF
  double MetalsHotGasSelf[NUM_METAL_CHANNELS];
#endif
  double MetalsICM[NUM_METAL_CHANNELS];

#ifdef TRACK_MASSGROWTH_CHANNELS
  double MassFromInSitu;
  double MassFromMergers;
  double MassFromBursts;
#endif
#ifdef TRACK_BURST
  double BurstMass;
#endif //TRACK_BURST


  /* misc */
  double PrimordialAccretionRate;
  double CoolingRate;
  double CoolingRate_beforeAGN;
  double CoolingRadius;
  double CoolingGas;
  double QuasarAccretionRate;
  double RadioAccretionRate;
  //double RadioMass;
  double AGNheatingFromCentral;
#ifdef H2_AND_RINGS
  double Sfr;
  double SfrRings[RNUM];
#else
  double Sfr;
#endif
  double SfrBulge;
  double StarMerge;
  double XrayLum;
  double BulgeSize;
  double DiskRadius;
  double ColdGasRadius;
  double StellarHalfMassRadius;
#ifdef GALAXYTREE
  int   DisruptOn;
#endif
  // double halfradius;
  //double periradius;
#ifdef COMPUTE_SPECPHOT_PROPERTIES
  double CosInclination; //angle between galaxy spin and the z-axis
#endif
#ifndef HT09_DISRUPTION
  double OriMergTime;
  double MergTime;
  double OriMvir;
  double OriRvir;
#else
  double OriMergRadius;
  double MergRadius;
  double OriMergmass;
#endif
  double MergeSat;
  double DistanceToCentralGal[3];
  int MergeOn;
#ifdef TRACK_SPLASHBACKS
  int flagSplashBack;
  double TimeSinceSplashBack;
#endif
#ifdef TRACK_NMERGERS
  double NMajorMergers;
  double NMinorMergers;
#endif

  /* luminosities in various bands */
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
#ifdef OUTPUT_REST_MAGS
  double Lum[NMAG][NOUT];
  double YLum[NMAG][NOUT];
  double LumBulge[NMAG][NOUT];
  double YLumBulge[NMAG][NOUT];
  double LumDust[NMAG][NOUT];
#ifdef ICL
  double ICLLum[NMAG][NOUT];
#endif
#endif //OUTPUT_REST_MAGS

#ifdef COMPUTE_OBS_MAGS
  double ObsLum[NMAG][NOUT];
  double ObsYLum[NMAG][NOUT];
  double ObsLumBulge[NMAG][NOUT];
  double ObsYLumBulge[NMAG][NOUT];
  double ObsLumDust[NMAG][NOUT];
#ifdef ICL
  double ObsICL[NMAG][NOUT];
#endif

#ifdef OUTPUT_MOMAF_INPUTS
  double dObsLum[NMAG][NOUT];
  double dObsYLum[NMAG][NOUT];
  double dObsLumBulge[NMAG][NOUT];
  double dObsYLumBulge[NMAG][NOUT];
  double dObsLumDust[NMAG][NOUT];
#ifdef ICL
  double dObsICL[NMAG][NOUT];
#endif
#endif
#endif //COMPUTE_OBS_MAGS

#endif //ndef POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

  double MassWeightAge[NOUT];
#ifdef STAR_FORMATION_HISTORY
  int sfh_ibin; //Index of highest bin currently in use
  double sfh_age; //Time in years of last call to sph_update_bins
  double sfh_dt[SFH_NBIN]; //Size of time interval in units of years
  double sfh_t[SFH_NBIN]; //Time at low-redshift edge of bin in same units
  int sfh_Nbins[SFH_NBIN]; //Number of bins on the time interval
  double sfh_DiskMass[SFH_NBIN]; //Stellar mass in disk, in bin in standard units
  double sfh_BulgeMass[SFH_NBIN]; //Stellar mass in bulge, in bin in standard units
#ifdef H2_AND_RINGS
  double sfh_DiskMassRings[RNUM][SFH_NBIN]; //Stellar mass in disk RINGS, in bin in standard units
#ifdef RINGS_IN_BULGES
  double sfh_BulgeMassRings[RNUM][SFH_NBIN]; //Stellar mass in disk RINGS, in bin in standard units
#endif
#endif
  double sfh_ICM[SFH_NBIN]; //Stellar mass in ICM, in bin in standard units
  double sfh_MetalsDiskMass[SFH_NBIN][NUM_METAL_CHANNELS]; //Metals locked up in stars in disk.
  double sfh_MetalsBulgeMass[SFH_NBIN][NUM_METAL_CHANNELS]; //Metals locked up in stars in bulge.
  double sfh_MetalsICM[SFH_NBIN][NUM_METAL_CHANNELS]; //Metals locked up in stars in ICM.
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
  double sfh_MassFromInSitu[SFH_NBIN]; //Stellar mass formed in situ, in standard units.
  double sfh_MassFromMergers[SFH_NBIN]; //Stellar mass accreted from mergers, in standard units.
  double sfh_MassFromBursts[SFH_NBIN]; //Stellar mass formed in bursts, in standard units.
#endif
#ifdef TRACK_BURST
  double sfh_BurstMass[SFH_NBIN]; //Stellar mass formed in bursts, in standard units.
#endif //TRACK_BURST

#endif //STAR_FORMATION_HISTORY

#ifdef INDIVIDUAL_ELEMENTS
  //All: [H][He][Cb][N][O][Ne][Mg][Si][S][Ca][Fe] or //Only [H][He][O][Mg][Fe]
#ifdef STAR_FORMATION_HISTORY
  double sfh_DiskMass_elements[SFH_NBIN][NUM_ELEMENTS];
  double sfh_BulgeMass_elements[SFH_NBIN][NUM_ELEMENTS];
  double sfh_ICM_elements[SFH_NBIN][NUM_ELEMENTS];
#endif
  double DiskMass_elements[NUM_ELEMENTS];
  double BulgeMass_elements[NUM_ELEMENTS];
  double ColdGas_elements[NUM_ELEMENTS];
  double HotGas_elements[NUM_ELEMENTS];
  //double ReheatedGas_elements[NUM_ELEMENTS];
  double ICM_elements[NUM_ELEMENTS];
  double EjectedMass_elements[NUM_ELEMENTS];
#ifdef EXCESS_MASS
  double ExcessMass_elements[NUM_ELEMENTS];
#endif
#ifdef H2_AND_RINGS
  double DiskMassRings_elements[RNUM][NUM_ELEMENTS];
#ifdef RINGS_IN_BULGES
  double BulgeMassRings_elements[RNUM][NUM_ELEMENTS];
#endif
  double ColdGasRings_elements[RNUM][NUM_ELEMENTS];
#endif
#endif //INDIVIDUAL_ELEMENTS
} *Gal, *HaloGal;

