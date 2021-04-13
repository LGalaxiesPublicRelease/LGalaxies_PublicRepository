/*
 * Galaxy structure for output.
 *
 * NOTE: due to the way that the HDF5 builder routines work, variables should be commented out
 * with slash-star ... star-slash and not //.  Otherwise they are included in the properties list.
 *
 */

#ifdef LIGHT_OUTPUT

struct GALAXY_OUTPUT {
    int   Type; // None // Galaxy type: 0 for central galaxies of a main halo, 1 for central galaxies in sub-halos, 2 for satellites without halo.
    int   SnapNum; // None // Number of snapshot in which this galaxy resides.
    float CentralMvir; // 10^10/h Msun // Virial mass of background (FOF) halo containing this galaxy
    float CentralRvir; // 10^10/h Msun // Virial mass of background (FOF) halo containing this galaxy
    float Pos[3]; // 1/h Mpc // Galaxy position (comoving)
    float Mvir; // 10^10/h Msun // Virial mass of the subhalo the galaxy is/was the centre of.
    float Rvir; // Mpc/h // Virial radius of the subhalo the galaxy is/was the centre of.
    float Vvir; // km/s // Virial velocity of the subhalo the galaxy is/was the centre of.
    float DistanceToCentralGal[3]; // Mpc/h // Distance to central galaxy (proper/comoving?)
    /* baryonic reservoirs */
    float ColdGas; // 10^10/h Msun // Mass of cold gas.
    float DiskMass; // 10^10/h Msun // Mass of stars in the disk
    float BulgeMass; // 10^10/h Msun // Mass of stars in the bulge
    float HotGas; // 10^10/h Msun // Mass of hot gas
    /* float ReheatedGas; // 10^10/h Msun // Mass of Reheated gas */
    float BlackHoleGas; // 10^10/h Msun // Mass of gas surrounding black hole
    float BlackHoleMass; // 10^10/h Msun // Mass of black hole
    /* magnitudes in various bands */
#ifdef   COMPUTE_SPECPHOT_PROPERTIES
#ifdef     OUTPUT_OBS_MAGS
    float ObsMagDust[NMAG]; // AB mag // Dust corrected, observed-frame absolute mags
#endif     //OUTPUT_OBS_MAGS
#ifdef     OUTPUT_REST_MAGS
    float MagDust[NMAG]; // AB mag // Dust corrected, rest-frame absolute mags
#endif     //OUTPUT_REST_MAGS
#endif   //COMPUTE_SPECPHOT_PROPERTIES
};

#else  //LIGHT_OUTPUT

#pragma pack(1)  //structure alignment for 1 Byte.
struct GALAXY_OUTPUT {
#ifdef   NO_PROPS_OUTPUTS
#ifdef     GALAXYTREE
    long long GalID; // None // ID of galaxy, unique within simulation and SAM run.
#endif     //GALAXYTREE
#else    //NO_PROPS_OUTPUTS
#ifdef     GALAXYTREE
    long long GalID; // None // ID of galaxy, unique within simulation and SAM run.
    long long HaloID; // None // Unique ID of MPA halo containing this galaxy
#endif     //GALAXYTREE
#ifdef     MBPID
    long long MostBoundID; // None // Most bound particle at centre of subhalo last associated with this galaxy.  
#endif     //MBPID
#ifdef     GALAXYTREE
    long long FirstProgGal; // None // Main progenitor of this galaxy. Also the first progenitor in a linked list representation of the merger tree.
    long long NextProgGal; // None // Next progenitor of this galaxy in linked list representation of merger tree
    long long LastProgGal; // None // Galaxies with id between this galaxyId and this lastProgenitorId form the merger tree rooted in this galaxy.
    long long FOFCentralGal; // None // The galaxy id of the central galaxy of the FOF group this galaxy is in.
#ifdef       NIFTY
    long long FOFCentralID; // None // The galaxy id of the central galaxy of the FOF group this galaxy is in?
#endif       //NIFTY
    long long FileTreeNr; // None // Number of the tree file
    long long DescendantGal; // None // Pointer to the descendant of this galaxy in its merger tree; -1 if there is no descendant
    long long MainLeafId; // None //galaxyId of the leaf on the main branch this galaxy is part of. Obtained by following firstProgenitorId as far as it goes.
    long long TreeRootId; // None //The galaxyId of the galaxy at the root of the merger tree this galaxy is in. Especially useful for speeding up queries for descendants for a given progenitor.
    long long SubID; // None //Id of the subhalo containing this galaxy as given by the column subhaloFileID in the MillenniumII..SubHalo miniMilII..SubHalo table (for MRII and mMRII) and by the column subhaloId in the MField.FOFSubHalo millimil..FOFSubHalo tables (for MR and mMR). Alternative to haloId.
    long long MMSubID; // None // fofId, the subhaloid of the subhalo at the center of the fof group
    int   PeanoKey; // None // Peano-Hilbert key, (bits=8), for position in 500/h Mpc box
    float Redshift; // None // redshift of the snapshot where this galaxy resides
#endif     //GALAXYTREE
    int   Type; // None //Galaxy type: 0 for central galaxies of a main halo, 1 for central galaxies in sub-halos, 2 for satellites without halo.
#ifdef SAVE_FOFHALO
    int   FileNr; // None // File number
    int   TreeNr; // None // Tree number
    int   FoFHaloNr; // None // Halo 
#endif
#ifndef    GALAXYTREE
    int   HaloIndex; // None // ?Unique ID of MPA halo containing this galaxy
#endif     //GALAXYTREE
#ifdef     HALOPROPERTIES
    float HaloM_Mean200; // 1e10 Msun/h // M200 cf mean last time this halo was a type 0
    float HaloM_Crit200; // 1e10 Msun/h // M200 cf critical last time this halo was a type 0
    float HaloM_TopHat; // 1e10 Msun/h // Virial mass last time this halo was a type 0
    float HaloPos[3]; // Mpc/h // Comoving position of halo.
    float HaloVel[3]; // km/s // Mean velocity of halo.
    float HaloVelDisp; // km/s // Velocity dispersion of halo.
    float HaloVmax; // km/s // Maximum circular velocity of halo.
    float HaloSpin[3]; // km/s Mpc/h // specific spin of the halo.
#endif     //HALOPROPERTIES
    int   SnapNum; // None //The snapshot number where this galaxy was identified.
    float LookBackTimeToSnap; // yr // The time from a given snapshot to z=0
    float CentralMvir; // 10^10/h Msun // virial mass of background (FOF) halo containing this galaxy
    float CentralRvir; // Mpc/h // Proper[?] R200 cf critical of background (FOF) halo containing this galaxy
    float DistanceToCentralGal[3];  // Mpc/h // Proper[?] components of the distance between this galaxy and the galaxy at the centre of the FoF group.
    float Pos[3]; // 1/h Mpc // Comoving galaxy/subhalo position
    float Vel[3]; // km/s // Galaxy/subhalo peculiar velocity
    int   Len; // None // Number of particles in the associated subhalo  
    /* properties of subhalo at the last time this galaxy was a central galaxy */
    float Mvir; // 10^10/h Msun // M200 cf critical of the halo last time galaxy was type 0
    float Rvir; // Mpc/h // R200 cf critical of the halo last time galaxy was type 0
    float Vvir; // km/s // Virial velocity of the halo last time galaxy was type 0
    float Vmax; // km/s //Maximum rotational velocity of the subhalo, or the last value for type 2's galaxies.
    float ColdGasSpin[3]; // Mpc/h km/s // The specific angular momentum of the cold gas disk
    float DiskSpin[3]; // Mpc/h km/s // The specific angular momentum of the stellar disk
    float InfallVmax; // km/s // Maximum rotational velocity of the host halo of this galaxy at infall (ie last time a type 0)
    float InfallVmaxPeak; // km/s // ? Peak Vmax along past history
    int InfallSnap; // None // Most recent (largest) snapnum at which this galaxy's type changed from 0 to 1 or 2
    float InfallHotGas; // 10^10 Msun/h // Mass in hot gas at the time of infall (same as hotGas for type 0 galaxies).
    float HotRadius; // Mpc/h // Proper[?] radius out to which hot gas extends: rvir for type 0; 0 for type 2; maximum radius out to which hot gas is not stripped for type 1.
    /*dynamical friction merger time*/
#ifndef    HT09_DISRUPTION
    float OriMergTime; // yr // Estimated dynamical friction time when the merger clock is set.
    float MergTime; //yr // Estimated remaining merging time. 
#else      //HT09_DISRUPTION
    float OriMergRadius; // Mpc/h // ? Proper[?] radius when merger clock is set
    float MergRadius; // Mpc/h // ? Current proper[?] radius
#endif     //HT09_DISRUPTION
#ifdef TRACK_SPLASHBACKS
    int flagSplashBack; // None // SplashBack flag
    float TimeSinceSplashBack; // yr // Time since SplashBack
#endif
#ifdef TRACK_NMERGERS
    float NMajorMergers; // None // Number of major mergers
    float NMinorMergers; // None // Number of minor mergers
#endif
    /* baryonic reservoirs */
    float ColdGas; // 10^10/h Msun // Mass in cold gas.
#ifdef OUTPUT_RINGS
    float ColdGasRings[RNUM]; // 1e10 Msun/h //  Mass of clod gas in each annulur ring.
    float H2fraction; // None //  Fraction of ColdGas in the form of H_2
    float H2fractionRings[RNUM]; // None //  H2 fraction within each annular ring.
#endif     //H2_AND_RINGS
    float StellarMass; // 10^10/h Msun // Total mass in stars in the disk and the bulge combined
    float DiskMass; // 10^10/h Msun // Mass of stars in the disk
    float BulgeMass; // 10^10/h Msun // Mass of stars in the bulge
#ifdef OUTPUT_RINGS
    float DiskMassRings[RNUM]; // 1e10 Msun/h // Mass of stars within each annular ring
    float BulgeMassRings[RNUM]; // 1e10 Msun/h // Mass of stars within each annular ring
#endif     //H2_AND_RINGS
    float HotGas; // 10^10/h Msun // Mass in hot gas
    /*float ReheatedGas; // 10^10/h Msun // Mass in reheated gas*/
    float EjectedMass; // 10^10/h Msun // Mass in ejected gas
#ifdef     EXCESS_MASS
    float ExcessMass; // 10^10/h Msun // Mass in excess of universal baryon fraction
#endif     //EXCESS_MASS
    float BlackHoleGas; // 10^10/h Msun // Mass of gas surrounding central black hole
    float BlackHoleMass; // 10^10/h Msun // Mass of central black hole
    /* ICL magnitude and mass*/
    float ICM; //10^10/h Msun //Total mass in metals in intra-cluster stars, for type 0,1
#ifdef TRACK_MASSGROWTH_CHANNELS
    float MassFromInSitu; // 1e10 Msun/h // Mass formed in situ.
    float MassFromMergers; // 1e10 Msun/h // Mass accreted from mergers.
    float MassFromBursts; // 1e10 Msun/h // Mass formed in starbursts
#endif
#ifdef TRACK_BURST
    float BurstMass; // 1e10 Msun/h // Mass formed in starbursts
#endif //TRACK_BURST
    float MetalsColdGas[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in cold gas.
#ifdef OUTPUT_RINGS
    float MetalsColdGasRings[RNUM][NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in cold gas in each annular ring
#endif       //H2_AND_RINGS
    float MetalsStellarMass[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in the disk
    float MetalsDiskMass[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in the disk
    float MetalsBulgeMass[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in the bulge
#ifdef OUTPUT_RINGS
    float MetalsDiskMassRings[RNUM][NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in stars in each annular ring
    float MetalsBulgeMassRings[RNUM][NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in stars in each annular ring
#endif       //H2_AND_RINGS
    float MetalsHotGas[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in the hot gas
    /* float MetalsReheatedGas[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in the Reheated gas */
    float MetalsEjectedMass[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in the ejected gas
#ifdef EXCESS_MASS
    float MetalsExcessMass[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals associated with ExcessMass component
#endif       //EXCESS_MASS
    float MetalsICM[NUM_METAL_CHANNELS]; // 10^10/h Msun // Mass in metals in intra-cluster stars, for type 0,1
#ifdef       METALS_SELF
    float MetalsHotGasSelf[NUM_METAL_CHANNELS]; // 10^10/h Msun // hot gas metals that come from self
#endif       //METALS_SELF

    /* misc */
    float PrimordialAccretionRate; // Msun/yr // Accretion rate of primordial gas.
    float CoolingRadius; // Mpc/h // The radius within which the cooling time scale is shorter than the dynamical timescale
    /* float CoolingGas; // 10^10/h Msun // Mass of cooling gas */
#ifdef COOLING_TEST
    float CoolingTest; // 10^10/hMsun? // Whatever is returned by CoolingTest
#ifdef BETAPROF
    float dt_ratio; // None // betaprof diagnostic
    float tau_ratio; // None // betaprof diagnostic
#endif
#endif
    float CoolingRate; // Msun/yr // Cooling rate of the hot gas
    float CoolingRate_beforeAGN; // Msun/yr // What the cooling rate of the hot gas would have been if there was no AGN feedback.
    float QuasarAccretionRate; // Msun/yr // Rate at which cold gas is accreted into the central black hole in the quasar mode.
    float RadioAccretionRate; // Msun/yr // Rate at which hot gas is accreted into the central black hole in the radio mode.
    float Sfr; // Msun/yr // Star formation rate
#ifdef OUTPUT_RINGS
    float SfrRings[RNUM]; // Msun/yr // Star formation rate within each annular ring
#endif     //H2_AND_RINGS
    float SfrBulge; // Msun/yr // Star formation rate in bulge.
    float XrayLum; // log10(erg/sec) // (log_10 of) X-Ray luminosity
    float BulgeSize; // Mpc/h // Half mass radius of bulge
    float DiskRadius; // Mpc/h // Size of the stellar disk, 3x the scale length.
    float ColdGasRadius; // Mpc/h // Size of the gas disk, 3x the scale length.
    float StellarHalfMassRadius; // Mpc/h // stellar Half mass radius
#ifdef COMPUTE_SPECPHOT_PROPERTIES
    float StellarHalfLightRadius; // Mpc/h // stellar Half light radius
    float CosInclination; // deg // Inclination of the galaxy. Derived from the angle between the stellar spins of the galaxy and the z-axis
#endif
    int   DisruptOn; // None // 0: galaxy merged onto merger center 1: galaxy was disrupted before merging onto its descendant, matter went into ICM of merger center;
#ifdef     MERGE01
    int   MergeOn; // None // 0: standard delucia-like merger behaviour for type 1 galaxy; 1: galaxy mass > halo mass, separate dynamical friction time calculated ....
#endif     //MERGE01
#endif   //NO_PROPS_OUTPUTS

       /* magnitudes in various bands */
#ifdef   COMPUTE_SPECPHOT_PROPERTIES
#ifdef     OUTPUT_REST_MAGS
    float MagDust[NMAG]; // AB mag // dust corrected, rest-frame absolute mags
    float Mag[NMAG]; // AB mag // rest-frame absolute mags
    float MagBulge[NMAG]; // AB mag // rest-frame absolute mags for the bulge
#ifdef       ICL
    float MagICL[NMAG]; // AB mag // rest-frame absolute mags of ICL
#endif       //ICL
#endif     //OUTPUT_REST_MAGS
#ifdef     OUTPUT_OBS_MAGS
    float ObsMagDust[NMAG]; // AB mag // dust-corrected, obs-frame absolute mags
    float ObsMag[NMAG]; // AB mag // obs-frame absolute mags
    float ObsMagBulge[NMAG]; // AB mag // obs-frame absolute mags for the bulge
#ifdef       ICL
    float ObsMagICL[NMAG]; // AB mag // observer-frame absolute mags for intra-cluster light
#endif       //ICL
#ifdef       OUTPUT_MOMAF_INPUTS
    /* define luminosities as if the galaxy were one snapshot earlier, i.e. higher redshift, than its actual snapshot */
    float ObsMagDust[NMAG]; // AB mag // Dust-corrected obs-frame absolute mags
    float ObsMag[NMAG]; // AB mag // Obs-frame absolute mags
    float ObsMagBulge[NMAG]; // AB mag // Obs-frame absolute mags for the bulge
#ifdef         ICL
    float dObsMagICL[NMAG];  // AB mag // Obs-frame absolute mags of ICL
#endif	       //ICL
#ifdef         KITZBICHLER
    /* define luminosities as if the galaxy were one snapshot later, i.e. lower redshift, than its actual snapshot */
    float dObsMagDust_forward[NMAG]; // AB mag // Dust-corrected, one snapshot later-frame absolute mags
    float dObsMag_forward[NMAG]; // AB mag // One snapshot later-frame absolute mags
    float dObsMagBulge_forward[NMAG]; // AB mag // One snapshot later-frame absolute mags for the bulge
#ifdef           ICL
    float dObsMagICL_forward[NMAG]; // AB mag // One snapshot later absolute mags for intra-cluster light 
#endif	         //ICL
#endif         //KITZBICHLER
#endif	     //OUTPUT_MOMAF_INPUTS
#endif	   //OUTPUT_OBS_MAGS
#endif   //COMPUTE_SPECPHOT_PROPERTIES

    float MassWeightAge; //10^9yr //The age of this galaxy weighted by mass of its components.
#ifdef   POST_PROCESS_MAGS
    float rbandWeightAge; // 10^9yr // The age of this galaxy weighted by mass of its components.
#endif   //POST_PROCESS_MAGS
#ifdef OUTPUT_SFH
    int sfh_ibin; // None // Index of highest star formation history bin currently in use
    int sfh_numbins; // None // Number of non empty star formation history bins
#ifndef    NORMALIZEDDB
    /* float sfh_time[SFH_NBIN]; // yr // lookback time to middle of star formation history bin. */
    /* float sfh_dt[SFH_NBIN]; // yr // Width of star formation history bin. */
    float sfh_DiskMass[SFH_NBIN]; // 10^10 Msun/h // Star formation history in the disk.
    float sfh_BulgeMass[SFH_NBIN]; // 10^10 Msun/h // Star formation history in the bulge.
#ifdef OUTPUT_RINGS
    float sfh_DiskMassRings[RNUM][SFH_NBIN]; // 10^10 Msun/h // Star formation history in the disk RINGS.
    float sfh_BulgeMassRings[RNUM][SFH_NBIN]; // 10^10 Msun/h // Star formation history in the bulge RINGS.
#endif
    float sfh_ICM[SFH_NBIN]; // 10^10 Msun/h // Star formation history in intra-cluster stars.
    float sfh_MetalsDiskMass[SFH_NBIN][NUM_METAL_CHANNELS]; // 10^10 Msun/h // Metal formation history in the disk.
    float sfh_MetalsBulgeMass[SFH_NBIN][NUM_METAL_CHANNELS]; // 10^10 Msun/h // Metal formation history in the bulge.
    float sfh_MetalsICM[SFH_NBIN][NUM_METAL_CHANNELS]; // 10^10 Msun/h // Metal formation history in the ICM.
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
    float sfh_MassFromInSitu[SFH_NBIN]; // 10^10 Msun/h // Star formation history of stars formed in situ.
    float sfh_MassFromMergers[SFH_NBIN]; // 10^10 Msun/h // Star formation history of stars accreted from mergers.
    float sfh_MassFromBursts[SFH_NBIN]; // 10^10 Msun/h // Star formation history of stars formed in starbursts.
#endif
#ifdef TRACK_BURST
    float sfh_BurstMass[SFH_NBIN]; // 10^10 Msun/h // Star formation history of stars formed in starbursts.
#endif //TRACK_BURST
#endif     //NORMALIZEDDB
#endif   //OUTPUT_SFH

#ifdef OUTPUT_ELEMENTS
    //All: [H][He][Cb][N][O][Ne][Mg][Si][S][Ca][Fe] or //Only [H][He][O][Mg][Fe]
#ifdef OUTPUT_SFH
    float sfh_DiskMass_elements[SFH_NBIN][NUM_ELEMENTS]; // Msun // History of mass of elements locked up in stars in disk.
    float sfh_BulgeMass_elements[SFH_NBIN][NUM_ELEMENTS]; // Msun // History of mass of elements locked up in stars in bulge.
    float sfh_ICM_elements[SFH_NBIN][NUM_ELEMENTS]; // Msun // History of mass of elements locked up in stars in the ICM.
#endif
    float DiskMass_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in stars in disk.
    float BulgeMass_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in stars in bulge.
#ifdef OUTPUT_RINGS
    float DiskMassRings_elements[RNUM][NUM_ELEMENTS]; // Msun // Mass of elements locked up in stars in each annular ring.
    float BulgeMassRings_elements[RNUM][NUM_ELEMENTS]; // Msun // Mass of elements locked up in stars in each annular ring.
#endif     //H2_AND_RINGS
    float ColdGas_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in stars in cold gas.
#ifdef OUTPUT_RINGS
    float ColdGasRings_elements[RNUM][NUM_ELEMENTS]; // Msun // Mass of elements locked up in cold gas in each annular ring.
#endif     //H2_AND_RINGS
    float HotGas_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in hot gas.
    /* float ReheatedGas_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in reheated gas. */
    float ICM_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in stars in the ICM
    float EjectedMass_elements[NUM_ELEMENTS]; // Msun // Mass of elements locked up in ejected gas.
#ifdef     EXCESS_MASS
    float ExcessMass_elements[NUM_ELEMENTS]; // Msun // Mass in elements associated with ExcessMass component
#endif     //EXCESS_MASS
#endif   //OUTPUT_ELEMENTS
};

// next only of interest to DB output, which generally requires complete tree
#ifdef OUTPUT_SFH
struct SFH_BIN {
    long long GalID; // None // ID of the galaxy
    short snapnum; // None // snapnum of the galaxy, repeated here for faster lookups of times etc
    short sfh_ibin; // None //Index of highest bin currently in use
    /* float sfh_time; // yr // Lookback time at the middle of bin. */
    /* float sfh_dt; // yr // time width of bin. */
    float sfh_DiskMass; // 1e10 Msun/h // SFH of disk
    float sfh_BulgeMass; // 1e10 Msun/h // SFH of bulge
#ifdef OUTPUT_RINGS
    float sfh_DiskMassRings[RNUM]; // 10^10 Msun/h // Star formation history in the disk RINGS.
    float sfh_BulgeMassRings[RNUM]; // 10^10 Msun/h // Star formation history in the bulge RINGS.
#endif
    float sfh_ICM; // 1e10 Msun/h // SFH of ICM
#ifdef TRACK_SFH_MASSGROWTH_CHANNELS
    float sfh_MassFromInSitu; // 10^10 Msun/h // Star formation history of stars formed in situ.
    float sfh_MassFromMergers; // 10^10 Msun/h // Star formation history of stars accreted from mergers.
    float sfh_MassFromBursts; // 10^10 Msun/h // Star formation history of stars formed in starbursts.
#endif
#ifdef TRACK_BURST
    float sfh_BurstMass; // 10^10 Msun/h // Star formation history of stars formed in starbursts.
#endif //TRACK_BURST
    float sfh_MetalsDiskMass[NUM_METAL_CHANNELS]; // 1e10 Msun/h // Metals locked up in stars in disk.
    float sfh_MetalsBulgeMass[NUM_METAL_CHANNELS]; // 1e10 Msun/h // Metals locked up in stars in bulge.
    float sfh_MetalsICM[NUM_METAL_CHANNELS]; // 1e10 Msun/h // Metals locked up in stars in ICM.
};

struct SFH_Time {
    int snapnum; // None // snapnum
    int bin; // None // index of current bin
    /* proposal: in output write the start of the bin and its end, rather than centre and dt */
    double lookbacktime; // yr // lookback time to centre of current bin
    double dt; // yr // width of the current bin
    int nbins; // None // number of highest resolution bins used to create current bin
};
#endif   //OUTPUT_SFH
#pragma pack()   //structure alignment ends.

#endif //LIGHT_OUTPUT
