//To understand the units in the code read through set_units in init.c!!!
#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10 // Dangerous to have such a short name
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24 // TODO this is read in from input.par
#define  PROTONMASS  1.6726e-24
#define  HUBBLE      3.2407789e-18   /* in h/sec */
#define  MUMH        0.59341*PROTONMASS // Could be a variable but we never change it!

//To understand the units in the code read through set_units in init.c!!!
#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#define MIN_ALLOC_NUMBER       1000
#define ALLOC_INCREASE_FACTOR  1.1
#define ALLOC_DECREASE_FACTOR  0.7
#define PRECISION_LIMIT 1.e-7
#define TINY_MASS 1.e-8
#define TINY_LENGTH 1.e-6

#ifdef GALAXYTREE
#undef  NOUT
#define NOUT MAXSNAPS
#endif

//WATCH OUT! In the case of MCMC running both MR and MRII the larger value is used to "allocate" all the arrays
//inside the code its LastDarkMatterSnapShot+1 that defines the extent of the loops
//(in MCMC MR_plus_MRII mode this are not always identical)
#ifdef MRII
#define  MAXSNAPS  68     /* Number of snapshots in the dark matter simulation */
#else

#ifdef PHOENIX
#define  MAXSNAPS  72
#else

#ifdef NIFTY
#define  MAXSNAPS  108
#else

#ifdef CATERPILLAR
#define  MAXSNAPS  320
#else

#define  MAXSNAPS  64  //NORMAL MILLENNIUM

#endif //CATERPILLAR
#endif //NIFTY
#endif //PHOENIX
#endif //else MRII

#define  MAXGALFAC 2.3 /*1.5/2.3 - maximum fraction of satellite without a halo (for memory allocation)*/


#ifdef FAST_TESTING_MODE
#define  STEPS 10
#else
#define  STEPS 20
#endif

#ifdef H2_AND_RINGS
#define RNUM 12          /* radially divide one disk into RNUM */
//#define RNUM 30          /* radially divide one disk into RNUM */
#define WARM_PHASE_FACTOR 1.3 // to use when deriving HI and H2 from total cold gas because 1/3 of it is ionized)
#define LENSIGMAH 40
#define LENZ 6
//#define LENSIGMAH 101
//#define LENZ 13
#endif

#ifdef METALRICHWIND
#ifdef GASDENSITYFWIND
#define NORMGASDENSITY 10. //Msun/pc^2 //ISM gas surface density to normalise to when calculating density-dependent direct ejection into HotGas (18-05-18)
#endif
#endif

				
#define  ALLOCPARAMETER 50.  /* new definition !!! THIS HAS TO BE 50 !!! DONT EVER EVER EVER CHANGE !!! */

#ifdef STAR_FORMATION_HISTORY
#define SFH_NMERGE 3  //  SFH_NMERGE=Nmax+1 (Nmax used in Shamshiri2014)
//#define SFH_NMERGE 10  //  SFH_NMERGE=Nmax+1 (Nmax used in Shamshiri2014)
#ifdef NIFTY
#define SFH_NBIN 22 //  NIFTY - 108 snapshots
#else
#ifdef CATERPILLAR
#define SFH_NBIN 24 //  CATERPILLAR - 320 snapshots
#else
#define SFH_NBIN 20
//#define SFH_NBIN 66
#endif //NIFTY
#endif //CATERPILLAR
#endif //STAR_FORMATION_HISTORY

