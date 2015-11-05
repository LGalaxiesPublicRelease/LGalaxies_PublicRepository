#include <stdio.h>
#include <stdlib.h>

#define NGAL 2
#define SFH_TIME_INTERVAL 1.e7
#define SFH_NBIN 20

struct GALAXY {
  float StellarMass;
  float MetalsStellarMass;
  int sfh_ibin; //Index of highest bin are currently in use
  double sfh_age; //Time of last call  to sph_update_bins
  int sfh_dt[SFH_NBIN]; //Size of time interval in units of SFH_TIME_INTERVAL
  int sfh_t[SFH_NBIN]; //Time at upper edge of bin in same units
  float sfh_StellarMass[SFH_NBIN]; //Stellar mass in bin in standard units
  float sfh_MetalsStellarMass[SFH_NBIN]; //Metals locked up in stars ditto
} Gal[NGAL];
