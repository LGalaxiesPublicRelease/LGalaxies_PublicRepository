#include "allvars.h"
#include "proto.h"
#include "hdf5.h"
#include "hdf5_hl.h"
 
#define NRECORDS (hsize_t) 0
#define HDF5_STRING_SIZE 2048
 
/* File to set all the HDF5 table properties */ 
int ifield;
int rowsize=0;
hsize_t chunk_size=CHUNK_SIZE;
int * fill_data=NULL;
hid_t file_id;
size_t * output_offsets;
hid_t * field_types;
size_t * output_sizes;
size_t output_size;
 
// Write the datatype for HDF5 to write out the data
 char types[]={
'i',
'i',
'i',
'f',
'f',
'f',
'3',
'3',
'3',
'i',
'f',
'f',
'f',
'f',
'3',
'3',
'f',
'f',
'i',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'i',
'i',
'o',
'o',
'o',
'f',
'f',
'i',
'i',
's',
's',
's',
's',
's',
's',
}; 
 
const char * field_names[]={
"Type",
"HaloIndex",
"SnapNum",
"LookBackTimeToSnap",
"CentralMvir",
"CentralRvir",
"DistanceToCentralGal",
"Pos",
"Vel",
"Len",
"Mvir",
"Rvir",
"Vvir",
"Vmax",
"ColdGasSpin",
"DiskSpin",
"InfallVmax",
"InfallVmaxPeak",
"InfallSnap",
"InfallHotGas",
"HotRadius",
"OriMergTime",
"MergTime",
"ColdGas",
"StellarMass",
"DiskMass",
"BulgeMass",
"HotGas",
"ReheatedGas",
"EjectedMass",
"BlackHoleMass",
"ICM",
"MassFromInSitu",
"MassFromMergers",
"MassFromBursts",
"MetalsColdGas",
"MetalsStellarMass",
"MetalsDiskMass",
"MetalsBulgeMass",
"MetalsHotGas",
"MetalsReheatedGas",
"MetalsEjectedMass",
"MetalsICM",
"PrimordialAccretionRate",
"CoolingRadius",
"CoolingRate",
"CoolingRate_beforeAGN",
"QuasarAccretionRate",
"RadioAccretionRate",
"Sfr",
"SfrBulge",
"XrayLum",
"BulgeSize",
"DiskRadius",
"ColdGasRadius",
"StellarHalfMassRadius",
"StellarHalfLightRadius",
"CosInclination",
"DisruptOn",
"MergeOn",
"MagDust",
"Mag",
"MagBulge",
"MassWeightAge",
"rbandWeightAge",
"sfh_ibin",
"sfh_numbins",
"sfh_DiskMass",
"sfh_BulgeMass",
"sfh_ICM",
"sfh_MetalsDiskMass",
"sfh_MetalsBulgeMass",
"sfh_MetalsICM",
};
 
//The number of fields in the data
int nfields=73;

// Define the dimensions for the HDF5 table
hsize_t    float_dims[1]={3};
#ifdef COMPUTE_SPECPHOT_PROPERTIES
hsize_t    nmag_dims[1]={NMAG};
#endif
#ifdef STAR_FORMATION_HISTORY
hsize_t    sfh_dims[1]={SFH_NBIN};
hsize_t    sfh_3_dims[1]={3*SFH_NBIN};
#ifdef INDIVIDUAL_ELEMENTS
hsize_t    numele_dims[1]={NUM_ELEMENTS};
hsize_t    numele_sfh_dims[2]={NUM_ELEMENTS,SFH_NBIN};
hsize_t    mag_sfh_dims[2]={3,SFH_NBIN};
#endif //INDIVIDUAL_ELEMENTS
#endif //STAR_FORMATION_HISTORY
