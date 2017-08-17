// TODO add documentation, also for all fields
#ifndef MCMC
extern struct  halo_ids_data
{
 long long HaloID;
 long long FileTreeNr;
 long long FirstProgenitor;
 long long LastProgenitor;
 long long NextProgenitor;
 long long Descendant;
 long long FirstHaloInFOFgroup;
 long long NextHaloInFOFgroup;
#ifdef MRII
 long long MainLeafID; 
#endif
 float    Redshift;
 int       PeanoKey;
 int       dummy;      /* need to use this padding for 64bit alignment */
} *HaloIDs, *HaloIDs_Data;
#else
extern struct  halo_ids_data
{
 long long FirstHaloInFOFgroup;
} *HaloIDs, *HaloIDs_Data;
#endif
