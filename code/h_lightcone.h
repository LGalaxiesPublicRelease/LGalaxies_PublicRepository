#ifdef ALL_SKY_LIGHTCONE
extern struct dist_table
{
  double log_time, logD, Redshift;
}
*TimeTable, *DistanceTable;

struct GALAXY_OUTPUT_LIGHTCONE
{
  int   Type; // Galaxy type: 0 for central galaxies of a main halo, 1 for central galaxies in sub-halos, 2 for satellites without halo.
  int   SnapNum; // The snapshot number where this galaxy was identified.
  float CentralMvir; // 10^10/h Msun virial mass of background (FOF) halo containing this galaxy
  float Pos[3]; // 1/h Mpc - Galaxy Positions
  float BulgeMass; // 10^10/h Msun - Mass in the bulge
  float DiskMass;
  float Distance;
  float Time;
  float Redshift;
  float i_mag;
  /* magnitudes in various bands */
#ifdef OUTPUT_REST_MAGS
  float MagDust[NMAG]; // rest-frame absolute mags
#endif
#ifdef OUTPUT_OBS_MAGS
  float ObsMagDust[NMAG]; // rest-frame absolute mags
#endif
};
#endif
