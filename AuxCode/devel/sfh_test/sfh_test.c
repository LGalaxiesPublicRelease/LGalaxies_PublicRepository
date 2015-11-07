/* Code to test the bin creation and merging in star_formation_history.c */

#include "allvars.h"
#include "proto.h"

int main(int argc, char **argv) {

  int i;

  //Initialise
  sfh_initialise(0);
  sfh_print(0);

  sfh_update_bins(0,0.3e7);
  Gal[0].sfh_StellarMass[Gal[0].sfh_ibin]+=1.;
  sfh_print(0);

  sfh_update_bins(0,1.3e7);
  Gal[0].sfh_StellarMass[Gal[0].sfh_ibin]+=2.;
  sfh_print(0);

  sfh_update_bins(0,5.0e7);
  Gal[0].sfh_StellarMass[Gal[0].sfh_ibin]+=3.;
  sfh_print(0);

  /* What happens if we go back in time? */
  sfh_update_bins(0,1.3e7);
  Gal[0].sfh_StellarMass[Gal[0].sfh_ibin]+=2.;
  sfh_print(0);
  /* Answer: correctly leaves bins unchanged,
   * but obviously adds mass into wrong bin.
   * I could add a function to add mass in at
   * the correct time, but this seems unnecessary,
   * so just add a usage comment into the
   * star_formation_history file instead.
   */

  /* What happens if we add in a very large number
   * of bins, exceeding our bin limit? */
  sfh_update_bins(0,32.1e7);
  Gal[0].sfh_StellarMass[Gal[0].sfh_ibin]+=4.;
  sfh_print(0);
  /* It works! */

  /* Create a new galaxy and add it to the first */
  sfh_initialise(1);
  sfh_update_bins(1,10.0e7);
  Gal[1].sfh_StellarMass[Gal[1].sfh_ibin]+=4.;
  sfh_update_bins(1,25.0e7);
  Gal[1].sfh_StellarMass[Gal[1].sfh_ibin]+=4.;
  sfh_update_bins(1,32.1e7);
  Gal[1].sfh_StellarMass[Gal[1].sfh_ibin]+=4.;
  sfh_print(1);
  sfh_merge(0,1);
  sfh_print(0);
  sfh_print(1);

  /* Test calculation of ages to write */
  for(i=Gal[0].sfh_ibin;i>0;i--)
    printf("%5d,%12f,%12f\n",Gal[0].sfh_t[i],
	   32.1e7-Gal[0].sfh_t[i-1]*SFH_TIME_INTERVAL,
	   Gal[0].sfh_StellarMass[i]);
  printf("%5d,%12f,%12f\n",Gal[0].sfh_t[0],
	 32.1e7,Gal[0].sfh_StellarMass[0]);

  /* And again, using sfh_age */
  for(i=Gal[0].sfh_ibin;i>0;i--)
    printf("%5d,%12f,%12f\n",Gal[0].sfh_t[i],
	   Gal[0].sfh_age-Gal[0].sfh_t[i-1]*SFH_TIME_INTERVAL,
	   Gal[0].sfh_StellarMass[i]);
  printf("%5d,%12f,%12f\n",Gal[0].sfh_t[0],
	 Gal[0].sfh_age,Gal[0].sfh_StellarMass[0]);

  return(0);

}

void sfh_print(int p) {
  /* prints out populated sfh_structure */
  int i;

  printf("For galaxy %d:\n",p);
  printf("sfh_init=%d\n",Gal[p].sfh_ibin);
  printf("  i    dt   t      Stars      Metals\n");
  for(i=0;i<SFH_NBIN;i++)
    if (Gal[p].sfh_dt[i]!=0)
      printf("%5d%5d%5d%12f%12f\n",i,Gal[p].sfh_dt[i],Gal[p].sfh_t[i],
	     Gal[p].sfh_StellarMass[i],Gal[p].sfh_MetalsStellarMass[i]);

}
