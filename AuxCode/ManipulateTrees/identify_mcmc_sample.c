//gcc -o identify_mcmc_sample.exe identify_mcmc_sample.c -lm

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#include "structures.h"

int main()
{
  int ii, jj, treenr; 
  int Ntrees, *TreeNHalos; 
  int totNHalos;
  int *Ntrees_PerSmallBox, *totNHalos_PerSmallBox, **TreeNHalos_SmallBox;

  int LastDarkMatterSnapShot=107;
  int Nfiles=64, CurrentFile; 
  float FullBoxSize=125., SmallBoxSize;
 
  char buf[1000];
  char SimulationDir[1000];
 
  FILE *f_fullbox;
  FILE *fdbids_fullbox;
  FILE **f_smallbox;
  FILE **fdbids_smallbox;
  int i1,i2,i3;


  ///home/nifty2014/SemiAnalytics/nIFTy2014/SAMs/LGALAXY/Codes/code/convertSUSSING2LGAL/treedata/

  sprintf(SimulationDir, "/galformod/scratch/bmh20/Workspace/MergerTrees/Nifty/convertSUSSING2LGAL/");
 
  Ntrees_PerSmallBox = malloc(sizeof(int) * Nfiles);
  totNHalos_PerSmallBox = malloc(sizeof(int) * Nfiles);
  TreeNHalos_SmallBox = malloc(sizeof(int*) * Nfiles);

  for(ii=0;ii<Nfiles;ii++)
    {
      Ntrees_PerSmallBox[ii] = 0;
      totNHalos_PerSmallBox[ii] = 0;  
    }

  //if there are 512 new small box, each side will be divided 512^1./3 times 
   SmallBoxSize=FullBoxSize/pow(Nfiles,1./3.);



  /* First read everything and figure out which trees will
   * belong to each file to create headers for individual files
   * Ntrees_SmallBox, totNHalos_SmallBox, TreeNHalos_SmallBox */

  sprintf(buf, "%s/single_tree//treedata/trees_%03d.0", SimulationDir, LastDarkMatterSnapShot); 
  if(!(f_fullbox = fopen(buf, "rb")))
    {
      printf("can't open file '%s' in io_tree.c %d\n", buf, __LINE__);
      exit(1);
    }

  fread(&Ntrees, sizeof(int), 1, f_fullbox);
  fread(&totNHalos, sizeof(int), 1, f_fullbox);
  TreeNHalos = malloc(sizeof(int) * Ntrees); 
  fread(TreeNHalos, sizeof(int), Ntrees, f_fullbox);
 
  for(ii=0;ii<Nfiles;ii++)
    TreeNHalos_SmallBox[ii] = malloc(sizeof(int) * Ntrees);

  printf("Ntrees=%d totNHalos=%d\n",Ntrees,totNHalos);

  for (treenr = 0; treenr < Ntrees; treenr++)
    {
 
      Halo = malloc(sizeof(struct halo_data) * TreeNHalos[treenr]);
      fread(Halo, sizeof(struct halo_data), TreeNHalos[treenr], f_fullbox);

      for(jj=0;jj<TreeNHalos[treenr];jj++)
	{
	if(Halo[jj].FirstHaloInFOFgroup==jj && Halo[jj].SnapNum==LastDarkMatterSnapShot)
	  {
	    CurrentFile = (float)((int)(Halo[jj].Pos[0]/SmallBoxSize)*pow(Nfiles,2./3.) +
	                          (int)(Halo[jj].Pos[1]/SmallBoxSize)*pow(Nfiles,1./3.) + 
				  (int)(Halo[jj].Pos[2])/SmallBoxSize);		  
	    Ntrees_PerSmallBox[CurrentFile]++;	  
	    totNHalos_PerSmallBox[CurrentFile]+=TreeNHalos[treenr];
	    TreeNHalos_SmallBox[CurrentFile][treenr]=TreeNHalos[treenr];
	    break;
	  }
	}
      free(Halo);
    }
  free(TreeNHalos);
  fclose(f_fullbox);



  /* At this point the Ntrees_SmallBox, totNHalos_SmallBox, TreeNHalos_SmallBox
   * are defined. Now the code opens all the new smaller files, writes the headers
   * and then writes one tree at the time (into the appropriate file) as they 
   * are read again */

  for(ii=0;ii<Nfiles;ii++)
    printf("file[%d] Ntrees=%ld totNHalos=%ld\n",ii,Ntrees_PerSmallBox[ii],totNHalos_PerSmallBox[ii]);
 

  //open files to contain small boxes
  f_smallbox = malloc (Nfiles * sizeof(FILE*));

  for (ii=0; ii<Nfiles; ++ii)
    {     
      sprintf(buf, "%s/treedata/trees_%03d.%d", SimulationDir, LastDarkMatterSnapShot,ii); 
      if(!(f_smallbox[ii] = fopen(buf, "wb")))
	{
	  printf("can't open file '%s' in io_tree.c %d\n", buf, __LINE__);
	  exit(1);
	}
      fwrite(&Ntrees_PerSmallBox[ii], sizeof(int), 1, f_smallbox[ii]);
      fwrite(&totNHalos_PerSmallBox[ii], sizeof(int), 1, f_smallbox[ii]);

      for(jj=0;jj<Ntrees_PerSmallBox[ii];jj++)
	fwrite(&TreeNHalos_SmallBox[ii][jj], sizeof(int), 1, f_smallbox[ii]);
    }
  
 //open files to contain small dbids boxes
  fdbids_smallbox = malloc (Nfiles * sizeof(FILE*));
  for (ii=0; ii<Nfiles; ++ii)
    {     
      sprintf(buf, "%s/treedata/tree_dbids_%03d.%d", SimulationDir, LastDarkMatterSnapShot,ii); 
      if(!(fdbids_smallbox[ii] = fopen(buf, "wb")))
	{
	  printf("can't open file '%s' in io_tree.c %d\n", buf, __LINE__);
	  exit(1);
	}  
    }
   
  printf("\nHeaders Written\n\n");


  //read trees again and write
  sprintf(buf, "%s/single_tree//treedata/trees_%03d.0", SimulationDir, LastDarkMatterSnapShot); 
  if(!(f_fullbox = fopen(buf, "rb")))
    {
      printf("can't open file '%s' in io_tree.c %d\n", buf, __LINE__);
      exit(1);
    }

  sprintf(buf, "%s/single_tree//treedata/tree_dbids_%03d.0", SimulationDir, LastDarkMatterSnapShot); 
  if(!(fdbids_fullbox = fopen(buf, "rb")))
    {
      printf("can't open file '%s' in io_tree.c %d\n", buf, __LINE__);
      exit(1);
    }

  fread(&Ntrees, sizeof(int), 1, f_fullbox);
  fread(&totNHalos, sizeof(int), 1, f_fullbox);
  TreeNHalos = malloc(sizeof(int) * Ntrees); 
  fread(TreeNHalos, sizeof(int), Ntrees, f_fullbox);
 
  for (treenr = 0; treenr < Ntrees; treenr++)
    {
 
      Halo = malloc(sizeof(struct halo_data) * TreeNHalos[treenr]);
      fread(Halo, sizeof(struct halo_data), TreeNHalos[treenr], f_fullbox);

      HaloIDs = malloc(sizeof(struct halo_ids_data) * TreeNHalos[treenr]);
      fread(HaloIDs, sizeof(struct halo_ids_data), TreeNHalos[treenr], fdbids_fullbox);

      for(jj=0;jj<TreeNHalos[treenr];jj++)
	{
	if(Halo[jj].FirstHaloInFOFgroup==jj && Halo[jj].SnapNum==LastDarkMatterSnapShot)
	  {
	    CurrentFile = (float)((int)(Halo[jj].Pos[0]/SmallBoxSize)*pow(Nfiles,2./3.) +
	                          (int)(Halo[jj].Pos[1]/SmallBoxSize)*pow(Nfiles,1./3.) + 
				  (int)(Halo[jj].Pos[2])/SmallBoxSize);	 
	    fwrite(Halo, sizeof(struct halo_data), TreeNHalos[treenr], f_smallbox[CurrentFile]);
	    fwrite(HaloIDs, sizeof(struct halo_ids_data), TreeNHalos[treenr], fdbids_smallbox[CurrentFile]);
	    break;
	  }
	}
      free(Halo);
    }
  free(TreeNHalos);
  fclose(f_fullbox);
  fclose(fdbids_fullbox);

  //close everything
  for (ii=0; ii<Nfiles; ++ii)
    {
      fclose(f_smallbox[ii]);
      fclose(fdbids_smallbox[ii]);
      free(TreeNHalos_SmallBox[ii]);
    }

  free(Ntrees_PerSmallBox); 
  free(totNHalos_PerSmallBox); 

  printf("done\n");

}
