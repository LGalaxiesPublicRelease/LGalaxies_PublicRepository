/*  Copyright (C) <2016>  <L-Galaxies>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"

#ifdef PARALLEL
#ifdef MCMC
#include <mpi.h>
#endif
#endif

/**@file io_tree.c
 * @brief Reads in the data from the dark matter simulation merger
 *        trees, creates output files and after calculations frees
 *        the allocated memory.
 *
 * There are three different input files: trees_** - normal tree
 * files; tree_dbids - file containing halo IDs (used if GALAXYTREE
 * ON); tree_aux - file containing the particle data (used if
 * UPDATETYPETWO ON).
 */

/**@brief Reads all the tree files if USE_MEMORY_TO_MINIMIZE_IO ON,
 *        otherwise reads in the headers for trees_** and trees_aux;
 *        Opens output files.
 *
 *  If USE_MEMORY_TO_MINIMIZE_IO ON the first step on this file is to
 *  call load_all_dbids(), load_all_auxdata(), load_all_treedata().
 *  These will read all the tree data for the current file into pointers
 *  instead of doing it once for every tree.
 *
 *  If USE_MEMORY_TO_MINIMIZE_IO OFF the trees are read independently
 *  from the files.
 *
 *  Modified versions of myfread/myfwrite/myfseek are called to either
 *  read individual trees from the files into structures or from the
 *  pointers with all the data for the current file into structures.
 *
 *  For each tree the code reads in the header in trees_**: Ntrees -
 *  number of trees in the current file (int); totNHalos - total number
 *  of halos summed over all the trees in current file (int);
 *  TreeNHalos - number of halos in each tree (Ntrees).
 *
 *  Then output files are opened SA_z**_** - for snapshot output;
 *  SA_galtree_** for GALAXYTREE option and SA_**.** for MOMAF.
 *
 *  If UPDATETYPETWO ON the header in tree_aux is also read: NtotHalos -
 *  total number of halos in the file (int); TotIds - total number of
 *  particle IDs (int); Ntrees - total number of trees (int); TotSnaps -
 *  total number of snapshots (int). Define some other quantities:
 *  CountIDs_snaptree - number of Ids for each tree at each snapshot
 *  (TotSnaps * Ntrees); OffsetIDs_snaptree (TotSnaps * Ntrees);
 *  CountIDs_halo - Number of Ids per halo (NtotHalos); OffsetIDs_halo
 *  (int). */
void load_tree_table(int filenr)
{
  int i, n, totNHalos, SnapShotInFileName;
  char buf[1000];


#ifdef  UPDATETYPETWO
  load_all_auxdata(filenr);
#endif

#ifdef PARALLEL
#ifndef MCMC
  printf("\nTask %d reading tree file %d\n", ThisTask, filenr);
#else
  //if def MCMC and PARALLEL only task 0 reads the representative treefile, then broadcasts
  if(ThisTask==0)
  {
  	 printf("\nTask %d reading trees_%d\n", ThisTask, filenr);
#endif
#endif

  SnapShotInFileName=LastDarkMatterSnapShot;

#ifdef MCMC
#ifdef MR_PLUS_MRII
  SnapShotInFileName=LastDarkMatterSnapShot_MRII;
#endif
#endif

#ifdef LOADIDS
#ifndef MRII
  sprintf(buf, "%s/treedata/tree_dbids_%03d.%d", SimulationDir, SnapShotInFileName, filenr);
#else
  sprintf(buf, "%s/treedata/tree_sf1_dbids_%03d.%d", SimulationDir, SnapShotInFileName, filenr);
#endif
  if(!(treedbids_file = fopen(buf, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "can't open file `%s'\n", buf);
      terminate(sbuf);
    }
#endif


#ifndef MRII
  sprintf(buf, "%s/treedata/trees_%03d.%d", SimulationDir, SnapShotInFileName, filenr);
#else
  sprintf(buf, "%s/treedata/trees_sf1_%03d.%d", SimulationDir, SnapShotInFileName, filenr);
#endif

  if(!(tree_file = fopen(buf, "r")))
    {
      char sbuf[1000];
      sprintf(sbuf, "can't open file place `%s'\n", buf);
      terminate(sbuf);
    }

  //read header on trees_** file
  myfread(&Ntrees, 1, sizeof(int), tree_file);
  myfread(&totNHalos, 1, sizeof(int), tree_file);

  TreeNHalos = mymalloc("TreeNHalos", sizeof(int) * Ntrees);
  TreeFirstHalo = mymalloc("TreeFirstHalo", sizeof(int) * Ntrees);
  TreeNgals[0] = mymalloc("TreeNgals[n]", NOUT * sizeof(int) * Ntrees);
  for(n = 1; n < NOUT; n++)
    TreeNgals[n] = TreeNgals[n - 1] + Ntrees;



  myfread(TreeNHalos, Ntrees, sizeof(int), tree_file);

  if(Ntrees)
    TreeFirstHalo[0] = 0;
  /*Define a variable containing the number you have to jump to
   * get from one firshalo to the next. */
  for(i = 1; i < Ntrees; i++)
    TreeFirstHalo[i] = TreeFirstHalo[i - 1] + TreeNHalos[i - 1];

#ifdef PRELOAD_TREES
  Halo_Data = mymalloc("Halo_Data", sizeof(struct halo_data) * totNHalos);
  myfseek(tree_file, sizeof(int) * (2 + Ntrees), SEEK_SET);
  myfread(Halo_Data, totNHalos, sizeof(struct halo_data), tree_file);
#ifdef PARALLEL
  printf("\nTask %d done loading trees_%d\n", ThisTask, filenr);
#endif

#ifdef LOADIDS
  HaloIDs_Data = mymalloc("HaloIDs_Data", sizeof(struct halo_ids_data) * totNHalos);
  myfseek(treedbids_file, 0, SEEK_SET);
  myfread(HaloIDs_Data, totNHalos, sizeof(struct halo_ids_data), treedbids_file);
#ifdef PARALLEL
  printf("\nTask %d done loading tree_dbids_%d\n", ThisTask, filenr);
#endif
#endif
#endif

  //if MCMC is turned only Task 0 reads the file and then broadcasts
#ifdef PARALLEL
#ifdef MCMC
  } // end if ThisTask==0

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Bcast(&Ntrees,sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&totNHalos,sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(ThisTask>0)
  	TreeNHalos = mymalloc("TreeNHalos", sizeof(int) * Ntrees);

  MPI_Bcast(TreeNHalos,sizeof(int)*Ntrees, MPI_BYTE, 0, MPI_COMM_WORLD);

  if(ThisTask>0)
  {
  	TreeNgals[0] = mymalloc("TreeNgals[n]", NOUT * sizeof(int) * Ntrees);
  	for(n = 1; n < NOUT; n++)
  		TreeNgals[n] = TreeNgals[n - 1] + Ntrees;

  	TreeFirstHalo = mymalloc("TreeFirstHalo", sizeof(int) * Ntrees);
  	if(Ntrees)
  		TreeFirstHalo[0] = 0;
    /*Define a variable containing the number you have to jump to
     * get from one firshalo to the next. */
  	for(i = 1; i < Ntrees; i++)
  		TreeFirstHalo[i] = TreeFirstHalo[i - 1] + TreeNHalos[i - 1];

  	Halo_Data = mymalloc("Halo_Data", sizeof(struct halo_data) * totNHalos);
  	HaloIDs_Data = mymalloc("HaloIDs_Data", sizeof(struct halo_ids_data) * totNHalos);
  }

  MPI_Bcast(HaloIDs_Data,totNHalos* sizeof(struct halo_ids_data), MPI_BYTE, 0, MPI_COMM_WORLD);

  size_t bytes=totNHalos* sizeof(struct halo_data);
  int ii, Nmessages=10000;
  int HaloChunks=1000000;

  //MPI_BCast has a limit of 2Gb so everything needs to be passed in smaller chunks
  for(ii=0;ii<Nmessages;ii++)
  {
  	//if next chunk is outside of array size, just pass whats left and then exit the loop
   	if((ii+1)*HaloChunks>totNHalos)
   	{
   		MPI_Bcast(&Halo_Data[ii*HaloChunks],bytes-ii*HaloChunks* sizeof(struct halo_data), MPI_BYTE, 0, MPI_COMM_WORLD);
   		break;
   	}
   	else
   		MPI_Bcast(&Halo_Data[ii*HaloChunks],HaloChunks* sizeof(struct halo_data), MPI_BYTE, 0, MPI_COMM_WORLD);
   }


  if(ThisTask==0)
  	printf("all tree data has now been broadcasted\n");
#endif
#endif

}




/**@brief Deallocates the arrays used to read in the headers of
 *        trees_** and tree_aux (the ones with more than one
 *        element); if PRELOAD_TREES ON, deallocates
 *        the pointers containing all the input and output data.*/
void free_tree_table(void)
{
#ifdef PRELOAD_TREES
#ifdef LOADIDS
  myfree(HaloIDs_Data);
#endif
  myfree(Halo_Data);
#endif

  myfree(TreeNgals[0]);

  //deallocates header from trees_**
  myfree(TreeFirstHalo);//derived from the header of trees_**
  myfree(TreeNHalos);

#ifdef UPDATETYPETWO
  myfree(TreeAuxData);
#endif

#ifdef LOADIDS
  fclose(treedbids_file);
#endif

  fclose(tree_file);
}


/**@brief Reads the actual trees into structures to be used in the
 *        code: Halo and HaloIDs; the galaxy structures (HaloGal
 *        and Gal) and the HaloAux are also allocated
 *
 *  If USE_MEMORY_TO_MINIMIZE_IO & NEW_IO are OFF, the trees_** files
 *  are opened every time a tree needs to be read in. Then the code's
 *  structure that will have the tree information is read: Halo =
 *  (sizeof(struct halo_data) * TreeNHalos[]) are read.
 *
 *  HaloAux structure is allocated =
 *  (sizeof(struct halo_aux_data) * TreeNHalos[])
 *
 *  Considering the number of halos allocated for the current tree
 *  the size of the structures containing the galaxies
 *  with a halo and the galaxies with and without a halo to be
 *  allocated is estimated:
 *
 *  For galaxies with a halo - MaxGals = MAXGALFAC * TreeNHalos[] and
 *  HaloGal = (sizeof(struct GALAXY) * MaxGals).
 *
 *  For all galaxies - FoF_MaxGals = 10000*15 and
 *  Gal = (sizeof(struct GALAXY) * FoF_MaxGals)
 *
 *  If GALAXYTREE ON, HaloIDs structure is read from tree_dbids =
 *  sizeof(struct halo_ids_data) * TreeNHalos[] */

void load_tree(int nr)
{
  int i;

#ifdef PRELOAD_TREES
  Halo = Halo_Data + TreeFirstHalo[nr];
  /*for(i=0;i<TreeNHalos[nr];i++)
  	printf("vel=%f\n",Halo[i].Vel[1]);*/
#ifdef LOADIDS
  HaloIDs = HaloIDs_Data + TreeFirstHalo[nr];
#endif
#else

  Halo = mymalloc("Halo", sizeof(struct halo_data) * TreeNHalos[nr]);
  myfseek(tree_file, sizeof(int) * (2 + Ntrees) + sizeof(struct halo_data) * TreeFirstHalo[nr], SEEK_SET);
  myfread(Halo, TreeNHalos[nr], sizeof(struct halo_data), tree_file);
#ifdef LOADIDS
  HaloIDs = mymalloc("HaloIDs", sizeof(struct halo_ids_data) * TreeNHalos[nr]);
  myfseek(treedbids_file, sizeof(struct halo_ids_data) * TreeFirstHalo[nr], SEEK_SET);
  myfread(HaloIDs, TreeNHalos[nr], sizeof(struct halo_ids_data), treedbids_file);
#endif

#endif

  //Allocate HaloAux and Galaxy structures.
  HaloAux = mymalloc("HaloAux", sizeof(struct halo_aux_data) * TreeNHalos[nr]);

  for(i = 0; i < TreeNHalos[nr]; i++)
    {
      HaloAux[i].DoneFlag = 0;
      HaloAux[i].HaloFlag = 0;
      HaloAux[i].NGalaxies = 0;
    }

  if(AllocValue_MaxHaloGal == 0)
    AllocValue_MaxHaloGal = 1 + TreeNHalos[nr] / (0.25 * (LastDarkMatterSnapShot+1));

  if(AllocValue_MaxGal == 0)
    AllocValue_MaxGal = 2000;

  MaxHaloGal = AllocValue_MaxHaloGal;
  NHaloGal = 0;
  HaloGal = mymalloc_movable(&HaloGal, "HaloGal", sizeof(struct GALAXY) * MaxHaloGal);
  HaloGalHeap = mymalloc_movable(&HaloGalHeap, "HaloGalHeap", sizeof(int) * MaxHaloGal);
  for(i = 0; i < MaxHaloGal; i++)
    HaloGalHeap[i] = i;

  MaxGal = AllocValue_MaxGal;
  Gal = mymalloc_movable(&Gal, "Gal", sizeof(struct GALAXY) * MaxGal);

#ifdef GALAXYTREE
  if(AllocValue_MaxGalTree == 0)
    AllocValue_MaxGalTree = 1.5 * TreeNHalos[nr];

  MaxGalTree = AllocValue_MaxGalTree;
  GalTree = mymalloc_movable(&GalTree, "GalTree", sizeof(struct galaxy_tree_data) * MaxGalTree);
#endif
}



/**@brief Frees all the Halo and Galaxy structures in the code. */
void free_galaxies_and_tree(void)
{
#ifdef GALAXYTREE
  myfree(GalTree);
#endif
  myfree(Gal);
  myfree(HaloGalHeap);
  myfree(HaloGal);
  myfree(HaloAux);

#ifndef PRELOAD_TREES
#ifdef LOADIDS
  myfree(HaloIDs);
#endif
  myfree(Halo);
#endif
}


/**@brief Reading routine, either from a file into a structure or
 *        from a pointer to a structure.
 *   */
size_t myfread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if(size * nmemb > 0)
    {
      if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
	{
	  if(feof(stream))
	    printf("I/O error (fread) has occured: end of file\n");
	  else
	    printf("I/O error (fread) has occured: %s\n", strerror(errno));
	  fflush(stdout);
	  terminate("read error");
	}
    }
  else
    nread = 0;

  return nread;
}

size_t myfwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if(size * nmemb > 0)
    {
      if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
	{
	  printf("I/O error (fwrite) has occured: %s\n", strerror(errno));
	  fflush(stdout);
	  terminate("write error");
	}
    }
  else
    nwritten = 0;

  return nwritten;
}

int myfseek(FILE * stream, long offset, int whence)
{
  if(fseek(stream, offset, whence))
    {
      printf("I/O error (fseek) has occured: %s\n", strerror(errno));
      fflush(stdout);
      terminate("write error");
    }

  return 0;
}
