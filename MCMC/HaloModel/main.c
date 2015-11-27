
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_deriv.h>

#include "main.h"
#include "allvars.h"
#include "halomodel.h"
#include "read_gals.h"
#include "read_fofs.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Status status;
  int initsize=1000;
  int stepsize=1000;
  double maxerr=0.1;
  char buf[500];
  FILE *fd;
  int i,j,k,maxfof,allowed,allowed_tot,bestk,diff,iter;
  float deltam,mingalmass,maxgalmass;
  double toterr,toterr_tmp,smallesterr;
  const gsl_rng_type *T;
  gsl_rng *random_generator;
  int* usedbymass;
  index_struct* indexing_array;
  double *r,*proj,*r_tmp,*proj_tmp;
  struct timespec start,end,tmp1,tmp2;
/* if (ThisTask==0 && NTasks%6!=0) {
    printf("NTasks must be a multiple of 6.\n");
    exit(0);
  } //if
  if (ThisTask==0) clock_gettime(CLOCK_MONOTONIC,&start);*/
  NR=60;
  massbins=60;
  minfofmass=10.;
  maxfofmass=16.;
  mingalmass=8.77+(ThisTask%6)*0.5;
  maxgalmass=8.77+(ThisTask%6+1)*0.5;
  sprintf(buf,"Galtrimmed.dat");
  read_gal(buf);
  //gal[i].
  sprintf(buf,"Foftrimmed.dat");
  read_fof(buf);
  //fof[i].
 // deltam=(maxfofmass-minfofmass)/(double)massbins;
 // numbymass=malloc(massbins*sizeof(int));
 /* for (i=0; i<massbins; ++i) numbymass[i]=0;
  for (i=0; i<Nfof; ++i) {
    j=floor((fof[i].fofmass-minfofmass)/deltam);
    if (j>=0 && j<massbins) numbymass[j]++;
  } //for
  indexbymass=malloc(massbins*sizeof(int*));
  for (i=0; i<massbins; ++i) {
    indexbymass[i]=malloc(numbymass[i]*sizeof(int));
    numbymass[i]=0;
  } //for
  for (i=0; i<Nfof; ++i) {
    j=floor((fof[i].fofmass-minfofmass)/deltam);
    if (j>=0 && j<massbins) {
      indexbymass[j][numbymass[j]]=i;
      numbymass[j]++;
    } //if
  } //for*/
 // maxfof=0;
 /* for (i=0; i<Nfof; ++i) if (fof[i].fofid>maxfof) maxfof=fof[i].fofid;
  FofidTable=malloc((maxfof+1)*sizeof(int));
  for (i=0; i<=maxfof; ++i) FofidTable[i]=-1;
  for (i=0; i<Ngal; i+=gal[i].ngal) FofidTable[gal[i].fofid]=i;
  gsl_rng_env_setup();
  T=gsl_rng_default;
  random_generator=gsl_rng_alloc(T);
  gsl_rng_set(random_generator,get_seed());
  for (i=0; i<massbins; ++i) {
    indexing_array=malloc(numbymass[i]*sizeof(index_struct));
    for (j=0; j<numbymass[i]; ++j) {
      indexing_array[j].my_index=indexbymass[i][j];
      indexing_array[j].random=gsl_rng_uniform(random_generator);
    } //for
    qsort(indexing_array,numbymass[i],sizeof(index_struct),compare_index);
    for (j=0; j<numbymass[i]; ++j) {
      indexbymass[i][j]=indexing_array[j].my_index;
    } //for
    free(indexing_array);
  } //for*/
  //usedbymass=malloc(massbins*sizeof(int));
  r=malloc(NR*sizeof(double));
  proj=malloc(NR*sizeof(double));
  //r_tmp=malloc(NR*sizeof(double));
  //proj_tmp=malloc(NR*sizeof(double));
  initialization();
  //offset=0;
  //clock_gettime(CLOCK_MONOTONIC,&tmp1);
  halomodel(r,proj,mingalmass,maxgalmass);
 /* clock_gettime(CLOCK_MONOTONIC,&tmp2);
  printf("Task %d, calculation 1 took %.2lf seconds.\n",ThisTask,(tmp2.tv_sec-tmp1.tv_sec)+(tmp2.tv_nsec-tmp1.tv_nsec)/1e9);
  if (argc>1) {
    sprintf(buf,argv[1]);
    fd=fopen(buf,"r");
    if (!fd) {
      printf("Input file could not be found!\n");
      exit(0);
    } //if
    for (i=0; i<massbins; ++i) {
      fscanf(fd,"%d",&initsize);
      usedbymass[i]=MIN(initsize,numbymass[i]);
    } //for
    fclose(fd);
  } //if
  else {
    for (i=0; i<massbins; ++i) usedbymass[i]=MIN(initsize,numbymass[i]);
  } //else
  //printf("Task %d going into sampled halomodel calculation...\n",ThisTask);
  MPI_Barrier(MPI_COMM_WORLD);
  clock_gettime(CLOCK_MONOTONIC,&tmp1);
  halomodel(r_tmp,proj_tmp,mingalmass,maxgalmass,usedbymass);
  clock_gettime(CLOCK_MONOTONIC,&tmp2);
  printf("Task %d, calculation 2 took %.2lf seconds.\n",ThisTask,(tmp2.tv_sec-tmp1.tv_sec)+(tmp2.tv_nsec-tmp1.tv_nsec)/1e9);
  //each processor now has the full projected correlation function for their mass bin stored
  //in proj, and a randomly sampled one in proj_tmp*/
  free(proj_tmp);
  free(r_tmp);
  free(proj);
  free(r);
  free(usedbymass);
  free(FofidTable);
  for (i=massbins-1; i>=0; --i) free(indexbymass[i]);
  free(indexbymass);
  free(numbymass);
  MPI_Finalize();
  return (1);
} //main

unsigned long int get_seed() {
  unsigned int seed;
  FILE *devrandom;
  devrandom=fopen("/dev/urandom","r");
  if (devrandom==NULL) {
    printf("Failed to open /dev/urandom.\n");
    exit(0);
  } //if
  fread(&seed,sizeof(seed),1,devrandom);
  fclose(devrandom);
  return(seed);
} //get_seed

int compare_index(const void *a, const void *b) {
  if (((index_struct *) a)->random < ((index_struct *) b)->random) return -1;
  if (((index_struct *) a)->random > ((index_struct *) b)->random) return +1;
  return 0;
} //compare_index
