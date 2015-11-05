#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "read_fofs.h"

long get_number_of_fofs(char *filename) {
  int np;
  char buf[255];
  FILE *fd;
  fd=fopen(filename,"r");
  if (!fd) {
    sprintf(buf,"%s.dat",filename);
    fd=fopen(buf,"r");
  } //if
  if (!fd) {
    printf("Problem opening %s.\n",filename);
    exit(0);
  } //if
  fread(&np,sizeof(int),1,fd);
  fclose(fd);
  return np;
} //get_number_of_fofs

void read_fof(char *filename) {
  FILE *fd;
  fd = fopen(filename,"r");
  if (!fd) {
    sprintf(buf,"%s.dat",filename);
    fd=fopen(buf,"r");
  } //if
  if (!fd) {
    printf("Problem opening %s.\n",filename);
    exit(0);
  } //if
  fread(&Nfof,sizeof(int),1,fd);
  fof=malloc(Nfof*sizeof(Fof));
  fread(fof,sizeof(Fof),Nfof,fd);
  fclose(fd);
} //read_gal
