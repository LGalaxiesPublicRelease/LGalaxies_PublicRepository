#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "read_gals.h"

long get_number_of_galaxies(char *filename) {
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
} //get_number_of_galaxies

void read_gal(char *filename) {
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
  fread(&Ngal,sizeof(int),1,fd);
  gal=malloc(Ngal*sizeof(Gal));
  fread(gal,sizeof(Gal),Ngal,fd);
  fclose(fd);
} //read_gal
