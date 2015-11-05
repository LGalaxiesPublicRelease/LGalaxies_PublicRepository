
#ifndef READGALS_H
#define READGALS_H

typedef struct Gal_ {
  long long fofid;
  float fofmass;
  float x;
  float y;
  float z;
  float galmass;
  int ngal;
} Gal;

Gal* gal;
int Ngal;

char buf[512];

long get_number_of_galaxies(char *filename);
void read_gal(char *filename);

#endif
