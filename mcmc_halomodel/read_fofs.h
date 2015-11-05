
#ifndef READFOFS_H
#define READFOFS_H

typedef struct Fof_ {
  long long fofid;
  float fofmass;
  int ngal;
} Fof;

Fof* fof;
int Nfof;

char buf[512];

long get_number_of_fofs(char *filename);
void read_fof(char *filename);

#endif
