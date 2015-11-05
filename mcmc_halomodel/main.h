
#ifndef MAIN_H
#define MAIN_H

typedef struct index_struct_ {
  int my_index;
  double random;
} index_struct;

unsigned long int get_seed();
int compare_index(const void *a, const void *b);

#endif
