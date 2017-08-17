#ifdef DETAILED_METALS_AND_MASS_RETURN
struct metals
{
  float type1a;
  float type2;
  float agb;
};

#ifdef INDIVIDUAL_ELEMENTS
//Number of chemical elements tracked:
#ifndef MAINELEMENTS
#define NUM_ELEMENTS 11 //All: [H][He][C][N][O][Ne][Mg][Si][S][Ca][Fe]
#else
#define NUM_ELEMENTS 5 //Only [H][He][O][Mg][Fe]
#endif
#endif //INDIVIDUAL_ELEMENTS
#endif //DETAILED_METALS_AND_MASS_RETURN
