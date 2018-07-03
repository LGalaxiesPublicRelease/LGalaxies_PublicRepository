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
struct elements_str
{
    float H;
    float He;
    float C;
    float N;
    float O;
    float Ne;
    float Mg;
    float Si;
    float S;
    float Ca;
    float Fe;
};
#else
#define NUM_ELEMENTS 5 //Only [H][He][O][Mg][Fe]
struct elements_str
{
    float H;
    float He;
    float O;
    float Mg;
    float Fe;
};
#endif  //MAINELEMENTS
/* Define two views into the same element structure/array.
   For example define:
      union elements DiskMassElements;
   Then access as either
      DiskMassElements.str.H
   or
      DiskMassElements.arr[0]
*/
union elements
{
    struct elements_str str;
    float arr[NUM_ELEMENTS];
};
#endif //INDIVIDUAL_ELEMENTS
#endif //DETAILED_METALS_AND_MASS_RETURN
