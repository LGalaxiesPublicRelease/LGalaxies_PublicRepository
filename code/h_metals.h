#ifdef DETAILED_METALS_AND_MASS_RETURN
/*struct elements_str
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
/*union elements
{
    struct elements_str str;
    float arr[NUM_ELEMENTS];
};*/
#define NUM_METAL_CHANNELS 3

#ifdef INDIVIDUAL_ELEMENTS
//Number of chemical elements tracked:
#ifndef MAINELEMENTS
#define NUM_ELEMENTS 11 //All: [H][He][C][N][O][Ne][Mg][Si][S][Ca][Fe]
#else
#define NUM_ELEMENTS 5 //Only [H][He][O][Mg][Fe]
#endif
#endif //INDIVIDUAL_ELEMENTS

#else
#define NUM_METAL_CHANNELS 1
#endif //DETAILED_METALS_AND_MASS_RETURN
