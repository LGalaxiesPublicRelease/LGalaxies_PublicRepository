#include "export.h" 
float p_qromo(float (*func)(float), float a, float b,
              float (*choose)(float(*)(float), float, float, int));
void p_polint(float xa[], float ya[], int n, float x, float *y, 
							float *dy);
float p_midpnt(float (*func)(float), float a, float b, int n);
float p_midinf(float (*func)(float), float a, float b, int n);
void p_free_vector(float *v, long nl, long nh);
float *p_vector(long nl, long nh);
void p_cisi(float x, float *ci, float *si);
int photfrac(int xnpix, int ynpix, float radius, float *frac, long xcen, 
             long ycen);
