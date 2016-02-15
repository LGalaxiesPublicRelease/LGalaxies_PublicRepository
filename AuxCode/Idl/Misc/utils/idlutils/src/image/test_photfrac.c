#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include "ph.h"

/*
 */

#define USAGE \
		{ fprintf(stderr,"Usage: cat <galaxy file> | fit_coeffs [--vfile <vfile> --lfile <lfile>\n"); \
		  fprintf(stderr,"            --ffile <ffile> ]\n"); }

int main(int argc,
				 char **argv)
{
  int i,j,xnpix,ynpix; 
  long c;
  float *frac,radius;

	/* read arguments */
  radius=3.75;
	i=0;
	while(1) {
		int option_index = 0;
 		static struct option long_options[] =
			{
				{"radius", 1, 0, 0}, 
			};
		static const char short_options[]="r:";
		static const char short_options_c[]="r";

		c=getopt_long(argc, argv, short_options, long_options, &option_index);
		if(c==-1) break;
		if(c==0) c=short_options_c[option_index];
		switch(c) {
		case 'r':
      radius=atof(optarg);
			break; 
		case '?':
			break;
		default: 
			printf("test_photfrac: getopt returned character code 0%o ??\n", 
             (unsigned int) c);
		}
		i++;
	}
	if(argc<0) {
    USAGE;
		exit(1);
	} /* end if */
#if 0
	/* read arguments */
	strcpy(vfile,"vmatrix.default.dat");
	strcpy(lfile,"lambda.default.dat");
	strcpy(ffile,"sdss_filters.dat");
	sprintf(path,"%s/data/templates",getenv("KCORRECT_DIR"));
	i=0;
	while(1) {
		int this_option_optind = optind ? optind : 1;
		int option_index = 0;
 		static struct option long_options[] =
			{
				{"vfile", 1, 0, 0}, 
				{"lfile", 1, 0, 0},
				{"path", 1, 0, 0},
				{"ffile", 1, 0, 0},
				{"help", 0, 0, 0}
			};
		static const char short_options[]="v:l:p:f:h";
		static const char short_options_c[]="vlpfh";

		c=getopt_long(argc, argv, short_options, long_options, &option_index);
		if(c==-1) break;
		if(c==0) c=short_options_c[option_index];
		switch(c) {
		case 'v':
			strcpy(vfile,optarg);
			break; 
		case 'l':
			strcpy(lfile,optarg);
			break; 
		case 'p':
			strcpy(path,optarg);
			break; 
		case 'f':
			strcpy(ffile,optarg);
			break; 
		case 'h':
      USAGE;
      exit(1);
			break; 
		case '?':
			break;
		default: 
			printf("fit_coeffs: getopt returned character code 0%o ??\n", c);
		}
		i++;
	}
	if(argc<0) {
    USAGE;
		exit(1);
	} /* end if */
#endif

  xnpix=9;
  ynpix=9;
  frac=(float *) malloc(xnpix*ynpix*sizeof(float));
  photfrac(xnpix,ynpix,radius,frac,0,0);
  for(i=0;i<xnpix;i++) 
    for(j=0;j<ynpix;j++) 
      printf("%d %d %e\n", i, j, frac[j*xnpix+i]);

	return(0);
} /* end main */
