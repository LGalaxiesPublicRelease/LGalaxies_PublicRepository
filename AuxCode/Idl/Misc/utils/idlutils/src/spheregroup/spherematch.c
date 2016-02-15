#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "chunks.h"

static IDL_LONG *nra=NULL, ndec;
static double **rabounds=NULL, *decbounds=NULL;
static double raoffset;
static double *xx2=NULL,*yy2=NULL,*zz2=NULL;
static IDL_VPTR vxx2,vyy2,vzz2;
static IDL_LONG **nchunk2=NULL, ***chunklist2=NULL; 

#define DEG2RAD .01745329251994

double separation(double xx1, double yy1, double zz1, double xx2, double yy2,
									double zz2);

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
	IDL_Deltmp(vxx2); xx2=NULL;
	IDL_Deltmp(vyy2); yy2=NULL;
	IDL_Deltmp(vzz2); zz2=NULL;
	if(nchunk2!=NULL) 
		unassignchunks(&nchunk2,&chunklist2,nra,ndec);
	if(rabounds!=NULL)
		unsetchunks(&rabounds,&decbounds,&nra,&ndec);
	nchunk2=NULL;
	chunklist2=NULL;
	rabounds=NULL;
	decbounds=NULL;
	nra=NULL;
}

/********************************************************************/
IDL_LONG spherematch
  (int      argc,
   void *   argv[])
{
   IDL_LONG    npoints1;
   double    *  ra1;
   double    *  dec1;
   IDL_LONG    npoints2;
   double    *  ra2;
   double    *  dec2;
   double    matchlength;
   double    minchunksize;
	 IDL_LONG *match1;
	 IDL_LONG *match2;
   double    *  distance12;
	 IDL_LONG *nmatch;

         double myx1,myy1,myz1;
	 IDL_LONG maxmatch, jmax;
	 double currra,sep;
	 IDL_LONG i,j,k,rachunk,decchunk;
	 IDL_LONG retval=1;

   /* 0. allocate pointers from IDL */
   npoints1 = *((IDL_LONG *)argv[0]);
   ra1 = (double *)argv[1];
   dec1 = (double *)argv[2];
   npoints2 = *((IDL_LONG *)argv[3]);
   ra2 = (double *)argv[4];
   dec2 = (double *)argv[5];
   matchlength = *(double *)argv[6];
   minchunksize = *(double *)argv[7];
   match1 = (IDL_LONG *)argv[8];
   match2 = (IDL_LONG *)argv[9];
   distance12 = (double *)argv[10];
   nmatch = (IDL_LONG *)argv[11];

	 /* 1. define chunks */
	 setchunks(ra1,dec1,npoints1,minchunksize,&rabounds,
						 &decbounds,&nra,&ndec,&raoffset);

	 /* 2. assign targets to chunks, with minFibreSpacing of leeway */
	 assignchunks(ra2,dec2,npoints2,raoffset,matchlength,minchunksize,&nchunk2,
		      &chunklist2,rabounds,decbounds,nra,ndec);

	 /* 3. make x, y, z coords -- compute (x,y,z)1 on the fly - DPF */
	 xx2=(double *) IDL_GetScratch(&vxx2,npoints2,sizeof(double));
	 yy2=(double *) IDL_GetScratch(&vyy2,npoints2,sizeof(double));
	 zz2=(double *) IDL_GetScratch(&vzz2,npoints2,sizeof(double));
	 for(i=0;i<npoints2;i++) {
		 xx2[i]=cos(DEG2RAD*ra2[i])*cos(DEG2RAD*dec2[i]);
		 yy2[i]=sin(DEG2RAD*ra2[i])*cos(DEG2RAD*dec2[i]);
		 zz2[i]=sin(DEG2RAD*dec2[i]);
	 } /* end for i */

	 /* 4. run matching */
	 maxmatch = (*nmatch);  /* if nmatch != 0 then fill arrays up to maxmatch */
	 (*nmatch)=0;
	 for(i=0;i<npoints1;i++) {
	   currra=fmod(ra1[i]+raoffset,360.);
	   getchunk(currra,dec1[i],&rachunk,&decchunk,rabounds,decbounds,nra,ndec);
	   jmax=nchunk2[decchunk][rachunk];
	   if(jmax>0) {
	     myx1=cos(DEG2RAD*ra1[i])*cos(DEG2RAD*dec1[i]);
	     myy1=sin(DEG2RAD*ra1[i])*cos(DEG2RAD*dec1[i]);
	     myz1=sin(DEG2RAD*dec1[i]);
	     for(j=0;j<jmax;j++) {
	       k=chunklist2[decchunk][rachunk][j];
	       sep=separation(myx1,myy1,myz1,xx2[k],yy2[k],zz2[k]);
	       if(sep<matchlength) {
		 if(maxmatch>(*nmatch)) {
		   match1[(*nmatch)]=i;
		   match2[(*nmatch)]=k;
		   distance12[(*nmatch)]=sep;
		 }
		 (*nmatch)++;
	       } /* end if */
	     } /* end for j */
	   } /* end if jmax>0 */
	 } /* end for i */
	 
	 /* 4. clean up after chunks */
	 unassignchunks(&nchunk2,&chunklist2,nra,ndec);
	 unsetchunks(&rabounds,&decbounds,&nra,&ndec);

	 /* 6. free memory */
	 free_memory();

   return retval;
}

/******************************************************************************/
