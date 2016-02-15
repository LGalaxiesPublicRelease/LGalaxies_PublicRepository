#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "export.h"
#include "chunks.h"

/*
 * This file contains routines to break up a sample into chunks in the
 * ra and dec directions. I allow you to set a minimize size (minSize)
 * of each chunk. It makes the length of each chunk in the dec
 * direction equal exactly to "minSize" and makes the length of each
 * chunk in the ra direction equal to "minSize" at the maximum
 * declination of each declination slice. This means that at different
 * declinations, there are a different number of chunks in the ra
 * direction.
 *
 * Doing this useful in order to do an efficient friends-of-friends (a
 * chunky one), and also for doing efficient matching.
 *
 * Now, setchunks feels free to redefine the zeropoint of the ra
 * coordinate system to simplify its work. It returns the offset it
 * applies in raOffset; this is not actually applied to the ra[]
 * array, it just needs to be senset to assignchunks to be used
 * there. If it finds that it needs to put a chunk within minSize of
 * the 0/360 boundary in a particular dec strip, it decides simply to
 * include the whole ra range (0. to 360.).  Likewise, if
 * getchunkbounds finds the ra range is such that an object on one
 * end of a declination strip can be linked to an object on the other
 * end, it allows raChunkMin and raChunkMax to have values of -1 and
 * nRa[i], if necessary. assignchunks then interprets these values
 * correctly.
 *
 * The spherical nature of the data makes this a somewhat complicated
 * piece of code, but I think I have managed to hide all of the
 * complications of doing efficient computation in a spherical geometry in
 * here, so they don't propagate to the rest of the code. In other words,
 * be happy that this works, and don't worry too much about the details.
 *
 * setchunks will set up the chunks according to where the targets
 * are, making sure that the chunks are always bigger than a given
 * size "minSize".  (unsetchunks cleans up memory from setchunks)
 * 
 * assignchunks will return a list of the points which are within
 * "marginSize" of each chunk; it actually includes more than this,
 * because it simply looks within marginSize arcsec in dec and in ra
 * (measured at the maximum declination of the chunk in question
 * (unassignchunks cleans up memory from assign_chunks)
 *
 * getchunkbounds will find the chunks which are within "marginSize" 
 * degrees of a given point
 *
 * getchunk will find the chunk which a given point is in (like calling
 * getchunkbounds with marginSize=0 but a bit faster)
 *
 * Michael Blanton
 * 4/2000 
 *
 * Fixed subtle 0/360 bug -- 12/2002
 */

#define NRA 6
#define DEG2RAD .01745329251994
#define EPS 1.e-5

/* set up the chunks based on where the objects are; this code
 * is allowed to redefine the ra array at its convenience, to avoid
 * stuff landing on the 0./360. boundary; it defines the arrays 
 * decBounds and raBounds which define the chunks; it leaves an 
 * extra 1 cell width boundary, as well; use unsetchunks
 * to clean up the memory */
CH_CODE 
setchunks(double ra[],         /* degrees */
					double dec[],         /* degrees */
					IDL_LONG nPoints,
					double minSize,       /* degrees */
					double ***raBounds,
					double **decBounds,
					IDL_LONG **nRa,
					IDL_LONG *nDec,
					double *raOffset)
{
	IDL_LONG i,j;
	double raMax,raMin,raRange;
	double raMaxTmp,raMinTmp,raRangeTmp;
	double decMax,decMin,decRange,cosDecMin;

	/* Make sure that minimum size makes sense */
	if(minSize<=0.) {
		fprintf(stderr,"minSize=%lf not positive in setchunks()\n",minSize);
		return(CH_ERROR);
	} /* end if */

	/* Make sure there are at least some points */
	if(nPoints<=0) {
		fprintf(stderr,"nPoints=%d not positive in setchunks()\n",(int)nPoints);
		return(CH_ERROR);
	} /* end if */

	/* Find maximum and minimum dec (in degrees) */
	decMax=-91.;
	decMin=91.;
	for(i=0;i<nPoints;i++) {
#if 0
		/* reassure ourselves it is in range */
		if(fabs(dec[i])>90.-DECLIMIT) {   
			fprintf(stderr,"dec=%lf for object %d too near pole in setchunks()\n", 
							dec[i],(int)i); 
			fprintf(stderr,"All objects must be > %lf degrees from either pole.\n",
							DECLIMIT);
			return(CH_ERROR);
		} /* end if */
#endif

		/* check minimum and maximum */
		decMin=(dec[i]<decMin) ? dec[i] : decMin;
		decMax=(dec[i]>decMax) ? dec[i] : decMax;
	} /* end for i */
	decRange=decMax-decMin;

	/* Find the declination boundaries; make them
	 * an integer multiple of minSize, with extra
	 * room (one cell) on the edges */
	(*nDec)=3+(IDL_LONG)floor(decRange/minSize);
	decRange=minSize*(double)(*nDec);
	decMin=decMin-0.5*(decRange-decMax+decMin);
	decMax=decMin+decRange;
	if(decMin<-90.+3.*minSize) decMin=-90.;
	if(decMax>90.-3.*minSize) decMax=90.;
	(*decBounds)=(double *) malloc(((*nDec)+1)*sizeof(double));
	for(i=0;i<=(*nDec);i++)
		(*decBounds)[i]=decMin+(decMax-decMin)*(double)i/(double)(*nDec);

	/* Find ra offset which minimizes the range
	 * in ra (this should take care of the case
	 * that ra crosses zero in some parts */
	cosDecMin=(fabs((*decBounds)[(*nDec)])>fabs((*decBounds)[0]))   
		? cos(DEG2RAD*(*decBounds)[(*nDec)]) : cos(DEG2RAD*(*decBounds)[0]);
	if(cosDecMin<=0.) {
		fprintf(stderr,"cosDecMin=%lf not positive in setchunks()\n",
						cosDecMin);
		return(CH_ERROR);
	} /* end if */
	rarange(ra,nPoints,minSize/cosDecMin,&raRange,raOffset);
	getraminmax(ra,(*raOffset),nPoints,&raMin,&raMax);
	raRange=raMax-raMin;
	
	/* For each declination slice, find the number of ra divisions
	 * necessary and set them */
	(*nRa)=(IDL_LONG *) malloc((*nDec)*sizeof(IDL_LONG));
	(*raBounds)=(double **) malloc((*nDec)*sizeof(double *));
	for(i=0;i<(*nDec);i++) {
		/* Get maximum declination and its cosine */
		cosDecMin=(fabs((*decBounds)[i])>fabs((*decBounds)[i+1]))   
			? cos(DEG2RAD*(*decBounds)[i]) : cos(DEG2RAD*(*decBounds)[i+1]);
		if(cosDecMin<=0.) {
			fprintf(stderr,"cosDecMin=%lf not positive in setchunks()\n",
							cosDecMin);
			return(CH_ERROR);
		} /* end if */
		
		/* get raBounds array for this declination array;
		 * again, leave an extra cell on each end */
		(*nRa)[i]=3+(IDL_LONG)floor(cosDecMin*raRange/minSize);
		raRangeTmp=minSize*(double)(*nRa)[i]/cosDecMin;
		raMinTmp=raMin-0.5*(raRangeTmp-raMax+raMin);
		raMaxTmp=raMinTmp+raRangeTmp;
		/* if we cannot avoid the 0/360 point, embrace it */
		if(raRangeTmp>=360. 
			 || raMinTmp<=minSize/cosDecMin
			 || raMaxTmp>=360.-minSize/cosDecMin
			 || (*decBounds)[i]==-90. || (*decBounds)[i]==90.) {
			raMinTmp=0.;
			raMaxTmp=360.;
			raRangeTmp=360.;
		} /* end if */
		if((*decBounds)[i]==-90. || (*decBounds)[i]==90.) 
			(*nRa)[i]=1;
		(*raBounds)[i]=(double *) malloc(((*nRa)[i]+1)*sizeof(double));
		for(j=0;j<=(*nRa)[i];j++)
			(*raBounds)[i][j]=raMinTmp+(raMaxTmp-raMinTmp)*(double)j
				/(double)(*nRa)[i];
	} /* end for i */
	
	return(CH_OK);
} /* setchunks */

/* cleans up the memory allocated in setchunks */
CH_CODE
unsetchunks(double ***raBounds,
						double **decBounds,
						IDL_LONG **nRa,
						IDL_LONG *nDec)
{
	IDL_LONG i;

	/* free up decBounds */
	if((*decBounds)!=NULL) {
		free((char *) (*decBounds));
		(*decBounds)=NULL;
	} /* end if */

	/* free up raBounds */
	if((*raBounds)!=NULL) {
		for(i=0;i<(*nDec);i++) {
			if((*raBounds)[i]!=NULL)
				free((char *) (*raBounds)[i]);
			(*raBounds)[i]=NULL;
		} 
		free((char *) (*raBounds));
		(*raBounds)=NULL;
	} /* end if */

	/* free up nRa and nDec */
	if((*nRa)!=NULL) {
		free((char *) (*nRa));
		(*nRa)=NULL;
	} /* end if */
	(*nDec)=0;
	
	return(CH_OK);
} /* unsetchunks */

/* take the objects and the chunks (already defined using setchunks)
 * and assign the objects to the appropriate chunks, with some leeway 
 * given by the parameter marginSize; use unassignchunks to get rid of 
 * the memory allocated here */
CH_CODE
assignchunks(double ra[],        /* degrees */
						 double dec[],       /* degrees */
						 IDL_LONG nPoints,
						 double raOffset,   /* ra offset to apply */
						 double marginSize,    /* degrees */
						 double minSize,       /* degrees */
						 IDL_LONG ***nChunk,     /* number of targets in each chunk */
						 IDL_LONG ****chunkList,   /* index of targets in each chunk */
						 double **raBounds,  /* ra divisions for each dec */
						 double *decBounds,     /* 1d array of declination divs */
						 IDL_LONG *nRa,          /* number of ra divs for each dec */
						 IDL_LONG nDec)          /* number of dec divisions */
{
	CH_CODE result;
	double currRa;
	IDL_LONG decChunk,raChunk,currRaChunk;
	IDL_LONG *raChunkMin, *raChunkMax, decChunkMin, decChunkMax;
	IDL_LONG i,j;

	/* initialize pointers */
	raChunkMin=raChunkMax=NULL;

	/* Check that marginSize is smaller than minSize */
	if(marginSize>=minSize) {
		fprintf(stderr,"marginSize>=minSize (%lf>=%lf) in assignchunks()\n",
						marginSize,minSize);
		return(CH_ERROR);
	} /* end if */

	/* Check number of points */
	if(nPoints<=0) {
		fprintf(stderr,"nPoints=%d not positive in assignchunks()\n",
						(int)nPoints);
		return(CH_ERROR);
	} /* end if */

	/* Check that setchunks has been called */
	if(raBounds==NULL || decBounds==NULL ||
		 nRa==NULL || nDec==0) {
		fprintf(stderr,"setchunks not called before assignchunks()?\n");
		return(CH_ERROR);
	} /* end if */

	/* Allocate memory for number in each chunk */
	(*nChunk)=(IDL_LONG **) malloc(nDec*sizeof(IDL_LONG *));
	for(i=0;i<nDec;i++) 
		(*nChunk)[i]=(IDL_LONG *) malloc(nRa[i]*sizeof(IDL_LONG));
	
	/* Reset nChunk */
	for(i=0;i<nDec;i++)
		for(j=0;j<nRa[i];j++)
			(*nChunk)[i][j]=0;

	/* Find number in each chunk */
	for(i=0;i<nPoints;i++) {
		/* get bounds */
		currRa=fmod(ra[i]+raOffset,360.);
		result=getchunkbounds(currRa,dec[i],marginSize,&raChunkMin,&raChunkMax,
													&decChunkMin,&decChunkMax,raBounds,decBounds,nRa,
													nDec);

		/* only perform assignment if we are in range */
		if(result!=CH_OK) {
			if(result!=CH_OUTOFRANGE) {
				fprintf(stderr,
								"getchunkbounds returned error %d in assignchunks()\n",
								(int)result);
				if(raChunkMin!=NULL) free((char *) raChunkMin);
				if(raChunkMax!=NULL) free((char *) raChunkMax);
				raChunkMin=raChunkMax=NULL;
				return(result);
			} /* end if */
		} else {

			/* check dec bounds */
			if(decChunkMin<0 || decChunkMin>=nDec ||
				 decChunkMax<0 || decChunkMax>=nDec || decChunkMin>decChunkMax) {
				fprintf(stderr,
								"decChunkMin=%d, decChunkMax=%d illegal in assignchunks()\n",
								(int)decChunkMin,(int)decChunkMax);
				if(raChunkMin!=NULL) free((char *) raChunkMin);
				if(raChunkMax!=NULL) free((char *) raChunkMax);
				raChunkMin=raChunkMax=NULL;
				return(CH_OUTOFRANGE);
			} /* end if */
			
			/* check ra arrays */
			if(raChunkMin==NULL || raChunkMax==NULL) {
				fprintf(stderr,
								"raChunkMin==NULL or raChunkMax==NULL in assignchunks()\n");
				if(raChunkMin!=NULL) free((char *) raChunkMin);
				if(raChunkMax!=NULL) free((char *) raChunkMax);
				raChunkMin=raChunkMax=NULL;
				return(CH_ERROR);
			} /* end if */
			
			/* go through each declination slice */
			for(decChunk=decChunkMin;decChunk<=decChunkMax;decChunk++) {
				
				/* check ra bounds */
				if(raChunkMin[decChunk-decChunkMin]<-1 || 
					 raChunkMin[decChunk-decChunkMin]>nRa[decChunk] ||
					 raChunkMax[decChunk-decChunkMin]<-1 || 
					 raChunkMax[decChunk-decChunkMin]>nRa[decChunk] || 
					 raChunkMin[decChunk-decChunkMin]>raChunkMax[decChunk-decChunkMin]) {
					fprintf(stderr,
									"raChunkMin=%d, raChunkMax=%d illegal in assignchunks()\n",
									(int)raChunkMin[decChunk-decChunkMin],
									(int)raChunkMax[decChunk-decChunkMin]);
					if(raChunkMin!=NULL) free((char *) raChunkMin);
					if(raChunkMax!=NULL) free((char *) raChunkMax);
					raChunkMin=raChunkMax=NULL;
					return(CH_OUTOFRANGE);
				} /* end if */
				
				/* go through each chunk in slice; be conservative and 
				 * include buffers when allocating memory */
				for(raChunk=raChunkMin[decChunk-decChunkMin]-1;
						raChunk<=raChunkMax[decChunk-decChunkMin]+1;
						raChunk++) {
					/* handle edge cases if necessary; setchunkbounds is
					 * supposed to decide if raChunk can be -1 or nra[decChunk] */
					if(raChunk<0) {
						currRaChunk=(raChunk+nRa[decChunk])%nRa[decChunk];
						if(currRaChunk>=0) 
							(*nChunk)[decChunk][currRaChunk]++;
					} else if (raChunk>nRa[decChunk]-1) {
						currRaChunk=(raChunk-nRa[decChunk])%nRa[decChunk];
						if(currRaChunk<=nRa[decChunk]-1)
							(*nChunk)[decChunk][currRaChunk]++;
					} else {
						(*nChunk)[decChunk][raChunk]++;
					} /* end if */
				} /* end for raChunk */
			} /* end for decChunk */
		} /* end if */

		/* free memory */
		if(raChunkMin!=NULL) free((char *) raChunkMin);
		if(raChunkMax!=NULL) free((char *) raChunkMax);
		raChunkMin=raChunkMax=NULL;
	} /* end for i */

	/* free memory */
	if(raChunkMin!=NULL) free((char *) raChunkMin);
	if(raChunkMax!=NULL) free((char *) raChunkMax);
	raChunkMin=raChunkMax=NULL;

	/* Allocate chunkList memory */
	(*chunkList)=(IDL_LONG ***) malloc(nDec*sizeof(IDL_LONG **));
	for(i=0;i<nDec;i++) {
		(*chunkList)[i]=(IDL_LONG **) malloc(nRa[i]*sizeof(IDL_LONG *));
		for(j=0;j<nRa[i];j++) {
			if((*nChunk)[i][j]>0)
				(*chunkList)[i][j]=(IDL_LONG *) 
					malloc((*nChunk)[i][j]*sizeof(IDL_LONG));
			else
				(*chunkList)[i][j]=NULL;
		} /* end for i */
	} /* end for i */
	
	/* Reset nChunk */
	for(i=0;i<nDec;i++)
		for(j=0;j<nRa[i];j++)
			(*nChunk)[i][j]=0;

	/* Construct chunkList */
	for(i=0;i<nPoints;i++) {

		/* get bounds */
		currRa=fmod(ra[i]+raOffset,360.);
		result=getchunkbounds(currRa,dec[i],marginSize,&raChunkMin,&raChunkMax,
													&decChunkMin,&decChunkMax,raBounds,decBounds,nRa,
													nDec);
		if(result!=CH_OK) {
			if(result!=CH_OUTOFRANGE) {
				fprintf(stderr,
								"getchunkbounds returned error %d in assignchunks()\n",
								(int)result);
				if(raChunkMin!=NULL) free((char *) raChunkMin);
				if(raChunkMax!=NULL) free((char *) raChunkMax);
				raChunkMin=raChunkMax=NULL;
				return(result);
			} /* end if */
		} else {

			/* check dec bounds */
			if(decChunkMin<0 || decChunkMin>=nDec ||
				 decChunkMax<0 || decChunkMax>=nDec || decChunkMin>decChunkMax) {
				fprintf(stderr,
								"decChunkMin=%d, decChunkMax=%d illegal in assignchunks()\n",
								(int)decChunkMin,(int)decChunkMax);
				if(raChunkMin!=NULL) free((char *) raChunkMin);
				if(raChunkMax!=NULL) free((char *) raChunkMax);
				raChunkMin=raChunkMax=NULL;
				return(CH_OUTOFRANGE);
			} /* end if */
			
			/* check ra arrays */
			if(raChunkMin==NULL || raChunkMax==NULL) {
				fprintf(stderr,
								"raChunkMin==NULL or raChunkMax==NULL in assignchunks()\n");
				if(raChunkMin!=NULL) free((char *) raChunkMin);
				if(raChunkMax!=NULL) free((char *) raChunkMax);
				raChunkMin=raChunkMax=NULL;
				return(CH_ERROR);
			} /* end if */
			
			/* go through each declination slice */
			for(decChunk=decChunkMin;decChunk<=decChunkMax;decChunk++) {
				
				/* check ra bounds */
				if(raChunkMin[decChunk-decChunkMin]<-1 || 
					 raChunkMin[decChunk-decChunkMin]>nRa[decChunk] ||
					 raChunkMax[decChunk-decChunkMin]<-1 || 
					 raChunkMax[decChunk-decChunkMin]>nRa[decChunk] || 
					 raChunkMin[decChunk-decChunkMin]>raChunkMax[decChunk-decChunkMin]) {
					fprintf(stderr,
									"raChunkMin=%d, raChunkMax=%d illegal in assignchunks()\n",
									(int)raChunkMin[decChunk-decChunkMin],
									(int)raChunkMax[decChunk-decChunkMin]);
					if(raChunkMin!=NULL) free((char *) raChunkMin);
					if(raChunkMax!=NULL) free((char *) raChunkMax);
					raChunkMin=raChunkMax=NULL;
					return(CH_OUTOFRANGE);
				} /* end if */
				
				/* go through each chunk in slice */
				for(raChunk=raChunkMin[decChunk-decChunkMin];
						raChunk<=raChunkMax[decChunk-decChunkMin];
						raChunk++) {
					
					/* handle edge cases if necessary; setchunkbounds is
					 * supposed to decide if raChunk can be -1 or nra[decChunk] */
					if(raChunk<0) {
						currRaChunk=(raChunk+nRa[decChunk])%nRa[decChunk];
						(*chunkList)[decChunk][currRaChunk]
							[(*nChunk)[decChunk][currRaChunk]]=i;
						(*nChunk)[decChunk][currRaChunk]++;
					} else if (raChunk>nRa[decChunk]-1) {
						currRaChunk=(raChunk-nRa[decChunk])%nRa[decChunk];
						(*chunkList)[decChunk][currRaChunk]
							[(*nChunk)[decChunk][currRaChunk]]=i;
						(*nChunk)[decChunk][currRaChunk]++;
					} else {
						(*chunkList)[decChunk][raChunk][(*nChunk)[decChunk][raChunk]]=i;
						(*nChunk)[decChunk][raChunk]++;
					} /* end if */
				} /* end for raChunk */
			} /* end for decChunk */
		} /* end if..else */

		/* free memory */
		if(raChunkMin!=NULL) free((char *) raChunkMin);
		if(raChunkMax!=NULL) free((char *) raChunkMax);
		raChunkMin=raChunkMax=NULL;
	} /* end for i */

	/* free memory */
	if(raChunkMin!=NULL) free((char *) raChunkMin);
	if(raChunkMax!=NULL) free((char *) raChunkMax);
	raChunkMin=raChunkMax=NULL;

	return(CH_OK);
} /* end assign_chunks */

/* clean up memory allocated in assign_chunks */
CH_CODE 
unassignchunks(IDL_LONG ***nChunk,     /* number of targets in each chunk */
							 IDL_LONG ****chunkList, /* index of targets in each chunk */
							 IDL_LONG *nRa,          /* number of ra divs for each dec */
							 IDL_LONG nDec)          /* number of dec divisions */
{
	IDL_LONG i,j;

	/* get rid of chunkList */
	if((*chunkList)!=NULL) {
		for(i=0;i<nDec;i++) {
			if((*chunkList)[i]!=NULL) {
				for(j=0;j<nRa[i];j++) {
					if((*chunkList)[i][j]!=NULL)
						free((char *) (*chunkList)[i][j]);
					(*chunkList)[i][j]=NULL;
				}
				free((char *) (*chunkList)[i]);
			} 
			(*chunkList)[i]=NULL;
		} /* end for i */
		free((char *) (*chunkList));
		(*chunkList)=NULL;
	} /* end if */

	/* get rid of nChunk */
	if((*nChunk)!=NULL) {
		for(i=0;i<nDec;i++) {
			if((*nChunk)[i]!=NULL) free((char *) (*nChunk)[i]);
			(*nChunk)[i]=NULL;
		}
		free((char *) (*nChunk));
		(*nChunk)=NULL;
	} /* end if */
	
	return(CH_OK);
} /* end unassignchunks */

/* utility to find the set of chunks which a given point belongs to;
 * if raBounds wraps around 0/360, it allows -1 and nRa[decChunk] to 
 * be used as raChunkMin or raChunkMax */
CH_CODE
getchunkbounds(double ra, 
							 double dec,
							 double marginSize,
							 IDL_LONG **raChunkMin,  
							 IDL_LONG **raChunkMax,
							 IDL_LONG *decChunkMin, 
							 IDL_LONG *decChunkMax,
							 double **raBounds,
							 double *decBounds,
							 IDL_LONG *nRa,
							 IDL_LONG nDec)
{
	double cosDecMin;
	IDL_LONG i;

	/* Check that setchunks has been called */
	if(raBounds==NULL || decBounds==NULL ||
		 nRa==NULL || nDec==0) {
		fprintf(stderr,"setchunks not called before getchunkbounds()?\n");
		return(CH_ERROR);
	} /* end if */

	/* find which declination slice we are in and make sure it is in 
	 * range */
	(*decChunkMin)=(*decChunkMax)=
		(IDL_LONG)floor((dec-decBounds[0])*(double)nDec/
										(decBounds[nDec]-decBounds[0]));
	if((*decChunkMin)<0 || (*decChunkMin)>nDec-1) 
		return(CH_OUTOFRANGE);
	
	/* set minimum and maximum bounds of dec; in fact, this 
	 * step depends on minSize>marginSize, since it only looks 
	 * at the nearest neighbor chunks  */
	if(dec-decBounds[(*decChunkMin)]<marginSize && (*decChunkMin)>0) 
		(*decChunkMin)--;
	if(decBounds[(*decChunkMax)+1]-dec<marginSize && (*decChunkMax)<nDec-1)
		(*decChunkMax)++;

	/* find ra chunk bounds for each dec chunk */
	(*raChunkMin)=
		(IDL_LONG *) malloc(((*decChunkMax)-(*decChunkMin)+1)*sizeof(IDL_LONG));
	(*raChunkMax)=
		(IDL_LONG *) malloc(((*decChunkMax)-(*decChunkMin)+1)*sizeof(IDL_LONG));
	for(i=(*decChunkMin);i<=(*decChunkMax);i++) {
		/* Get maximum declination and its cosine */
		cosDecMin=(fabs(decBounds[i])>fabs(decBounds[i+1]))   
			? cos(DEG2RAD*decBounds[i]) : cos(DEG2RAD*decBounds[i+1]);

		/* find which declination slice we are in and make sure it is in
		 * range */
		(*raChunkMin)[i-(*decChunkMin)]=(*raChunkMax)[i-(*decChunkMin)]=
			(IDL_LONG)floor((ra-raBounds[i][0])*(double)nRa[i]
											/(raBounds[i][nRa[i]]-raBounds[i][0]));
		if((*raChunkMin)[i-(*decChunkMin)]<0 || 
			 (*raChunkMin)[i-(*decChunkMin)]>nRa[i]-1) {
			if((*raChunkMin)!=NULL) free((char *) (*raChunkMin));
			if((*raChunkMax)!=NULL) free((char *) (*raChunkMax));
			(*raChunkMin)=(*raChunkMax)=NULL;
			return(CH_OUTOFRANGE);
		} /* end if */

		/* set minimum and maximum bounds of ra; in fact, this 
		 * step depends on minSize>marginSize, since it only looks 
		 * at the nearest neighbor chunks  */
		if((ra-raBounds[i][(*raChunkMin)[i-(*decChunkMin)]])*cosDecMin<marginSize)
			(*raChunkMin)[i-(*decChunkMin)]--;
		if((raBounds[i][(*raChunkMax)[i-(*decChunkMin)]+1]-ra)*cosDecMin
			 <marginSize)
			(*raChunkMax)[i-(*decChunkMin)]++;
		
#if 0
		/* if the ra range is such that objects can link over the 
		 * 0/360 border, allow -1 and nRa[i] to be included; if not,
		 * fix it here. */
		if(raBounds[nRa[i]]-raBounds[0]>360.-marginSize/cosDecMin) { 
			if((*raChunkMin)[i-(*decChunkMin)]<0)
				(*raChunkMin)[i-(*decChunkMin)]=0;
			if((*raChunkMax)[i-(*decChunkMin)]>nRa[i]-1)
				(*raChunkMax)[i-(*decChunkMin)]=nRa[i]-1;
		} /* end if */
#endif
	} /* end for i */

	return(CH_OK);
} /* end getchunkbounds */ 

/* utility to find the chunk which a given point belongs to */
CH_CODE
getchunk(double ra, 
				 double dec,
				 IDL_LONG *raChunk,  
				 IDL_LONG *decChunk, 
				 double **raBounds,
				 double *decBounds,
				 IDL_LONG *nRa,
				 IDL_LONG nDec)
{
	/* Check that setchunks has been called */
	if(raBounds==NULL || decBounds==NULL ||
		 nRa==NULL || nDec==0) {
		fprintf(stderr,"setchunks not called before getchunk()?\n");
		return(CH_ERROR);
	} /* end if */

	/* find dec chunk */
	(*decChunk)=
		(IDL_LONG)floor((dec-decBounds[0])*(double)nDec/
										(decBounds[nDec]-decBounds[0]));
	if((*decChunk)<0 || (*decChunk)>nDec-1) return(CH_OUTOFRANGE);
	
	/* find ra chunk */
	if((*decChunk)<nDec && (*decChunk)>=0) {
		(*raChunk)=
			(IDL_LONG)floor((ra-raBounds[(*decChunk)][0])*(double)nRa[(*decChunk)]
											/(raBounds[(*decChunk)][nRa[(*decChunk)]]
												-raBounds[(*decChunk)][0]));
		if((*raChunk)<0 || 
			 (*raChunk)>nRa[(*decChunk)]-1) 
			return(CH_OUTOFRANGE);
	} else {
		(*raChunk)=-1;
	} /* end if */

	return(CH_OK);
} /* end getchunk */ 


