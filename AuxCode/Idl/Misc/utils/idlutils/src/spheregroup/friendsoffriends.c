#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "export.h"
#include "friendsoffriends.h"

/*
 * Friends-of-friends using chunked data. That is, we run the friends
 * of friends on each chunk, which should have been designed to include
 * an extra layer of thickness linkSep, to guarantee that we get the 
 * groups linked over chunk borders.
 *
 * Then, we tranfer the group information from each chunk array to
 * the full array, equating groups which overlap.
 *
 * MB 5/2000
 */

static IDL_LONG *chunkFirstGroup=NULL; /* first member of group i in chunk */
static IDL_LONG *chunkMultGroup=NULL;  /* multiplicity of group i in chunk */
static IDL_LONG *chunkNextGroup=NULL; /* next member of group of element i
                                    in this chunk */
static IDL_LONG *chunkInGroup=NULL;     /* group of element i in this chunk */
static IDL_LONG chunkNGroups;      /* number of groups in this chunk */
static IDL_LONG *mapGroups=NULL;  /* equivalency mapping; mapGroup[igroup]=
															* index of an earlier group which is equivalent
															* to group igroup */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
	FREEVEC(chunkFirstGroup);
	FREEVEC(chunkMultGroup);
	FREEVEC(chunkNextGroup);
	FREEVEC(chunkInGroup);
	FREEVEC(mapGroups);
}

IDL_LONG 
friendsoffriends(double x[],
								 double y[],
								 double z[],
								 IDL_LONG nPoints,
								 double linkSep,
								 IDL_LONG **nChunk,
								 IDL_LONG ***chunkList,
								 IDL_LONG *nRa,
								 IDL_LONG nDec,
								 IDL_LONG firstGroup[],
								 IDL_LONG multGroup[],
								 IDL_LONG nextGroup[],
								 IDL_LONG inGroup[],
								 IDL_LONG *nGroups)
{
	IDL_LONG result;
	IDL_LONG i,j,k,l;
	IDL_LONG minEarly,checkEarly,tmpEarly,nMapGroups;
	IDL_LONG nChunkMax;
	
	/* 
	 * Make sure chunks have been properly set
	 */
	if(nChunk==NULL || chunkList==NULL || nRa==NULL || nDec==0) {
		fprintf(stderr,
						"Chunk lists not properly assigned in friendsoffriends().\n");
		free_memory();
		return(0);
	} /* end if */

	/*
	 * Find maximum number of galaxies in a chunk 
	 */
	nChunkMax=0;
	for(i=0;i<nDec;i++)
		for(j=0;j<nRa[i];j++)
			nChunkMax=(nChunk[i][j]>nChunkMax) ? nChunk[i][j] : nChunkMax;
	
	/* 
	 * Allocate memory using the maximum number in each chunk, so I don't
	 * constantly allocate and free
	 */
	chunkFirstGroup=(IDL_LONG *) malloc(nChunkMax*sizeof(IDL_LONG));
	chunkMultGroup=(IDL_LONG *) malloc(nChunkMax*sizeof(IDL_LONG));
	chunkNextGroup=(IDL_LONG *) malloc(nChunkMax*sizeof(IDL_LONG));
	chunkInGroup=(IDL_LONG *) malloc(nChunkMax*sizeof(IDL_LONG));

	/*
	 * mapGroups contains an equivalency mapping of groups. mapGroup[i]=j
	 * means i and j are actually the same group. j<=i always, by design.
	 * The greatest number groups you can get (assuming chunks have 
	 * linkSep<marginSize<minSize) is 9 times the number of targets.
	 */
	mapGroups=(IDL_LONG *) malloc(9*nPoints*sizeof(IDL_LONG)); 

	/*
	 * Initialize inGroup and mapGroups values for main list 
	 */
	for(i=0;i<nPoints;i++) 
		inGroup[i]=-1;
	for(i=0;i<9*nPoints;i++) 
		mapGroups[i]=-1;

	/*
	 * Run fof for each chunk and then merge it into the main 
	 * list of objects 
	 */
	nMapGroups=0;
	for(i=0;i<nDec;i++)
		for(j=0;j<nRa[i];j++) {

			/* Run friends of friends for each chunk */
			result=chunkfriendsoffriends(x, y, z, chunkList[i][j], nChunk[i][j], 
																	 linkSep, chunkFirstGroup, chunkMultGroup, 
																	 chunkNextGroup, chunkInGroup, 
																	 &chunkNGroups);
			if(result!=1) {
				fprintf(stderr,
								"chunkfriendsoffriends error %d in friendsoffriends()\n",
								(int) result);
				free_memory();
				return(result);
			} /* end if */

			/* see which clumps include objects in one or more
			 * earlier clumps, and find the earliest by following
			 * the map_clump links down all the way */
			for(k=0;k<chunkNGroups;k++) {

				/* make sure the group is real */
				if(chunkMultGroup[k]<=0) {
					fprintf(stderr,
									"chunkMultGroup[%d]=%d in friendsoffriends()\n",
									(int) k,(int) chunkMultGroup[k]);
					free_memory();
					return(0);
				} /* end if */

				/* search for links in group with earlier chunks and find earliest */
				minEarly=9*nPoints;
				for(l=chunkFirstGroup[k];l!=-1;l=chunkNextGroup[l]) {
					/* has this member been previously assigned? */
					if(inGroup[chunkList[i][j][l]]!=-1) {  
						/* previously assigned, set minEarly to earliest assignment */
						checkEarly=inGroup[chunkList[i][j][l]];
						while(mapGroups[checkEarly]!=checkEarly)
							checkEarly=mapGroups[checkEarly];
						minEarly=(minEarly<checkEarly) ? minEarly : checkEarly;
					} else {                          
						/* not previously assigned, assign to latest group
						 * as a placekeeper for these targets */
						inGroup[chunkList[i][j][l]]=nMapGroups;
					} /* end if */
				} /* end for l */

				/* go to each earlier group which any object in the current group
				 * has been assigned to, and map that group to the earliest */
				if(minEarly==9*nPoints) {   /* all are new, map group to itself */
					mapGroups[nMapGroups]=nMapGroups;
				} else {                /* at least one is old */
					mapGroups[nMapGroups]=minEarly; /* map current group, so 
																					 * new members will be reassigned */
					for(l=chunkFirstGroup[k];l!=-1;l=chunkNextGroup[l]) {
						/* find and reassign all earlier maps to earliest */
						checkEarly=inGroup[chunkList[i][j][l]];
						while(mapGroups[checkEarly]!=checkEarly) {
							tmpEarly=mapGroups[checkEarly];
							mapGroups[checkEarly]=minEarly;
							checkEarly=tmpEarly;
						} /* end while */
						mapGroups[checkEarly]=minEarly;
					} /* end for l */
				} /* end if..else */
				nMapGroups++;

				if(nMapGroups>=9*nPoints) {
					fprintf(stderr,
									"nMapGroups=%d has reached limit in friendsoffriends()\n",
									(int) nMapGroups);
					free_memory();
					return(0);
				} /* end if */
			} /* end for k */
		} /* end for i j */

	/* now all clumps which are mapped to themselves are 
	 * the "real" clumps; make sure the mappings are set
	 * up to go all the way down */
	(*nGroups)=0;
	for(i=0;i<nMapGroups;i++) 
		if(mapGroups[i]!=-1) {
			if(mapGroups[i]==i) {
				mapGroups[i]=(*nGroups);
				(*nGroups)++;
			} else {
				mapGroups[i]=mapGroups[mapGroups[i]];
			} 
		} else {
			fprintf(stderr,"mapGroups[%d]=%d in friendsoffriends()\n",
							(int) i,(int) mapGroups[i]);
			free_memory();
			return(0);
		}/* end if */

	/* reassign the inclump values to account for the mapping */
	for(i=0;i<nPoints;i++)
		inGroup[i]=mapGroups[inGroup[i]];
	
	/* Now set up clumps and llclumps based on inclump[] */
	for(i=0;i<nPoints;i++) 
		firstGroup[i]=-1;
	for(i=nPoints-1;i>=0;i--) {
		nextGroup[i]=firstGroup[inGroup[i]];
		firstGroup[inGroup[i]]=i;
	} /* end for i */
	
	/* Finally, return multiplicity of each group */
	for(i=0;i<(*nGroups);i++) {
		multGroup[i]=0;
		for(j=firstGroup[i];j!=-1;j=nextGroup[j]) {
			multGroup[i]++;
		} /* end for j */
		if(multGroup[i]<=0) {
			fprintf(stderr,"multGroup[%d]=%d in friendsoffriends()\n",
							(int) i,(int) multGroup[i]);
			free_memory();
			return(0);
		} /* end if */
	} /* end for i */
	
	/* 
	 * Free memory 
	 */
	chunkNGroups=0;
	free_memory();

	return(1);
} /* end fof */
