#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "export.h" 

/*
 * Does friends of friends on a sample within x, y, z which is
 * defined by the index list chunkList[], with linking length 
 * linkSep.
 *
 * Returns:
 *  firstGroup[]   (first member of group i)
 *  multGroup[]   (number of members in group i)
 *  nextGroup[]   (next members of group which element i is in)
 *  inGroup[]   (group which element i is in)
 *  nGroups    (number of groups)
 *
 */

#define RAD2DEG 57.29577951

static IDL_LONG *renumberedCFOF=NULL;

double separation(double x1, double y1, double z1, double x2, double y2,
									double z2);

IDL_LONG 
chunkfriendsoffriends(double x[],
												double y[],
												double z[],
												IDL_LONG chunkList[],
												IDL_LONG nTargets,
												double linkSep,
												IDL_LONG firstGroup[],
												IDL_LONG multGroup[],
												IDL_LONG nextGroup[],
												IDL_LONG inGroup[],
												IDL_LONG *nGroups)
{
	IDL_LONG i,j,k,minGroup,nTmp;
	double sep;

	/* initialization */
	(*nGroups)=0;
	for(i=0;i<nTargets;i++) {
		firstGroup[i]=-1;
		inGroup[i]=nTargets;
		nextGroup[i]=-1;
	} /* end for i */

	/* Find all of the other targets associated with each target;
	 * use multGroup as temporary storage here; find minimum group
	 * which any of the targets are in 
	 */
	(*nGroups)=0;
	for(i=0;i<nTargets;i++) {
		nTmp=0;
		minGroup=(*nGroups);
		for(j=0;j<nTargets;j++) {
			sep=separation(x[chunkList[i]],y[chunkList[i]],z[chunkList[i]],
										 x[chunkList[j]],y[chunkList[j]],z[chunkList[j]]);
			if(sep<=linkSep) {
				multGroup[nTmp]=j;
				minGroup=(minGroup>inGroup[j]) ? inGroup[j] : minGroup;
				nTmp++;
			} /* end if */
		} /* end for j */
			
		/* Use this minimum for all, including me! Note that
		 * if inGroup[multGroup[j]]<nTargets but is not the minimum,
		 * a group number has been eliminated, and we will
		 * have to go back later and renumber the groups */
		for(j=0;j<nTmp;j++) {
			if(inGroup[multGroup[j]]<nTargets) 
				for(k=firstGroup[inGroup[multGroup[j]]];k!=-1;k=nextGroup[k])
					inGroup[k]=minGroup;
			inGroup[multGroup[j]]=minGroup;
		} /* end for j */

		/* If it is a new group (no earlier groups), increment nGroups */
		if(minGroup==(*nGroups))
			(*nGroups)++;

		/* Now set up clumps and llclumps based on inclump() */
		for(j=0;j<=i;j++)
			firstGroup[j]=-1;
		for(j=i;j>=0;j--) {
			nextGroup[j]=firstGroup[inGroup[j]];
			firstGroup[inGroup[j]]=j;
		} /* end for i */
	} /* end for i */

	/* renumber the clumps to get rid of the 
	 * clump numbers which were skipped before */
	renumberedCFOF=(IDL_LONG *) malloc(nTargets*sizeof(IDL_LONG));
	for(i=0;i<nTargets;i++) 
		renumberedCFOF[i]=0;
	nTmp=(*nGroups);
	(*nGroups)=0;
	for(i=0;i<nTargets;i++) {
		if(!renumberedCFOF[i]) {
			for(j=firstGroup[inGroup[i]];j!=-1;j=nextGroup[j]) {
				inGroup[j]=(*nGroups);
				renumberedCFOF[j]=1;
			} /* end for j */
			(*nGroups)++;
		} /* end if */
	} /* end for i */
	free((char *) renumberedCFOF);
	renumberedCFOF=NULL;

	/* Now set up clumps and llclumps based on inGroup() */
	for(i=0;i<nTargets;i++)
		firstGroup[i]=-1;
	for(i=nTargets-1;i>=0;i--) {
		nextGroup[i]=firstGroup[inGroup[i]];
		firstGroup[inGroup[i]]=i;
	} /* end for i */

	/* Finally, return multiplicity of each group */
	for(i=0;i<(*nGroups);i++) {
		multGroup[i]=0;
		for(j=firstGroup[i];j!=-1;j=nextGroup[j]) 
			multGroup[i]++;
		if(multGroup[i]<=0) {
			fprintf(stderr,"multGroup[%d]=%d in chunkfriendsoffriends()\n",
							(int)i,(int)multGroup[i]);
			return(0);
		} /* end if */
	} /* end for i */

	return(1);
	
} /* end chunkfriendsoffriends */
