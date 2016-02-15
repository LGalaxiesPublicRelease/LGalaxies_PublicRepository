#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * dfind.c
 *
 * Find non-zero objects in a binary image.
 *
 * Mike Blanton
 * 1/2006 */

#define PI 3.14159265358979

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

static int *matches=NULL;
static int *nmatches=NULL;
static int *mapgroup=NULL;

int dfind(int *image, 
          int nx, 
          int ny,
          int *object)
{
  int i,ip,j,jp,k,kp,l,ist,ind,jst,jnd,igroup,minearly,checkearly,tmpearly;
  int ngroups;

  mapgroup=(int *) malloc((size_t) nx*ny*sizeof(int));
  matches=(int *) malloc((size_t) nx*ny*9*sizeof(int));
  nmatches=(int *) malloc((size_t) nx*ny*sizeof(int));

  for(k=0;k<nx*ny;k++)
    object[k]=-1;
  for(k=0;k<nx*ny;k++)
    mapgroup[k]=-1;
  for(k=0;k<nx*ny;k++)
    nmatches[k]=0;
  for(k=0;k<nx*ny*9;k++)
    matches[k]=-1;
  
  /* find matches */
  for(j=0;j<ny;j++) {
    jst=j-1;
    jnd=j+1;
    if(jst<0) jst=0;
    if(jnd>ny-1) jnd=ny-1;
    for(i=0;i<nx;i++) {
      ist=i-1;
      ind=i+1;
      if(ist<0) ist=0;
      if(ind>nx-1) ind=nx-1;
      k=i+j*nx;
      if(image[k]) {
        for(jp=jst;jp<=jnd;jp++) 
          for(ip=ist;ip<=ind;ip++) {
            kp=ip+jp*nx;
            if(image[kp]) {
              matches[9*k+nmatches[k]]=kp;
              nmatches[k]++;
            }
          }
      } /* end if */
    }
  }

  /* group pixels on matches */
  igroup=0;
  for(k=0;k<nx*ny;k++) {
    if(image[k]) {
      minearly=igroup;
      for(l=0;l<nmatches[k];l++) {
        kp=matches[9*k+l];
        checkearly=object[kp];
        if(checkearly>=0) {
          while(mapgroup[checkearly]!=checkearly) {
            checkearly=mapgroup[checkearly];
          }
          if(checkearly<minearly) minearly=checkearly;
        }
      }
      
      if(minearly==igroup) {
        mapgroup[igroup]=igroup;
        for(l=0;l<nmatches[k];l++) {
          kp=matches[9*k+l];
          object[kp]=igroup;
        }
        igroup++;
      } else {
        for(l=0;l<nmatches[k];l++) {
          kp=matches[9*k+l];
          checkearly=object[kp];
          if(checkearly>=0) {
            while(mapgroup[checkearly]!=checkearly) {
              tmpearly=mapgroup[checkearly];
              mapgroup[checkearly]=minearly;
              checkearly=tmpearly;
            }
            mapgroup[checkearly]=minearly;
          }
        }
        for(l=0;l<nmatches[k];l++) {
          kp=matches[9*k+l];
          object[kp]=minearly;
        }
      }
    }
  }

  ngroups=0;
  for(i=0;i<nx*ny;i++) {
    if(mapgroup[i]>=0) {
      if(mapgroup[i]==i) {
        mapgroup[i]=ngroups;
        ngroups++;
      } else {
        mapgroup[i]=mapgroup[mapgroup[i]];
      }
    }
  }

  for(i=0;i<nx*ny;i++) 
    object[i]=mapgroup[object[i]];
  
  for(i=0;i<nx*ny;i++) 
    mapgroup[i]=-1;
  igroup=0;
  for(k=0;k<nx*ny;k++) {
    if(image[k]>0 && mapgroup[object[k]]==-1) {
      mapgroup[object[k]]=igroup;
      igroup++;
    }
  }

  for(i=0;i<nx*ny;i++) 
    if(image[i]>0)
      object[i]=mapgroup[object[i]];
    else 
      object[i]=-1;
  
  FREEVEC(matches);
  FREEVEC(nmatches);
  FREEVEC(mapgroup);
  
	return(1);
} /* end dfind */
