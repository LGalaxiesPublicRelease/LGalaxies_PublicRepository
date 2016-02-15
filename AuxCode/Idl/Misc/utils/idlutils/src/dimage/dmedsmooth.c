#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * dmedsmooth.c
 *
 * Median smooth an image
 *
 * Mike Blanton
 * 1/2006 */

#define PI 3.14159265358979

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

static float *arr=NULL;
static int *xgrid=NULL;
static int *ygrid=NULL;
static float *grid=NULL;
static int *xlo=NULL;
static int *xhi=NULL;
static int *ylo=NULL;
static int *yhi=NULL;

float dselip(unsigned long k, unsigned long n, float *arr);

int dmedsmooth(float *image, 
               float *invvar,
							 int nx, 
							 int ny,
							 int box,
							 float *smooth)
{
  int i,j,ip,jp,ist,jst,nxt,nyt,nb,ind,jnd,sp;
  int xoff, yoff, nm, nxgrid, nygrid;
  int ypsize, ymsize, xpsize, xmsize;
  float dx,dy,xkernel,ykernel;

	/* guarantee odd */
	sp=box;

	/* get grids */
  nxgrid=nx/sp+2;
	xlo=(int *) malloc(nxgrid*sizeof(int));
	xhi=(int *) malloc(nxgrid*sizeof(int));
	xgrid=(int *) malloc(nxgrid*sizeof(int));
	xoff=(nx-1-(nxgrid-3)*sp)/2;
	for(i=1;i<nxgrid-1;i++)
		xgrid[i]=(i-1)*sp+xoff;
	xgrid[0]=xgrid[1]-sp;
	xgrid[nxgrid-1]=xgrid[nxgrid-2]+sp;
	for(i=0;i<nxgrid;i++) {
    xlo[i]=xgrid[i]-sp;
    if(xlo[i]<0) xlo[i]=0;
    xhi[i]=xgrid[i]+sp;
    if(xhi[i]>nx-1) xhi[i]=nx-1;
  }

  nygrid=ny/sp+2;
	ylo=(int *) malloc(nygrid*sizeof(int));
	yhi=(int *) malloc(nygrid*sizeof(int));
	ygrid=(int *) malloc(nygrid*sizeof(int));
	yoff=(ny-1-(nygrid-3)*sp)/2;
	for(i=1;i<nygrid-1;i++)
		ygrid[i]=(i-1)*sp+yoff;
	ygrid[0]=ygrid[1]-sp;
	ygrid[nygrid-1]=ygrid[nygrid-2]+sp;
  
	for(i=0;i<nygrid;i++) {
    ylo[i]=ygrid[i]-sp;
    if(ylo[i]<0) ylo[i]=0;
    yhi[i]=ygrid[i]+sp;
    if(yhi[i]>ny-1) yhi[i]=ny-1;
  }

	grid=(float *) malloc(nxgrid*nygrid*sizeof(float));

  arr=(float *) malloc((sp*2+5)*(sp*2+5)*sizeof(float));

  for(j=0;j<nygrid;j++) 
    for(i=0;i<nxgrid;i++) 
      grid[i+j*nxgrid]=0.;

  for(j=0;j<nygrid;j++) {
    jst=ylo[j];
    jnd=yhi[j];
    nyt=jnd-jst+1;
    for(i=0;i<nxgrid;i++) {
      ist=xlo[i];
      ind=xhi[i];
      nxt=ind-ist+1;
      nb=0;
      for(jp=jst;jp<=jnd;jp++) 
        for(ip=ist;ip<=ind;ip++) {
          if(invvar[ip+jp*nx]>0.) {
            arr[nb]=image[ip+jp*nx];
            nb++;
          }
        }
      if(nb>1) {
        nm=nb/2;
        grid[i+j*nxgrid]=dselip(nm,nb,arr);
      } else {
        grid[i+j*nxgrid]=image[(long) xlo[i]+((long) ylo[j])*nx];
      }
    }
  }

  for(j=0;j<ny;j++) 
    for(i=0;i<nx;i++) 
      smooth[i+j*nx]=0.;
  for(j=0;j<nygrid;j++) {
    jst=(long) ( (float) ygrid[j] - sp*1.5);
    jnd=(long) ( (float) ygrid[j] + sp*1.5);
    if(jst<0) jst=0;
    if(jnd>ny-1) jnd=ny-1;
    ypsize=sp;
    ymsize=sp;
    if(j==0) ypsize=ygrid[1]-ygrid[0];
    if(j==1) ymsize=ygrid[1]-ygrid[0];
    if(j==nygrid-2) ypsize=ygrid[nygrid-1]-ygrid[nygrid-2];
    if(j==nygrid-1) ymsize=ygrid[nygrid-1]-ygrid[nygrid-2];
    for(i=0;i<nxgrid;i++) {
      ist=(long) ( (float) xgrid[i] - sp*1.5);
      ind=(long) ( (float) xgrid[i] + sp*1.5);
      if(ist<0) ist=0;
      if(ind>nx-1) ind=nx-1;
      xpsize=sp;
      xmsize=sp;
      if(i==0) xpsize=xgrid[1]-xgrid[0];
      if(i==1) xmsize=xgrid[1]-xgrid[0];
      if(i==nxgrid-2) xpsize=xgrid[nxgrid-1]-xgrid[nxgrid-2];
      if(i==nxgrid-1) xmsize=xgrid[nxgrid-1]-xgrid[nxgrid-2];

      for(jp=jst;jp<=jnd;jp++) {
        dy=((float) jp-ygrid[j]);
        ykernel=0.;
        if(dy>-1.5*ymsize && dy<=-0.5*ymsize) 
          ykernel=0.5*(dy/ymsize+1.5)*(dy/ymsize+1.5);
        else if(dy>-0.5*ymsize && dy<0.) 
          ykernel=-(dy*dy/ymsize/ymsize-0.75);
        else if(dy<0.5*ypsize && dy>=0.) 
          ykernel=-(dy*dy/ypsize/ypsize-0.75);
        else if(dy>=0.5*ypsize && dy <1.5*ypsize) 
          ykernel=0.5*(dy/ypsize-1.5)*(dy/ypsize-1.5);
        for(ip=ist;ip<=ind;ip++) {
          dx=((float) ip-xgrid[i]);
          xkernel=0.;
          if(dx>-1.5*xmsize && dx<=-0.5*xmsize) 
            xkernel=0.5*(dx/xmsize+1.5)*(dx/xmsize+1.5);
          else if(dx>-0.5*xmsize && dx<0.) 
            xkernel=-(dx*dx/xmsize/xmsize-0.75);
          else if(dx<0.5*xpsize && dx>=0.) 
            xkernel=-(dx*dx/xpsize/xpsize-0.75);
          else if(dx>=0.5*xpsize && dx <1.5*xpsize) 
            xkernel=0.5*(dx/xpsize-1.5)*(dx/xpsize-1.5);
          smooth[ip+jp*nx]+=xkernel*ykernel*grid[i+j*nxgrid];
        }
      }
    }
  }
  
#if 0
  for(j=0;j<nygrid;j++) 
    for(i=0;i<nxgrid;i++) 
      smooth[i+j*nx]=grid[i+j*nxgrid];
#endif
    
  FREEVEC(arr);
  FREEVEC(grid);
  FREEVEC(xgrid);
  FREEVEC(ygrid);
  FREEVEC(xlo);
  FREEVEC(ylo);
  FREEVEC(xhi);
  FREEVEC(yhi);
  
	return(1);
} /* end dmedsmooth */
