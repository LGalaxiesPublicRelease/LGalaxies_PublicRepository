#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ph.h"

/*
 * reject_cr_psf.c
 *
 * Based on an estimate of the PSF, reject pixels in an image 
 * which appear to have violated the PSF bounds, and are thus likely
 * to be PSFs. Does not alter original image, just returns the indices
 * of the rejected pixels. 
 *
 * Uses RHL's prescription for this a la the PHOTO paper. 
 *
 * Does not check edge pixels at ALL. 
 *
 * Mike Blanton
 * 10/2003 */

#define PI 3.14159265358979
#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

int reject_cr_psf(float *image, 
                  float *image_ivar, 
                  int xnpix, 
                  int ynpix,
                  float nsig, /* number of sigma above background required */
                  float *psfvals, /* psf value at radii of 1 pix and 
                                     sqrt(2) pix */
                  float cfudge, /* number of sigma inconsistent with
                                   PSF required */
                  float c2fudge, /* fudge factor applied to PSF */
                  int *rejected, 
                  int *ignoremask)
{
  float *sigmaback,*back,imcurr,ival,lval,invsigma,goodback[3][3];
  int i,j,ip,jp;

  sigmaback=(float *) malloc(4*sizeof(float));
  back=(float *) malloc(4*sizeof(float));
  
  for(j=1;j<ynpix-1;j++) {
    for(i=1;i<xnpix-1;i++) {
      rejected[j*xnpix+i]=0;
      if(image_ivar[j*xnpix+i]>0. && ignoremask[j*xnpix+i]==0) {
        invsigma=sqrt(image_ivar[j*xnpix+i]);
        imcurr=image[j*xnpix+i];
        
        /* check if it exceeds background for ALL four pairs */
        ival=invsigma*imcurr;
        for(ip=-1;ip<=1;ip++)
          for(jp=-1;jp<=1;jp++)
            goodback[jp+1][ip+1]=(float) (image_ivar[(j+jp)*xnpix+(i+ip)]>0.);
        if((goodback[1][0]+goodback[1][2])>0) {
          back[0]=(image[j*xnpix+(i-1)]*goodback[1][0]
                   +image[j*xnpix+(i+1)]*goodback[1][2])/
            (goodback[1][0]+goodback[1][2]);
          if(ival<back[0]*invsigma+nsig) continue;
        } /* end if */
        if((goodback[0][1]+goodback[2][1])>0) {
          back[1]=(image[(j-1)*xnpix+i]*goodback[0][1]+
                   image[(j+1)*xnpix+i]*goodback[2][1])/
            (goodback[0][1]+goodback[2][1]);
          if(ival<back[1]*invsigma+nsig) continue;
        } /* end if */
        if((goodback[0][0]+goodback[2][2])>0) {
          back[2]=(image[(j-1)*xnpix+(i-1)]*goodback[0][0]+
                   image[(j+1)*xnpix+(i+1)]*goodback[2][2])/
            (goodback[0][0]+goodback[2][2]);
          if(ival<back[2]*invsigma+nsig) continue;
        } /* end if */
        if((goodback[2][0]+goodback[0][2])>0) {
          back[3]=(image[(j+1)*xnpix+(i-1)]*goodback[2][0]+
                   image[(j-1)*xnpix+(i+1)]*goodback[0][2])/
            (goodback[2][0]+goodback[0][2]);
          if(ival<back[3]*invsigma+nsig) continue;
          ival=invsigma*imcurr;
        } /* end if */

        /* if it does, now check if ANY pair violates PSF conditions */
        sigmaback[0]=
          sqrt((goodback[1][0]/(image_ivar[j*xnpix+(i-1)]+1.-goodback[1][0])+ 
                goodback[1][2]/(image_ivar[j*xnpix+(i+1)]+1.-goodback[1][2]))/
               (goodback[1][0]+goodback[1][2]));
        lval=(ival-cfudge)*c2fudge*psfvals[0];
        if(lval>invsigma*(back[0]+cfudge*sigmaback[0])) {
          rejected[j*xnpix+i]=1;
          image[j*xnpix+i]=back[0];
          continue;
        }
        sigmaback[1]=
          sqrt((goodback[0][1]/(image_ivar[(j-1)*xnpix+i]+1.-goodback[0][1])+ 
                goodback[2][1]/(image_ivar[(j+1)*xnpix+i]+1.-goodback[2][1]))/
               (goodback[0][1]+goodback[2][1]));
        if(lval>invsigma*(back[1]+cfudge*sigmaback[1])) {
          rejected[j*xnpix+i]=1;
          image[j*xnpix+i]=back[1];
          continue;
        }
        sigmaback[2]=
          sqrt((goodback[0][0]/(image_ivar[(j-1)*xnpix+(i-1)]+1.- 
                                goodback[0][0])+ 
                goodback[2][2]/(image_ivar[(j+1)*xnpix+(i+1)]+1.- 
                                goodback[2][2]))/
               (goodback[0][0]+goodback[2][2]));
        lval=(ival-cfudge)*c2fudge*psfvals[1];
        if(lval>invsigma*(back[2]+cfudge*sigmaback[2])) {
          rejected[j*xnpix+i]=1;
          image[j*xnpix+i]=back[2];
          continue;
        }
        sigmaback[3]=
          sqrt((goodback[2][0]/(image_ivar[(j+1)*xnpix+(i-1)]+1.- 
                                goodback[2][0])+ 
                goodback[0][2]/(image_ivar[(j-1)*xnpix+(i+1)]+1.- 
                                goodback[0][2]))/
               (goodback[2][0]+goodback[0][2]));
        if(lval>invsigma*(back[3]+cfudge*sigmaback[3])) {
          rejected[j*xnpix+i]=1;
          image[j*xnpix+i]=back[3];
          continue;
        }
        
      } /* end if */
    } /* end for j */
  } /* end for i */
  
  FREEVEC(sigmaback);
  FREEVEC(back);
  return(0);
  
} /* end photfrac */
