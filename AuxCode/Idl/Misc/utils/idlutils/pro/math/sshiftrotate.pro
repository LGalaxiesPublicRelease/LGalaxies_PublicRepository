;+
; NAME:
;   sshiftrotate
;
; PURPOSE:
;   Routine to sinc-shift and rotate a 2-D image, but preserving scale.
;
; CALLING SEQUENCE:
;   newimg =  sshiftrotate( image, [ theta, xshift=, yshift=, xcen=, ycen=, $
;    /bigger, xoffset=, yoffset= ]
;
; INPUTS:
;   image      - Image (2-dimensional)
;
; OPTIONAL KEYWORDS:
;   theta      - Rotate image counter-clockwise this angle [degrees] about
;                the 0-indexed point XCEN,YCEN; default to 0 degrees
;   xshift    - Shift in X direction
;   yshift    - Shift in Y direction
;   xcen       - Center X position for rotation; default to the center of
;                the image
;   ycen       - Center Y position for rotation; default to the center of
;                the image
;   bigger     - If set, then keep the bigger image necessary for containing
;                the shifted + rotated image.
;
; OUTPUTS:
;   newimg     - Rotated and shifted image
;
; OPTIONAL OUTPUTS:
;   xoffset    - If /BIGGER is set, then this contains the integer pixel
;                offset in the X direction of the enlarged image.
;   yoffset    - If /BIGGER is set, then this contains the integer pixel
;                offset in the Y direction of the enlarged image.
;
; COMMENTS:
;   When both a shifT (XSHIFT,YSHIFT) and a rotation (THETA) are specified,
;   the resulting image is as if the shift is performed first, and the
;   rotation second.
;
; EXAMPLES:
;   Generate a random image and rotate by 30 degrees:
;     IDL> image = smooth(randomu(1234,200,200),5)
;     IDL> newimg = sshiftrotate(image,30)
;
; BUGS:
;   The sinc shifts need not do all pixels in each row each time, only
;     the "active" area!???  This will just be for a speed improvement.
;   Special-case rotations of 0,90,180,270 !???
;   Optionally return a mask of the illuminated region???
;   Optionally fill in missing regions with some value???
;   Allow double-precision for the image, or for the arithmatic???
;
; PROCEDURES CALLED:
;   cirrange
;   sshift()
;
; REVISION HISTORY:
;   18-Sep-2002  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function sshiftrotate, image, theta1, xshift=xshift1, yshift=yshift1, $
 xcen=xcen1, ycen=ycen1, bigger=bigger1, xoffset=xoffset, yoffset=yoffset

   if (n_params() LT 1) then $
    message, 'Incorrect number of arguments'
   if (keyword_set(xshift1)) then xshift = xshift1 $
    else xshift = 0
   if (keyword_set(yshift1)) then yshift = yshift1 $
    else yshift = 0
   dims = size(image, /dimens)
   nx = dims[0]
   ny = dims[1]
   if (n_elements(xcen1) NE 0) then xcen = xcen1 $
    else xcen = 0.5 * nx + 0.5
   if (n_elements(ycen1) NE 0) then ycen = ycen1 $
    else  ycen = 0.5 * ny + 0.5
   if (keyword_set(theta1)) then theta = theta1 $
    else theta = 0
   ; Put this angle in the range [-45,+45).
   theta = theta - 360 * floor((theta+45)/360)
   if (n_elements(bigger1) NE 0) then bigger = bigger1 $
    else bigger = 0

   ;----------
   ; If the rotation angle is not in the range [-45,+45) degrees, then rotate
   ; the image first and call this routine recursively.

   irot = long( (theta+45) / 90 )
   if (irot NE 0) then begin
      thisang = theta - irot * 90

      ; Rotate the image by a multiple of 90 degrees.
      thisimg = rotate(image, irot)

      xycen = [xcen, ny-ycen, nx-xcen, ycen]
      xnewcen = xycen[irot]
      ynewcen = xycen[irot-1]
      xoffset = xnewcen - xcen
      yoffset = ynewcen - ycen

      xyshift = [xshift, -yshift, -xshift, yshift]
      xnewshift = xyshift[irot]
      ynewshift = xyshift[irot-1]

      thisimg = sshiftrotate(thisimg, thisang, $
       xshift=xnewshift, yshift=ynewshift, $
       xcen=xnewcen, ycen=ynewcen, bigger=1, $
       xoffset=xoff1, yoffset=yoff1)

      xoffset = xoffset + xoff1 + xnewshift - xshift
      yoffset = yoffset + yoff1 + ynewshift - yshift

      ; Trim image to original size if /BIGGER not set.
      if (NOT keyword_set(bigger)) then begin
         thisdim = size(thisimg, /dimens)
         trimimg = make_array(size=size(image)) ; Same var type as for IMAGE
         ix = -xoffset * (xoffset LT 0)
         jx = xoffset * (xoffset GT 0)
         iy = -yoffset * (yoffset LT 0)
         jy = yoffset * (yoffset GT 0)
         xsz = (nx - ix) < (thisdim[0] - jx)
         ysz = (ny - iy) < (thisdim[1] - jy)
         if (xsz GT 0 AND ysz GT 0) then $
          trimimg[ix:ix+xsz-1,iy:iy+ysz-1] = thisimg[jx:jx+xsz-1,jy:jy+ysz-1]
         xoffset = 0
         yoffset = 0
         return, trimimg
      endif

      return, thisimg
   endif

   ;----------
   ; Compute the offset functions for each of the 3 sinc shifts

   sint = sin(theta/!radeg)
   cost = cos(theta/!radeg)
   aslope = -sint / (1. + cost)
   bslope = sint
   cslope = aslope
   ssx = xshift - cslope * yshift
   ssy = -bslope * xshift + (1+bslope*cslope) * yshift

   ;----------
   ; Compute the size of the super-image for containing the shifts

   ypad1 = -min( bslope * ([0,(nx-1.)] - xcen), max=ypad2)
   xpad1 = -min( aslope * ([0,(ny-1.)] - ycen) $
               + aslope * ([-ypad1,(ny-1.)+ypad2] - ycen), max=xpad2)

   xpad1 = xpad1 + ((-xshift)>0)
   xpad2 = xpad2 + ((xshift)>0)
   ypad1 = ypad1 + ((-yshift)>0)
   ypad2 = ypad2 + ((yshift)>0)

   xpad1 = ceil(xpad1) > 0L
   xpad2 = ceil(xpad2) > 0L
   ypad1 = ceil(ypad1) > 0L
   ypad2 = ceil(ypad2) > 0L

   ;----------
   ; Construct the output image, padding by the necessary amount in
   ; each dimension

   nbigx = nx + xpad1 + xpad2
   nbigy = ny + ypad1 + ypad2
;   print, 'Resize image from ', nx, ny, ' to ', nbigx, nbigy

   newimg = fltarr(nbigx, nbigy) ; What if this is double???
   newimg[xpad1:xpad1+nx-1,ypad1:ypad1+ny-1] = image

   xvec = findgen(nbigx) - xcen - xpad1
   yvec = findgen(nbigy) - ycen - ypad1

   ;----------
   ; Do the sinc shifts

   for iy=0,nbigy-1 do $
    if (abs(yvec[iy]*aslope+ssx) LT 0.95*nbigx) then $
     newimg[*,iy] = sshift(newimg[*,iy],yvec[iy]*aslope+ssx)

   for ix=0,nbigx-1 do $
    if (abs(xvec[ix]*bslope+ssy) LT 0.95*nbigy) then $
     newimg[ix,*] = sshift(newimg[ix,*],xvec[ix]*bslope+ssy)

   for iy=0,nbigy-1 do $
    if (abs(yvec[iy]*cslope) LT 0.95*nbigx) then $
     newimg[*,iy] = sshift(newimg[*,iy],yvec[iy]*cslope)

   ;----------
   ; Trim the output image

   if (keyword_set(bigger)) then begin
      xoffset = xpad1
      yoffset = ypad1
   endif else begin
      xoffset = 0
      yoffset = 0
      newimg = newimg[xpad1:xpad1+nx-1,ypad1:ypad1+ny-1]
   endelse

   return, newimg
end
;------------------------------------------------------------------------------
