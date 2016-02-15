;------------------------------------------------------------------------------
;+
; NAME:
;   exact_photfrac
;
; PURPOSE:
;   Create list of contribution of pixels to circular or annular aperture
;
; CALLING SEQUENCE:
;   exact_photfrac, xcen, ycen, radius [, fracs=, xdimen=, ydimen=, ]
;           pixnum=, xpixnum=, ypixnum= ]
;
; INPUTS:
;   xcen - X center(s) (LONG) 
;   ycen - Y center(s) (LONG) 
;   radius - radius of aperture (if 2-element array, inner and outer
;            radii of annulus) 
;
; OPTIONAL INPUTS:
;   xdimen - number of pixels upon a side to output [default - radius+1]
;   ydimen - number of pixels upon a side to output [default - radius+1]
;   safefactor - we set strictly to zero all pixels outside
;                max(radius)*safefactor [default 1.2]
;
; OUTPUTS:
;   fracs- contribution of each pixel to image
;   pixnum - Pixel number, 0-indexed, for referencing array using one index.
;   xPixNum - Pixel number in X, 0-indexed.
;   yPixNum - Pixel number in Y, 0-indexed.
;
; COMMENTS:
;   Uses Robert Lupton's Aperture Photometry scheme to measure seeing-
;   and pixel-convolved aperture photometry in band-limited images.
;
;   xcen and ycen MUST be integers. This simplifies the caching 
;   of the weights considerably. Note that for band-limited images 
;   (the only kind that this code works for) you can always sshift 
;   the image to get the center of the object at the center of a 
;   pixel (ie. an integer pixel number). 
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Started - 22-Aug-2003 M. Blanton (NYU)
;-
;------------------------------------------------------------------------------
; utility to cache exact_photfrac values used many times
function get_exact_photfrac, nx, ny, radius, xc, yc

common com_exact_photfrac, exact_fracs, soname

; Set source object name
if(NOT keyword_set(soname)) then $
  soname=filepath('libimage.'+idlutils_so_ext(), $
                  root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')

for i=0L, n_elements(exact_fracs)-1L do $
  if(nx eq exact_fracs[i].nx and $
     ny eq exact_fracs[i].ny and $
     radius eq exact_fracs[i].radius and $
     xc eq exact_fracs[i].xc and $
     yc eq exact_fracs[i].yc) then $
  return, *exact_fracs[i].fracs

fracs=fltarr(nx,ny)
retval=call_external(soname, 'idl_photfrac', long(nx), long(ny), $
                     float(radius),float(fracs),long(xc),long(yc))

new_exact_fracs={nx:nx, $
                 ny:ny, $
                 radius:radius, $
                 fracs:ptr_new(fracs), $
                 xc:xc, $
                 yc:yc}
if(n_tags(exact_fracs) gt 0) then $
  exact_fracs=[exact_fracs,new_exact_fracs] $
else $
  exact_fracs=new_exact_fracs

return, fracs

end
;
pro exact_photfrac, xcen, ycen, radius, fracs=fracs, ydimen=ydimen, $
                    xdimen=xdimen, pixnum=pixnum, xpixnum=xpixnum, $
                    ypixnum=ypixnum, safefactor=safefactor

if(n_params() lt 2) then begin
    print, 'Syntax - exact_photfrac, xcen, ycen, radius [, fracs=, xdimen=, ydimen=, '
    print, '            pixnum=, xpixnum=, ypixnum= ] '
    return
endif

if(n_elements(xdimen) eq 0) then xdimen=long(radius)+2L
if(n_elements(ydimen) eq 0) then ydimen=xdimen
if(n_elements(safefactor) eq 0) then safefactor=2.
if(n_elements(xcen) eq 0) then xcen=0L
if(n_elements(ycen) eq 0) then ycen=0L

ini=where(long(xcen) ne xcen OR long(ycen) ne ycen, nni)
if(nni gt 0) then begin
    message, 'xcen and ycen MUST be integers' 
endif

; define region to cut out
if(keyword_set(safefactor)) then begin
    safedistance=(long(max(radius)*safefactor))>10L
    xstart=(xcen-safedistance)>0L
    xend=(xcen+safedistance)<(xdimen-1L)
    ystart=(ycen-safedistance)>0L
    yend=(ycen+safedistance)<(ydimen-1L)
endif else begin
    xstart=0
    xend=xdimen-1L
    ystart=0
    yend=ydimen-1L
endelse 
nx=xend-xstart+1L
ny=yend-ystart+1L
xc=xcen-xstart
yc=ycen-ystart

fracs=get_exact_photfrac(nx,ny,radius[0],xc,yc)
if(n_elements(radius) gt 1) then begin
    fracs1=get_exact_photfrac(nx,ny,radius[1],xc,yc)
    fracs=fracs1-fracs
endif

if(arg_present(xpixnum) OR arg_present(pixnum)) then $
  xpixnum=(xstart+lindgen(nx))#replicate(1L,ny)
if(arg_present(ypixnum) OR arg_present(pixnum)) then $
  ypixnum=replicate(1L,nx)#(ystart+lindgen(ny))
if(arg_present(pixnum)) then $
  pixnum= ypixnum*xdimen+xpixnum

end
;------------------------------------------------------------------------------
