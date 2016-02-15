;------------------------------------------------------------------------------
;+
; NAME:
;   quick_photfrac
;
; PURPOSE:
;   Create a list of pixel numbers and their fractional contribution to
;   an annular region.
;
; CALLING SEQUENCE:
;   quick_photfrac, xcen, ycen, Rvec, xdimen=, ydimen=, /ragged, $
;    [ xPixNum=, yPixNum=, pixnum=, fracs=, fillfrac= ]
;
; INPUTS:
;   xcen:       X center(s)
;   ycen:       Y center(s)
;   Rvec:       Either a 2-element array with two radii to define an annulus,
;               or a scalar to define a circular aperature.
;
; OPTIONAL INPUTS:
;   xdimen:     Number of X pixels.
;   ydimen:     Number of Y pixels.
;   ragged:     Use ragged edges (weights either 0 or 1) - faster
;
; OUTPUTS:
;   pixnum:     Pixel number, 0-indexed, for referencing array using one index.
;   xPixNum:    Pixel number in X, 0-indexed.
;   yPixNum:    Pixel number in Y, 0-indexed.
;   fracs:      Return value of covering fraction of the annulus
;               over the pixel number.
;   fillfrac:   Ratio of returned pixel areas to the annulus area;
;               this ratio should equal 1.0 if the aperature falls completely
;               within the image boundaries
;
; COMMENTS:
;   The total counts within this region is given by
;     totcounts = total( pData(pixnum) * fracs )
;   The area within this region is given by
;     area = total(fracs)
;   The average counts is given by
;     totcounts = total( pData(pixnum) * fracs ) / total(fracs)
;
;   If no pixels within the given annulus are found, then return pixnum=-1.
;
; BUGS:
;   The area can be calculated with TOTAL(FRACS), and will differ
;   slightly from the analytic area within a circle.
;
; PROCEDURES CALLED:
;   djs_ceil()
;   djs_floor()
;
; REVISION HISTORY:
;   Written by D. Finkbeiner, 2000-Nov-02
;  derived from djs_photfrac
;   Written D. Schlegel, 27 November 1996, Durham
;   
;-
pro quick_photfrac, xcen, ycen, Rvec, xdimen=xdimen, ydimen=ydimen, $
 ragged=ragged, xPixNum=xPixNum, yPixNum=yPixNum, pixnum=pixnum, $
 fracs=fracs, fillfrac=fillfrac


; Set return values in the event of an error
  pixnum = -1L
  fracs = 0
  
; Need 2 parameters
  if N_params() LT 2 then begin
     print, 'Syntax - quick_photfrac, xcen, ycen, Rvec, $'
     print, ' xdimen=, ydimen=, xPixNum=, yPixNum=, pixnum=, $'
     print, ' fracs=, fillfrac='
     return
  endif

; If Rvec contains one element, then use the annulus [0,Rvec],
; otherwise use the annulus [Rvec[0],Rvec[1]].
  
  nr = N_elements(Rvec)
  rad = (nr EQ 1) ? [0, rvec] : rvec
   
  iStart = long(0 > djs_floor(xcen + 0.5 - rad[1]))
  jStart = long(0 > djs_floor(ycen + 0.5 - rad[1]))
  
  sz = long(Rad[1]*2+2)

  if (NOT keyword_set(xbox)) OR (NOT keyword_set(ybox)) then begin 
     xbox = lindgen(sz, sz) mod sz
     ybox = lindgen(sz, sz)  /  sz
  endif
  
  xoff = xcen-istart
  yoff = ycen-jstart
  r2 = (xbox-xoff)^2 + (ybox-yoff)^2
  if keyword_set(ragged) then begin 
     if nr EQ 1 then mask = (r2 LT rad[1]^2) ELSE  $
       mask = (r2 LT rad[1]^2) AND (r2 GT rad[0]^2)
     pix1 = where(mask, ct)
     if ct NE 0 then begin 
        xPixNum = xbox[pix1] + istart
        yPixNum = ybox[pix1] + jstart
        fracs   = fltarr(ct)+1.
     endif else return
  endif else begin 
     r = sqrt(r2)
     an = ((Rad[1]+0.5)-r) < 1.
     
     if nr EQ 2 then an = an-(((Rad[0]+0.5)-r) > 0.)
     
; Set the return values
  
; Limit the return values to only those pixels with non-zero contributions
     pix1    = where(an GT 0.0, ct)
     if ct NE 0 then begin 
        xPixNum = xbox[pix1] + istart
        yPixNum = ybox[pix1] + jstart
        fracs   = an[pix1]
     endif else return
  endelse 

  if keyword_set(xdimen) AND keyword_set(ydimen) then begin 
     pix2 = where((xPixNum LT xdimen) AND (yPixNum LT ydimen), ct)
     if ct NE n_elements(xPixNum) then begin 
        if ct NE 0 then begin 
           xPixNum = xPixNum[pix2]
           yPixNum = yPixNum[pix2]
           fracs   = fracs[pix2]
        endif else return
     endif 
     if arg_present(pixnum) then pixnum  = 0L + xPixNum + xdimen * yPixNum
  endif 
; Test to see if aperature exceeds image boundary by computing the 
; ratio of the filled pixels to the area of the annulus.  If all
; pixels are within the image boundary, then fillfrac=1.0.
  if arg_present(fillfrac) then $
    fillfrac = total(fracs) / (!pi * (rad[1]^2 - rad[0]^2))
  
  return
end

