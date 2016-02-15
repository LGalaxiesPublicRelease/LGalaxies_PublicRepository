;------------------------------------------------------------------------------
;+
; NAME:
;   squarephot
;
; PURPOSE:
;   Simple (square) aperature photometry.
;
; COMMENTS:
;   This routine has all the same calling arguments as DJS_PHOT(), but
;   does not use most of them.  The sky is computed as the median of the
;   entire image.  The flux is computed from square aperatures.  Nothing
;   else is computed.  This proc is basically for testing purposes.
;
; CALLING SEQUENCE:
;   03-Jun-2000  Written by D. Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function squarephot, xcen, ycen, objrad, skyrad, image, $
 calg=calg, cbox=cbox, cmaxiter=cmaxiter, cmaxshift=cmaxshift, $
 fwhm=fwhm, fixfw=fixfw, ceps=ceps, $
 salg=salg, srejalg=srejalg, smaxiter=smaxiter, $
 lorej=lorej, hirej=hirej, $
 flerr=flerr, skyval=skyval, skyrms=skyrms, skyerr=skyerr

   nobj = n_elements(xcen)

   dims = size(image, /dimens)
   xdimen = dims[0]
   ydimen = dims[1]
   nrad = n_elements(objrad)

   ; Allocate memory for return values
   flux = fltarr(nobj,nrad)
   flerr = fltarr(nobj,nrad)
   skyval = fltarr(nobj)
   skyrms = fltarr(nobj)
   skyerr = fltarr(nobj)

   ; Compute the sky from the entire image
   skyval[*] = median(image)

   ;----- LOOP THROUGH EACH OBJECT -----
   for iobj=0L, nobj-1 do begin

      x1 = round(xcen[iobj] - objrad) > 0
      x2 = round(xcen[iobj] + objrad) < xdimen-1
      y1 = round(ycen[iobj] - objrad) > 0
      y2 = round(ycen[iobj] + objrad) < ydimen-1
      area = (x2-x1) * (y2-y1)
      for irad=0, n_elements(objrad)-1 do begin
         flux[iobj,irad] = total( image[x1[irad]:x2[irad],y1[irad]:y2[irad]] ) $
          - skyval[iobj] * area[irad]
      endfor
   endfor

   return, flux
end
;------------------------------------------------------------------------------
