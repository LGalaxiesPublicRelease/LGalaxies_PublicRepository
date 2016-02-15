;------------------------------------------------------------------------------
;+
; NAME:
;   djs_neff
;
; PURPOSE:
;   Return the neff statistic (effective number of pixels)
;     neff = [ SUM(image) ]^2 / SUM(image^2)
;
; CALLING SEQUENCE:
;   neff = djs_neff(image, [sigimg, nerr=nerr] )
;
; INPUTS:
;   image:      An image of any number of dimensions
;
; OPTIONAL INPUTS:
;   sigimg:     Image of errors for computing nerr
;
; OUTPUTS:
;   neff:       Return value of neff
;
; OPTIONAL OUTPUTS:
;   nerr:       Error in neff, if sigimg is specified
;
; PROCEDURES CALLED:
;
; COMMENTS:
;   If computing this statistic for an object on an image, the background
;   (sky) level should first be removed.  However, the image need not be
;   renormalized.  For an object of constant surface brightness, Neff equals
;   the number of pixels subtended by the object.  For a Gaussian profile
;   that is well-sampled, Neff = 4 * pi * sigma^2.
;
; REVISION HISTORY:
;   10-Sep-1998  Written by D Schlegel, Princeton
;-
;------------------------------------------------------------------------------
function djs_neff, image, sigimg, nerr=nerr
 
   ; Need 1 parameter
   if N_params() LT 1 then begin
      print, 'Syntax - neff = djs_neff( image, [ sigimg, nerr=nerr ] )'
      return, -1
   endif

   sumx = total(image, /double)
   sumxx = total(image^2, /double)
   neff = sumx^2 / sumxx

   if (keyword_set(sigimg)) then begin
      sumss = total(sigimg^2, /double)
      sumxss = total(image*sigimg^2, /double)
      sumxxss = total(image^2*sigimg^2, /double)
      nerr = sqrt( 2*neff * $
       ( sumss/(sumx^2) - 2*sumxss/(sumx*sumxx) + sumxxss/(sumxx^2) ) )
   endif

   return, neff
end
;------------------------------------------------------------------------------
