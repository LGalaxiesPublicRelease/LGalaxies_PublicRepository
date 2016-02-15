;+
; Returns the number of healpix pixels
; for a given resolution.
;
; npix = healnpix(res,[/nside])
;
; If nside is set, then it returns nside
;
; Nikhil Padmanabhan, Princeton, 
; July 29,2003
;-

function healnpix, res, nside=nside
   if (n_elements(res) EQ 0) then $
	message,'ERROR : Set resolution'

   side = 2L^res
  
   if (keyword_set(nside)) then $
	return, side $
   else begin
     npix = 12L*side*side
     return, npix
   endelse

end