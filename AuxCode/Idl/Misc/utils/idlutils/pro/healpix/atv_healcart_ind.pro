;+
; NAME:
;   atv_healcart_ind
;
; PURPOSE:
;   Generate and cache healcart pixel index array and header for atv
;
; CALLING SEQUENCE:
;   ind = atv_healcart_ind(image, nest=, header=)
;
; INPUTS:
;   image   - healpix array
;
; KEYWORDS:
;   nest    - indicate nested ordering
;
; OUTPUTS:
;   ind     - index array
;
; OPTIONAL OUTPUTS:
;   header  - mock FITS header for nonstandard HCT projection
;
; COMMENTS:
;   Used by atv. 
;   IDLUTILS version of wcsxy2sph.pro recognizes this projection, even
;    though it is NOT STANDARD FITS!
;   
; REVISION HISTORY:
;   2003-May-10  Written by Douglas Finkbeiner, Princeton
;----------------------------------------------------------------------

function atv_healcart_ind, image, nest=nest, header=header

; -------- save index array in common block  
  common atv_healpix, nside_save, ind_save, nest_save

  if NOT keyword_set(nest) then nest = 0B

; -------- check size of inputs
  if NOT keyword_set(image) then message, 'must set image'
  npix = n_elements(image) 
  nside = long(sqrt(npix/12))

; -------- see if we have index array already
  if keyword_set(ind_save) then begin 
     if (nside eq nside_save) and (nest_save eq nest) then $
       ind = ind_save
  endif 

  if NOT keyword_set(ind) then begin  ; otherwise generate index array
     ind = healcart_ind(image, nest=nest)
     ind_save = ind
     nside_save = nside
     nest_save = nest
  endif

; -------- generate header
  if arg_present(header) then begin 
     header = healcart_header(ind)
  endif 
  
  return, ind
end
