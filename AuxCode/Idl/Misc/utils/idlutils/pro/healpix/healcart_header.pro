;+
; NAME:
;   healcart_header
;
; PURPOSE:
;   Generate mock header for healcart image
;
; CALLING SEQUENCE:
;   header=healcart_header([ind, nside=nside])
;
; OPTIONAL INPUTS:
;   ind      - healcart index array, used only for size information
;
; KEYWORDS:
;   nside    - used for size information
;
; OUTPUTS:
;   header   - FITS header appropriate for atv
;
; RESTRICTIONS:
;   Must be used with idlutils version of atv and xy2ad.
;
; EXAMPLES:
;   
; COMMENTS:
;   Defines header for the approximately Cartesian reprojection of the 
;    healpix map defined by healcart_ind.pro.
;    atv knows how to use this header, even though it is not standard
;    FITS. 
;
;   Either nside or ind MUST be set. 
;
; REVISION HISTORY:
;   2003-May-10  Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function healcart_header, ind, nside=nside

  if NOT keyword_set(nest) then nest = 0B

; -------- check size of inputs
  if n_elements(ind) EQ 0 then begin 
     if keyword_set(nside) then begin 
        sz = [8L*nside, 4L*nside-1]
        mkhdr, head, fltarr(sz[0], sz[1], /nozero)
     endif else begin 
        print, 'HEALCART_HEADER: Must set set either ind or Nside'
        return, 0
     endelse
  endif else begin 
     mkhdr, head, ind
     sz = size(ind, /dimen)
  endelse 

  nlon = sz[0]
  nlat = sz[1]

  cdelt = 360.d/nlon

; X-axis
  sxaddpar, head, 'CTYPE1', 'GLON-HCT', /pdu
  sxaddpar, head, 'CRPIX1', nlon/2.+0.5, /pdu
  sxaddpar, head, 'CRVAL1', 0.0, /pdu
  
; Y-axis
  sxaddpar, head, 'CTYPE2', 'GLAT-HCT', /pdu
  sxaddpar, head, 'CRPIX2', nlat/2.+0.5, /pdu
  sxaddpar, head, 'CRVAL2', 0.0, /pdu
  
; Rotation/scaling matrix (CD)
  sxaddpar, head, 'CD1_1', -1*cdelt, /pdu 
  sxaddpar, head, 'CD1_2', 0.d, /pdu
  sxaddpar, head, 'CD2_1', 0.d, /pdu
  sxaddpar, head, 'CD2_2', +1*cdelt, /pdu
  
  sxaddpar, head, 'LONPOLE', 180, /pdu
  sxaddpar, head, 'PROJP1', 0, /pdu
  sxaddpar, head, 'PROJP2', 0, /pdu
  
  sxaddpar, head, 'AUTHOR', 'Doug Finkbeiner', /pdu
  sxaddpar, head, 'COMMENT', ' Mock header for HEALCart projection --', /pdu
  sxaddpar, head, 'COMMENT', ' This header uses the non-standard HCT projection.', /pdu
  
  return, head
end
