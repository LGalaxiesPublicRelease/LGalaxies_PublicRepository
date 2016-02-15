;+
; NAME:
;   healpix_ring_weight
;
; PURPOSE:
;   read healpix ring weight files
;
; CALLING SEQUENCE:
;   wt = healpix_ring_weight(nside, iring=iring)
;
; INPUTS:
;   nside   - healpix nside
;   
; OUTPUTS:
;   wt = the ring weights (4*nside-1 element array)
;
; OPTIONAL OUTPUTS:
;   iring = ring number index array for each pixel (12*nside^2 array)
;   
; RESTRICTIONS:
;   breaks at nside > 8192
;
; BUGS:
;   The method of computing iring is really dumb and slow, but works. 
;   
; COMMENTS:
;   Gorski stores weight-1 as a float; we return weight as a double.
;   Because the weight array is symmetric, Gorski only stores half; we 
;     return the whole array for simplicity. 
;   The new mrdfits() (in v5_0_2b) replaces a " " in binary fits table
;   field names with "_" instead of removing it.  Current version of
;   this routine works with both new and old mrdfits().
;   
; REVISION HISTORY:
;   2003-Mar-11  Written by Douglas Finkbeiner, Princeton
;   2003-Nov-12  Cache outputs - DPF & NP
;   2004-Aug-09  Fix fatal bug with new mrdfits().
;
;----------------------------------------------------------------------
function healpix_ring_weight, nside, iring=iring

  common healpix_ring_weight_common, nside_save, iring_save, wt_save

  if keyword_set(nside_save) then begin 
     if nside_save eq nside then begin 
        if arg_present(iring) then iring = iring_save
        if NOT keyword_set(wt_save) then message, 'I am a bad person'
        return, wt_save+1.d
     endif 
  endif 

  healpix_dir  = getenv('HEALPIX_DIR')+'/'
  wt_ring_file = string(healpix_dir+'data/weight_ring_n', nside, '.fits', $
                        format='(A,I5.5,A)')

  if file_test(wt_ring_file) ne 1 then message, 'file not found' 
  wt_str = mrdfits(wt_ring_file, 1)
  tind = where(strpos(tag_names(wt_str),'TEMP') EQ 0, nind)
  if nind NE 1 then message, 'structure does not have a field with TEMPERATURE WEIGHTS'
  twt = wt_str.(tind)
  twt = reform(twt, n_elements(twt)) ; handle 1024 case
  wt = [twt, reverse(twt[0:n_elements(twt)-2] )]

  healgen, nside, theta, phi
  qind = uniq(theta)
  nq   =  n_elements(qind)
  npix = 12L*nside*nside
  iring = lindgen(npix)
  for q=0, nq-1 do begin 
     q0    = (q eq 0) ? 0 : qind[q-1]+1
     q1    = qind[q]
     iring[q0:q1] = q
  endfor

; -------- Cache outputs
  nside_save = nside
  iring_save = iring
  wt_save    = wt


  return, wt+1.d
end
