;+
; NAME:
;   read_ds_cooling
;
; PURPOSE:
;   Read in Dopita & Sutherland 1993 cooling function
;
; CALLING SEQUENCE:
;   read_ds_cooling, fname, logT, loglambda
;
; INPUTS:
;   fname - one of
;            m-00.cie
;            m-05.cie
;            m+05.cie
;            m-10.cie
;            m-15.cie
;            m-20.cie
;            m-30.cie
;            mzero.cie
;   logT  - log10 T in Kelvin
;   loglambda - log of cooling function [erg/s cm^3]
; 
; OPTIONAL OUTPUTS:
;   logT  - if not passed, then logT is set to values in file
;
; RESTRICTIONS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   
; REVISION HISTORY:
;   2007-Mar-19  Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
pro read_ds_cooling, fname, logT, loglambda

  idlutils_dir = getenv('IDLUTILS_DIR')
  path = concat_dir(idlutils_dir, 'data/cooling')
  if NOT file_test(path, /directory) then message, 'Cannot find '+path

  fullname = concat_dir(path, fname)
  if NOT file_test(fullname) then message, 'Cannot find '+fullname

  readcol, fullname, logTin, nelec, nH, nt, loglambdanet, loglambdanorm, $
    skip=4

  if keyword_set(logT) then begin 
     if max(logT) GT max(logTin) then message, 'logT values too large!', /info
     if min(logT) LT min(logTin) then message, 'logT values too small!', /info
     loglambda = interpol(loglambdanet, logTin, logT)
  endif else begin 
     logT = logTin
     loglambda = loglambdanet
  endelse
  
  return
end
