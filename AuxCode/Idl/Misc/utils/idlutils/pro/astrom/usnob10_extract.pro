;+
; NAME:
;   usnob10_extract
;
; PURPOSE:
;   Extract information from USNO-B1.0 binary data files
;
; CALLING SEQUENCE:
;   a = usnob10_extract(data)
;   
; INPUTS:
;   data    - ulong array [20, Nstar] read from binary file
;
; OUTPUTS:
;   a       - IDL data structure with positions, proper motions, and 
;                up to five magnitudes. 
;   
; COMMENTS:
;   See http://ftp.nofs.navy.mil/data/fchpix/usnob_format.html
;     for descriptions of the fields
;
; REVISION HISTORY:
;   2002-Nov-26   Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function create_usnob10struct, n, simple=simple

; Similar to the ASCII tables returned by database queries

  if keyword_set(simple) then begin 

     ftemp = create_struct( $
                            'ra',  0.0D, $   ; [deg] (J2000)
                            'dec', 0.0D, $
                            'sra', 0, $      ; [mas]
                            'sde', 0, $
                            'epoch', 0.0, $  ; [yr]
                            'mura', 0, $     ; [mas/yr]
                            'mudec', 0, $
                            'muprob', 0B, $
                            'muflag', 0B, $
                            'smura', 0, $    ; [mas/yr]
                            'smude', 0, $
                            'sFitRA', 0, $   ; [mas]
                            'sFitDE', 0, $
                            'NFitPt', 0B, $
                            'B1', 0.0, $     ; [mag]
                            'R1', 0.0, $
                            'B2', 0.0, $
                            'R2', 0.0, $
                            'I2', 0.0)
  endif else begin 
     
     ftemp = create_struct( $
                            'ra',  0.0D, $   ; [deg] (J2000)
                            'dec', 0.0D, $
                            'sra', 0, $      ; [mas]
                            'sde', 0, $
                            'epoch', 0.0, $  ; [yr]
                            'mura', 0, $     ; [mas/yr]
                            'mudec', 0, $
                            'muprob', 0B, $
                            'muflag', 0B, $
                            'smura', 0, $    ; [mas/yr]
                            'smude', 0, $
                            'sFitRA', 0, $   ; [mas]
                            'sFitDE', 0, $
                            'NFitPt', 0B, $
                            'mag', fltarr(5), $     ; [mag] B1,R1,B2,R2,I2
                            'FldID', intarr(5), $   ; 4-digit field ID
                            'SG', bytarr(5), $      ; star/galaxy sep (0-11)
                            'Fldepoch', fltarr(5) ) ; field epoch
  endelse 

  usnostruct = replicate(ftemp, n)
  
  return, usnostruct
end



pro usnob10_extract_mag, data, ind, str

  filter = ['B1', 'R1', 'B2', 'R2', 'I2']
  arr = transpose(data[5+ind, *])
  mag   = (arr mod 10000L) * 0.01
  magind = where(tag_names(str) eq filter[ind], simple)

  if simple then begin 
     str.(magind) = mag
     return
  endif 

  str.mag[ind]   = mag
  str.fldid[ind] = (arr mod 100000000L) / 10000L
  str.sg[ind]    = (arr / 100000000L)
  return
end



function usnob10_extract, data, simple=simple

; -------- create empty structure  
  if size(data, /n_dim) EQ 1 then n = 1 else $
    n = (size(data,/dimens))[1]
  a = create_usnob10struct(n, simple=simple)

; -------- fill structure

  a.ra  = transpose(data[0, *]) /3.6d5
  a.dec = transpose(data[1, *]) /3.6d5 - 90.d
  
  mu = transpose(data[2, *])
  a.mura   = ((mu mod 10000L) - 5000) * 2
  a.mudec  = ((mu mod 100000000L)/10000L - 5000) * 2
  a.muprob = (mu mod 1000000000L) / 100000000L
  a.muflag = mu / 1000000000L

  sig = transpose(data[3, *])
  a.smura  = (sig mod 1000L)
  a.smude  = (sig mod 1000000L)/1000L
  a.sfitra = (sig mod 10000000L)/1000000L *100
  a.sfitde = (sig mod 100000000L)/10000000L *100
  a.nfitpt = (sig mod 1000000000L)/100000000L

  sig = transpose(data[4, *])
  a.sra    = (sig mod 1000L)
  a.sde    = (sig mod 1000000L)/1000L
  a.epoch  = (sig mod 1000000000L)/1000000L *0.1 + 1950

  for ind=0, 4 do usnob10_extract_mag, data, ind, a

  return, a
end
