;+
; NAME:
;   iau_to_radec
; PURPOSE:
;   convert an IAU name to RA and DEC
; CALLING SEQUENCE:
;   iau_to_radec, name, ra, dec
; INPUTS:
;   name - IAU name
; OUTPUTS:
;   ra, dec - coords
; REVISION HISTORY:
;   2005-10-24  Blanton (NYU)
;-
pro iau_to_radec, name, ra, dec 

numbers=(stregex(name, '.*J([0-9|\.]*)([+|-][0-9|\.]*).*', /extr, /sub))

rahr=strmid(numbers[1,*], 0, 2)
ramin=strmid(numbers[1,*], 2, 2)
rasec=strmid(numbers[1,*], 4)

decdeg=strmid(numbers[2,*], 0, 3)
decmin=strmid(numbers[2,*], 3, 2)
decsec=strmid(numbers[2,*], 5)

string2radec, rahr, ramin, rasec, decdeg, decmin, decsec, ra, dec

ra=reform(ra, n_elements(ra))
dec=reform(dec, n_elements(dec))

end
