;+
; NAME:
;   hogg_iau_name
; PURPOSE:
;   properly format astronomical source names to the IAU convention
; INPUTS:
;   ra         - J2000 deg
;   dec        - J2000 deg
;   prefix     - survey prefix
; OPTIONAL INPUTS:
;   precision  - how many decimal places to put on dec; ra gets one more
; COMMENTS:
;   Defaults to SDSS name conventions.
;   Don't EVER use goddard/pro/astro/adstring.pro.
; BUGS:
;   Enforces J2000 "J".
;   Doesn't deal properly with precision<0 names, which *can* exist.
; REVISION HISTORY:
;   2004-10-24  started by Hogg (NYU)
;-
function hogg_iau_name, ra,dec,prefix,precision=precision
if (n_elements(prefix) eq 0) then prefix= 'SDSS'
if (n_elements(precision) EQ 0) then precision= 1

rah= floor(ra/15.0)
ram= floor(((ra/15.0)-double(rah))*60.0)
ras= (((ra/15.0)-double(rah))*60.0-double(ram))*60.0
ras= floor(ras*10.0^(precision+1))/10.0^(precision+1)
rasformat= '(F0'+strtrim(string(precision+4),2)+ $
  '.'+strtrim(string(precision+1),2)+')'

desgn= replicate('+',n_elements(dec))
negative= where(dec LT 0,nneg)
if (nneg GT 0) then desgn[negative]= '-'
adec= abs(dec)
ded= floor(adec)
dem= floor((adec-double(ded))*60.0)
des= ((adec-double(ded))*60.0-double(dem))*60.0
des= floor(des*10.0^precision)/10.0^precision
desformat= '(F0'+strtrim(string(precision+3),2)+ $
  '.'+strtrim(string(precision),2)+')'
if (precision EQ 0) then desformat= '(I2.2)'

adstr= string(rah,format='(I2.2)') $
  +string(ram,format='(I2.2)') $
  +string(ras,format=rasformat) $
  +desgn $
  +string(ded,format='(I2.2)') $
  +string(dem,format='(I2.2)') $
  +string(des,format=desformat)
for jj=0,1 do begin
    for kk=0,n_elements(adstr)-1 do begin
        tmp= adstr[kk]
        pos= strpos(tmp,' ')
        if (pos GE 0) then strput, tmp,'0',pos
        adstr[kk]= tmp
    endfor
endfor

if (prefix EQ '') then jstr= 'J' else jstr= ' J'
return, prefix+jstr+adstr
end
