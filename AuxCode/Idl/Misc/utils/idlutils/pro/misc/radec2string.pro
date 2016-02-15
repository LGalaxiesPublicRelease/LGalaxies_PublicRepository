;+
; NAME:
;   radec2string
; PURPOSE:
;   converts RA and DEC to string versions
; CALLING SEQUENCE:
;   radec2string, ra, dec, outstring [, rahr=, ramin=, rasec=, $
;     decdeg=, decmin=, decsec=, rastr=, decstr=, precision= ]
; INPUTS:
;   ra - deg
;   dec - deg
; OPTIONAL INPUTS:
;   precision  - how many decimal places to put on dec; ra gets one more
; OUTPUTS:
;   outstring - output string of form xx:xx:xx.xxx[+|-]xx:xx:xx.xx
; OPTIONAL OUTPUTS:
;   rahr, ramin, rasec - components of first half of string
;   decdeg, decmin, decsec - components of second half of string
;   rastr, decstr - first and second halves separately
; COMMENTS:
;   Defaults to SDSS name conventions.
;   Don't EVER use goddard/pro/astro/adstring.pro.
;   Hacked from hogg_iau_name
; BUGS:
;   Doesn't deal properly with precision<0 
; REVISION HISTORY:
;   2005-10-24  Blanton (NYU)
;-
pro radec2string, ra,dec, outstring, precision=precision, $
                  rahr=rahr, ramin=ramin, rasec=rasec, $
                  decdeg=decdeg, decmin=decmin, decsec=decsec, $
                  rastr=rastr, decstr=decstr 

if (n_elements(prefix) eq 0) then prefix= 'SDSS'
if (n_elements(precision) EQ 0) then precision= 1

rah= floor(ra/15.0)
ram= floor(((ra/15.0)-double(rah))*60.0)
ras= (((ra/15.0)-double(rah))*60.0-double(ram))*60.0
ras= round(ras*10.0^(precision+1))/10.0^(precision+1)
rasformat= '(F0'+strtrim(string(precision+4),2)+ $
  '.'+strtrim(string(precision+1),2)+')'

desgn= replicate('+',n_elements(dec))
negative= where(dec LT 0,nneg)
if (nneg GT 0) then desgn[negative]= '-'
adec= abs(dec)
ded= floor(adec)
dem= floor((adec-double(ded))*60.0)
des= ((adec-double(ded))*60.0-double(dem))*60.0
des= round(des*10.0^precision)/10.0^precision
desformat= '(F0'+strtrim(string(precision+3),2)+ $
  '.'+strtrim(string(precision),2)+')'
if (precision EQ 0) then desformat= '(I2.2)'

rahr= string(rah,format='(I2.2)') 
ramin= string(ram,format='(I2.2)') 
rasec= string(ras,format=rasformat) 
decdeg= desgn+string(ded,format='(I2.2)') 
decmin=string(dem,format='(I2.2)') 
decsec= string(des,format=desformat)

rastr=rahr+':'+ramin+':'+rasec
decstr=decdeg+':'+decmin+':'+decsec

outstring=rastr+decstr

end
