;+
; NAME:
;   local_extragalactic
; PURPOSE:
;   returns list of local extragalactic locations and names to plot
; CALLING SEQUENCE:
;   local_extragalactic, ra, dec, cz_lg, names, list=list
; OPTIONAL INPUTS:
;   list - list of names to match to (default: 
;       ['VIRGO_CLUSTER', 'URSA_MAJOR_CLUSTER', $
;        'PISCES_CLUSTER', 'COMA_CLUSTER' ] )
; OUTPUTS:
;   ra - [N] ra (J2000 deg)
;   dec - [N] dec (J2000 deg)
;   cz_lg - [N] local group frame cz
;   names - [N] object name
;-
pro local_extragalactic, out_ra, out_dec, out_cz_lg, out_names, list=list

if(NOT keyword_set(list)) then $
  list= ['VIRGO_CLUSTER', 'URSA_MAJOR_CLUSTER', 'PISCES_CLUSTER', $
         'COMA_CLUSTER']

infile=getenv('IDLUTILS_DIR')+'/data/objects/ned-query-clusters-lt0.03.txt'
nclusters=numlines(infile)
ra=dblarr(nclusters)
dec=dblarr(nclusters)
cz_lg=fltarr(nclusters)
names=strarr(nclusters)
openr,unit,/get_lun,infile
line=' ' 
for i=0L, nclusters-1L do begin
    readf, unit, line
    words=strsplit(line, ' ', /extr)
    names[i]=words[0]
    string2radec, words[1], words[2], words[3], words[4], words[5], words[6], $
      tmp_ra, tmp_dec
    ra[i]=tmp_ra
    dec[i]=tmp_dec
    cz_lg[i]=helio_to_lg(float(words[7])/2.99792e+5, ra[i], dec[i])*2.99792e+5
endfor
free_lun, unit

use=bytarr(n_elements(names))
for i=0L, n_elements(list)-1L do begin
    iuse=where(strmatch(names, list[i]), nuse)
    if(nuse gt 0) then use[iuse]=1
endfor

iuse=where(use, nuse)
if(nuse gt 0) then begin
    out_ra=ra[iuse]
    out_dec=dec[iuse]
    out_cz_lg=cz_lg[iuse]
    out_names=names[iuse]
endif


end
