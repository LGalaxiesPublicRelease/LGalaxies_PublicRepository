;+
; NAME:
;   hogg_ned_name2radec
; PURPOSE:
;   Takes a name, throws it at NED, retreives RA, Dec (J2000)
;   USE AT YOUR OWN RISK; NO IMPLIED WARRANTY
; INPUTS:
;   name       - name
; OUTPUTS:
;   ra         - J2000 deg
;   dec        - J2000 deg
; OPTIONAL INPUTS:
; COMMENTS:
; BUGS:
;   Does nasty things to the name to format it for URL-ifying.
;   USE AT YOUR OWN RISK; NO IMPLIED WARRANTY.
;   Relies on NED not changing in any substantial way!
; DEPENDENCIES:
;   wget
; REVISION HISTORY:
;   2005-02-01  started by you know who
;-
pro hogg_ned_name2radec, name,ra,dec
ra= -1
dec= -1
splog, name
urlname= name
urlname= strjoin(strsplit(urlname,'+',/extract),'%2B')
urlname= strjoin(strsplit(urlname,' ',/extract),'+')
url= 'http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?objname=' $
  +urlname+'&extend=no'
;;cmd= 'wget -O - "'+url+'" | grep "Equatorial" | grep "(J2000.0)"'
cmd= 'wget -O - "'+url+'" | grep "Equatorial (J2000.0)"'
splog, cmd
spawn, cmd,output,stderr
splog, output
output= strsplit(output[0],/extract)
if (n_elements(output) GT 3) then begin
    ra= double(output[2])
    dec= double(output[3])
endif
splog, ra,dec
return
end
