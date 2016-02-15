;+
; NAME:
;   hogg_name2radec
; PURPOSE:
;   Takes properly-formatted iau name in "J" format and converts it into ra,dec
;   USE AT YOUR OWN RISK; NO IMPLIED WARRANTY
; EXAMPLES:
;   IDL> hogg_name2radec, 'SDSS_J013713.84-092820.4',ra,dec
;   IDL> print, ra,dec
;        24.307671      -9.4723336
;   IDL> print, hogg_iau_name(ra,dec)
;   SDSS J013713.84-092820.4
; INPUTS:
;   name       - name
; OUTPUTS:
;   ra         - J2000 deg
;   dec        - J2000 deg
; OPTIONAL INPUTS:
; COMMENTS:
;   Doesn't require anything particular *before* the J, but it assumes
;     (dumbly) that the RA begins at the last J.
; BUGS:
;   USE AT YOUR OWN RISK; NO IMPLIED WARRANTY
;   Requires J2000 "J".
;   Ultra-brittle.
; REVISION HISTORY:
;   2004-10-22  started by you know who
;-
pro hogg_name2radec, name,ra,dec
ra= -1
dec= -1
splog, name
rastr= strcompress(name,/remove_all)
jindx= strpos(rastr,'J',/reverse_search)
if (jindx EQ -1) then begin
    splog, 'ERROR: no "J" found'
    return
endif
rastr= strmid(rastr,jindx+1)
splog, rastr
signindx= stregex(rastr,'[+-]')
if (signindx EQ -1) then begin
    splog, 'ERROR: no sign found'
    return
endif
decstr= (stregex(strmid(rastr,signindx),'([+-0123456789\.]+)',/subexpr,/extract))[0]
rastr= strmid(rastr,0,signindx)
splog, rastr
splog, decstr
rahstr= strmid(rastr,0,2)
ramstr= strmid(rastr,2,2)
rasstr= strmid(rastr,4)
desgnstr= strmid(decstr,0,1)+'1'
dedstr= strmid(decstr,1,2)
demstr= strmid(decstr,3,2)
desstr= strmid(decstr,5)
help, rahstr,ramstr,rasstr,dedstr,demstr,desstr
tiny= 1D-3 ; added so the output truncates back to the same name!
ra= (double(rahstr)+double(ramstr)/6D1 $
     +(double(rasstr)+tiny)/3.6D3)*1.5D1
dec= (double(dedstr)+double(demstr)/6D1 $
      +(double(desstr)+tiny)/3.6D3)*double(desgnstr)
return
end
