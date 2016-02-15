;+
; NAME:
;   string2radec
;
; PURPOSE:
;   Convert hours, min, sec + deg, min, sec to ra dec in floating degrees
;
; CALLING SEQUENCE:
;   string2radec, rahour, ramin, rasec, decdeg, decmin, decsec, ra, dec, $
;    [ rastr=, decstr= ]
; INPUTS:
;   rahour - hours in ra (string)
;   ramin - minutes in ra (string)
;   rasec - seconds in ra (string)
;   decdeg - degrees in dec (string) 
;   decmin - arcminutes in dec (string)
;   decsec - arcseconds in dec (string)
; OPTIONAL KEYWORDS:
;   rastr - If set, then override RAHOUR,RAMIN,RASEC with a single string
;           with those values separated by colons
;   decstr- If set, then override DECDEG,DECMIN,DECSEC with a single string
;           with those values separated by colons
; OUTPUTS:
;   ra - ra in degrees
;   dec - dec in degrees
; COMMENTS:
;   The sign in dec is set by the sign in "decdeg", so that you need
;   to propagate the negative in "-00" to the whole expression
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   18-Nov-2003  Michael Blanton, NYU
;-
;------------------------------------------------------------------------------
pro string2radec, rahour, ramin, rasec, decdeg, decmin, decsec, ra, dec, $
 rastr=rastr, decstr=decstr

if(n_params() lt 6) then begin
    print, 'Syntax - string2radec, rahour, ramin, rasec, decdeg, decmin, decsec, ra, dec'
    return
endif

if(keyword_set(rastr)) then begin
    rahour=strarr(n_elements(rastr))
    ramin=strarr(n_elements(rastr))
    rasec=strarr(n_elements(rastr))
    decdeg=strarr(n_elements(rastr))
    decmin=strarr(n_elements(rastr))
    decsec=strarr(n_elements(rastr))
    for i=0L, n_elements(rastr)-1L do begin
        words=strsplit(rastr[i], ':', /extr)
        rahour[i]=words[0]
        ramin[i]=words[1]
        rasec[i]=words[2]
        words=strsplit(decstr[i], ':', /extr)
        decdeg[i]=words[0]
        decmin[i]=words[1]
        decsec[i]=words[2]
    endfor
endif

ra=(double(rahour)+double(ramin)/60.+double(rasec)/3600.)*360./24
dec=(abs(double(decdeg))+double(decmin)/60.+double(decsec)/3600.)
dec=(1.-2.*(stregex(decdeg,'-') ge 0))*dec

end
;------------------------------------------------------------------------------
