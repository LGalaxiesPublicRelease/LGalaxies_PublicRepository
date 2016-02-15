;+
; NAME:
;   hms2dec
;
; PURPOSE:
;   convert base-sixty string to decimal
;
; CALLING SEQUENCE:
;   decimal = hms2dec(hmsstring)
;
; INPUTS:
;   hmsstring - e.g. 'HH:MM:SS.S' or 'HHhMMmSS.S' or 'HH MM SS.S' etc. 
;
; OUTPUTS:
;   decimal  - double-precision decimal number
;
; RESTRICTIONS:
;   
; EXAMPLES:
;   ra = hms2dec(ra_string)*15
;
; COMMENTS:
;   The strpos() loop below looks awkward but is probably a lot 
;    faster than repstr()
;   I left most of the 20th century code intact, so I wouldn't
;    introduce any change in behavior. 
;
; REVISION HISTORY:
;   Written by Douglas Finkbeiner in ancient times
;   2000-Nov-29 fixed to handle "d" properly! - DPF
;   2005-Sep-16 call self recursively for array - DPF
;
;----------------------------------------------------------------------
function hms2dec,hmsstring

; -------- if hmsstring is an array, call self recursively
  n_angle = n_elements(hmsstring)
  if n_angle GT 1 then begin
     angle = dblarr(n_angle)
     for i=0L, n_angle-1 do begin
        angle[i] = hms2dec(hmsstring[i])
     endfor
     return, angle
  endif

; -------- Replace any h, m, or s in the string with a space
  hms=strlowcase(hmsstring+' ')
  delim = ['h', 'd', 'm', 's', ':', ':']
  for i=0L, n_elements(delim)-1 do begin 
     hpos = strpos(hms, delim[i])
     if (hpos ge 0) then strput, hms, ' ', hpos
  endfor

; -------- Now continue with the ancient code (it works...)
;              (stregex() would be nice here...)

; strip leading spaces
while (strpos(hms,' ') eq 0) do hms=strmid(hms,1,strlen(hms)-1)
neg = (strmid(hms,0,1) eq '-') 
if neg then strput,hms,' ',0
while (strpos(hms,' ') eq 0) do hms=strmid(hms,1,strlen(hms)-1)

cut=strpos(hms,' ')
if (cut gt 0) then begin
  hstr=strmid(hms,0,cut)
  hms=strmid(hms,cut+1,strlen(hms))
  reads,hstr,h
endif else h=0
while (strpos(hms,' ') eq 0) do hms=strmid(hms,1,strlen(hms)-1)
cut=strpos(hms,' ')
if (cut gt 0) then begin
  reads,strmid(hms,0,cut),m
  hms=strmid(hms,cut+1,strlen(hms))
endif else m=0

while (strpos(hms,' ') eq 0) do hms=strmid(hms,1,strlen(hms)-1)
cut=strpos(hms,' ')
if (cut gt 0) then begin
  reads,strmid(hms,0,cut),s
  hms=strmid(hms,cut+1,strlen(hms))
endif else s=0
dec=h+m/60.0d + s/3600.0d
if neg then dec=-dec
return,dec
end
