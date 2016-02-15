;+
; NAME:
;   dec2hms
;
; PURPOSE:
;   convert decimal number to HH:MM:SS (base sixty)
;
; CALLING SEQUENCE:
;   result = dec2hms(angle_in, double=double)
;
; INPUTS:
;   angle_in  - angle [hours for RA, degrees for declination]
;
; KEYWORDS:
;   double    - double precision math
;
; OUTPUTS:
;   result    - string containing HH:MM:SS or DD:MM:SS
;
; EXAMPLES:
;   ra_string  = dec2hms(ra_degree/15)
;   dec_string = dec2hms(dec_degree)
;
; COMMENTS:
;   This function does not convert from hours to degrees!
;   Pass type double, or suffer 0.003 arcsec error
;   There is some tom-foolery to prevent roundoff problems
;     (like 1 -> '00 59 60.0') but if you aren't using /double
;     you have no right to expect high precision anyway. 
;
; REVISION HISTORY:
;   Written by Douglas Finkbeiner in ancient times
;   2000-Nov-29 - double keyword added            - DPF
;   2005-Sep-16 - call self recursively for array - DPF
;
;----------------------------------------------------------------------
function dec2hms,angle_in, double=double  

; -------- if angle_in is an array, call self recursively
  n_angle = n_elements(angle_in) 
  if n_angle GT 1 then begin 
     hmsarr = strarr(n_angle)
     for i=0L, n_angle-1 do begin 
        hmsarr[i] = dec2hms(angle_in[i], double=double)
     endfor 
     return, hmsarr
  endif

; -------- if there is only one angle, carry on with old routine
  eps = (machar(double=double)).eps ; machine precision

  angle = double(angle_in)
  neg   = angle LT 0.0d
  angle = abs(angle)
  h     = floor(angle+eps)
  angle = (angle-h)*60.0d
  m     = floor(angle+eps*60)
  s     = ((angle-m)*60.0d) > 0 ; can be slightly LT 0 due to roundoff
  
  if keyword_set(double) then begin 
     if ((s ge 59.999949d) and (s le 60.00d)) then s=59.999949d
     hms=string(h,m,s,format='(I3.2,I3.2,F8.4)') 
  endif else begin 
     if ((s ge 59.949d) and (s le 60.00d)) then s=59.949d
     hms=string(h,m,s,format='(I3.2,I3.2,F5.1)')
  endelse 

  if (strmid(hms,7,1) eq ' ') then strput,hms,'0',7
  if neg then strput,hms,'-',0

  return, hms
end
