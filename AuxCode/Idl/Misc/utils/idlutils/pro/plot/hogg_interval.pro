;+
; NAME:
;   hogg_interval
; PURPOSE:
;   find best interval for tickmarks on a plot
; INPUTS:
;   range     - the range for the axis
; OPTIONAL INPUTS:
;   nticks    - the approximate, desired number of ticks; default 5
; BUGS:
;   Goes insane if range contains NaN.
; REVISION HISTORY:
;   2002-03-25  written - Hogg
;-
function hogg_interval, range,nticks=nticks
if abs(range[1]-range[0]) EQ 0.0 then begin
    return, 0.0
endif else begin
    if NOT keyword_set(nticks) then nticks= 5.0
    loginterval= alog10(abs(range[1]-range[0])/double(nticks))
    interval= 1d1^round(loginterval)
    dloginterval= loginterval-alog10(interval)
    if dloginterval LT (-alog10(2d0)/2d0) then interval= interval/2d0
    if dloginterval GT (alog10(2d0)/2d0) then interval= interval*2d0
    return, interval
endelse
end
