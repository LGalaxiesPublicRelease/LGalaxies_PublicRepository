;+
; NAME:  
;       fillspan
;
; PURPOSE: 
;       return an array evenly sampling a given range. The endpoints
;       are always included.
;
; CALLING SEQUENCE:
;       span = fillspan(lo, hi, [spacing=spacing], [cnt=cnt])
;
; INPUTS:
;       lo, hi: the endpoints of the span. 
;
;   one of the following must be specified:
;
;       spacing=spacing: the desired spacing between points. Will be
;       rounded to allow the closest even spacing, but must be <=
;       hi-lo
;
;       cnt=cnt: how many elements to get. Must be >= 2, to allow the
;       endpoints to fit.
;-
function fillspan,lo,hi,spacing=spacing,cnt=cnt
    if (keyword_set(spacing) + keyword_set(cnt)) NE 1 then $
       message, 'either one of cnt or spacing must be specified'
    if keyword_set(spacing) then $
       cnt = 1 + fix((hi-lo)/spacing)
    if (cnt LT 2) then $
       message, 'the requested spacing is wider than the bounds.' 

    return,lo + (findgen(cnt)/(cnt-1)*(hi-lo))
 end
