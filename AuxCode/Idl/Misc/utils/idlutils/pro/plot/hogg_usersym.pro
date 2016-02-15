;+
; NAME:
;   hogg_usersym
; PURPOSE:
;   make an n-sided plotting point
; USAGE:
;   hogg_usersym, N [,scale=scale,/diagonal]
;   plot, x,y,psym=8
; INPUTS:
;   N           - number of sides on the polygon
; OPTIONAL INPUTS:
;   scale       - linear size
;   _extra      - keywords for usersym (see usersym help page)
;                 eg, /fill or thick=thick
; KEYWORDS
;   diagonal    - rotate symbol through 1/2 of 1/N turns
;   stellar     - make a stellar symbol rather than convex polygon
;                 (only works for odd values of N!)
; REVISION HISTORY:
;   2002-04-09  written - Hogg
;-
pro hogg_usersym, N,diagonal=diagonal,scale=scale,stellar=stellar, $
                  _EXTRA=KeywordsForUsersym
  if keyword_set(diagonal) then delta= 0D0 else delta= 5D-1
  if NOT keyword_set(scale) then scale= 1D0
  if keyword_set(stellar) then tscale= 1.75*scale else $
    tscale= scale/(1.0-1.0/double(N))
  if keyword_set(stellar) then $
    theta= 2D0*!PI*(dindgen(N+1)+delta)*(double(N-1)/(2*double(N))) $
  else theta= 2D0*!PI*(dindgen(N+1)+delta)/double(N)
  usersym, tscale*sin(theta),-tscale*cos(theta),_EXTRA=KeywordsForUsersym
end
