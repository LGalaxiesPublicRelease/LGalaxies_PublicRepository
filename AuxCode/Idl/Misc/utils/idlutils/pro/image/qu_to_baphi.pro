;+
; NAME:
;   qu_to_baphi
; PURPOSE:
;   transform q and u values to b/a and pih
; CALLING SEQUENCE:
;   qu_to_baphi, q, u, ba, phi
; INPUTS:
;   q - Stokes q parameter
;   u - Stokes u parameter
; OPTIONAL INPUTS:
; OUTPUTS:
;   ba - b/a
;   phi - angle phi
; OPTIONAL INPUT/OUTPUTS:
; DATA DEPENDENCIES:
; COMMENTS:
; EXAMPLES:v
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   10-Aug-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro qu_to_baphi, q, u, ba, phi

tanphi=u/q
phi=atan(tanphi)*180./!DPI
yy=q^2+u^2
ba=(1.+yy-2.*sqrt(yy))/(1-yy)

end

