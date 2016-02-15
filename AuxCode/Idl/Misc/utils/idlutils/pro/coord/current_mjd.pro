;+
; NAME:
;   current_mjd
; PURPOSE:
;   Transform systime() to current MJD.
; CALLING SEQUENCE:
;   mjd= current_mjd()
; INPUTS:
; OPTIONAL KEYWORDS:
; OUTPUTS:
; OPTIONAL OUTPUTS:
; COMMENTS:
;   Your clock better be right when you ask for UT!
; EXAMPLES:
; BUGS:
;   Relies on IDL 5.4 feature /utc for systime().
; PROCEDURES CALLED:
; REVISION HISTORY:
;   2001-Feb-07  written by Hogg, NYU
;-
;------------------------------------------------------------------------------
function current_mjd
   return, (systime(/julian,/utc)-2400000.5D)
end
