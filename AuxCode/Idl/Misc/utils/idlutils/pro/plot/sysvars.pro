;+
; NAME:
;   sysvars
;
; PURPOSE:
;   store IDL system variables in a structure for later restoration 
;
; CALLING SEQUENCE:
;   state = sysvars(print=print, X=X)
;
; KEYWORDS:
;   print  - set default Postscript settings
;   X      - set default X settings
;
; OUTPUTS:
;   state  - structure containing previous values of !p,!d,!x,!y,!z
;
; EXAMPLES:
;   state = sysvars(/print)
;     <code to print stuff>
;   restore_sysvars, state
;
; COMMENTS:
;   Use restore_sysvars to restore !{p,x,y,z}.  !d is a readonly, but
;   is carried along for the ride. 
;
; REVISION HISTORY:
;   2001-Aug-06  Written by Douglas Finkbeiner, Princeton
;----------------------------------------------------------------------
function sysvars, print=print, X=X

  state = {p:!p, $
           d:!d, $
           x:!x, $
           y:!y, $
           z:!z}
  
  if keyword_set(print) then begin ; set up postscript defaults
     !p.font      = 0
     !p.thick     = 2
     !x.thick     = 2
     !y.thick     = 2
     !p.charsize  = 1.5
     !p.charthick = 2
  endif 

  if keyword_set(X) then begin ; set up X defaults
     !p.font      = -1
     !p.thick     = 1
     !x.thick     = 1
     !y.thick     = 1
     !p.charsize  = 1.5
     !p.charthick = 1
  endif 

  return, state
end

