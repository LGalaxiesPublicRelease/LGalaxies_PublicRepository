;+
; NAME:
;   restore_sysvars
;
; PURPOSE:
;   restore plot system variabls from structure created by sysvars
;
; CALLING SEQUENCE:
;   restore_sysvars, state
;
; INPUTS:
;   state  - structure containing previous values of !p,!d,!x,!y,!z
;   
; EXAMPLES:
;   state = sysvars(/print)
;     <code to print stuff>
;   restore_sysvars, state
;   
; REVISION HISTORY:
;   2001-Aug-06  Written by Douglas Finkbeiner, Princeton
;----------------------------------------------------------------------
pro restore_sysvars, state

  !p = state.p
  !x = state.x
  !y = state.y
  !z = state.z

  return
end
