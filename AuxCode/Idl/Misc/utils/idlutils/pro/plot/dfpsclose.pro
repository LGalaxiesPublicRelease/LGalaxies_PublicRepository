;+
; NAME:
;   dfpsclose
;
; PURPOSE:
;   Finkbeiner's routine to close a PostScript file previously opened
;   with DFPSOPEN, and revert to sending plots to the X-display.
;
; CALLING SEQUENCE:
;   dfpsclose
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   The-Beginning-of-Time  Written by Doug Finkbeiner, Berkeley.
;   05-Sep-1999  Modified and commented by David Schlegel, Princeton.
;   06-Aug-2001  Check if X device is set already - DPF, Princeton
;-
;------------------------------------------------------------------------------
pro dfpsclose
  
  if !d.name EQ 'X' then begin 
     print, 'Postscript device not open!'
  endif else begin 
     device, /close
     set_plot, 'X'
  endelse 

  return
end
;------------------------------------------------------------------------------
