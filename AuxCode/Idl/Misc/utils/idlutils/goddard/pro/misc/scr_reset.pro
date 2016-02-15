PRO SCR_RESET, dum
;+
; NAME:
;	SCR_RESET
; PURPOSE:
;	To reset the terminal.
; CALLING SEQUENCE:
;	scr_reset
; INPUTS:
;	None.
; OUTPUTS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	This procedure will only work with DEC compatible equipment (or
;	terminal emulators).
; PROCEDURE:
;	A string containing the appropriate DEC terminal command is put 
;	together and printed.  NOTE:  In general, the DEC commands correspond
;	to the ANSI escape sequences.
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, STX, May 1990.
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;			Set up the command string.
;
scmd = strtrim(27B,2) + 'c'
;
;			Issue the command.
;
fmt = "(A" + strtrim(strlen(scmd),2) + ",$)"
print, format=fmt, scmd
;
RETURN
END
