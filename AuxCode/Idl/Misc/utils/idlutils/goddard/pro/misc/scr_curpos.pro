PRO SCR_CURPOS, lin, col
;+
; NAME:
;	SCR_CURPOS
; PURPOSE:
;	To position the cursor at the specified screen location.  
; EXPLANATION:
;	Unspecified coordinates are set to one.  Please note that the ESCAPE 
;	sequence expects the coordinates to be counted from (1,1).
; CALLING SEQUENCE:
;	scr_curpos [, lin, col]
; INPUTS:
;	lin  --  The screen line coordinate.
;	col  --  The screen column coordinate.
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
;			Check arguments.
;
IF n_params(0) LT 1 THEN lin = 1
IF n_params(0) LT 2 THEN col = 1
;
;			Set up the command string.
;
scmd = strtrim(27B,2) + "[" + strtrim(lin,2) + ";" + strtrim(col,2) + "H"
;
;			Issue the command.
;
fmt = "(A" + strtrim(strlen(scmd),2) + ",$)"
print, format=fmt, scmd
;
RETURN
END
