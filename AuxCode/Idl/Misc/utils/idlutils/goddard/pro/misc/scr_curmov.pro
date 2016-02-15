PRO SCR_CURMOV, cmd, n
;+
; NAME:
;	SCR_CURMOV
; PURPOSE:
;	To mov the cursor around the screen relative to its original position.
; CALLING SEQUENCE:
;	scr_curmov [, cmd, n]
; INPUTS:
;	cmd  --  An integer indicating the direction in which to move the curs.
;	              0 : Up
;	              1 : Down  (Default)
;	              2 : Left
;	              3 : Right
;	n    --  The number of spaces to move the cursor.  If not specified
;	         (or if less than or equal to zero), this is set to one.
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
IF n_params(0) LT 1 THEN cmd = 1
IF n_params(0) LT 2 THEN n = 1
n = n > 1
;
;			Set up the command string.
;
CASE cmd OF
	   0 : scmd = strtrim(27B,2) + "[" + strtrim(n,2) + "A"	; Up
	   2 : scmd = strtrim(27B,2) + "[" + strtrim(n,2) + "D" ; Left
	   3 : scmd = strtrim(27B,2) + "[" + strtrim(n,2) + "C" ; Right
	ELSE : scmd = strtrim(27B,2) + "[" + strtrim(n,2) + "B" ; Down
ENDCASE
;
;			Issue the command.
;
fmt = "(A" + strtrim(strlen(scmd),2) + ",$)"
print, format=fmt, scmd
;
RETURN
END
