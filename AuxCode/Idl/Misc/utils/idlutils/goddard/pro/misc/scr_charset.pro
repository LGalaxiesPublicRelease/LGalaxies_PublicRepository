PRO SCR_CHARSET, g, cset
;+
; NAME:
;	SCR_CHARSET
; PURPOSE:
;	To change the character sets.
; CALLING SEQUENCE:
;	scr_charset [, g, cset]
; INPUTS:
;	g     --  The terminal character set to change (either 0, for the
;	          G0 designator, or 1, for the G1 designator).  0 = default.
;	cset  --  The character set to use:
;	               0 : United Kingdom.
;	               1 : United States (USASCII)  --  default.
;	               2 : Special graphics characters and line drawing set.
;	               3 : Alternate character ROM.
;	               4 : Alternate character ROM special graphics chars.
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
IF n_params(0) LT 1 THEN g = 0
IF n_params(0) LT 2 THEN cset = 1
;
;			Set up the command string.
;
IF g EQ 1 THEN mid = ')' ELSE mid = '('
CASE cset OF
	   0 : scmd = strtrim(27B,2) + '[' + mid + 'A'	; Up
	   2 : scmd = strtrim(27B,2) + '[' + mid + '0' ; Left
	   3 : scmd = strtrim(27B,2) + '[' + mid + '1' ; Right
	   4 : scmd = strtrim(27B,2) + '[' + mid + '2' ; Right
	ELSE : scmd = strtrim(27B,2) + '[' + mid + 'B' ; Down
ENDCASE
;
;			Issue the command.
;
fmt = "(A" + strtrim(strlen(scmd),2) + ",$)"
print, format=fmt, scmd
;
RETURN
END
