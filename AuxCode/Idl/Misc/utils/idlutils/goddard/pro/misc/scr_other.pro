PRO SCR_OTHER, str
;+
; NAME:
;	SCR_OTHER
; PURPOSE:
;	To allow the user to issue any ESCAPE sequence.
; CALLING SEQUENCE:
;	scr_other, str
; INPUTS:
;	str  --  A string containing the escape sequence.  The initial ESCAPE
;	         should not be included; this will be added by this procedure.
;	         This parameter is NOT optional; if not available, the 
;	         procedure will return without doing anything.
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
;			Check argument.
;
IF n_params(0) GE 1 THEN BEGIN
;
;			Set up the command string.
;
	scmd = strtrim(27B,2) + str
;
;			Issue the command.
;
	fmt = "(A" + strtrim(strlen(scmd),2) + ",$)"
	print, format=fmt, scmd
;
ENDIF
;
RETURN
END
