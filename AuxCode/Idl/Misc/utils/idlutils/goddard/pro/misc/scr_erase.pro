PRO SCR_ERASE, cmd
;+
; NAME:
;	SCR_ERASE
; PURPOSE:
;	To erase portions of the terminal screen.
; CALLING SEQUENCE:
;	scr_erase [, cmd]
; INPUTS:
;	cmd  --  An integer telling the procedure what part of the screen to
;	         erase.  If not specified, it is set to 5.  Key:
;	                 0 : From cursor to end-of-line.
;	                 1 : From beginning-of-line to cursor.
;	                 2 : Entire line containing cursor.
;	                 3 : From cursor to end-of-screen.
;	                 4 : from beginning-of-screen to cursor.
;	              ELSE : Entire screen.
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
IF n_params(0) LT 1 THEN cmd = 5
;
;			Set up the command string.
;
CASE cmd OF
	   0 : scmd = strtrim(27B,2) + "[0K"	; From cursor to end-of-line.
	   1 : scmd = strtrim(27B,2) + "[1K"	; From beg.-of-line to cursor.
	   2 : scmd = strtrim(27B,2) + "[2K"	; Entire line containing cursor.
	   3 : scmd = strtrim(27B,2) + "[0J"	; From cursor to end-of-screen.
	   4 : scmd = strtrim(27B,2) + "[1J"	; from beg.-of-screen to cursor.
	ELSE : scmd = strtrim(27B,2) + "[2J"	; Entire screen.
ENDCASE
;
;			Issue the command.
;
fmt = "(A" + strtrim(strlen(scmd),2) + ",$)"
print, format=fmt, scmd
;
RETURN
END
