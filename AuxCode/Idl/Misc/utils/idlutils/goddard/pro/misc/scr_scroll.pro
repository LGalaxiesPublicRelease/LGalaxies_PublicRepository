PRO SCR_SCROLL, top, bot
;+
; NAME:
;	SCR_SCROLL
; PURPOSE:
;	Define the scrolling area on the screen.
; EXPLANATION:  
;	Please note that the line coordinates should be counted from 1.
; CALLING SEQUENCE:
;	scr_scroll [, top, bot]
; INPUTS:
;	top  --  The line to be the top of the scrolling area.
;	         The default value is 1 and the maximum value is 23.
;	bot  --  The line to be the bottom of the scrolling area.
;	         The default value is 24 and the minimum value is 2.
; OUTPUTS:
;	None.
; SIDE EFFECTS:
;	NOTE:  The screen coordinate system is NOT effected.  (1,1) is not
;	       the top of the scrolling area but the top of the screen.
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
IF n_params(0) LT 1 THEN top = 1
top = ((top > 1) < 23)
IF n_params(0) LT 2 THEN bot = 24
bot = ((bot < 24) > 2)
;
;			Set up the command string.
;
scmd = strtrim(27B,2) + "[" + strtrim(top,2) + ";" + strtrim(bot,2) + "r"
;
;			Issue the command.
;
fmt = "(A" + strtrim(strlen(scmd),2) + ",$)"
print, format=fmt, scmd
;
RETURN
END
