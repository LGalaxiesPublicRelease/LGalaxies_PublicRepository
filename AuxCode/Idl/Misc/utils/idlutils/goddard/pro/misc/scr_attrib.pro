PRO SCR_ATTRIB, a1, a2, a3, a4, a5
;+
; NAME:
;	SCR_ATTRIB
; PURPOSE:
;	To set the screen attribute to those given, in the given order.
; CALLING SEQUENCE:
;	scr_attrib [, a1, a2, a3, a4, a5]
; INPUTS:
;	a1 - a5  --  The attribute codes.  The attributes are set in the
;	             command string in the given order.  Thus, if a1 turns
;	             the attributes off and a2 sets reverse video, the final
;	             attribute will reset and then set to reverse video.  If
;	             the order were reversed, then the current attribute 
;	             would have reverse video added to it, and then would be
;	             reset, leaving the terminal with all attributes off.  Up
;	             to five attribute codes may be specified.  The codes are:
;	                  0 : all attributes off  (default)
;	                  1 : bold on
;	                  2 : underscore on
;	                  3 : blink on
;	                  4 : reverse video on
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
;			Check arguments.  Put the attributes into an array.
;
n = n_params(0) < 5
attcod = ['0;', '1;', '4;', '5;', '7;']
attrib = replicate(attcod[0], (n > 1))
IF n GE 1 THEN BEGIN
	IF (a1 LT 0) OR (a1 GT 4) THEN a1 = 0
	attrib[0] = attcod[a1]
ENDIF
IF n GE 2 THEN BEGIN
	IF (a2 LT 0) OR (a2 GT 4) THEN a2 = 0
	attrib[1] = attcod[a2]
ENDIF
IF n GE 3 THEN BEGIN
	IF (a3 LT 0) OR (a3 GT 4) THEN a3 = 0
	attrib[2] = attcod[a3]
ENDIF
IF n GE 4 THEN BEGIN
	IF (a4 LT 0) OR (a4 GT 4) THEN a4 = 0
	attrib[3] = attcod[a4]
ENDIF
IF n GE 5 THEN BEGIN
	IF (a5 LT 0) OR (a5 GT 4) THEN a5 = 0
	attrib[4] = attcod[a5]
ENDIF
;
;			Set up the command string.
;
scmd = strtrim(27B,2) + '[' + attrib[0]
IF (n GT 1) THEN BEGIN
	FOR i = 1, (n-1) DO scmd = scmd + attrib[i]
ENDIF
n = strlen(scmd)
strput, scmd, 'm', (n - 1)
;
;			Issue the command.
;
fmt = "(A" + strtrim(strlen(scmd),2) + ",$)"
print, format=fmt, scmd
;
RETURN
END
