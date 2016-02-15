function sunsymbol
;+
; NAME:
;	SUNSYMBOL
; PURPOSE:
;	Return the Sun symbol as a subscripted postscript character string
; EXPLANATION:
;	Returns the Sun symbol (circle with a dot in the middle) as a 
;	(subscripted) postscript character string.    Needed because although 
;	the Sun symbol	is available using the vector fonts as the string 
;	'!9n', it is not in the standard postscript set.   
;
; CALLING SEQUENCE:
;	result = SUNSYMBOL()
;
; INPUTS:
;	None
;
; OUTPUTS:
;	result - a scalar string representing the Sun symbol.   A different
;		string is output depending (1) the device is postscript and
;		hardware fonts are used (!P.FONT=0), (2) vector fonts are used,
;		or (3) hardware fonts are used on a non-postscript device.
;		For case (3), SUNSYMBOL simply outputs the 3 character string
;		'Sun'
;
; EXAMPLE:
;	To make the X-axis of a plot read  M/M_Sun
;	IDL>  plot,indgen(10),xtit = 'M / M' + sunsymbol()
;
; RESTRICTIONS:
;	(1) The postscript output does not have the dot perfectly centered in 
;		the circle
;	(2) SUNSYMBOL() includes subscript output positioning commands in the 
;		output string.
;       (3) True type fonts (!p.font = 1) are not supported
; REVISION HISTORY:
;	Written,  W. Landsman,    HSTX          April, 1997
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
 On_error,2

 if (!P.FONT EQ -1) then return,'!D!9n!N!X' else $
 if (!D.NAME NE 'PS') then return,'!DSun!N' else begin

;Want to use /AVANTGARDE,/BOOK which is the default font 17, but to make sure
;that ISOLATIN encoding is turned off, we'll define our own font.

   device,/AVANTGARDE,/BOOK,ISOLATIN=0,FONT_INDEX = 20

   return, '!20!S!DO!R!I ' + string(183b) + '!X!N'
 endelse
 end
