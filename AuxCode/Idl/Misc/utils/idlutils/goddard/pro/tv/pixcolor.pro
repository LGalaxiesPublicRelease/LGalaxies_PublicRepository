pro pixcolor, pix_value, color
;+
; NAME:
;	PIXCOLOR
; PURPOSE:
;	Assign colors to specified pixel values in a color lookup table
;
; CALLING SEQUENCE:
;      	PIXCOLOR, pixvalue, color         ;Set color at specified pixel values
;
; OPTIONAL INPUT PARMETERS:
;	pixvalue - value or range of pixel value whose color will be modified.
;		A single pixel value may be specified by an integer
;		If a range of values is specified, then it must be written
;		as a string, with a colon denoting the range (e.g.'102:123')
;		If omitted, program will prompt for this parameter.
;
;	color -    single character string giving specified color values.
;		Available options are 'R' (red), 'B' (blue), 'G' (green)
;		'Y' (yellow), 'T' (turquoise), 'V' (violet), 'W' (white)
;		or 'D' (dark).  If omitted, program will prompt for this 
;		parameter.
;
; OUTPUTS:
;	None
; PROCEDURE:
;	TVLCT is used in RGB mode to load the specified pixel values.
;
; EXAMPLE:
;	Set pixel values of 245 to a color of red
;
;	IDL> pixcolor,245,'R'
;
; REVISION HISTORY:
;	Written, W. Landsman ST Systems Corp.		February, 1987
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
 On_error,2

 npar = N_params()

 if ( N_elements(pix_value) EQ 0) then begin
	pix_value = ''
	print,'Enter pixel value(s) to be assigned a color value'
	print,'Value may be either number or a range (e.g. 102:123)'
	read,'Pixel Value(s): ',pix_value
 endif

 type = size(pix_value)
 if ( type[1] EQ 7 ) then begin
	pixmin = fix(gettok(pix_value,':')) >0
	if strlen(pix_value) eq 0 then pixmax = fix(pixmin)  $
		else pixmax = fix(pix_value) > pixmin < 255
 endif else begin                                               
	pixmin = fix(pix_value)>0<255
	pixmax = pixmin
 endelse 
 npts = pixmax - pixmin + 1

GETCOL: if ( npar LT 2 ) then begin
	color = ''
	print,'Enter first letter of color which pixel(s) will be asssigned'
	print,'Available options are '
	print,'Red (R), Blue (B), Green (G), Yellow (Y), Turquoise (T),
	print,'Violet (V), White (W), or Dark (D)
        read,color
 endif

 red = 0B & green = 0B & blue = 0B

 case strupcase(color) of
	'R': red = 1
	'G': green = 1
	'B': blue = 1
	'Y': begin & red=1 & green =1 & end
	'T': begin & blue=1 & green=1 & end
	'V': begin & red =1 & blue =1 & end
	'W': begin & red=1 & blue=1 & green=1 & end
	'D':
	else: begin
	message, 'Color of '+ $
           string(color)+ ' is not an available options',/continue
	npar =1
	goto,GETCOL 
	end
 endcase
 r = replicate(255*red,npts)      ;Multiply three color vectors by 
 g = replicate(255*green,npts)    ;color gun masks
 b = replicate(255*blue,npts)

 tvlct,r,g,b,pixmin

 return
 end
