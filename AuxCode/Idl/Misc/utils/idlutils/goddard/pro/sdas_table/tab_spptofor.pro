pro tab_spptofor,sppformat,format,width
;+
; NAME:
;	TAB_SPPTOFOR  
; PURPOSE:
;	This procedure converts an spp format specification to a normal
;	Fortran format specification.
;
; CALLING SEQUENCE:
;	tab_spptofor, sppformat, format, width
;
; INPUTS:
;	sppformat - spp format specification (without preceeding %)
;
; OUTPUTS:
;	forformat - fortran format specification (string)
;	width - field width (integer)
;
; HISTORY:
;	version 1  D. Lindler  Jan 1989
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;--------------------------------------------------------------------------
;
; initialization
;
	zero=byte('0') & zero=zero[0]
	nine=byte('9') & nine=nine[0]
	period=byte('.') & period=period[0]
	minus=byte('-') & minus=minus[0]
	O=byte('O') & O=O[0]
;
; determine if leading character is - or O. if so ignore it
;
	b=byte(sppformat) & nb=n_elements(b)
	if (b[0] eq O) or (b[0] eq minus) then begin
		nb=nb-1
		b=b[1:nb]
	endif
;
; get field width by searching for consecutive digits
;
	w = 0
	while (b[0] ge zero) and (b[0] le nine) do begin
		w = w*10 + b[0]-zero
		nb=nb-1
		b=b[1:nb]
	end
;
; if next character is a period search for number of digits of
; precision
;
	d = 0
	if b[0] eq period then begin
		nb=nb-1
		b=b[1:nb]
		while (b[0] ge zero) and (b[0] le nine) do begin
			d = d*10 + b[0]-zero
			nb=nb-1
			b=b[1:nb]
		end
	end
;
; convert to strings for use in fortran format specification
;
	w=w>d
	w=w>1
	width=w
	sw=strtrim(w,2)
	sd=strtrim(d,2)
;
; next character should be the format character
;
	case string(b[0]) of
	   'b' : format='I'+sw
	   'c' : format='A'+sw
	   'd' : format='I'+sw
	   'e' : format='E'+sw+'.'+sd
	   'f' : format='F'+sw+'.'+sd
	   'g' : format='G'+sw+'.'+sd
	   'h' : format='F'+sw+'.'+sd
	   'm' : format='F'+sw+'.'+sd
	   'o' : format='O'+sw
	   'r' : format='I'+sw
	   's' : format='A'+sw
	   'u' : format='I'+sw
	   'x' : format='Z'+sw
	   else: begin
		print,'TAB_SPPTOFOR-- '+sppformat+' spp format not supported'
		retall
		end
	endcase
	return
end
