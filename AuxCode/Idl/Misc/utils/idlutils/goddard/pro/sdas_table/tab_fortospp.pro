pro tab_fortospp,format,sppformat
;+
; NAME:
;	TAB_FORTOSPP
; PURPOSE:
;	Procedure to convert a FORTRAN format to an SPP format specfication.
;
; CALLING SEQUENCE:
;	sppformat, format, sppformat
;
; INPUTS:
;	format - fortran format specification
;
; OUTPUTS:
;	sppformat - sppformat specification
;
; HISTORY:
;	version 1  D. Lindler   Jan, 1989
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;------------------------------------------------------------------------
;
; decode fortran specification
;
	f=strtrim(strupcase(format),2)
	type=strmid(f,0,1)
	f=strmid(f,1,strlen(f)-1)
	w=gettok(f,'.')
	case type of
		'I' : sppformat = w+'d'
		'Z' : sppformat = w+'x'
		'O' : sppformat = w+'o'
		'F' : sppformat = w+'.'+f+'f'
		'E' : sppformat = w+'.'+f+'e'
		'D' : sppformat = w+'.'+f+'e'
		'G' : sppformat = w+'.'+f+'g'
		'A' : sppformat = w+'s'
		else: begin
			print,'TAB_FORTOSPP -- format '+format+' not supported'
			retall
			end
	endcase
	return
end
