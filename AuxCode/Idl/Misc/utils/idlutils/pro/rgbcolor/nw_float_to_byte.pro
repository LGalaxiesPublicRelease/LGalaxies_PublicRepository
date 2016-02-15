;+
;NAME:
;  nw_float_to_byte
;PURPOSE:
;  Converts an array of floats in [0.0,1.0] to bytes in [0,255].
;INPUTS:
;  image       - image array of any dimensions; data range should be
;                0.0 to 1.0
;OPTIONAL INPUTS:
;  bits        - number of bits per element; default 8; must be LE 8
;KEYWORDS:
;  none
;OUTPUTS:
;  The float-value image
;REVISION HISTORY:
;  2003-03-10  written - wherry
;  2004-03-20  bits keyword - Hogg
;-
FUNCTION nw_float_to_byte,image,bits=bits
if (NOT keyword_set(bits)) then bits= 8
if (bits GT 8) then begin
    splog, 'warning: bits > 8 not allowed; setting bits=8'
    bits= 8
endif
fmax= 2.0^bits
bmax= byte(2L^bits-1L)
RETURN, byte((floor(image * fmax) > 0) < bmax)
END
