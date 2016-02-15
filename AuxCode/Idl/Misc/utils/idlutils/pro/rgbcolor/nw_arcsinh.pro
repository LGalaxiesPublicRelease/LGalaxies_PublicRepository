;+
;NAME:
;  nw_arcsinh
;PURPOSE:
;  scales the fits image by a degree of nonlinearity specified by user
;INPUTS:
;  colors      - (NXxNYx3) array that contains the R/G/B images
;OPTIONAL INPUTS:
;  nonlinearity- 'b'
;              - b=0 for linear fit, b=Inf for logarithmic
;              - default is 3
;KEYWORDS:
;
;OUTPUTS:
;  The image
;COMMENTS:
;  The input image must be background-subtracted (ie, have zero background).
;BUGS:
;  
;REVISION HISTORY:
;  11/07/03 written - wherry
;  11/12/03 changed radius - wherry
;-
FUNCTION nw_arcsinh,colors,nonlinearity=nonlinearity, radius=radius, $
           inplace=inplace
;set default nonlinearity
IF NOT n_elements(nonlinearity) THEN nonlinearity=3
IF (nonlinearity eq 0.) THEN BEGIN 
    return, colors
ENDIF

dim = size(colors,/dimensions)
NX = LONG(dim[0])
NY = LONG(dim[1])

radius = total(colors,3)
radius = radius+(radius eq 0)
; -------- asinh returns type double by default
fac = float(asinh(radius*nonlinearity)/nonlinearity)/radius

if keyword_set(inplace) then begin 
   for bb= 0,2 do colors[*,*,bb] = colors[*,*,bb]*fac
   return, colors
endif else begin 
   fitted_colors = fltarr(NX, NY, 3, /nozero)
   for bb= 0,2 do fitted_colors[*,*,bb] = colors[*,*,bb]*fac
endelse
RETURN,fitted_colors
END
