;+
; NAME:
;   labelloc
; PURPOSE:
;   convert fractional positions in X and Y to those appropriate for xyouts
; CALLING SEQUENCE:
;   labelloc, xfrac, yfrac, xloc, yloc
; INPUTS:
;   xfrac,yfrac - fractional distances from axes in x and y
; OUTPUTS:
;   xloc, yloc - units appropriate for xyouts
; REVISION HISTORY:
;   2003-11-21  started - Blanton
;-
pro labelloc, xfrac, yfrac, xloc, yloc

xloc=!X.CRANGE[0]+(!X.CRANGE[1]-!X.CRANGE[0])*xfrac
yloc=!Y.CRANGE[0]+(!Y.CRANGE[1]-!Y.CRANGE[0])*yfrac

end
