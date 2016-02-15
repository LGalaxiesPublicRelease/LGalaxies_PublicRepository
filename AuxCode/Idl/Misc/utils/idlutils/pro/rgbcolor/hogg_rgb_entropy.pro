;+
; NAME:
;   hogg_rgb_entropy
; PURPOSE:
;   Compute "entropy" of an RGB representation of an image
; INPUTS:
;   colors    - [nx,ny,3] image of values in range [0.0,1.0]
;   bits      - number of bits to assign per pixel per color (default 4)
; OUTPUT:
;   entropy   - entropy of the representation
; BUGS:
; REVISION HISTORY:
;   2004-03-20  started - Hogg
;-
function hogg_rgb_entropy, colors,bits
if (NOT keyword_set(bits)) then bits= 4
bcolors= nw_float_to_byte(colors,bits=bits)
dim= size(colors,/dimens)
lcolors= lonarr(dim[0],dim[1])
lcolors[*,*]= bcolors[*,*,0]+2L^bits*(bcolors[*,*,1]+2L^bits*bcolors[*,*,2])
bcolors= 0
ncolor= 2L^(3L*long(bits))
npixel= long(dim[0]*dim[1])
entropy= 0D0
for cc=0L,ncolor-1L do begin
    foo= where(lcolors EQ cc,ncc)
    if (ncc GT 0) then begin
        entropy= entropy-double(ncc)*alog10(double(ncc)/double(npixel))
    endif
endfor
return, entropy
end
