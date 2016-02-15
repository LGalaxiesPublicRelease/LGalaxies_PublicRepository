;+
; NAME:
;   pixfont
;
; PURPOSE:
;   pixelize a font for inclusion in a pixelized image
;
; CALLING SEQUENCE:
;   im = pixfont(str, [ xsize=, ysize=, charsize=, align= ])
;
; INPUTS:
;   str      - String(s) to pixelize; if this is an array
;   
; KEYWORDS:
;   xsize    - xsize of pixelized area; default to 1024 pixels
;   ysize    - ysize of pixelized area; default to 200 pixels
;   charsize - character size (default ysize/10)
;   align    - Alignment, 0=left, 0.5=center, 1=right; default to 0 (left)
;
; OUTPUTS:
;   im       - Byte image of dimensions XSIZE,YSIZE 
;
; RESTRICTIONS:
;   If you make the font too big, it will just run off the edge. 
; 
; EXAMPLES:
;   One string, default size:
;     im = pixfont('Hello, world!')
;   same with smaller letters:
;     im = pixfont('Hello, world!',charsize=10)
;   Two lines of text:
;     im = pixfont(['Hello, world!','- Anonymous programmer'], $
;                     charsize=10,align=1)
;
; REVISION HISTORY:
;   2003-May-23  Written by Douglas Finkbeiner, Princeton
;   2007-May-15  moved from photoop to idlutils - DPF
;
;----------------------------------------------------------------------
function pixfont, str, xsize=xsize, ysize=ysize, charsize=charsize, $
                  align=align

  if NOT keyword_set(align) then align = 0
  if NOT keyword_set(xsize) then xsize = 1024
  if NOT keyword_set(ysize) then ysize = 200
  
  !p.font = 1
  if NOT keyword_set(charsize) then charsize = ysize/10
  window, 32, /free, /pixmap, xsize=xsize, ysize=ysize
  nstr = n_elements(str)
  for i=0, nstr-1 do begin 
     xpos = 0.01 > align < 0.99 ; avoid edges
     xyouts, xpos, float(i+0.25)/nstr, str[nstr-i-1], charsize=charsize, $
       align=align, charthick=3
  endfor 
  im = tvrd()

  return, im
end
