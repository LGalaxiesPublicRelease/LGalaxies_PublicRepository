;-----------------------------------------------------------------------
;+
; NAME:
;   djs_icolor
;
; PURPOSE:
;   Internal routine for converting a color name to an index.
;
; CALLING SEQUENCE:
;   icolor = djs_icolor(color)
;
; INPUT:
;   color:
;
; OUTPUTS:
;   icolor
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Written by D. Schlegel, 27 September 1997, Durham
;-
;-----------------------------------------------------------------------
function djs_icolor, color

   if (NOT keyword_set(color)) then color = !p.color

   ncolor = N_elements(color)

   ; If COLOR is a string or array of strings, then convert color names
   ; to integer values
   if (size(color,/tname) EQ 'STRING') then begin ; Test if COLOR is a string

      ; Detemine the default color for the current device
      if (!d.name EQ 'X') then defcolor = 27 $ ; white for X-windows
       else defcolor = 0 ; black otherwise

      ; Load a simple color table with the basic 27 colors
;      red   = long(255.0*(indgen(8) mod 2)/1.0)
;      green = long(255.0*(fix(findgen(8)/2) mod 2)/1.0)
;      blue  = long(255.0*fix(findgen(8)/4)/1.0)
;print, red, green, blue
      red   = long(255.0*(indgen(27) mod 3)/2.0)
      green = long(255.0*(fix(findgen(27)/3) mod 3)/2.0)
      blue  = long(255.0*fix(findgen(27)/9)/2.0)
      red = [red, 70, 200]
      green = [green, 70, 200]
      blue = [blue, 70, 200]

;for ii=0,26 do print, ii,red[ii],green[ii],blue[ii]

      tvlct, red, green, blue
      icolor = 0 * (color EQ 'black') $
             + 1 * (color EQ 'dark red') $
             + 1 * (color EQ 'brown') $
             + 2 * (color EQ 'red') $
             + 3 * (color EQ 'dark green') $
             + 4 * (color EQ 'dark yellow') $
             + 5 * (color EQ 'yellow red') $
             + 5 * (color EQ 'orange') $
             + 6 * (color EQ 'green') $
             + 7 * (color EQ 'yellow green') $
             + 8 * (color EQ 'yellow') $
             + 9 * (color EQ 'dark blue') $
             + 9 * (color EQ 'navy') $
             + 10 * (color EQ 'dark magenta') $
             + 10 * (color EQ 'purple') $ 
             + 11 * (color EQ 'magenta red') $
             + 12 * (color EQ 'dark cyan') $
             + 13 * (color EQ 'grey') $
             + 13 * (color EQ 'gray') $
             + 14 * (color EQ 'light red') $
             + 14 * (color EQ 'pink') $
             + 15 * (color EQ 'cyan green') $
             + 16 * (color EQ 'light green') $
             + 17 * (color EQ 'light yellow') $
             + 18 * (color EQ 'blue') $
             + 19 * (color EQ 'magenta blue') $
             + 20 * (color EQ 'magenta') $
             + 21 * (color EQ 'cyan blue') $
             + 22 * (color EQ 'light blue') $
             + 23 * (color EQ 'light magenta') $
             + 24 * (color EQ 'cyan') $
             + 25 * (color EQ 'light cyan') $
             + 26 * (color EQ 'white') $
             + 27 * (color EQ 'dark grey') $
             + 27 * (color EQ 'dark gray') $
             + 28 * (color EQ 'light grey') $
             + 28 * (color EQ 'light gray') $
             + defcolor * (color EQ 'default')

      if(strupcase(!D.NAME) eq 'X' OR strupcase(!D.NAME) eq 'WIN') then $
     device, get_decompose=de_compose $
        else $
        de_compose=1
      if((de_compose EQ 1) AND (!d.N_colors EQ 16777216)) then begin
        colors = red + ishft(green,8) + ishft(blue,16)
        icolor = colors[icolor]
      endif

   endif else begin
      icolor = color
   endelse

   return, icolor
end 
;-----------------------------------------------------------------------
