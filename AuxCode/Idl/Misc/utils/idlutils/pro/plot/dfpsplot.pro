;+
; NAME:
;   dfpsplot
;
; PURPOSE:
;   Finkbeiner's routine to open a PostScript file for plotting commands.
;   Close with DFPSCLOSE.
;
; CALLING SEQUENCE:
;   dfpsplot, filename, [/square, /landscape, ysize=ysize, $
;    /encap, /color, _EXTRA=KeywordsForDevice ]
;
; INPUTS:
;   filename   - File name to open
;
; OPTIONAL INPUTS:
;   square     - Make the plotting area square.
;   landscape  - Use landscape paper; default is to use portrait.
;   ysize      - For portrait mode, the YSIZE can be changed from its
;                default of 8-inches.  For landscape mode, the value
;                of YSIZE is ignored.
;   encap      - Force non-encapsulated file unless this keyword is set.
;   color      - Force non-color file unless this keyword is set.
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   The-Beginning-of-Time  Written by Doug Finkbeiner, Berkeley.
;   05-Sep-1999  Modified and commented by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro dfpsplot, filename, square=square, landscape=landscape, ysize=ysize, $
 xsize=xsize, encap=encap, color=color, _EXTRA=KeywordsForDevice

   if (N_params() LT 1) then begin
      print, 'Syntax - dfpsplot, filename, [/square, /landscape, ysize=ysize, $,'
      print, ' /encap, /color, _EXTRA=KeywordsForDevice ]'
      return
   end

   if (N_elements(encap) EQ 0) then encap=0
   if (N_elements(color) EQ 0) then color=0

   if (keyword_set(landscape)) then begin
      xs = 10.5                          ; full page
      IF keyword_set(square) THEN xs = 8 ; square page
      IF keyword_set(xsize) THEN xs = xsize
      set_plot, 'PS'
      device, file=filename, /landscape, xsize=xs, ysize=8, $
       xoff=0, yoff=11.5+(10-xs), /inch, encap=encap, color=color, $
       _EXTRA=KeywordsForDevice

   endif else begin
      xs = 8.0
      ys = 10.0                          ;full page
      IF keyword_set(square) THEN ys = 8 ; square page
      IF keyword_set(ysize) THEN ys = ysize
      IF keyword_set(xsize) THEN xs = xsize
      set_plot, 'PS'
      device, file=filename, /portrait, xsize=xs, ysize=ys, $
       xoff=0, yoff=0.5+(10-ys), /inch, encap=encap, color=color, $
       _EXTRA=KeywordsForDevice

   endelse

   return
end
;------------------------------------------------------------------------------
