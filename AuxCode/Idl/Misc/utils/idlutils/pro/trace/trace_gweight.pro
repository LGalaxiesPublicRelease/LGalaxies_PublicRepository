;+
; NAME:
;   trace_gweight
;
; PURPOSE:
;   Recenter a trace using gaussian-weighted centers.
;
; CALLING SEQUENCE:
;   xnew = trace_fweight( fimage, xcen, ycen, [sigma=sigma, xerr=xerr, 
;    invvar=invvar] )
;
; INPUTS:
;   fimage     - Image
;   xcen       - Initial guesses for X centers
;   ycen       - Y positions corresponding to "xcen" (expected as integers)
;
; OPTIONAL KEYWORDS:
;   radius     - Sigma in pixels; default to 1.0
;   invvar     - Inverse variance of image used only in computing errors XERR.
;                If not set, then INVVAR=1 is used.
;
; OUTPUTS:
;   xnew       - New X centers
;
; OPTIONAL OUTPUTS:
;   xerr       - Formal errors for XNEW; set equal to 999.0 if there are any
;                masked pixels in a centroiding region (e.g., if INVVAR=0)
;
; COMMENTS:
;
; EXAMPLES:
;
; REVISION HISTORY:
;   17-Jan-2000  Written by Scott Burles, Chicago
;-
;------------------------------------------------------------------------------
function trace_gweight, fimage, xcen, ycen, sigma=sigma, xerr=xerr, $
 invvar=invvar, idl=idl

   ; Need 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - xnew = trace_gweight( fimage, xcen, ycen, $'
      print, '       [sigma=sigma, xerr=xerr, invvar=invvar] )'
      return, -1
   endif
   if (NOT keyword_set(sigma)) then sigma = 1.0

   nx = (size(fimage))[1]
   ny = (size(fimage))[2]
   ncen = N_elements(xcen)
   xnew = float(xcen)
   xerr = 0.0 * xnew ; Allocate memory

   if (NOT keyword_set(invvar)) then invvar = 0.0 * fimage + 1.0

   if keyword_set(idl) then begin

     lower = xcen - 3.0*sigma
     upper = xcen + 3.0*sigma
     x_int = round(xcen)
     nstep = 2*long(3.0*sigma) - 1
     i1 = x_int - nstep/2
     i2 = i1 + nstep - 1

     weight = xcen *0.
     numer  = xcen *0
     meanvar = xcen *0
     bad = lonarr(ncen)
     xtemp = (findgen(nstep)-nstep/2)/5.  * sigma

     for i=0, nstep-1 do begin
       xh = x_int - nstep/2 + i
       xtemp = (xh - xcen - 0.5)/sigma/sqrt(2.0)
       g_int = (errorf(xtemp+1./sigma/sqrt(2.0)) - errorf(xtemp))/2.
       xs = (xh > 0) < (nx -1)
       cur_weight = fimage[xs, ycen] * (invvar[xs, ycen] GT 0) * g_int * $
                                       (xh GE 0 AND xh LT nx) 
       weight = weight + cur_weight 
       numer = numer + cur_weight * xh
       meanvar = meanvar + cur_weight * cur_weight * (xcen-xh)^2 / $
                            (invvar[xs, ycen] + (invvar[xs, ycen] EQ 0))
       bad = bad OR (xh LT 0 OR xh GE nx)
     endfor

     xerr[*] = 999.0
     good = where(bad EQ 0 AND weight GT 0)
     if good[0] NE -1 then begin
        xnew[good] = numer[good]/weight[good]
        xerr[good] = sqrt(meanvar[good])/weight[good]
     endif
     
   endif else begin

     soname = filepath('libtrace.'+idlutils_so_ext(), $
      root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')
     result = call_external(soname, 'trace_gweight', $
      nx, ny, float(fimage), float(invvar), $
      float(sigma), ncen, xnew, long(ycen), xerr)
   endelse

   return, xnew
end
;------------------------------------------------------------------------------
