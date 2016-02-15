;+
; NAME:
;   trace_crude
;
; PURPOSE:
;   Create a crude trace set given one position (eg, a center) in each trace.
;
; CALLING SEQUENCE:
;   xset = trace_crude( fimage, [ invvar, xstart=, ystart=, radius=, yset=, $
;    nave=, nmed=, thresh=, maxerr=, maxshifte=, maxshift0=, xerr=, /double ] )
;
; INPUTS:
;   fimage     - Image
;
; OPTIONAL INPUTS:
;   invvar     - Inverse variance (weight) image
;   xstart     - Initial guesses for X centers (one for each trace).
;                If not set, then this code searches for all peaks at YSTART.
;   ystart     - Y positions corresponding to "xstart" (expected as integers).
;                There are three options for this parameter:
;                (1) One element of YSTART for each value of XSTART,
;                (2) A scalar value that is used for every XSTART, or
;                (3) Not set, in which case the center row is used.
;   radius     - Radius for centroiding; default to 3.0
;   nmed       - Median filtering size down columns before performing trace;
;                default to 1
;   nave       - Averaging size down columns before performing trace;
;                default to 5
;   thresh     - Threshold for initial peak finding; if not set, then use
;                1.0 times the median of the row(s) used for the initial peaks.
;   maxerr     - Maximum error in centroid allowed for valid recentering;
;                default to 0.2
;   maxshifte  - Maximum shift in centroid allowed for valid recentering;
;                default to 0.1
;   maxshift0  - Maximum shift in centroid allowed for initial row;
;                default to 0.5
;   double     - If set, then return values are double-precision; values are
;                already double-precision if FIMAGE or XSTART already are
;
; OUTPUTS:
;   xset       - X centers for all traces
;
; OPTIONAL OUTPUTS:
;   yset       - Y centers for all traces
;   xerr       - Errors for XSET
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   This linked C code is always single-precision, even if /DOUBLE is set.
;
; PROCEDURES CALLED:
;   djs_laxisgen()
;   splog
;
;   Dynamic link to trace_crude.c
;
; REVISION HISTORY:
;   14-May-1999  Written by David Schlegel, Princeton.
;   12-Jul-1999  Added optional output YSET (DJS).
;   06-Aug-1999  Added optional outpust XERR (DJS).
;-
;------------------------------------------------------------------------------
function trace_crude, fimage, invvar, xstart=xstart1, ystart=ystart1, $
 radius=radius1, yset=yset, nave=nave1, nmed=nmed1, thresh=thresh, $
 maxerr=maxerr1, maxshifte=maxshift_in, maxshift0=maxshift0_in, xerr=xerr, $
 double=double1, idl=idl

   ; Need 1 parameter
   if (N_params() LT 1) then begin
      print, 'Syntax - xset = trace_crude( fimage, [ invvar, xstart=, ystart=, $'
      print, ' radius=, nave=, nmed=, maxerr=, maxshifte=, maxshift0=, xerr= ] )'
      return, -1
   endif

   ndim = size(fimage, /n_dimen)
   if (ndim NE 2) then message, 'FIMAGE must be 2-dimensional image'
   dims = size(fimage, /dimens)
   nx = dims[0]
   ny = dims[1]
   if (keyword_set(ystart1)) then ystart = ystart1 $
    else ystart = long(ny/2)
   if (keyword_set(radius1)) then radius = radius1 $
    else radius = 3.0
   if (keyword_set(nmed1)) then nmed = nmed1 > 1 $
    else nmed = 1
   if (keyword_set(nave1)) then nave = nave1 < ny $
    else nave = 5 < ny
   if (keyword_set(maxerr1)) then maxerr = maxerr1 $
    else maxerr = 0.2
   if (keyword_set(maxshift_in)) then maxshift = maxshift_in $
    else maxshift = 0.1
   if (keyword_set(maxshift0_in)) then maxshift0 = maxshift0_in $
    else maxshift0 = 0.5
   if (size(fimage,/tname) EQ 'DOUBLE' OR size(xstart,/tname) EQ 'DOUBLE') $
    OR keyword_set(double1) then double = 1B

   ; Make a copy of the image and error map
   imgtemp = float(fimage)
   if (keyword_set(invvar)) then begin
      invtemp = float(invvar)
   endif else begin
      invtemp = 1.0 / (imgtemp > 1)
   endelse

   ; Median filter the entire image along columns by NMED rows
   if (nmed GT 1) then $   
    for ix=1, nx-1 do imgtemp[ix,*] = median(transpose(imgtemp[ix,*]), nmed)

   ; Boxcar-sum the entire image along columns by NAVE rows
   if (nave GT 1) then begin
      kernel = transpose(intarr(nave) + 1.0)
      if (nx EQ 1) then $
       imgconv = reform(convol(reform(imgtemp*invtemp,ny), $
        reform(kernel,nave), /edge_truncate), 1, ny) $
      else $
       imgconv = convol(imgtemp*invtemp, kernel, /edge_truncate)

      ; Add the variances
      if (nx EQ 1) then $
       invtemp = reform(convol(reform(invtemp,ny), $
        reform(kernel,nave), /edge_truncate), 1, ny) $
      else $
       invtemp = convol(invtemp, kernel, /edge_truncate)

      ; Look for pixels with infinite errors - replace with original values
      ibad = where(invtemp EQ 0, nbad)
      if (nbad GT 0) then begin
         invtemp[ibad] = 1
         imgconv[ibad] = imgtemp[ibad]
      endif

      ; Renormalize the summed image by the weights
      imgtemp = imgconv / invtemp
   endif

   if (keyword_set(xstart1)) then begin
      xstart = xstart1
   endif else begin
      ; Automatically find peaks for XSTART

      ; Extract NSUM rows from the image at YSTART
      nsum = 1
      imrow = $
       imgtemp[*,long(ystart[0]-0.5*(nsum-1)):long(ystart[0]+0.5*(nsum-1))]
      imrow = rebin(imrow, nx, 1)

      if (keyword_set(thresh)) then mthresh = thresh $
       else mthresh = median(imrow)

      ; Boxcar smooth along X
      imrow = smooth(imrow,radius/2 > 3)

      ; Find all local peaks that are also above THESH times the median
      if (nx GT 1) then begin
         rderiv = imrow[1:nx-1] - imrow[0:nx-2]
         izero = where( rderiv[0:nx-3] GT 0 AND rderiv[1:nx-2] LE 0 $
          AND imrow[1:nx-2] GT mthresh, nzero)
      endif else begin
         nzero = 0
      endelse
      if (nzero EQ 0) then begin
         yset = 0
         xset = 0
         xerr = 0
         splog, 'Warning: No peaks found'
         return, 0
      endif

      xstart = izero + 0.5 + rderiv[izero] / (rderiv[izero] - rderiv[izero+1])
   endelse

   if (N_elements(ystart) EQ 1) then $
    ypass = replicate(long(ystart), N_elements(xstart)) $
    else ypass = long(ystart)

   if (N_elements(xstart) NE N_elements(ypass)) then $
    message, 'Wrong number of elements for YSTART'

   ntrace = N_elements(xstart)
   xset = fltarr(ny, ntrace)
   xerr = fltarr(ny, ntrace)

   if keyword_set(idl) then begin
    ;
    ;  Need xset and xerr
    ; 

     trace_crude_idl, imgtemp, invtemp, radius, xstart, ypass, $
                           xset, xerr, maxerr, maxshift, maxshift0

   endif else begin
     soname = filepath('libtrace.'+idlutils_so_ext(), $
      root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')

     result = call_external(soname, 'trace_crude', $
      nx, ny, imgtemp, invtemp, float(radius), ntrace, float(xstart), ypass, $
      xset, xerr, float(maxerr), float(maxshift), float(maxshift0))
   endelse

   if (keyword_set(double)) then begin
      xset = double(xset)
      xerr = double(xerr)
   endif

   yset = djs_laxisgen([ny,nTrace], iaxis=0)

   return, xset
end
;------------------------------------------------------------------------------
