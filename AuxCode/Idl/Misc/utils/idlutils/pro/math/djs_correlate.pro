;+
; NAME:
;   djs_correlate
;
; PURPOSE:
;   Compute a cross-correlation function using weights (or masks).
;
; CALLING SEQUENCE:
;   result = djs_correlate( x, y, [ lags, xweight=, yweight= ] )
;
; INPUTS:
;   x          - Vector
;   y          - Vector, which may have a different number of elements from X
;
; OPTIONAL INPUTS:
;   lags       - A scalar or integer vector specifying the lags at which
;                to compute the cross-correlation; default to one lag at 0.
;   xweight    - Vector of weights for X; default to 1 for all points
;   yweight    - Vector of weights for Y; default to 1 for all points
;
; OUTPUTS:
;   result     - The output vector, with one result per LAG value.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine is similar to the IDL routine C_CORRELATE(), but with a
;   few notable differences.  The X and Y vectors do not wrap when they are
;   shifted, but rather only overlapping elements are compared at each lag.
;   Because of this, X and Y do not have to have the same number of dimensions.
;   A weight (or mask) can be assigned to each element of X,Y using the
;   XWEIGHT,YWEIGHT keywords.  These weights can effectively be used to
;   mask out regions of each vector by setting the weight to 1 for good
;   pixels and 0 for bad ones.
;
;   Each pixel of both X and Y are effectively weighted by XWEIGHT*YWEIGHT
;   appropriately shifted before the multiplication.
;
; EXAMPLES:
;
; BUGS:
;   The C routine only supports type FLOAT.
;
; PROCEDURES CALLED:
;   Dynamic link to ccorrelate.c
;
; REVISION HISTORY:
;   07-Jul-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function djs_correlate, x, y, lags, xweight=xweight, yweight=yweight

   ; Need at least 2 parameters
   if (N_params() LT 2) then begin
      print, 'Syntax - result = djs_correlate( x, y, [ lags, xweight=, yweight= ] )'
      return, -1
   endif

   if (n_elements(lags) EQ 0) then lags = 0

   nlag = n_elements(lags)
   nx = n_elements(x)
   ny = n_elements(y)

   if (NOT keyword_set(xweight)) then xweight = fltarr(nx) + 1 $
    else if (n_elements(xweight) NE nx) then $
    message, 'X and XWEIGHT have inconsistent dimensions'
   if (NOT keyword_set(yweight)) then yweight = fltarr(ny) + 1 $
    else if (n_elements(yweight) NE ny) then $
    message, 'Y and YWEIGHT have inconsistent dimensions'

   ; Allocate memory for the output array
   result = fltarr(n_elements(lags))

   soname = filepath('libmath.'+idlutils_so_ext(), $
    root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')
   retval = call_external(soname, 'ccorrelate', $
    nx, float(x), float(xweight), ny, float(y), float(yweight), $
    nlag, long(lags), result)

   return, result
end
;------------------------------------------------------------------------------
