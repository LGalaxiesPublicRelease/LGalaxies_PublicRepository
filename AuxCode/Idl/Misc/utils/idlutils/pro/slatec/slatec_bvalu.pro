;+
; NAME:
;   slatec_bvalu
;
; PURPOSE:
;   Evaluate a bspline 
;
; CALLING SEQUENCE:
;   
;    y = slatec_bvalu(x, fullbkpt, coeff, ideriv=ideriv)
;
; INPUTS:
;   x          - Vector of positions to evaluate
;   fullbkpt   - Breakpoint vector returned by SLATEC_EFC()
;   coeff      - B-spline coefficients calculated by SLATEC_EFC()
;
; OPTIONAL KEYWORDS:
;   ideriv     - Derivative to evaluate at x; default to 0
;
; OUTPUTS:
;   y          - Evaluations corresponding to x positions
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   Dynamic link to bvalu_idl in slatec/src/idlwrapper.c,
;   which calls bvalu.f in the library "libslatecidl.so".
;
; REVISION HISTORY:
;   15-Oct-1999  Written by Scott Burles, Chicago
;-
;------------------------------------------------------------------------------

function slatec_bvalu, x, fullbkpt, coeff, ideriv=ideriv

   if (NOT keyword_set(ideriv)) then ideriv = 0L

   ideriv = long(ideriv)

   nbkpt = n_elements(fullbkpt)
   ncoeff = n_elements(coeff)

   k = nbkpt - ncoeff
   n = nbkpt - k
   if (ideriv LT 0 OR ideriv GE k) then begin
      print, 'IDERIV is ', ideriv
      message, 'IDERIV must be >= 0 and < NORD'
   endif

   if (k LT 1) then message, 'NORD must be > 0'
   if (n LT k) then message, 'NBKPT must be >= 2*NORD'

   minbkpt = fullbkpt[k-1]
   maxbkpt = fullbkpt[n]

;   if (min(x) LT minbkpt) then $
;       message, 'x values fall below the first breakpoint'

;   if (max(x) GT maxbkpt) then $
;       message, 'x values fall above the last breakpoint'

   work = fltarr(3*k)
   y = float(x)
   nx = n_elements(x)

   inbv = 1L
   xtemp = float(x < maxbkpt)
   xtemp = (xtemp > minbkpt)

   soname = filepath('libslatec.'+idlutils_so_ext(), $
    root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')

   rr = call_external(soname, 'bvalu_idl', $
    float(fullbkpt), float(coeff), n, k, ideriv, xtemp, nx, inbv, work, y)

   return, y
end
