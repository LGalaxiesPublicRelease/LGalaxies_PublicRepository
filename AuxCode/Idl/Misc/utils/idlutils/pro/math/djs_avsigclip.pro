;+
; NAME:
;   djs_avsigclip
;
; PURPOSE:
;   Average multiple images with sigma-rejection.
;
; CALLING SEQUENCE:
;   result = djs_avsigclip( array, [ dimension, sigrej=, maxiter=, $
;    inmask=, outmask= ] )
;
; INPUTS:
;   array      - N-dimensional array
;
; OPTIONAL INPUTS:
;   dimension  - The dimension over which to collapse the data.  If ommitted,
;                assume that the last dimension is the one to collapse.
;   sigrej     - Sigma for rejection; default to 3.0.
;   maxiter    - Maximum number of sigma rejection iterations.  One iteration
;                does no sigma rejection; default to 10 iterations.
;   inmask     - Input mask, setting =0 for good elements
;
; OUTPUTS:
;   result     - The output array.
;   outmask    - Output mask, setting =0 for good elements, =1 for bad.
;                Any pixels masked in INMASK are also masked in OUTMASK.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The DIMENSION input is analogous to that used by the IDL built-in
;   function TOTAL.
;
; EXAMPLES:
;   Create a data cube of 10 random-valued 100x200 images.  At each pixel in
;   the image, compute the average of the 10 values, but rejecting 3-sigma
;   outliers:
;   > array = randomu(123,100,200,10)
;   > ave = djs_avsigclip(array, sigrej=3)
;
;
;   If all points are masked in any given vector or array, a mean and
;   dispersion are computed for all the points.  Is this the behaviour we want?
;   If you want to replace those values with zeros instead, look at OUTMASK:
;   > array = randomu(123,100,200)
;   > inmask = bytarr(100,200)
;   > inmask[*,8] = 1 ; mask all of row #8
;   > ave = djs_avsigclip(array, 1, inmask=inmask, outmask=outmask)
;   > ibad = where( total(1-outmask, 1) EQ 0)
;   > if (ibad[0] NE -1) then ave[ibad] = 0 ; zero-out bad rows
;
; BUGS:
;
; PROCEDURES CALLED:
;   Dynamic link to arravsigclip.c
;
; REVISION HISTORY:
;   07-Jul-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function djs_avsigclip, array, dim, sigrej=sigrej, maxiter=maxiter, $
 inmask=inmask, outmask=outmask

   ; Need at least 1 parameter
   if (N_params() LT 1) then begin
      print, 'Syntax - result = djs_avsigclip( array, [ dimension, sigrej=, maxiter=, $'
      print, ' inmask=, outmask= ] ); '
      return, -1
   endif

   if (NOT keyword_set(dim)) then dim = size(array, /n_dim)
   if (NOT keyword_set(sigrej)) then sigrej = 3.0
   if (N_elements(maxiter) EQ 0) then maxiter = 10

   sz = N_elements(array)
   dimvec = size(array, /dimensions)
   ndim = N_elements(dimvec)

   if (dim GT ndim OR dim LT 1) then begin
      message, 'DIM must be between 1 and '+string(ndim)+' inclusive'
   endif

   ; Allocate memory for the output array
   if (ndim GT 1) then $
    newdimvec = dimvec[ where(lindgen(ndim)+1 NE dim) ] $
   else $
    newdimvec = [1]
   newsize = N_elements(array) / dimvec[dim-1]
   avearr = reform(fltarr(newsize), newdimvec)

   soname = filepath('libmath.'+idlutils_so_ext(), $
    root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')

   if (arg_present(inmask) OR arg_present(outmask)) then begin

      if (NOT keyword_set(inmask)) then inmask = reform(bytarr(sz), dimvec)
      outmask = reform(bytarr(sz), dimvec)

      retval = call_external(soname, 'arravsigclip_mask', $
       ndim, dimvec, float(array), long(dim), float(sigrej), float(sigrej), $
       long(maxiter), avearr, byte(inmask), outmask)

   endif else begin

      retval = call_external(soname, 'arravsigclip', $
       ndim, dimvec, float(array), long(dim), float(sigrej), float(sigrej), $
       long(maxiter), avearr)

   endelse

   return, avearr
end
;------------------------------------------------------------------------------
