;+
; NAME:
;   djs_median
;
; PURPOSE:
;   Return the median of an image either with a filtering box or by collapsing
;   the image along one of its dimensions.
;
; CALLING SEQUENCE:
;   result = djs_median( array, [ dimension, width=, boundary=, /idl ] )
;
; INPUTS:
;   array      - N-dimensional array
;
; OPTIONAL INPUTS:
;   dimension  - The dimension over which to compute the median, starting
;                at one.  If this argument is not set, the median of all array
;                elements (or all elements within the median window described
;                by WIDTH) are medianed.
;   width      - Width of median window; scalar value.
;                It is invalid to specify both DIMENSION and WIDTH.
;   boundary   - Boundary condition:
;                'none': Do not median filter within WIDTH/2 pixels of
;                        the edge; this is the default for both this
;                        routine and MEDIAN().
;                'nearest': Use the value of the nearest boundary pixel.
;                        NOT IMPLEMENTED
;                'reflect': Reflect pixel values around the boundary.
;                'wrap': Wrap pixel values around the boundary.
;                        NOT IMPLEMENTED
;                These boundary conditions only take effect if WIDTH is set,
;                and if ARRAY is either 1-dimensional or 2-dimensional.
;
; OUTPUTS:
;   result     - The output array.  If neither DIMENSION nor WIDTH are set,
;                then RESULT is a scalar.  If DIMENSION is not set and WIDTH
;                is set, then RESULT has the same dimensions as ARRAY.
;                If DIMENSION is set and WIDTH is not
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The DIMENSION input is analogous to that used by the IDL built-in
;   function TOTAL.
;
;   I should like to add the functionality of having WIDTH be an N-dimensional
;   smoothing box.  For example, one should be able to median a 2-D image
;   with a 3x5 filtering box.
;
; EXAMPLES:
;   Create a 2-D image and compute the median of the entire image:
;   > array = findgen(100,200)
;   > print, djs_median(array)
;
;   Create a data cube of 3 random-valued 100x200 images.  At each pixel in
;   the image, compute the median of the 3:
;   > array = randomu(123,100,200,3)
;   > medarr = djs_median(array,3)
;
;   Create a random-valued 2-D image and median-filter with a 9x9 filtering box:
;   > array = randomu(123,100,200)
;   > medarr = djs_median(array,9)
;
; BUGS:
;   The C routine only supports type FLOAT.
;
; PROCEDURES CALLED:
;   Dynamic link to arrmedian.c
;
; REVISION HISTORY:
;   06-Jul-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function djs_median, array, dim, width=width, boundary=boundary, idl=idl

   ; Need at least 1 parameter
   if (N_params() LT 1) then begin
      print, 'Syntax - result = djs_median( array, [ dimension, width= ] )'
      return, -1
   endif

   if (NOT keyword_set(boundary)) then boundary = 'none'

   dimvec = size(array, /dimensions)
   ndim = N_elements(dimvec)

   if (NOT keyword_set(dim) AND NOT keyword_set(width)) then begin

      if (n_elements(array) EQ 1) then medarr = array[0] $
       else medarr = median(array, /even)

   endif else if (NOT keyword_set(dim)) then begin

      if (boundary EQ 'none') then begin
         npix = n_elements(array)
         if (npix EQ 1) then medarr = array[0] $
          else if (width EQ 1) then medarr = array $
          else medarr = median(array, (width < npix), /even)
      endif else begin
         padsize = ceil(width/2.)  ; This padding will be at least 1 pixel
         zero = array[0] - array[0] ; Zero in the type of ARRAY
         if (ndim EQ 1) then begin
            bigarr = bytarr(dimvec[0]+2*padsize) + zero
            bigarr[padsize:padsize+dimvec[0]-1] = array
         endif else if (ndim EQ 2) then begin
            bigarr = bytarr(dimvec[0]+2*padsize, dimvec[1]+2*padsize) + zero
            bigarr[padsize:padsize+dimvec[0]-1, padsize:padsize+dimvec[1]-1] $
             = array
         endif else begin
            message, 'Unsupported number of dimensions with this b.c.'
         endelse

         if (ndim EQ 1) then begin
            bigarr[0:padsize-1] = reverse(array[0:padsize-1])
            bigarr[padsize+dimvec[0]:padsize*2+dimvec[0]-1] = $
             array[dimvec[0]-padsize:dimvec[0]-1]

            if (width GT 1) then $
             bigarr = temporary( median(bigarr, width, /even) )
            medarr = bigarr[padsize:padsize+dimvec[0]-1]

         endif else begin

            case boundary of
            'nearest': begin
               message, 'This boundary condition not implemented.'
            end
            'reflect': begin

               ; Copy into left + right
               bigarr[0:padsize-1,padsize:dimvec[1]+padsize-1] = $
                reverse(array[0:padsize-1,*],1)
               bigarr[padsize+dimvec[0]:padsize*2+dimvec[0]-1, $
                 padsize:dimvec[1]+padsize-1] = $
                reverse(array[dimvec[0]-padsize:dimvec[0]-1,*],1)

               ; Copy into bottom + top
               bigarr[padsize:dimvec[0]+padsize-1,0:padsize-1] = $
                reverse(array[*,0:padsize-1],2)
               bigarr[padsize:padsize+dimvec[0]-1, $
                padsize+dimvec[1]:padsize*2+dimvec[1]-1] = $
                reverse(array[*,dimvec[1]-padsize:dimvec[1]-1],2)

               ; Copy into lower left
               bigarr[0:padsize-1,0:padsize-1] = $
                reverse(reverse(array[0:padsize-1,0:padsize-1],1),2)

               ; Copy into lower right
               bigarr[padsize+dimvec[0]:padsize*2+dimvec[0]-1,0:padsize-1] = $
                reverse(array[dimvec[0]-padsize:dimvec[0]-1,0:padsize-1],2)

               ; Copy into upper left
               bigarr[0:padsize-1,padsize+dimvec[1]:padsize*2+dimvec[1]-1] = $
                reverse(reverse(array[0:padsize-1, $
                dimvec[1]-padsize:dimvec[1]-1],1),2)

               ; Copy into upper right
               bigarr[padsize+dimvec[0]:padsize*2+dimvec[0]-1, $
                padsize+dimvec[1]:padsize*2+dimvec[1]-1] = $
                reverse(reverse(array[dimvec[0]-padsize:dimvec[0]-1, $
                dimvec[1]-padsize:dimvec[1]-1],1),2)

               if (width GT 1) then $
                bigarr = temporary( median(bigarr, width, /even) )
               medarr = bigarr[padsize:padsize+dimvec[0]-1, $
                padsize:padsize+dimvec[1]-1]

            end
            'wrap': begin
               message, 'This boundary condition not implemented.'
            end
            endcase
         endelse
      endelse

   endif else if (NOT keyword_set(width)) then begin

;      if (ndim LE 1) then begin
;         message, 'ARRAY must be multi-dimensional if DIM is specified'
;      endif
      if (dim GT ndim OR dim LT 1) then begin
         message, 'DIM must be between 1 and '+string(ndim)+' inclusive'
      endif

      if (ndim EQ 1) then begin
         medarr = median(array)
      endif else begin
         ; Allocate memory for the output array
         newdimvec = dimvec[ where(lindgen(ndim)+1 NE dim) ]
         newsize = N_elements(array) / dimvec[dim-1]
         medarr = reform(fltarr(newsize), newdimvec)

         if keyword_set(idl) then begin
           medarr = median(array, dimension=dim, /even)
         endif else begin
           soname = filepath('libmath.'+idlutils_so_ext(), $
            root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')
           retval = call_external(soname, 'arrmedian', $
            ndim, dimvec, float(array), long(dim), medarr)
         endelse
         
      endelse

   endif else begin
      message, 'Invalid to specify both DIMENSION and WIDTH'
   endelse

   return, medarr
end
;------------------------------------------------------------------------------
