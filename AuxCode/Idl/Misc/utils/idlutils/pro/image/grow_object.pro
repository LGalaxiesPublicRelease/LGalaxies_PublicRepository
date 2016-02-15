;+
; NAME:
;   grow_object
;
; PURPOSE:
;   Identify objects as the contiguous non-zero pixels in an image.
;
; CALLING SEQUENCE:
;   mask = grow_object( image, [ xstart=, ystart=, putval=, /diagonal, nadd= ]
;
; INPUTS:
;   image      - Integer-valued image vector or array, where non-zero pixel
;                values indicate that an object touches that pixel.
;
; OPTIONAL INPUTS:
;   xstart     - Starting X position(s) for assembling the object; default to
;                settting all pixels where IMAGE != 0.
;   ystart     - Starting Y position(s) for assembling the object; default to
;                settting all pixels where IMAGE != 0.
;   putval     - Object ID(s) to put in MASK as positive-valued long integer;
;                default to a unique integer (starting at 1) for each object.
;                This can either be a scalar, or a vector with one element
;                per XSTART,YSTART position.
;   diagonal   - If set, then consider diagonally-offset pixels as contigous
;                as well as pixels simply to the left, right, down, or up.
;
; OUTPUTS:
;   mask       - Mask with object IDs; zeros indicate that there is no object
;                in that pixel, and positive values are used as object IDs.
;                Negative values are not allowed.
;
; OPTIONAL OUTPUTS:
;   nadd       - Number of pixels added to all objects
;
; COMMENTS:
;   Find the pixels that make up an "object" as the contiguous non-zero
;   pixels in IMAGE that touch the pixel XSTART,YSTART.  All such pixels
;   have MASK set to PUTVAL.
;
;   If XSTART,YSTART,PUTVAL are not specified, then all objects are found
;   in the image and assigned unique object IDs in MASK starting at 1.
;   Note that in this case, max(MASK) is the number of objects.
;
;   The memory usage is 9*(nx+2)*(ny+2) bytes in addition to the input
;   image, where [nx,ny] are the dimensions of the input image.
;
; EXAMPLES:
;   Create a random image of 0s and 1s, and identify all contiguous pixels
;   as objects:
;   IDL> image=smooth(randomu(123,100,100),5) GT 0.55 & mask = 0
;   IDL> mask = grow_object(image)
;
; BUGS:
;
; PROCEDURES CALLED:
;   Dynamic link to grow_obj.c
;
; REVISION HISTORY:
;   20-May-2003  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
function grow_obj1, image, mask, iloc, putval=putval, $
 diagonal=diagonal, nx=nx, ny=ny, workarray=workarray

   ; Don't bother calling the C code if we know it won't do anything
   if (image[iloc] EQ 0) then return, 0

   if (putval LT 0) then message, 'PUTVAL cannot be negative'

   soname = filepath('libimage.'+idlutils_so_ext(), $
    root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')

   ; Note that NY is not actually used by the C code, because the code
   ; quits when it no longer finds any masked pixels.  Because the image
   ; has been padded around the edges with zeros, we know that condition
   ; will be met without ever having to explicitly test that we've hit
   ; the end of the image.
   qdiag = long(keyword_set(diagonal))
   nadd = call_external(soname, 'grow_obj', $
    nx, ny, image, mask, iloc, putval, qdiag, workarray)

   return, nadd
end
;------------------------------------------------------------------------------
function grow_object, image, xstart=xstart1, ystart=ystart1, putval=putval1, $
 diagonal=diagonal, nadd=nadd

   ndim = size(image, /n_dimen)
   dims = size(image, /dimens)
   nx = dims[0]
   if (ndim EQ 1) then ny = 1L $
    else ny = dims[1]
   nxcen = n_elements(xstart1)
   nycen = n_elements(ystart1)
   if (nxcen NE nycen) then $
    message, 'Number of elements in XSTART,YSTART must agree'

   ; Set default return values
   nadd = 0L

   ; Pad everything by 1, which was necessary to make the C code fast.
   image_pad = bytarr(nx+2,ny+2)
   image_pad[1:nx,1:ny] = image
   mask_pad = lonarr(nx+2,ny+2)
   workarray = lonarr(nx+2,ny+2)

   if (nxcen GE 1) then begin
      if (keyword_set(putval1)) then $
       objid = long( putval1[i<(n_elements(putval1)-1)] ) $
      else $
       objid = 1L
      for i=0L, nxcen-1 do begin
         iloc = long( (xstart1[i]+1) + (ystart1[i]+1) * (nx+2) )
         nadd1 = grow_obj1(image_pad, mask_pad, iloc, $
          putval=objid, diagonal=diagonal, $
          nx=nx+2, ny=ny+2, workarray=workarray)
         if (nadd1 GT 0) then begin
            nadd = nadd + nadd1
            if (NOT keyword_set(putval1)) then objid = objid + 1L
         endif
      endfor
   endif else begin
      indx = [where(image_pad NE 0 AND mask_pad EQ 0, ct), 0]
      jj = 0L
      if (keyword_set(putval1)) then objid = long(putval1[0]) $
       else objid = 1L
      while (jj LT ct) do begin
         nadd1 = grow_obj1(image_pad, mask_pad, indx[jj], $
          putval=objid, diagonal=diagonal, $
          nx=nx+2, ny=ny+2, workarray=workarray)
         if (nadd1 GT 0) then begin
            nadd = nadd + nadd1
            if (NOT keyword_set(putval1)) then objid = objid + 1L
         endif
         ; Pick the next pixel that can be assigned...
         while ((image_pad[indx[jj]] EQ 0 OR mask_pad[indx[jj]] NE 0) $
          AND jj LT ct) do jj = jj + 1
      endwhile
   endelse

   workarray = 0 ; Clear memory

   ; Un-pad the output image by 1
   mask = mask_pad[1:nx,1:ny]

   return, mask
end
;------------------------------------------------------------------------------
