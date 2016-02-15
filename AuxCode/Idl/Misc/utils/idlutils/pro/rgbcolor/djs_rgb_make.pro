;+
;NAME:
;  djs_rgb_make
;
; PURPOSE:
;   Creates JPEG (or TIFF) from three images or FITS files
;
; CALLING SEQUENCE:
;   djs_rgb_make, rimage, gimage, bimage, [ name=, origin=, scales=, $
;    nonlinearity=, satvalue=, rebinfactor=, overlay=, quality=, $
;    /tiff, dpitiff=, bits_per_channel ]
;
; INPUTS:
;   rimage,gimage,bimage - Input 2-dimensional images or names of FITS files;
;                 the dimensions of all images must agree
;
; OPTIONAL KEYWORDS:
;   name        - Name of the output JPEG file; default to 'test.jpg',
;                 or to 'test.tiff' if /TIFF is set
;   origin      - Subtract these zero-point values from the input images
;                 before any other scalings; default to [0,0,0]
;   scales      - Multiplicative scaling for each image; default to [1,1,1]
;   nonlinearity- Nonlinearity constant, =0 for linear, =Inf for logarithmic;
;                 default to 3
;   satvalue    - Saturation value (before applying rescaling with SCALES, but
;                 after applying ORIGIN) for which we should group pixels
;                 into saturated "objects"; default to 100;
;                 set to zero to disable
;   rebinfactor - Integer by which to rebin pixels in the x and y
;                 directions; eg, a rebinfactor of 2 halves the number
;                 of pixels in each direction and quarters the total
;                 number of pixels in the image
;   overlay     - Optional overlay image, which must be dimensionsed as
;                 [NX/REBINFACTOR,NY/REBINFACTOR,3]
;   quality     - Quality for WRITE_JPEG; default to 75 per cent
;   tiff        - Set to make TIFF instead of JPEG if either this keyword
;                 or DPITIFF is set
;   dpitiff     - Set TIFF "dots per inch" resolution, and force /TIFF option
;
; OUTPUTS:
;
; COMMENTS:
;   This routine is based upon Nick Wherry's code NW_RGB_MAKE.
;   The main difference is that saturated pixels are grouped into
;   contiguous "objects", which are then assigned a color based upon
;   the sum of all the pixels in that object.
;
;   The nonlinearity function applied is
;     RIMAGE = RIMAGE * asinh(b*r)/(b*r)
;     GIMAGE = GIMAGE * asinh(b*r)/(b*r)
;     BIMAGE = BIMAGE * asinh(b*r)/(b*r)
;   where "b" is the input NONLINEARITY parameter and we define at each pixel
;     r = (RIMAGE + GIMAGE + BIMAGE)
;
;   Note that there are two types of saturation.  The first is that objects
;   can be considered saturated if they exceed SATVALUE in any of the input
;   images.  For such objects, contiguous saturated pixels are combined into
;   one "object" with the mean color of all included pixels.
;   The second type of saturation is of the RGB image.  This saturation is
;   dealt with at the pixel level, and each pixel rescaled in all three images
;   such that the brightest color hits the JPEG limit (of 255), but the
;   colors (ratios between the RGB images) are preserved.
;
; EXAMPLES:
;
; BUGS:
;
; REVISION HISTORY:
;   10-May-2004 - Written by D. Schlegel, Princeton;
;                 based upon Nick Wherry's code NW_RGB_MAKE
;   12-Dec-2005 - Change saturated pixel code; add TIFF option - DPF
;   28-Apr-2006 - bits_per_channel keyword added for Warner Bros - DPF
;-
;------------------------------------------------------------------------------
pro djs_rgb_make, rimage, gimage, bimage, name=name, $
 origin=origin1, scales=scales1, nonlinearity=nonlinearity1, $
 satvalue=satvalue1, rebinfactor=rebinfactor1, overlay=overlay, $
 quality=quality, tiff=tiff1, dpitiff=dpitiff, $
 bits_per_channel=bits_per_channel

   t0 = systime(1)
   thismem = float(ulong(memory()))

   ;----------
   ; Set defaults

   tiff = keyword_set(tiff1)
   if (keyword_set(dpitiff)) then tiff = 1B ; assume tiff if dpitiff set
   if (NOT keyword_set(bits_per_channel)) then bits_per_channel = 8
   if (bits_per_channel NE 8 AND keyword_set(tiff) EQ 0) then $
    message, 'BITS_PER_CHANNEL must be 8 for JPEGs'

   suffix = keyword_set(tiff) ? '.tif' : '.jpg'
   if (NOT keyword_set(name)) THEN name = 'test'
; -------- force name suffix to agree with file type
   if (strpos(name, '.jpg') EQ -1) and (strpos(name, '.tif') EQ -1) then begin 
      name = name+suffix
   endif else begin 
      extensions = ['.tiff', '.TIFF', '.tif', '.TIF', '.jpeg', '.jpg', '.JPEG', '.JPG']
      for iext=0L, n_elements(extensions)-1 do $
        name = repstr(name, extensions[iext], suffix)
   endelse
   
   if (n_elements(origin1) EQ 0) THEN origin = [0,0,0] $
    else origin = float(origin1)
   if (n_elements(scales1) EQ 0) THEN scales = [1,1,1] $
    else scales = float(scales1)
   if (NOT keyword_set(rebinfactor1)) THEN rebinfactor = 1 $
    else rebinfactor = float(rebinfactor1)
   if (NOT keyword_set(quality)) then quality = 75
   if (n_elements(nonlinearity1) EQ 0) then nonlinearity = 3. $
    else nonlinearity = float(nonlinearity1[0])
   if (n_elements(satvalue1) EQ 0) then satvalue = 100. $
    else satvalue = float(satvalue1[0])

   if (n_elements(origin) NE 3) then $
    message, 'ORIGIN must have 3 elements'
   if (n_elements(scales) NE 3) then $
    message, 'SCALES must have 3 elements'

   ;----------
   ; Read the 3 images, and sanity-check that they are the same dimensions

   if (size(rimage,/tname) EQ 'STRING') then begin
      rimg = mrdfits(rimage+'*', /silent)
      gimg = mrdfits(gimage+'*', /silent)
      bimg = mrdfits(bimage+'*', /silent)
   endif else begin
      rimg = float(rimage)
      gimg = float(gimage)
      bimg = float(bimage)
   endelse

   dims = size(rimg, /dimens)
   if (size(rimg, /n_dimen) NE 2) then begin
      print, 'Images must be 2-dimensional arrays!'
      return
   endif
   if (total(size(gimg, /dimens) NE dims) NE 0 $
    OR total(size(bimg, /dimens) NE dims) NE 0) then begin
      print, 'Dimensions of all 3 images must agree!'
      return
   endif

   ;----------
   ; Optionally rebin the images

   if (rebinfactor NE 1) then begin
      dims = round(dims / rebinfactor)
      if (rebinfactor EQ round(rebinfactor)) then begin
         rimg = rebin(rimg, dims[0], dims[1], /sample)
         gimg = rebin(gimg, dims[0], dims[1], /sample)
         bimg = rebin(bimg, dims[0], dims[1], /sample)
      endif else begin
         rimg = congrid(rimg, dims[0], dims[1], /interp)
         gimg = congrid(gimg, dims[0], dims[1], /interp)
         bimg = congrid(bimg, dims[0], dims[1], /interp)
      endelse
   endif

   ;----------
   ; Apply optional zero-point offsets or rescalings

   if (scales[0] NE 1 OR origin[0] NE 0) then $
    rimg = (scales[0] * (rimg - origin[0])) > 0 $
   else $
    rimg = rimg > 0
   if (scales[1] NE 1 OR origin[1] NE 0) then $
    gimg = (scales[1] * (gimg - origin[1])) > 0 $
   else $
    gimg = gimg > 0
   if (scales[2] NE 1 OR origin[2] NE 0) then $
    bimg = (scales[2] * (bimg - origin[2])) > 0 $
   else $
    bimg = bimg > 0

   ;----------
   ; Loop through each saturated object, and replace all pixels in each object
   ; with the mean color of that object

   if (keyword_set(satvalue)) then begin
      ; Determine where in the image we are saturating any of the 3 colors
      satmask = (rimg GT satvalue*scales[0]) $
       OR (gimg GT satvalue*scales[1]) $
       OR (bimg GT satvalue*scales[2])
      isat = where(satmask, nsat)

      ; This function groups all contiguous saturated pixels into one object.
      ; (This uses quite a bit of memory, though)
      if (nsat GT 0) then $
       objmask = grow_object(satmask)
   endif else begin
      nsat = 0
   endelse

   if (nsat GT 0) then begin
      nobj = max(objmask[isat])

      ; Produce the list of all saturated object pixels, sorted by object ID
      if (nobj EQ 1) then begin
         i1 = 0L
         i2 = nsat - 1
      endif else begin
         sortmask = objmask[isat]
         isort = sort(sortmask)
         isat = isat[isort]
         sortmask = sortmask[isort]
         isort = 0 ; clear memory
         i1 = where(sortmask NE shift(sortmask,1))
         sortmask = 0 ; clear memory
         i2 = [i1[1:nobj-1]-1, nsat-1]
      endelse
      for iobj=0L, nobj-1 do begin
         indx = isat[i1[iobj]:i2[iobj]]

; -------- DPF - get neighbors and compute fill value from them
         ix = indx mod dims[0]
         iy = indx  /  dims[0]
         
         ix0 = (ix-1) > 0
         ix1 = (ix+1) < (dims[0]-1)
         iy0 = (iy-1) > 0
         iy1 = (iy+1) < (dims[1]-1)
         jx = [ix, ix, ix0, ix1, ix] ; non-diagonal neighbors
         jy = [iy0, iy1, iy, iy, iy]
         indj = jx+jy*dims[0]
         indregion = indj[uniq(indj, sort(indj))]
         
         n_neighbor = n_elements(indregion)-n_elements(indx) 

         rimg[indx] = (total(rimg[indregion])-total(rimg[indx]))/n_neighbor
         gimg[indx] = (total(gimg[indregion])-total(gimg[indx]))/n_neighbor
         bimg[indx] = (total(bimg[indregion])-total(bimg[indx]))/n_neighbor
; -------- DPF - end change
      endfor
   endif

   objmask = 0 ; clear memory
   isat = 0
   i1 = 0
   i2 = 0

   ;----------
   ; Compute the nonlinear mapping

   radius = rimg + gimg + bimg
   radius = nonlinearity * radius
   radius = radius + (radius LE 0)
   nonlinfac = asinh(radius) / radius
   radius = 0 ; clear memory

   ;----------
   ; Determine where in the image we are saturating any of the 3 colors

   maxval = rimg > gimg > bimg
   satmask = (rimg * nonlinfac GT 1) $
    OR (gimg * nonlinfac GT 1) $
    OR (bimg * nonlinfac GT 1)
   isat = where(satmask, nsat)
   satmask = 0 ; clear memory
   if (nsat GT 0) then nonlinfac[isat] = 1. / maxval[isat]
   maxval = 0 ; clear memory
   isat = 0 ; clear memory

   ;----------
   ; Apply the nonlinearity corrections

   rimg = rimg * nonlinfac
   gimg = gimg * nonlinfac
   bimg = bimg * nonlinfac
   nonlinfac = 0 ; clear memory

   ;----------
   ; Optionally add the overlay images

   if (keyword_set(overlay)) then begin
      if (total(size(overlay, /dimens) NE [dims[0],dims[1],3]) NE 0) then begin
         splog, 'Ignoring overlay since dimensions disagree!'
      endif else begin
         rimg = rimg + overlay[*,*,0]
         gimg = gimg + overlay[*,*,1]
         bimg = bimg + overlay[*,*,2]
      endelse
   endif

   ;----------
   ; Convert from a floating-point to byte-scaled image

   case bits_per_channel of 
      8: begin ; Good enough for mortals.
         byteimg = bytarr(dims[0], dims[1], 3)
         byteimg[*,*,0] = byte((floor(temporary(rimg) * 256) > 0) < 255)
         byteimg[*,*,1] = byte((floor(temporary(gimg) * 256) > 0) < 255)
         byteimg[*,*,2] = byte((floor(temporary(bimg) * 256) > 0) < 255)
      endcase
      16: begin ; Hollywood likes this.
         byteimg = uintarr(dims[0], dims[1], 3)
         byteimg[*,*,0] = uint((floor(temporary(rimg) * 65536.0) > 0) < 65535)
         byteimg[*,*,1] = uint((floor(temporary(gimg) * 65536.0) > 0) < 65535)
         byteimg[*,*,2] = uint((floor(temporary(bimg) * 65536.0) > 0) < 65535)
         qshort = 1B
      endcase
      32: begin ; You would have to be insane...
         byteimg = ulonarr(dims[0], dims[1], 3)
         byteimg[*,*,0] = ulong((floor(temporary(rimg) * 4294967296.0) > 0) < 4294967295)
         byteimg[*,*,1] = ulong((floor(temporary(gimg) * 4294967296.0) > 0) < 4294967295)
         byteimg[*,*,2] = ulong((floor(temporary(bimg) * 4294967296.0) > 0) < 4294967295)
         qlong = 1B
      endcase
      else: message, 'value of bits_per_channel not supported!'
   endcase         
   ;----------
   ; Generate the JPEG (or TIFF) image

   IF keyword_set(tiff) THEN BEGIN
      colors = reverse(temporary(byteimg),2)
      splog, 'writing tiff file: ', name
      write_tiff, name, planarconfig=2, red=colors[*,*,0], $
        green=colors[*,*,1], blue=colors[*,*,2], xresol=dpitiff, yresol=dpitiff, $
        short=qshort, long=qlong
   ENDIF ELSE BEGIN
      splog, 'writing jpeg file: ', name
      write_jpeg, name, byteimg, true=3, quality=quality
   ENDELSE


   thismem = float(ulong(memory()))
   splog, 'Max memory usage = ', thismem[3]/1.d6, ' MB'
   splog, 'Elapsed time = ', systime(1)-t0, ' sec'

   return
end
;------------------------------------------------------------------------------
