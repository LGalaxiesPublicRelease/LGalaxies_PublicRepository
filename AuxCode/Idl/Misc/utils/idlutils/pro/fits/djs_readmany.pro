;-----------------------------------------------------------------------
;+
; NAME:
;   djs_readmany
;
; PURPOSE:
;   Read many FITS 2-D images into a data cube.
;
; CALLING SEQUENCE:
;   cube = djs_readmany( files, [hdr=] )
;
; INPUTS:
;   files:      FITS file names (array of strings)
;
; OUTPUTS:
;   cube:       Data cube with dimensions [NAXIS1, NAXIS2, nfile]
;
; OPTIONAL OUTPUTS:
;   hdr:        Header for first image
;
; COMMENTS:
;   Additional keywords are passed to the function READFITS().
;   At present, the output image is always typed FLOAT.
;
; PROCEDURES CALLED:
;   readfits()
;
; REVISION HISTORY:
;   07-Jul-1999  Written by David Schlegel, Princeton.
;-
;-----------------------------------------------------------------------
function djs_readmany, files, hdr=hdr, _EXTRA=KeywordsForReadfits

   ; Need at least 1 parameter
   if N_params() LT 1 then begin
      print, 'Syntax - cube = djs_readmany( files )'
      return, -1
   endif

   ; Read the first image
   image = readfits(files[0], hdr, _EXTRA=KeywordsForReadfits)
   if (size(image, /n_dimen) NE 2) then begin
      message, 'Image not 2-dimensional'
   endif

   ; Allocate the data cube, assuming all images are 2-D and the same size
   sz = size(image, /dimen)
   naxis1 = sz[0]
   naxis2 = sz[1]
   nfile = size(files, /n_elements)
   cube = make_array(naxis1, naxis2, nfile, /float)
   cube[*,*,0] = image

   ; Read additional images into the data cube
   for ifile=1, nfile-1 do begin
      image = readfits(files[ifile], _EXTRA=KeywordsForReadfits)
      szt = size(image, /dimen)
      if (size(image, /n_dimen) NE 2) then begin
         message, 'Image not 2-dimensional'
      endif else if (sz[0] NE szt[0] OR sz[1] NE szt[1]) then begin
         message, 'Images not same dimensions'
      endif
      cube[*,*,ifile] = image
   endfor

   return, cube
end 
;-----------------------------------------------------------------------
