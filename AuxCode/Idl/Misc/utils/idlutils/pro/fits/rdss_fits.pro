;+
; NAME:
;   rdss_fits
;
; PURPOSE:
;   Read a FITS file into IDL data and header variables
;
; CALLING SEQUENCE:
;   image = rdss_fits( filename, [ hdr, /nofloat, /silent ] )
;
; INPUTS:
;   filename   - Scalar string containing the name of the FITS file  
;                (including extension) to be read.   If the filename has
;                a *.gz extension, it will be treated as a gzip compressed
;                file.   If it has a .Z extension, it will be treated as a
;                Unix compressed file.
;
; OPTIONAL KEYWORDS:
;   nofloat    - If set, then keep data as unsigned integers.
;   silent     - suppress informational messages
;
; OUTPUTS:
;   image      - FITS data array constructed from designated record.
;                If the specified file was not found, then return -1.
;
; OPTIONAL OUTPUTS:
;   hdr        - String array containing the header from the FITS file.
;
; COMMENTS:
;   This routine will read a simple FITS image, or convert a non-standard
;   SDSS image that uses unsigned 16-bit integers.  One can pass any other
;   keywords expected by READFITS().
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   mrdfits()
;   sxaddpar
;   sxdelpar
;   sxpar()
;
; REVISION HISTORY:
;   13-May-1999  Written by David Schlegel, Princeton.
;   07-Jan-2001  Finkbeiner - added /silent because of U16 message
;-
;------------------------------------------------------------------------------
function rdss_fits, filename, hdr, nofloat=nofloat, _EXTRA=KeywordsForReadfits, silent=silent

   ; Need at least 1 parameter
   if (N_params() LT 1) then begin
      print, 'Syntax - image = rdss_fits( filename, [ hdr ] )'
      return, -1
   endif

   ;----------
   ; Read the image and header

   image = mrdfits(filename, 0, hdr, _EXTRA=KeywordsForReadfits, silent=silent)

   ;----------
   ; Remove extraneous or invalid FITS cards

   sxdelpar, hdr, 'SDSS'     ; This is an invalid FITS card, with no "=" sign
   sxdelpar, hdr, ' '        ; Remove blank header cards

   ;----------
   ; Test to see if the image is stored in the non-standard unsigned integer
   ; format of SDSS.  Convert to floating-point values unless /NOFLOAT
   ; is specified.

   qsimple = strpos( (str_sep(hdr[0],'/'))[0], 'F', 8) ; GE 0 if SIMPLE=F
   result = sxpar(hdr, 'UNSIGNED', count=ct1)
   if (qsimple GE 0 AND ct1 NE 0) then begin

      bitpix = sxpar(hdr, 'BITPIX')

      ; Convert from unsigned 16-bit integers to floats
      if (NOT keyword_set(nofloat)) then begin
         if NOT keyword_set(silent) then $
           print, 'Converting from U16...'
         sxaddpar, hdr, 'SIMPLE', 'T', 'FITS STANDARD'
         sxdelpar, hdr, 'UNSIGNED' ; Remove UNSIGNED keyword
         ineg = where(image LT 0)
         image = temporary(image) + 0.0
         if (ineg[0] NE -1) then image[ineg] = 2.0^bitpix + image[ineg]
      endif else begin
         image = temporary( uint(image) )
      endelse

   endif

   return, image
end
