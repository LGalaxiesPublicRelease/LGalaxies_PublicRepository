;------------------------------------------------------------------------------
;+
; NAME:
;   first_coverage
;
; PURPOSE:
;   Read the FIRST survey coverage maps.
;
; CALLING SEQUENCE:
;   skyrms = first_coverage(ra, dec)
;
; OPTIONAL INPUTS:
;   ra:          Right ascensions [J2000 deg]
;   dec:         Declinations [J2000 deg]
;
; OUTPUTS:
;   skyrms:      Returned RMS noise in mJy/beam
;
; COMMENTS:
;   These images give the RMS in mJy/beam tabulated on a ~3 arcmin grid
;   in RA and DEC.  However, the WCS headers in these FITS files are invalid.
;
; BUGS:
;   The coordinates in these FITS headers have been interpreted to be
;   for the *center* of each pixel, though there is no documentation as
;   to whether this is the correct interpretation.
;
; DATA FILES:
;   The following data files can be copied from:
;     http://sundog.stsci.edu/first/catalogs/
;   and should be put in a directory pointed to by the environment
;   variable $FIRST_DIR:
;     $FIRST_DIR/coverage-north-3arcmin-03apr11.fits
;     $FIRST_DIR/coverage-south-3arcmin-01oct15.fits
;
; PROCEDURES CALLED:
;   headfits()
;   mrdfits()
;
; REVISION HISTORY:
;   Written D. Schlegel, 18 July 2003, Princeton
;    31 July 2003 - /silent keyword added to read - DPF
;-
;------------------------------------------------------------------------------
function first_coverage, ravec1, decvec

   common com_first_coverage, ramin, ramax, decmin, decmax, dra, ddec

   npts = n_elements(ravec1)
   if (npts EQ 0) then return, 0
   ravec = ravec1
   cirrange, ravec ; Change to the range [0,360] degrees

   first_dir = getenv('FIRST_DIR')
   if (NOT keyword_set(first_dir)) then begin
      splog, 'WARNING: FIRST_DIR environment variable not set.'
      splog, 'This must point to a data directory with files downloade from:'
      splog, '  http://sundog.stsci.edu/first/catalogs/'
      return, 0
   endif

   nfile = 2
   filename = filepath([ $
    'coverage-north-3arcmin-03apr11.fits', $
    'coverage-south-3arcmin-01oct15.fits'], root_dir=first_dir)

   ;----------
   ; Start by reading (and caching) coordinates from the two file headers
   ; if not done so already.

   if (NOT keyword_set(rabin)) then begin
      ramax = fltarr(nfile)
      ramin = fltarr(nfile)
      ramax = fltarr(nfile)
      decmin = fltarr(nfile)
      decmax = fltarr(nfile)
      dra = fltarr(nfile)
      ddec = fltarr(nfile)
      for ifile=0L, nfile-1 do begin
         hdr = headfits(filename[ifile])
         if (NOT keyword_set(hdr)) then begin
            print, 'ERROR: File not found: ' + filename[ifile]
            ramin = 0
         endif

         dra[ifile] = sxpar(hdr,'CDELT1')
         ramin[ifile] = sxpar(hdr,'CRVAL1') - 0.5 * dra[ifile]
         ramax[ifile] = ramin[ifile] + dra[ifile] * (sxpar(hdr,'NAXIS1')+1)

         ddec[ifile] = sxpar(hdr,'CDELT2')
         decmin[ifile] = sxpar(hdr,'CRVAL2') - 0.5 * ddec[ifile]
         decmax[ifile] = decmin[ifile] + ddec[ifile] * (sxpar(hdr,'NAXIS2')+1)

      endfor
   endif

   ;----------
   ; Loop over each file, and read the relevant portion of each image
   ; for the requested positions.

   skyrms = fltarr(npts)
   for ifile=0L, nfile-1 do begin
      ; Deal with the case for the southern file where RAMIN=-48 deg.
      indx = where(((ravec GE ramin[ifile] AND ravec LT ramax[ifile]) $
       OR (ravec GE ramin[ifile]+360 AND ravec LT ramax[ifile]+360)) $
       AND decvec GE decmin[ifile] AND decvec LT decmax[ifile], ct)
      if (ct GT 0) then begin
         xbin = floor( (ravec[indx] - 360*(ravec[indx] GT ramax[ifile]) $
          - ramin[ifile]) / dra[ifile] )
         ybin = floor( (decvec[indx] - decmin[ifile]) / ddec[ifile] )
         ymin = min(ybin, max=ymax)
         img = mrdfits(filename[ifile], range=[ymin,ymax], /silent)
         skyrms[indx] = img[xbin,ybin-ymin]
      endif
   endfor

   return, skyrms
end
;------------------------------------------------------------------------------
