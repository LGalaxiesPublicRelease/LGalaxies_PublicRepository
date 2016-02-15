;------------------------------------------------------------------------------
;+
; NAME:
;   wcs_getval
;
; PURPOSE:
;   Read value(s) from NGP+SGP polar projections.
;
; CALLING SEQUENCE:
;   value = wcs_getval(files, lonvec, latvec, [ path=path, $
;    /interp, /noloop, /verbose ] )
;
; INPUTS:
;   files:      File name(s); if Lambert or ZEA projection, one can pass
;               two file names where the first is used for northern points
;               and the second for southern points
;   lonvec:     Longitude(s) [degrees]
;   latvec:     Latitude(s) [degrees]
;
; OPTIONAL KEYWORDS:
;   path:       File name path; default to ''
;   interp:     Set this flag to return a linearly interpolated value
;               from the 4 nearest pixels
;   noloop:     Set this flag to read all values at once without a FOR loop.
;               This is a faster option for reading a large number of values,
;               but requires reading an entire FITS image into memory.
;               (Actually, the smallest possible sub-image is read.)
;   verbose:    Set this flag for verbose output, printing pixel coordinates
;               and map values.  Setting NOLOOP disables this option.
;
; OUTPUTS:
;   value:      Value(s) from maps
;
; PROCEDURES CALLED:
;   fxread
;   headfits()
;   sxpar()
;   wcs_coord2pix
;
; REVISION HISTORY:
;   31-Mar-1999  Modified from lambert_getval() - D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
function wcs_getval, files, lonvec, latvec, path=path, $
 interp=interp, noloop=noloop, verbose=verbose

   ; Need 3 parameters
   if N_params() LT 3 then begin
      print, 'Syntax - value = wcs_getval(files, lonvec, latvec, $'
      print, ' [ path=path, /noloop, /interp ] )'
      return, -1
   endif

   if (NOT keyword_set(path)) then path = ''

   ; Allocate output data array as floating-point values
   if (N_elements(lonvec) EQ 1) then value = 0.0 $
    else value = float(0*lonvec)

   ; Loop through each file
   Nfile = N_elements(files)
   for ifile=0, Nfile-1 do begin

      fullname = path + files[ifile]

      ; If more than 1 file, assume the 1st file is for northern points,
      ; and the 2nd file is for southern
      if (Nfile EQ 1) then begin
         indx = lindgen(latvec)
      endif else begin
         if (ifile EQ 0) then indx = where(latvec GE 0) $
          else indx = where(latvec LT 0)
      endelse

      if (indx[0] NE -1) then begin

         ; Read the header only
         hdr = headfits(fullname)

         naxis1 = sxpar(hdr, 'NAXIS1')
         naxis2 = sxpar(hdr, 'NAXIS2')

         if (NOT keyword_set(interp)) then begin ; NEAREST PIXELS

            ; Determine the nearest pixel coordinates
            wcs_coord2pix, lonvec[indx], latvec[indx], hdr, xpix, ypix

;            xmax = max(xpix, min=xmin)
;            ymax = max(ypix, min=ymin)
;            print, 'xmin, xmax, ymin, ymax', xmin, xmax, ymin, ymax

            ; Force pixel locations to fall within the image bounds
            xpix = (temporary(xpix) > 0) < (naxis1-1)
            ypix = (temporary(ypix) > 0) < (naxis2-1)
            xmax = max(xpix, min=xmin)
            ymax = max(ypix, min=ymin)

            if (keyword_set(noloop)) then begin ; READ FULL IMAGE

               fxread, fullname, subimg, hdr, $
                xmin, xmax, ymin, ymax
               value[indx] = subimg[temporary(xpix)-xmin,temporary(ypix)-ymin]

            endif else begin ; READ ONE VALUE AT A TIME

               ; Read one pixel value at a time from data file
               for ii=0L, N_elements(indx)-1 do begin
                  fxread, fullname, onedat, hdr, $
                   xpix[ii], xpix[ii], ypix[ii], ypix[ii]
                  value[indx[ii]] = onedat
                  if (keyword_set(verbose)) then $
                   print, format='(f8.3,f8.3,i2,i9,i9,e13.5)', $
                    lonvec[indx], latvec[indx], $
                    ifile, xpix[ii], ypix[ii], value[indx[ii]]
               endfor

            endelse

         endif else begin ; INTERPOLATE

           ; Determine the pixel coordinates for this projection

            wcs_coord2pix, lonvec[indx], latvec[indx], hdr, xr, yr, $
             /fractional
             xpix1 = fix(xr)
             ypix1 = fix(yr)
             dx = xpix1 - float(xr) + 1.0
             dy = ypix1 - float(yr) + 1.0

            ; Force pixel values to fall within the image boundaries.
            ; Any pixels outside the image are changed to the boundary pixels.
            ibad = where(xpix1 LT 0)
            if (ibad[0] NE -1) then begin
               xpix1[ibad] = 0
               dx[ibad] = 1.0
            endif
            ibad = where(ypix1 LT 0)
            if (ibad[0] NE -1) then begin
               ypix1[ibad] = 0
               dy[ibad] = 1.0
            endif
            ibad = where(xpix1 GE naxis1-1)
            if (ibad[0] NE -1) then begin
               xpix1[ibad] = naxis1-2
               dx[ibad] = 0.0
            endif
            ibad = where(ypix1 GE naxis2-1)
            if (ibad[0] NE -1) then begin
               ypix1[ibad] = naxis2-2
               dy[ibad] = 0.0
            endif

            ; Create Nx4 arry of interpolation weights
            weight = [ [    dx  *    dy  ] , $
                       [ (1-dx) *    dy  ] , $
                       [    dx  * (1-dy) ] , $
                       [ (1-dx) * (1-dy) ] ]

; Clear memory
            dx = 0
            dy = 0

            if (keyword_set(noloop)) then begin ; READ FULL IMAGE

               xmax = max(xpix1, min=xmin)+1.
               ymax = max(ypix1, min=ymin)+1.
               fxread, fullname, subimg, hdr, $
                xmin, xmax, ymin, ymax
               value[indx] = $
                subimg[xpix1-xmin    ,ypix1-ymin    ] * weight[*,0] + $
                subimg[xpix1-(xmin-1),ypix1-ymin    ] * weight[*,1] + $
                subimg[xpix1-xmin    ,ypix1-(ymin-1)] * weight[*,2] + $
                subimg[xpix1-(xmin-1),ypix1-(ymin-1)] * weight[*,3]
; Clear memory
               weight = 0 & subimg = 0
               xpix1  = 0 & ypix1  = 0
            endif else begin ; READ ONE VALUE AT A TIME

               ; Read 2x2 array at a time from data file
               for ii=0L, N_elements(indx)-1 do begin
;                  weight = [ [    dx[ii]  *    dy[ii]  ,  $
;                               (1-dx[ii]) *    dy[ii]  ], $
;                             [    dx[ii]  * (1-dy[ii])  , $
;                               (1-dx[ii]) * (1-dy[ii]) ]  ]
                  fxread, fullname, onedat, hdr, $
                   xpix1[ii], xpix1[ii]+1, ypix1[ii], ypix1[ii]+1
                  value[indx[ii]] = total(onedat[*] * weight[ii,*])

                  if (keyword_set(verbose)) then $
                   print, format='(f8.3,f8.3,i2,f9.3,f9.3,e13.5)', $
                    lonvec[indx], latvec[indx], $
                    ifile, xr[ii], yr[ii], value[indx[ii]]
               endfor

            endelse

         endelse

      endif

   endfor

   return, value
end
;------------------------------------------------------------------------------
