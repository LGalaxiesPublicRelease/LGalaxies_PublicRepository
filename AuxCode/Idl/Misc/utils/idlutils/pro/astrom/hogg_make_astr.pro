;+
; NAME:
;   hogg_make_astr
; PURPOSE:
;   Generate the astrometric header for a particular pointing and
;   orientation.
; COMMENTS:
;   Adds NAXIS to the astrom structure.
; CALLING SEQUENCE:
;   astr= hogg_make_astr(racen,deccen,dra,ddec $
;                        [,pixscale=pixscale] [etc])
; INPUTS:
;   racen    - Central RA [deg]
;   deccen   - Central DEC [deg]
;   dra      - Size in the X dimension; default to 0.5 deg
;   ddec     - Size in the Y dimension; default to 0.5 deg
;
; OPTIONAL INPUTS:
;   pixscale    - Pixel scale (deg); default 1.0/3600
;   orientation - angle (deg) for the north vector to make relative to
;                 straight up (y direction), CCW being positive
;   npixround   - round array dimensions (naxis1 and naxis2) to nearest
;                 factor of npixround; default to 8.
; KEYWORDS:
;   orthographic  - Make orthographic "-SIN" header instead of
;                   gnomonic "-TAN" header.
; OUTPUTS:
;   astr     - Astrometry structure with NAXIS keyword added
; OPTIONAL OUTPUTS:
;   pixscale  - Set, if not input
; EXAMPLES:
; BUGS:
; REVISION HISTORY:
;   2003-Nov-21   Written by D. Schlegel (Princeton) & D. Hogg (NYU)
;   2003-Dec-02   Modified to produce simpler headers - Hogg
;   2004-Apr-18   Use make_astr to define astrom structure - DPF
;   2004-Jun-17   Modified to make "orthographic" headers on request - Hogg
;   2005-Aug-31   Changed name and put into idlutils - Hogg
;------------------------------------------------------------------------------
function hogg_make_astr, racen,deccen,dra1,ddec1,pixscale=pixscale, $
                         orientation=orientation,orthographic=orthographic, $
                         npixround=npixround
   if (n_params() LT 2) then begin
      print, 'Must specify RACEN, DECCEN'
      return, 0
   endif
   if (keyword_set(dra1)) then dra = dra1 else dra = 0.5
   if (keyword_set(ddec1)) then ddec = ddec1 else ddec = 0.5
   if (NOT keyword_set(pixscale)) then pixscale = 1d0 / 3600d0
   if (NOT keyword_set(npixround)) then npixround = 8
   naxis1= round(dra/pixscale/float(npixround))*npixround
   naxis2= round(ddec/pixscale/float(npixround))*npixround
   if keyword_set(orthographic) then ctype= ['RA---SIN','DEC--SIN'] else $
     ctype= ['RA---TAN','DEC--TAN']
   if (NOT keyword_set(orientation)) then orientation=0D0
   theta= orientation*!DPI/180D0
   ct= cos(theta)
   st= sin(theta)
   make_astr, bigast1, $
     CD       = double([[-pixscale*ct,-pixscale*st], $
                        [-pixscale*st, pixscale*ct]]), $
     DELT     = double([1.0,1.0]), $
     CRPIX    = double([0.5,0.5]+0.5*[naxis1,naxis2]), $  ; NB: FITS CONVENTION
     CRVAL    = double([racen,deccen]), $
     CTYPE    = ctype, $
     LONGPOLE = 1.8D2
   bigast= struct_addtags(bigast1,{naxis: [naxis1,naxis2]})
   return, bigast
end
;------------------------------------------------------------------------------
