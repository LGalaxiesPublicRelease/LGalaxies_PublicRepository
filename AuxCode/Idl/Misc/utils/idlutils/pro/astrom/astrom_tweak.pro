;+
; NAME:
;   astrom_tweak
;
; PURPOSE:
;   Tweak astrometric solution, given a good initial guess
;
; CALLING SEQUENCE:
;   gsa_out = astrom_tweak( gsa_in, catra, catdec, imx, imy, $
;    [ dtheta=, errflag=, nmatch=, catind=, obsind=, /verbose ] )
;
; INPUTS:
;   gsa_in     - Initial guess for astrometric solution (struct)
;   cat        - Structure (with fields .ra, .dec) of catalogue positions
;   im         - Structure (with fields .x, .y) of image star positions
;
; OPTIONAL KEYWORDS:
;   dtheta     - Match distance between catalog and image stars;
;                default to 5 arcsec
;   verbose    - If set, then print sizes of offsets
;  
; OUTPUTS:
;   gsa_out    - returned guess for astrometric solution (struct);
;                0 if solution failed
;
; OUTPUT OUTPUTS:
;   errflag    - Set to 1 if fatal error occurs, 0 otherwise
;   nmatch     - Number of matched stars
;   catind     - Indices of CATLON,CATLAT for matched objects.
;   obsind     - Indices of XPOS,YPOS for matched objects.
;
; COMMENTS:
;   Uses preliminary solution given in astr structure to match image
;   and catalogue stars within maxsep pixels of each other.  These
;   are then used by astrom_warp to determine a new solution, returned
;   in astr.
; 
; BUGS:
;   The terms AMDX[4],AMDY[4] in the GSSS structure should actually *not*
;    be fit for since they are higher-order, but this was the easiest
;    way to implement this code???
;
; PROCEDURES CALLED:
;   djs_angle_match()
;   gsssadxy
;   gsssxyad
;
; REVISION HISTORY:
;   02-Feb-2003  Written by D. Schlegel and D. Hogg, APO
;-
;-----------------------------------------------------------------------------
function astrom_tweak, gsa_in, catra, catdec, imx, imy, dtheta=dtheta, $
 errflag=errflag, nmatch=nmatch, catind=catind, obsind=obsind, $
 verbose=verbose

   errflag = 0
   if (NOT keyword_set(dtheta)) then dtheta = 5./3600.

   ;----------
   ; Find matches between catalog (RA,DEC) and image X,Y

   gsssxyad, gsa_in, imx, imy, imra, imdec
   nmatch = djs_angle_match(imra, imdec, catra, catdec, dtheta=dtheta, $
    mcount=mcount, mindx=mindx)
   if (nmatch LT 5) then begin
      splog, 'Too few matches; returning original GSA'
      errflag = 1
      catind = -1L
      obsind = -1L
      return, gsa_in
   endif
   obsind = where(mcount GT 0)
   catind = mindx[obsind]

   ;----------
   ; Compute GSA internal coordinate system chi,eta,obx,oby
   ; This code is replicated from GSSSADXY

   radeg = 180.0d/!DPI
   arcsec_per_radian= 3600.0d*radeg
   dec_rad = catdec[catind] / radeg
   ra_rad = catra[catind] / radeg
   pltra = gsa_in.crval[0] / radeg
   pltdec = gsa_in.crval[1] / radeg

   cosd = cos(dec_rad)
   sind = sin(dec_rad)
   ra_dif = ra_rad - pltra

   div = ( sind*sin(pltdec) + cosd*cos(pltdec) * cos(ra_dif))
   xi = cosd*sin(ra_dif) * arcsec_per_radian / div
   eta = (sind*cos(pltdec) - cosd*sin(pltdec) * cos(ra_dif)) $
    * arcsec_per_radian / div

   obx = ( gsa_in.ppo3 - (gsa_in.xll + imx[obsind] + 0.5d) * gsa_in.xsz) / 1000.d
   oby = (-gsa_in.ppo6 + (gsa_in.yll + imy[obsind] + 0.5d) * gsa_in.ysz) / 1000.d

   ;----------
   ; Compute bilinear transformation between obx,oby <-> imx,imy
   ; The terms kx[1,1] and ky[1,1] should actually *not* be fit for since they
   ; are higher-order, but the POLYWARP procedure does not give you that option.
   ; This should probably be fixed!???

   polywarp, xi, eta, obx, oby, 1, kx, ky

   gsa_out = gsa_in
   gsa_out.amdx[2] = kx[0,0]
   gsa_out.amdy[2] = ky[0,0]
   gsa_out.amdx[0] = kx[0,1]
   gsa_out.amdy[0] = ky[1,0]
   gsa_out.amdx[1] = kx[1,0]
   gsa_out.amdy[1] = ky[0,1]
   gsa_out.amdx[4] = kx[1,1]
   gsa_out.amdy[4] = ky[1,1]

   ;----------
   ; Report the RMS differences between catalog and CCD positions

   if (keyword_set(verbose)) then begin
      gsssadxy, gsa_in, catra[catind], catdec[catind], catx, caty
      xdiff = imx[obsind] - catx
      ydiff = imy[obsind] - caty
      splog, 'Input mean/stdev offset in X = ', mean(xdiff), stddev(xdiff)
      splog, 'Input mean/stdev offset in Y = ', mean(ydiff), stddev(ydiff)

      gsssadxy, gsa_out, catra[catind], catdec[catind], catx, caty
      xdiff = imx[obsind] - catx
      ydiff = imy[obsind] - caty
      splog, 'Output mean/stdev offset in X = ', mean(xdiff), stddev(xdiff)
      splog, 'Output mean/stdev offset in Y = ', mean(ydiff), stddev(ydiff)
   endif

;gsssadxy,gsa_in,catra[catind],catdec[catind],catx,caty
;splot, imx,imy,ps=4
;soplot,catx,caty,ps=5,color='red'
;gsssadxy,gsa_out,catra[catind],catdec[catind],outx,outy
;soplot,outx,outy,ps=4,color='green'

   return, gsa_out
end 
;-----------------------------------------------------------------------------
