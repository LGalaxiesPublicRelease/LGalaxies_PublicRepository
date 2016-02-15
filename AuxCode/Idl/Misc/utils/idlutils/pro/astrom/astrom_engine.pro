;+
; NAME:
;   astrom_engine
;
; PURPOSE:
;   Compute astrometric solution for a list of stars & catalogue stars
;
; CALLING SEQUENCE:
;   gsa_out = astrom_engine( xpos, ypos, catlon, catlat, gsa_in, $
;    [ search_rad=, search_scale=, search_angle=, $
;    poserr=, nmatch=, catind=, obsind=, /radial, /verbose ] )
;
; INPUTS:
;   xpos       - X positions in CCD coordinates
;   ypos       - Y positions in CCD coordinates
;   catlon     - Catalog star longitudes in the same coordinate system as GSA_IN
;   catlat     - Catalog star latitutes in the same coordinate system as GSA_IN
;   gsa_in     - Input GSSS structure with initial guess for astrometric
;                solution
;   radial     - (Not used.)
;
; OPTIONAL INPUTS:
;   search_rad   - Unused ???
;   search_scale - Unused ???
;   search_angle - If set, then search for rotations offset by +/-SEARCH_ANGLE
;                  relative to the input astrometric guess.
;   poserr       - Maximum position error in CCD coordinates; default to 1.
;                  No stars will be matched at distances further than this.
;   verbose      - If set, then be verbose.
;
; OUTPUTS:
;   gsa_out    - Output GSSS structure with astrometric solution;
;                return 0 if astrometry failed
;
; OPTIONAL OUTPUTS:
;   nmatch     - Number of matched objects with the input catalog.
;   catind     - Indices of CATLON,CATLAT for matched objects.
;   obsind     - Indices of XPOS,YPOS for matched objects.
;
; COMMENTS:
;   We assume that we know the scale and rotation well enough, then solve
;   for the X,Y offsets by correlating with catalog stars.
;
; BUGS:
;   The match distances are **hard-wired** to 6 arcsec on the first iteration,
;   and 3 arcsec on the 2nd iteration???
;
; PROCEDURES CALLED:
;   angle_from_pairs()
;   astrom_tweak
;   gsssadxy
;   gsssxyad
;   offset_from_pairs
;
; REVISION HISTORY:
;   10-Jun-2002  Written by D. Schlegel & D. Finkbeiner, Princeton.
;-
;------------------------------------------------------------------------------
function astrom_engine, xpos, ypos, catlon, catlat, gsa_in, $
 search_rad=search_rad, $
 search_scale=search_scale, search_angle=search_angle, $
 poserr=poserr, nmatch=nmatch, catind=catind, $
 obsind=obsind, radial=radial, verbose=verbose, start_angle=start_angle

   if (n_params() LT 5) then begin
      doc_library, 'astrom_engine'
      return, 0
   endif
   if (size(gsa_in,/tname) NE 'STRUCT') then $
    message, 'Must pass gsa structure with initial guess'

   if (NOT keyword_set(poserr)) then poserr = 1.0

   ;----------
   ; Set default return values

   catind = -1L
   obsind = -1L
   nmatch = 0L

   ;----------
   ; Set the bin size used in the ANGLE_FROM_PAIRS code,
   ; based upon the value of POSERR.
   ; The factor 0.6 deg is pulled out of thin air, but should be based
   ; upon the step size used to search for the position angle error.

   xsep = (max(xpos) - min(xpos)) > poserr
   ysep = (max(ypos) - min(ypos)) > poserr
   binsz = 0.5 * (poserr + 0.6 * sqrt(xsep^2 + ysep^2) / !radeg)
   dmax = xsep > ysep

   gsa1 = gsa_in
   gsssadxy, gsa1, catlon, catlat, catx, caty

   ;----------
   ; Search for the angle between the catalogue stars and image stars

   if (search_angle GT 1. OR keyword_set(start_angle) gt 0) then begin 
      if(search_angle gt 1.) then begin
	      if(NOT keyword_set(start_angle)) then start_angle=0.
        ang = angle_from_pairs(catx, caty, xpos, ypos, $
         dmax=dmax, binsz=binsz, bestsig=bestsig, $
         angrange=[-1.,1.]*search_angle+start_angle, verbose=verbose)
      endif else begin
        ang=start_angle
        bestsig=1000
      endelse

      if (bestsig gt 12) then begin 
         if (keyword_set(verbose)) then $
          splog, 'Best angle: ', ang, '  sigma: ', bestsig
         cd = dblarr(2, 2)
         cd[0,0] = gsa1.amdx[0]
         cd[0,1] = gsa1.amdx[1]
         cd[1,1] = gsa1.amdy[0]
         cd[1,0] = gsa1.amdy[1]
         angrad = ang * !dpi / 180.d0
         mm = [[cos(angrad), sin(angrad)], [-sin(angrad), cos(angrad)]]
         cd = cd # mm
         gsa1.amdx[0] = cd[0,0]
         gsa1.amdx[1] = cd[0,1]  
         gsa1.amdy[0] = cd[1,1]
         gsa1.amdy[1] = cd[1,0]
         gsssadxy, gsa1, catlon, catlat, catx, caty
      endif else begin
         if (keyword_set(verbose)) then $
          splog, 'Warning: I think I am lost, but I will try anyway...'
      endelse
   endif 
  
   ;----------
   ; Search for the X,Y offset between the catalogue stars and image stars

   xyshift = offset_from_pairs(catx, caty, xpos, ypos, $
    dmax=dmax, binsz=binsz, errflag=errflag, bestsig=bestsig, verbose=verbose)

   if (errflag NE 0) then begin
      splog, 'XY shift FAILED'
      return, 0
   endif

   splog, 'XYSHIFT: ', xyshift, ' pix'
  
   xcen = gsa1.ppo3 / gsa1.xsz - 0.5d
   ycen = gsa1.ppo6 / gsa1.ysz - 0.5d

   ; NOTE: FITS crpix is 1-indexed but argument of xy2ad is 0-indexed
   refpix = [xcen, ycen] - xyshift
   gsssxyad, gsa1, refpix[0], refpix[1], racen, deccen

   ; Update astrometry structure with new CRVALs
   gsa1.crval = [racen, deccen]

   ;----------
   ; Tweak astrometry structure with cat (ra,dec) and im (x,y) comparison

   ; Call twice, since more stars may be matched on the 2nd call.
   ; The match distances are **hard-wired** to 6 arcsec on the first
   ; iteration, and 3 arcsec on the 2nd iteration???
   gsa1 = astrom_tweak(gsa1, catlon, catlat, xpos, ypos, errflag=errflag, $
    nmatch=nmatch, catind=catind, obsind=obsind, dtheta=6./3600., $
    verbose=verbose)
   gsa1 = astrom_tweak(gsa1, catlon, catlat, xpos, ypos, errflag=errflag, $
    nmatch=nmatch, catind=catind, obsind=obsind, dtheta=3./3600., $
    verbose=verbose)

   splog, 'Number of catalog matches = ', nmatch
   if (errflag NE 0) then begin
      splog, 'WARNING: XY shift FAILED in astrom_tweak'
      return, 0
   endif

   ; Compute rotation angle w.r.t. initial guess
   xvec0 = gsa_in.amdx
   xvec1 = gsa1.amdx
   angerr = acos((transpose(xvec0[0:1]) # xvec1[0:1]) $
    / sqrt(total(xvec0[0:1]^2)*total(xvec1[0:1]^2))) / !dtor
   splog, 'Initial guess rotated by: ', angerr, ' deg'

   ; Compute the offset of the corner pixel (0,0)
   gsssxyad, gsa_in, 0.0, 0.0, ra_in, dec_in
   gsssxyad, gsa1, 0.0, 0.0, ra_out, dec_out
   gsssadxy, gsa1, ra_in, dec_in, xoff, yoff
   splog, 'RA offset = ', (ra_out - ra_in) * 3600. / cos(dec_out/!radeg), ' arcsec'
   splog, 'DEC offset = ', (dec_out - dec_in) * 3600., ' arcsec'
   splog, 'X offset = ', xoff, ' pix'
   splog, 'Y offset = ', yoff, ' pix'

   return, gsa1
end
;------------------------------------------------------------------------------
