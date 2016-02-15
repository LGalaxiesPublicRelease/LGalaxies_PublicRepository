;------------------------------------------------------------------------------
;+
; NAME:
;   hst_guidetest
;
; PURPOSE:
;   Decide if a list of potential HST guide stars are good by looking
;   at the 2MASS and UCAC catalogs
;
; CALLING SEQUENCE:
;   hst_guidetest, [ ra, dec, name=, filename=, searchrad=, select_tags=, $
;    goodrad= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   ra         - RA for each star [degrees]
;   dec        - Declination for each star [degrees]
;   name       - Name for each star; default to STAR-1, STAR-2, etc.
;   filename   - Input ASCII file name with at least 3 columns:
;                name, RA, Dec (where the coordinates are in degrees);
;                if set, then this overrides the inputs RA, DEC, NAME
;   searchrad  - Search radius between input coorinates and 2MASS and UCAC;
;                default to 5./3600 degrees
;   goodrad    - Maximum distance for a "good" HST guide star; default
;                to 1.5/3600 arcsec
;   select_tags- Names of tags to print
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The default is to return the following:
;     NAME - The object names either from the first column of FILENAME,
;            or from the input NAME
;     DIST - Distance in arcsec to the nearest 2MASS star
;     BL_FLG - Blend flags in J,H,K-bands, should be "111"
;     CC_FLG - Contamination flags in J,H,K-bands, should be "000"
;     GAL_CONTAM - Extended source contamination flag, should be "0"
;     MP_FLG - Minor planet flag, should be "0"
;     UCAC_DIST - Distance in arcsec to the nearest UCAC star
;     GOOD - Set to "BAD" if there is not a 2MASS and UCAC isolated star
;            within GOODRAD arcsec, or if any of the other, above flags have
;            suspicious values; note that it will always say "BAD" for
;            stars at high declination where UCAC does not have coverage
;
;   Note that the UCAC astrometric catalog includes stars in the
;   magnitude range R = [8,16] (approximately) from declination -90
;   to between Dec +40 and +52.  It discards stars within 3 arcsec
;   of any other star.
;
;   Refer to the following 2MASS documentation for a description
;   of the 2MASS data structure:
;     http://tdc-www.harvard.edu/software/catalogs/tmc.format.html
;
;   Refer to the following web site for the UCAC astrometric catalog:
;     http://ad.usno.navy.mil/ucac/
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_diff_angle()
;   struct_addtags()
;   struct_print
;   struct_selecttags()
;   tmass_read()
;   ucac_read()
;
; REVISION HISTORY:
;   09-Jul-2005  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro hst_guidetest, radeg, decdeg, name_in, filename=filename, $
 searchrad=searchrad1, select_tags=select_tags, goodrad=goodrad1

   if (keyword_set(searchrad1)) then searchrad = searchrad1[0] $
    else searchrad = 5./3600
   if (keyword_set(goodrad1)) then goodrad = goodrad1[0] $
    else goodrad = 1.5/3600
   if (NOT keyword_set(select_tags)) then $
    select_tags = ['NAME', 'TMASS_DIST', $
     'TMASS_BL_FLG', 'TMASS_CC_FLG', 'TMASS_GAL_CONTAM', 'TMASS_MP_FLG', $
     'UCAC_DIST','GOOD']

   if (keyword_set(filename)) then begin
      readcol, filename, name, radeg, decdeg, format='(a,d,d)', /silent
      nstar = n_elements(radeg)
      if (nstar EQ 0) then begin
         splog, 'No objects in file ', filename
         return
      endif
   endif else begin
      nstar = n_elements(radeg)
      if (nstar EQ 0 OR n_elements(decdeg) NE nstar) then begin
         splog, 'Wrong number of elements for RA,DEC'
         return
      endif
      if (n_elements(name_in) EQ nstar) then name = name_in $
       else name = 'STAR-' + strtrim(string(lindgen(nstar)+1),2)
   endelse

   for istar=0L, nstar-1L do begin
      obj1 = tmass_read(radeg[istar], decdeg[istar], searchrad)
      if (keyword_set(objs) EQ 0 AND keyword_set(obj1)) then begin
         blankobj = obj1
         struct_assign, {junk:0}, blankobj
         objs = replicate(blankobj, nstar)
      endif
      if (keyword_set(obj1)) then objs[istar] = obj1[0]
   endfor

   ; Prepend NAME,TMASS_DIST to the structure
   addtags = replicate(create_struct('NAME', '', 'TMASS_DIST', 0.), nstar)
   addtags.name = name
   if (keyword_set(objs)) then objs = struct_addtags(addtags, objs) $
    else objs = addtags
   if (tag_exist(objs, 'TMASS_RA')) then begin
      adiff = djs_diff_angle(radeg, decdeg, $
       objs.tmass_ra, objs.tmass_dec) * 3600.
      objs.tmass_dist = adiff
   endif

   ;----------
   ; Decide if there is a UCAC match

   addtags = replicate(create_struct('UCAC_DIST', 0.), nstar)
   objs = struct_addtags(objs, addtags)
   for istar=0L, nstar-1L do begin
      obj2 = ucac_read(racen=radeg[istar], deccen=decdeg[istar], $
       radius=searchrad)
      if (keyword_set(obj2)) then $
       objs[istar].ucac_dist = djs_diff_angle(radeg[istar], decdeg[istar], $
        obj2.ramdeg, obj2.demdeg) * 3600.
   endfor

   ;----------
   ; Decide if each star is a good HST guide star

   addtags = replicate(create_struct('GOOD', ''), nstar)
   objs = struct_addtags(objs, addtags)
   qbad = objs.tmass_dist EQ 0 OR objs.tmass_dist GT goodrad*3600. $
    OR objs.ucac_dist EQ 0 OR objs.ucac_dist GT goodrad*3600.
   if (tag_exist(objs,'TMASS_BL_FLG')) then $
    qbad = qbad $
     OR objs.tmass_bl_flg NE 111 $
     OR objs.tmass_cc_flg NE '000' $
     OR objs.tmass_gal_contam NE 0 $
     OR objs.tmass_mp_flg NE 0 $
   else $
    qbad[*] = 1B
   ibad = where(qbad, nbad)
   if (nbad GT 0) then objs[ibad].good = 'BAD'

   ;----------
   ; Print results

   printobj = struct_selecttags(objs, select_tags=select_tags)

   ; Strip '2MASS_' from the beginning of the tag names for printing purposes
   alias = strarr(2,n_elements(select_tags))
   alias[0,*] = select_tags
   alias[1,*] = repstr(select_tags, 'TMASS_', '')
   struct_print, printobj, alias=alias

   return
end
;------------------------------------------------------------------------------
