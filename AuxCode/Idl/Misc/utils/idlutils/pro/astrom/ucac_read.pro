;+
; NAME:
;   ucac_read()
;
; PURPOSE:
;   Read the UCAC catalog
;
; CALLING SEQUENCE:
;   outdat = ucac_read( [ racen=, deccen=, radius=, node=, incl=, hwidth=, $
;    decrange= ] )
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   racen:       Central RA for selecting a region of stars [J2000 deg]
;   deccen:      Central DEC for selecting a region of stars [J2000 deg]
;   radius:      Radius for selecting a region of stars [deg]
;   node:        Node of great circle for selecting a stripe of stars [deg]
;   incl:        Inclination of great circle for selecting a stripe [deg]
;   hwidth:      Half-width of great circle for selecting a stripe [deg]
;   decrange   - Declination range for data; default to [-90,90] degrees
;
; OUTPUT:
;   outdat     - Structure with UCAC data in its raw catalog format;
;                return 0 if no stars found
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Either RACEN, DECCEN, RADIUS must be set, or NODE, INCL, HWIDTH.
;   If all keywords are set, the the catalog is trimmed with both sets
;   of conditions.
;
; EXAMPLES:
;   a=ucac_read(racen=180.,deccen=0.,radius=0.1)
;
; BUGS:
;
; PROCEDURES CALLED:
;   cirrange
;   djs_diff_angle()
;   radec_to_munu
;   ucac_readindex()
;   ucac_readzone()
;
; REVISION HISTORY:
;   16-Aug-2005  Written by D. Schlegel (LBL)
;-
;------------------------------------------------------------------------------
function ucac_intersect, dec0, node, incl, nu

   DRADEG = 180.d0 / !dpi
   sini = sin(incl/DRADEG)
   cosi = cos(incl/DRADEG)
   sinnu = sin(nu/DRADEG)
   cosnu = cos(nu/DRADEG)
   sinmu = (sin(dec0/DRADEG) - sinnu * cosi) / (cosnu * sini)
   sinmu = ((sinmu) < 1.d0) > (-1.d0) ; Sanity-check this for round-off errors
   cosmu = sqrt(1 - sinmu^2)
   cosmu = [-cosmu,cosmu]
   yy = sinmu * cosnu * cosi - sinnu * sini
   xx = cosmu * cosnu
   ra = node + atan(yy,xx) * DRADEG
   cirrange, ra

   return, ra
end
;------------------------------------------------------------------------------
function ucac_readgc_add, outdat, thiszone, ravec

   if (ravec[0] GT 360) then ravec = ravec - 360
   newdat = ucac_readzone(thiszone, ravec[0], ravec[1])
   if (NOT keyword_set(newdat)) then return, outdat
   if (NOT keyword_set(outdat)) then return, newdat
   return, [outdat, newdat]
end
;------------------------------------------------------------------------------
; Note that DECMIN,DECMAX is the declination range of zone number THISZONE.
function ucac_readcircle1, racen, deccen, radius, thiszone, decmin, decmax

   if ((decmax LT deccen-radius) OR (decmin GT deccen+radius)) then return, 0

   DRADEG = 180.d0 / !dpi
   if (deccen GT decmax) then begin
      ravec = racen + [-1,1] * radius / cos(decmax/DRADEG)
   endif else if (deccen LT decmin) then begin
      ravec = racen + [-1,1] * radius / cos(decmin/DRADEG)
   endif else begin
      ravec = racen + [-1,1] * radius / cos(deccen/DRADEG)
   endelse

   ; Now read the data, but watch for wrapping at RA=360 degrees!
   outdat = 0
   if (ravec[1] - ravec[0] GT 360) then begin
      outdat = ucac_readgc_add(outdat, thiszone, [0,360])
   endif else if (ravec[0] LT 0) then begin
      outdat = ucac_readgc_add(outdat, thiszone, [360-ravec[0],360])
      outdat = ucac_readgc_add(outdat, thiszone, [0,ravec[1]])
   endif else if (ravec[1] GT 360) then begin
      outdat = ucac_readgc_add(outdat, thiszone, [ravec[0],360])
      outdat = ucac_readgc_add(outdat, thiszone, [0,ravec[1]-360])
   endif else begin
      outdat = ucac_readgc_add(outdat, thiszone, ravec)
   endelse

   return, outdat
end
;------------------------------------------------------------------------------
; Note that DECMIN,DECMAX is the declination range of zone number THISZONE.
function ucac_readgc1, node, incl, hwidth, thiszone, decmin, decmax

   ; These are the declinations of the great circle
   dec1max = incl - hwidth
   dec2max = incl + hwidth
   dec1min = -incl - hwidth ; = -dec2max
   dec2min = -incl + hwidth ; = -dec1max

   ; CASE: Zone is completely above the great circle
   if (dec2max LE decmin) then return, 0

   ; CASE: Zone is completely below the great circle
   if (decmax LE dec1min) then return, 0

   ; CASE: Zone overlaps for entire 360 degrees
   if (dec2min GE decmin AND dec1max LE decmax) then begin
      ravec = [0.d0, 360.d0]

   ; CASE: Take a cap at the top of the great circle
   endif else if (decmax GE dec1max) then begin
      nu = incl - dec1max
      ravec = ucac_intersect(decmin, node, incl, nu)
      if (nu GT 0) then ravec = reverse(ravec)

   ; CASE: Take a cap at the bottom of the great circle
   endif else if (decmin LE dec2min) then begin
      nu = incl - dec2max
      ravec = ucac_intersect(decmax, node, incl, nu)
      if (nu GT 0) then ravec = reverse(ravec)

   ; CASE: Take two intersections...
   endif else begin
      ravec1 = ucac_intersect(decmax, node, incl, dec1max-incl)
      ravec2 = ucac_intersect(decmin, node, incl, dec2max-incl)
      ravec = [ravec1[0], ravec2[0], ravec2[1], ravec1[1]]
   endelse

   ; Now read the data, but watch for wrapping at RA=360 degrees!
   outdat = 0
   if (ravec[0] LT ravec[1]) then begin
      outdat = ucac_readgc_add(outdat, thiszone, ravec[0:1])
   endif else begin
      outdat = ucac_readgc_add(outdat, thiszone, [ravec[0], 360])
      outdat = ucac_readgc_add(outdat, thiszone, [0, ravec[1]])
   endelse
   if (n_elements(ravec) GT 2) then begin
      if (ravec[2] LT ravec[3]) then begin
         outdat = ucac_readgc_add(outdat, thiszone, ravec[2:3])
      endif else begin
         outdat = ucac_readgc_add(outdat, thiszone, [ravec[2], 360])
         outdat = ucac_readgc_add(outdat, thiszone, [0, ravec[3]])
      endelse
   endif

   return, outdat
end
;------------------------------------------------------------------------------
function ucac_read, racen=racen, deccen=deccen, radius=radius, $
 node=node1, incl=incl1, hwidth=hwidth, decrange=decrange1

   common com_ucac, uindex

   outdat = 0

   ;----------
   ; Check inputs

   if (n_elements(decrange1) EQ 2) then decrange = decrange1 $
    else decrange = [-90,90]

   qreadcircle = (n_elements(racen) EQ 1 AND n_elements(deccen) EQ 1 $
    AND n_elements(radius) EQ 1)
   qreadstripe = (n_elements(node1) EQ 1 AND n_elements(incl1) EQ 1 $
    AND n_elements(hwidth) EQ 1)
   if (qreadcircle EQ 0 AND qreadstripe EQ 0) then return, 0

   ; For all our case statements to work, we need the node to be
   ; in the range [0,+90] degrees, so just transform the higher-inclination
   ; great circles.  This means that the mu values along the great circle
   ; are incorrect, but we don't care since we're returning the full circle.
   if (qreadstripe) then begin
      if (incl1 LE 90) then begin
         node = node1
         incl = incl1
      endif else begin
         node = node1 + 180
         cirrange, node
         incl = 180 - incl1
      endelse
   endif

   ;----------
   ; Read the index file

   uindex = ucac_readindex()

   ;----------
   ; Loop over all zones

   for thiszone=min(uindex.zn), max(uindex.zn) do begin
      jj = (where(uindex.zn EQ thiszone, ct))[0]
      if (ct GT 0) then begin
         decmax = uindex[jj].dcmax
         decmin = decmax - 0.5d0
         if (decmax GE decrange[0] AND decmin LE decrange[1]) then begin
            if (qreadstripe) then begin
               moredat = ucac_readgc1(node, incl, hwidth, thiszone, $
                decmin, decmax)
            endif else begin
               moredat = ucac_readcircle1(racen, deccen, radius, thiszone, $
                decmin, decmax)
            endelse
            if (keyword_set(moredat)) then $
             outdat = keyword_set(outdat) ? [outdat,moredat] : moredat
         endif
      endif
   endfor

   ;----------
   ; Now trim to objects exactly within the great circle bounds
   ; or to the requested circle on the sky

   if (keyword_set(outdat)) then begin
      qkeep = bytarr(n_elements(outdat)) + 1B
      if (qreadcircle) then begin
         adiff = djs_diff_angle(outdat.ramdeg, outdat.demdeg, racen, deccen)
         qkeep = qkeep AND (adiff LE radius)
      endif
      if (qreadstripe) then begin
         radec_to_munu, outdat.ramdeg, outdat.demdeg, node=node, incl=incl, $
          mu, nu
         qkeep = qkeep AND (abs(nu) LE hwidth)
      endif
      ikeep = where(qkeep, nkeep)
      if (nkeep EQ 0) then outdat = 0 $
       else outdat = outdat[ikeep]
   endif

   ;----------
   ; Now trim to objects exactly within the declination range

   if (keyword_set(decrange1) AND keyword_set(outdat)) then begin
      ikeep = where(outdat.demdeg GE decrange[0] $
       AND outdat.demdeg LE decrange[1], nkeep)
      if (nkeep EQ 0) then outdat = 0 $
       else outdat = outdat[ikeep]
   endif

   return, outdat
end
;------------------------------------------------------------------------------
