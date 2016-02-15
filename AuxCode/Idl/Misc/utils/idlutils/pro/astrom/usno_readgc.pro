;+
; NAME:
;   usno_readgc()
;
; PURPOSE:
;   Read the UCAC data files for a great circle on the sky.
;
; CALLING SEQUENCE:
;   outdat = usno_readgc(node=node, incl=incl, hwidth=, [ decrange= ])
;
; INPUTS:
;   node       - Node of great circle [degrees]
;   incl       - Inclination of great circle [degrees]
;   hwidth     - Half-width of great circle for selecting a stripe [deg]
;
; OPTIONAL INPUTS:
;   decrange   - Declination range for data; default to [-90,90] degrees
;
; OUTPUT:
;   outdat     - Structure with UCAC data in its raw catalog format;
;                return 0 if no stars found
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   a=usno_readgc(node=95.,incl=40.,hwidth=2.5)
;   a=usno_readgc(node=95.,incl=10.,hwidth=1.0)
;
; BUGS:
;
; PROCEDURES CALLED:
;   radec_to_munu
;   usno_readindex()
;   usno_readzone()
;
; REVISION HISTORY:
;   28-May-2003  Written by D. Schlegel and N. Padmanabhan, Princeton.
;-
;------------------------------------------------------------------------------
function usno_intersect, dec0, node, incl, nu

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
function usno_readgc_add, outdat, thiszone, ravec

   common com_usno_readgc, catpath ; Cache between calls

   if (NOT keyword_set(catpath)) then catpath = getenv('USNO_DIR')

   if (ravec[0] GT 360) then ravec = ravec - 360

   thispath = concat_dir(catpath, string(thiszone/10,format='(i3.3)'))
   usno_readzone, thispath, thiszone, ravec[0], ravec[1], 80L, 'b', $
    newdat, /swap_if_big_endian

   if (NOT keyword_set(newdat)) then return, outdat
   if (NOT keyword_set(outdat)) then return, newdat
   return, [[outdat], [newdat]]
end
;------------------------------------------------------------------------------
function usno_readgc1, node, incl, hwidth, thiszone, decmin, decmax

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
      ravec = usno_intersect(decmin, node, incl, nu)
      if (nu GT 0) then ravec = reverse(ravec)

   ; CASE: Take a cap at the bottom of the great circle
   endif else if (decmin LE dec2min) then begin
      nu = incl - dec2max
      ravec = usno_intersect(decmax, node, incl, nu)
      if (nu GT 0) then ravec = reverse(ravec)

   ; CASE: Take two intersections...
   endif else begin
      ravec1 = usno_intersect(decmax, node, incl, dec1max-incl)
      ravec2 = usno_intersect(decmin, node, incl, dec2max-incl)
      ravec = [ravec1[0], ravec2[0], ravec2[1], ravec1[1]]
   endelse

   ; Now read the data, but watch for wrapping at RA=360 degrees!
   outdat = 0
   if (ravec[0] LT ravec[1]) then begin
      outdat = usno_readgc_add(outdat, thiszone, ravec[0:1])
   endif else begin
      outdat = usno_readgc_add(outdat, thiszone, [ravec[0], 360])
      outdat = usno_readgc_add(outdat, thiszone, [0, ravec[1]])
   endelse
   if (n_elements(ravec) GT 2) then begin
      if (ravec[2] LT ravec[3]) then begin
         outdat = usno_readgc_add(outdat, thiszone, ravec[2:3])
      endif else begin
         outdat = usno_readgc_add(outdat, thiszone, [ravec[2], 360])
         outdat = usno_readgc_add(outdat, thiszone, [0, ravec[3]])
      endelse
   endif

   return, outdat
end
;------------------------------------------------------------------------------
function usno_readgc, node=node1, incl=incl1, hwidth=hwidth, decrange=decrange1

   outdat = 0

   ;----------
   ; Check inputs

   if (n_elements(node1) NE 1 OR n_elements(incl1) NE 1 $
    OR n_elements(hwidth) NE 1) then begin
      print, 'Wrong number of parameters!'
      return, 0
   endif
   if (keyword_set(decrange1)) then decrange = decrange1 $
    else decrange = [-90,90]

   ; For all our case statements to work, we need the node to be
   ; in the range [0,+90] degrees, so just transform the higher-inclination
   ; great circles.  This means that the mu values along the great circle
   ; are incorrect, but we don't care since we're returning the full circle.
   if (incl1 LE 90) then begin
      node = node1
      incl = incl1
   endif else begin
      node = node1 + 180
      cirrange, node
      incl = 180 - incl1
   endelse

   ;----------
   ; Loop over all zones

   zonewidth = 0.1
   zone0 = floor((90+decrange[0])/zonewidth) < long(180.0/zonewidth - 1)
   zone1 = floor((90+decrange[1])/zonewidth) < long(180.0/zonewidth - 1)

   for thiszone=zone0, zone1 do begin
      decmin = -90 + thiszone * zonewidth
      decmax = decmin + zonewidth
      moredat = usno_readgc1(node, incl, hwidth, thiszone, decmin, decmax)
      if (keyword_set(moredat)) then $
       outdat = keyword_set(outdat) ? [[outdat],[moredat]] : moredat
   endfor

   ;----------
   ; Convert from binary format data to something somewhat rational.

   outdat = usnob10_extract(outdat)

   ;----------
   ; Now trim to objects exactly within the great circle bounds

   if (keyword_set(outdat)) then begin
      radec_to_munu, outdat.ra, outdat.dec, mu, nu, node=node, incl=incl
      ikeep = where(abs(nu) LE hwidth, nkeep)
      if (nkeep EQ 0) then outdat = 0 $
       else outdat = outdat[ikeep]
   endif

   ;----------
   ; Now trim to objects exactly within the declination range

   if (keyword_set(decrange1) AND keyword_set(outdat)) then begin
      ikeep = where(outdat.dec GE decrange[0] AND outdat.dec LE decrange[1], $
       nkeep)
      if (nkeep EQ 0) then outdat = 0 $
       else outdat = outdat[ikeep]
   endif

   return, outdat
end
;------------------------------------------------------------------------------
