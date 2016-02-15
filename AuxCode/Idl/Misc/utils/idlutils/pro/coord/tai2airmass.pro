;+
; NAME:
;   tai2airmass
;
; PURPOSE:
;   Compute airmass.
;
; CALLING SEQUENCE:
;   airmass = tai2airmass( ra, dec, [ equinox, jd=, tai=, mjd=, $
;    longitude=, latitude=, altitude=, ha=, ipa= ] )
;
; INPUTS:
;   ra             - Right ascension [degrees]
;   dec            - Declination [degrees]
;   equinox        - Equinox of observation for RA, DEC; default to 2000.
;
; OPTIONAL KEYWORDS:
;   jd             - Decimal Julian date.  Note this should probably be
;                    type DOUBLE.
;   tai            - Number of seconds since Nov 17 1858
;                    Note this should probably either be type DOUBLE or LONG64.
;   mjd            - Modified Julian date.
;   longitude      - Longitude of observatory;
;                    default to (360-105.820417) deg for APO
;   latitute       - Latitude of observatory; default to 32.780361 deg for APO
;   altitude       - Altitude of observatory; default to 2788 m for APO
;   node           - Node of great circle on the sky (degrees); required
;                    if returning IPA.
;   incl           - Inclination of great circle on the sky (degrees); required
;                    if returning IPA.
;
; OUTPUTS:
;   airmass        - Airmass; 1.0 for zenith
;
; OPTIONAL OUTPUTS:
;   ha             - Hour angle (degrees)
;   ipa            - Position angle for image rotator (degrees)
;
; COMMENTS:
;   TAI, JD, or MJD must be specified.
;
;   This routine only returns sec(z) for the airmass.
;   Formula from Smart, Spherical Astronomy.
;
; EXAMPLES:
;
; BUGS:
;   Outputs SLIGHTLY different airmasses from those computed by the PT
;     system.  We think that they may be going to second order.
;   EQUINOX does nothing except for the IPA calculation!
;   ALTITUDE is unused!
;
; PROCEDURES CALLED:
;   ct2lst
;   ll2uv()
;   precess
;
; REVISION HISTORY:
;   10-May-2000  Written by D. Schlegel, Princeton, & D. Hogg, IAS
;   02-Jun-2000  Fixed minor bugs, Schlegel
;   05-Nov-2000  Added HA keyword
;-
;------------------------------------------------------------------------------
function tai2air_crossprod, aa, bb

  if (n_elements(A) GT 3) or (n_elements(B) GT 3) then begin
     print, 'Sorry - only 3-vectors, one at a time'
     return, 0
  endif

  A = reform(aa, 3) & B=reform(bb, 3)
  C = (A-A)*1.                       ; zero vector of type float or double
  C[0] = determ([[1, 0, 0], [A], [B]])
  C[1] = determ([[0, 1, 0], [A], [B]])
  C[2] = determ([[0, 0, 1], [A], [B]])

  return, C
end
;------------------------------------------------------------------------------
; angle (Degrees) between any vectors A, B (not necessarily unit vectors)

function tai2air_ang, a, b
  c = (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
  c = ((c/sqrt(total(a^2)*total(b^2))) < 1) > (-1)
  theta = acos(c)*180./!dpi

  return, theta
end
;------------------------------------------------------------------------------
function tai2airmass, ra, dec, equinox1, jd=jd, tai=tai, mjd=mjd, $
 longitude=longitude, latitude=latitude, altitude=altitude, $
 node=node, incl=incl, ha=hadeg, ipa=ipa

   ; Default to location of Apache Point Observatory
   if (n_elements(equinox1) NE 0) then equinox = equinox1 $
    else equinox = 2000.d0
   if (NOT keyword_set(longitude)) then longitude = 360. - 105.820417d0
   if (NOT keyword_set(latitude)) then latitude = 32.780361d0
   if (NOT keyword_set(altitude)) then altitude = 2788.

   if (NOT keyword_set(jd)) then begin
      if (keyword_set(tai)) then begin
         jd = 2400000.5D + tai / (24.D*3600.D)
      endif else if (keyword_set(mjd)) then begin
         jd = 2400000.5D + mjd
      endif
   endif

   if (NOT keyword_set(jd)) then begin
      message, 'Must specify TAI, JD or MJD', /cont
      return, 0
   endif

   DRADEG = 180.d0 / !DPI

   ;----------
   ; Compute the hour angle, HA, in degrees

   ct2lst, lst, longitude, junk, jd
   lst = 15. * lst ; convert from hours to degrees
   hadeg = lst - ra

   decrad = dec / DRADEG
   harad = hadeg / DRADEG
   latrad = latitude / DRADEG

   ;----------
   ; Compute airmass with spherical trig

   cosz = sin(decrad)*sin(latrad)+ cos(decrad)*cos(harad)*cos(latrad)
   airmass = 1.d0 / cosz

   ;----------
   ; Optionally compute IPA

   nra = n_elements(ra)
   if (arg_present(ipa) AND n_elements(node) GT 0 AND n_elements(incl) GT 0) $
    then begin
      ipa = dblarr(nra)
   endif

   for i=0L, n_elements(ipa)-1 do begin
      thisnode = node[i < (n_elements(node)-1)]
      thisincl = incl[i < (n_elements(incl)-1)]
      thisequnx = equinox[i < (n_elements(equinox)-1)]

      ; Precess to date of observation
      day = jd[i<(nra-1)] - 2451544.5d0
      yr = 2000.d0 + day / 365.2425d0
      radate = ra[i]
      decdate = dec[i]
      precess, radate, decdate, thisequnx, yr

      ; Rotation matrix
      xrot = [[1, 0, 0], $
              [0, cos(thisincl/DRADEG), -sin(thisincl/DRADEG)], $
              [0, sin(thisincl/DRADEG), cos(thisincl/DRADEG)]]

      ; Unit vectors in (HA, dec) coords
      bvec   = ll2uv([[hadeg[i]], [decdate]], /double)
      zenith = ll2uv([[0.], [latitude]], /double)
      pole   = ll2uv([[0.], [90.000000d]], /double)

      ; We want the normal vector of the great circle (of the stripe)
      ; in HA,DEC coordinates!
      node2 = hadeg[i] + radate - thisnode
      zrot2 = [[cos(node2/DRADEG), -sin(node2/DRADEG), 0], $
       [sin(node2/DRADEG), cos(node2/DRADEG), 0], $
       [0, 0, 1]]
      gcnorm = zrot2 ## (invert(xrot) ## (invert(zrot2) ## pole))

      if (total(finite(bvec)) NE 3) then message, 'Floating exception'
      plane1 = tai2air_crossprod(zenith, bvec); plane containing zenith and bvec
      plane2 = tai2air_crossprod(pole, bvec)  ; plane containing pole and bvec
      plane3 = tai2air_crossprod(gcnorm, bvec); plane containing GC norm and b

      theta  = tai2air_ang(plane1, plane2)
      pasky  = tai2air_ang(plane3, plane2)

      ; Fix sign of theta if at positive hour angle...
      if (abs( tai2air_ang(tai2air_crossprod(plane1, plane2), bvec) ) LT 90) $
       then theta = -theta

      ; Fix sign of pasky
      if (abs( tai2air_ang(tai2air_crossprod(plane3, plane2), bvec) ) LT 90) $
       then pasky = -pasky

      ; theta plus rotator angle give image rotation in focal plane
      ipa[i] = (pasky + 90 - theta + 720) MOD 360
   endfor

   return, airmass
end
;------------------------------------------------------------------------------
