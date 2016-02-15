;+
; NAME:
;   heliocentric
;
; PURPOSE:
;   Compute correction term to add to velocities to convert to heliocentric.
;
; CALLING SEQUENCE:
;   vcorr = heliocentric( ra, dec, [ epoch, jd=, tai=, $
;    longitude=, latitude=, altitude= ] )
;
; INPUTS:
;   ra             - Right ascension [degrees]
;   dec            - Declination [degrees]
;   epoch          - Epoch of observation for RA, DEC; default to 2000.
;
; OPTIONAL KEYWORDS:
;   jd             - Decimal Julian date.  Note this should probably be
;                    type DOUBLE.
;   tai            - Number of seconds since Nov 17 1858; either JD or TAI
;                    must be specified.  Note this should probably either
;                    be type DOUBLE or LONG64.
;   longitude      - Longitude of observatory;
;                    default to (360-105.820417) deg for APO
;   latitute       - Latitude of observatory; default to 32.780361 deg for APO
;   altitude       - Altitude of observatory in meters;
;                    default to 2788 m for APO
;
; OUTPUTS:
;   vcorr          - Velocity correction term, in km/s, to add to measured
;                    radial velocity to convert it to the heliocentric frame.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   baryvel
;   ct2lst
;
; REVISION HISTORY:
;   09-May-2000  Written by S. Burles & D. Schlegel
;-
;------------------------------------------------------------------------------

function heliocentric, ra, dec, epoch, jd=jd, tai=tai, $
 longitude=longitude, latitude=latitude, altitude=altitude 

   if (NOT keyword_set(epoch)) then epoch = 2000.0

   ; Default to location of Apache Point Observatory
   if (NOT keyword_set(longitude)) then longitude = 360. - 105.820417
   if (NOT keyword_set(latitude)) then latitude = 32.780361
   if (NOT keyword_set(altitude)) then altitude = 2788.

   if (NOT keyword_set(jd)) then begin
      if (keyword_set(tai)) then begin
         jd = 2400000.5D + tai / (24.D*3600.D)
      endif else begin
         message, 'Must specify either JD or TAI', /cont
         return, 0
      endelse
   endif

   DRADEG = 180.d0 / !DPI

   ;----------
   ; Compute baryocentric velocity

   baryvel, jd, epoch, dvelh, dvelb

   ; Project velocity toward star
   vbarycen = dvelb[0]*cos(dec/DRADEG)*cos(ra/DRADEG) + $
            dvelb[1]*cos(dec/DRADEG)*sin(ra/DRADEG) + dvelb[2]*sin(dec/DRADEG) 

   ;----------
   ; Compute rotational velocity of observer on the Earth

   ; LAT is the latitude in radians.
   latrad = latitude / DRADEG

   ; Reduction of geodetic latitude to geocentric latitude (radians).
   ; DLAT is in arcseconds.

   dlat = -(11.d0 * 60.d0 + 32.743000d0) * sin(2.d0 * latrad) + $
            1.163300d0 * sin(4.d0 * latrad) -0.002600d0 * sin(6.d0 * latrad)
   latrad  = latrad + (dlat / 3600.d0) / DRADEG

   ; R is the radius vector from the Earth's center to the observer (meters).
   ; VC is the corresponding circular velocity
   ; (meters/sidereal day converted to km / sec).
   ; (sidereal day = 23.934469591229 hours (1986))

   r = 6378160.0d0 * (0.998327073d0 + 0.00167643800d0 * cos(2.d0 * latrad) - $
       0.00000351d0 * cos(4.d0 * latrad) + 0.000000008d0 * cos(6.d0 * latrad)) $
       + altitude
   vc = 2.d0 * !DPI * (r / 1000.d0)  / (23.934469591229d0 * 3600.d0)

   ; Compute the hour angle, HA, in degrees
   ct2lst, LST, longitude, junk, jd
   LST = 15. * LST ; convert from hours to degrees
   HA = LST - ra

   ; Project the velocity onto the line of sight to the star.
   vrotate = vc * cos(latrad) * cos(dec/DRADEG) * sin(HA/DRADEG)

   return, (-vbarycen + vrotate)
end

