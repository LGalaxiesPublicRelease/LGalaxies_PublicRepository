;-----------------------------------------------------------------------
;+
; NAME:
;   dotproduct_sphere
; PURPOSE:
;   Compute the cosign of the angle between two unit vectors on the
;   sphere.  This formula is from Jackson, pg 101.  (Or see my notes
;   of 6 Dec 92).  The angles must be in the following ranges:
;     0 <= phi < 360
;     0 <= theta <= 180
;   where theta=0 corresponds to the N pole, and theta=180 is the S pole.
;   If you want the dot product between RA and DEC coordinates,
;   pass the following arguments (in radians):
;     RA1, DEC1+90, RA2, DEC2+90
;
; CALLING SEQUENCE:
;   dotproduct_sphere( phi1, theta1, phi2, theta2, [/degrees, /hrdeg] )
;
; INPUTS:
;   phi1:       RA of first point(s) in radians
;   theta1:     DEC of first point(s) in radians
;   phi2:       RA of second point(s) in radians
;   theta2:     DEC of second point(s) in radians
;
; OPTIONAL INPUTS:
;   degrees:    If set, then all angles are in degrees
;   hrdeg:      If se, then RA angles in hours and DEC angles in degrees
;
; OUTPUTS:
;   cosgamma:   Cosine of the angle between the two positions
;-
;-----------------------------------------------------------------------
 
function dotproduct_sphere, phi1, theta1, phi2, theta2, $
 degrees=degrees, hrdeg=hrdeg
 
   ; Need 4 parameters
   if N_params() LT 4 then begin
      print, 'Syntax - dotproduct_sphere( phi1, theta1, phi2, theta2, [/degrees, /hrdeg] )'
      return, -1
   endif
 
   convRA = 1.d0
   convDEC = 1.d0
   if (keyword_set(degrees)) then begin
      convRA = !dpi / 180.d0
      convDEC = !dpi / 180.d0
   endif
   if (keyword_set(hrdeg)) then begin
      convRA = !dpi / 12.d0
      convDEC = !dpi / 180.d0
   endif

   cosgamma= sin(theta1*convDEC) * sin(theta2*convDEC) $
    * cos((phi1-phi2)*convRA) $
    + cos(theta1*convDEC) * cos(theta2*convDEC)
 
   return, cosgamma
end
;-----------------------------------------------------------------------
