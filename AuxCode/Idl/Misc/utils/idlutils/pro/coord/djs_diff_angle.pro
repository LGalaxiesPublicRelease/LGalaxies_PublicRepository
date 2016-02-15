;+
; NAME:
;   djs_diff_angle
;
; PURPOSE:
;   Compute the angular distance between two points on a sphere.
;
; CALLING SEQUENCE:
;   adist = djs_diff_angle( ra, dec, ra0, dec0, [ units=units ] )
;
; INPUTS:
;   ra1:        RA of first point(s) in radians/degrees/hours
;   dec1:       DEC of first point(s) in radians/degrees
;   ra2:        RA of second point(s) in radians/degrees/hours
;   dec2:       DEC of second point(s) in radians/degrees
;
; OPTIONAL INPUTS:
;   units:      Set to
;                  degrees - All angles in degrees
;                  hrdeg - RA angles in hours, DEC angles and output in degrees
;                  radians - All angles in radians
;               Default to "degrees".
;
; OUTPUTS:
;   adist:      Angular distance(s) in radians if UNITS is set to 'radians',
;               or in degrees otherwise
;
; COMMENTS:
;   Note that either (ra1,dec1) or (rap,decp) must be scalars or 1-element
;   arrays, or all must be arrays of the same length.
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   14-May-1997  Written by D. Schlegel, Durham
;-
;------------------------------------------------------------------------------
function djs_diff_angle, ra1, dec1, ra2, dec2, units=units1

   DPIBY2 = 0.5d0 * !dpi

   ; Need 4 parameters
   if N_params() LT 4 then begin
      print, 'Syntax - adist = djs_diff_angle(ra1, dec1, ra2, dec2, [units= ] )'
      return, -1
   endif

   num1 = n_elements(ra1)
   num2 = n_elements(ra2)
   if (num1 NE n_elements(dec1) OR num2 NE n_elements(dec2)) then $
    message, 'Dimensions of inputs are incompatible'
   if (num1 NE 1 AND num2 NE 1 AND num1 NE num2) then $
    message, 'Dimensions of inputs are incompatible'

   if (keyword_set(units1)) then units = units1 $
    else units = 'degrees'

   case units of
      'hrdeg' : begin
         convRA = !dpi / 12.d0
         convDEC = !dpi / 180.d0
      end
      'radians' : begin
         convRA = 1.d0
         convDEC = 1.d0
      end
      'degrees' : begin
         convRA = !dpi / 180.d0
         convDEC = !dpi / 180.d0
      end
      else : message, 'Unknown UNITS='+string(units)
   endcase

   ; The following allows the inputs to be 1-element arrays rather than
   ; scalars, by recasting those 1-element arrays as scalars.
   if (num1 EQ 1) then begin
      theta1 = dec1[0] * convDEC + DPIBY2
      theta2 = dec2 * convDEC + DPIBY2
      cosgamma= sin(theta1) * sin(theta2) $
       * cos((ra1[0] - ra2) * convRA)  + cos(theta1) * cos(theta2)
   endif else if (num2 EQ 1) then begin
      theta1 = dec1 * convDEC + DPIBY2
      theta2 = dec2[0] * convDEC + DPIBY2
      cosgamma= sin(theta1) * sin(theta2) $
       * cos((ra1 - ra2[0]) * convRA)  + cos(theta1) * cos(theta2)
   endif else begin
      theta1 = dec1 * convDEC + DPIBY2
      theta2 = dec2 * convDEC + DPIBY2
      cosgamma= sin(theta1) * sin(theta2) $
       * cos((ra1 - ra2) * convRA)  + cos(theta1) * cos(theta2)
   endelse

   return, acos(cosgamma < 1.d0) / convDEC
end
;------------------------------------------------------------------------------
