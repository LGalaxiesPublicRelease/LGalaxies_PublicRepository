;------------------------------------------------------------------------------
;+
; NAME:
;   djs_angle_nmatch
;
; PURPOSE:
;   Given two lists of coordinates on a sphere, find matches within an
;   angular distance.  For each entry in list B, return the number of
;   matches in list A that lie within an angular distance dtheta.
;
;   This function loops through list A, so it is very slow if that list is long.
;
;   The angle dtheta can be the same for each object in A, or may be set
;   to a vector of length A.
;
;   A list of indices where B has a match in A is where(nmatch GT 0).
;
; CALLING SEQUENCE:
;   nmatch = djs_angle_nmatch( raA, decA, raB, decB, dtheta, $
;    [ units=units ]
;
; INPUTS:
;   raA:        RA of first point(s) in radians/degrees/hours
;   decA:       DEC of first point(s) in radians/degrees
;   raB:        RA of second point(s) in radians/degrees/hours
;   decB:       DEC of second point(s) in radians/degrees
;   dtheta:     Maximum angular distance for points to be considered matches
;
; OPTIONAL INPUTS:
;   units:      Set to
;                  degrees - All angles in degrees
;                  hrdeg - RA angles in hours, DEC angles and output in degrees
;                  radians - All angles in radians
;               Default to "degrees".
;
; OUTPUTS:
;   nmatch:     For each B, number of matches with A
;
; PROCEDURES CALLED:
;   djs_diff_angle()
;
; REVISION HISTORY:
;   18-Jul-1997  Written by David Schlegel, Durham
;                Modified from djs_angle_match().
;   24-Feb-1999  Converted to IDL 5 (DJS)
;-
;------------------------------------------------------------------------------
function djs_angle_nmatch, raA, decA, raB, decB, dtheta, units=units

   ; Need 5 parameters
   if N_params() LT 5 then begin
      print, 'Syntax - ntot = djs_angle_nmatch( raA, decA, raB, decB, dtheta, $'
      print, ' [ units=units ]'
      return, -1
   endif

   if (NOT keyword_set(units)) then units="degrees"

   case units of
      "hrdeg" : begin
         convRA = !dpi / 12.d0
         convDEC = !dpi / 180.d0
         allRA = 24.d0
         highDEC = 90.d0
      end
      "radians" : begin
         convRA = 1.d0
         convDEC = 1.d0
         allRA = 2. * !dpi
         highDEC = 0.5 * !dpi
      end
      else : begin
         convRA = !dpi / 180.d0
         convDEC = !dpi / 180.d0
         allRA = 360.d0
         highDEC = 90.d0
      end
   endcase

   ; Allocate arrays
   numA = N_elements(raA)
   numB = N_elements(raB)
   nmatch = intarr(numB)

   if (N_elements(dtheta) EQ numA) then vtheta = dtheta $
    else vtheta = fltarr(numA) + dtheta

print, 'Looking for matches'
   for iA=0L, numA-1 do begin

      maxdec = abs( decA[iA] ) + vtheta[iA]
      if ( maxdec GE highDEC ) then dalpha = allRA + 1 $
       else dalpha = (convRA/convDEC) * vtheta[iA] / cos(maxdec*convDEC)

      ixbox = where( ( abs(decA[iA]-decB) LE vtheta[iA] ) AND $
       ( abs(raA[iA]-raB      ) LE dalpha OR $
         abs(raA[iA]-raB+allRA) LE dalpha OR $
         abs(raA[iA]-raB-allRA) LE dalpha ) )

      if (ixbox[0] NE -1) then begin
         adist = djs_diff_angle( raB[ixbox], decB[ixbox], raA[iA], decA[iA], $
          units=units )
         ii = where( adist LT vtheta[iA] )
         if (ii[0] NE -1) then begin
            inear = ixbox[ii]
            nmatch[inear] = nmatch[inear] + 1
         endif
      endif

   endfor

   return, nmatch
end
;------------------------------------------------------------------------------
