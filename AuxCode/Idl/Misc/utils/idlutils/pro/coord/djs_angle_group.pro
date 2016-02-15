;------------------------------------------------------------------------------
;+
; NAME:
;   djs_angle_group
; PURPOSE:
;   Group objects using their coordinates on the sphere.
;
;   Any coordinates within dtheta of one another are put in the same group.
;   Note that if there is a string of successive objects on the sky, each
;   separated by less than dtheta, then all of these objects are assigned
;   to the same group.  This is incorrect in the sense that the first and
;   last objects in the string may have a large separation; however, this
;   is the only unambigious answer to the problem.
;
; CALLING SEQUENCE:
;   ngroup = djs_angle_group( ra, dec, dtheta, $
;    [gstart=gstart, gcount=gcount, gindx=gindx, units=units] )
;
; INPUTS:
;   ra:         RA of point(s) in radians/degrees/hours
;   dec:        DEC of point(s) in radians/degrees
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
;   ngroup:     Total number of groups.  If no matches are found, then this
;               equals the number of objects.
;
; OPTIONAL OUTPUTS:
;   gstart:     Vector of length "ngroup" with the starting index of each group.
;   gcount:     Vector of length "ngroup" with the number of objects in each
;               group.
;   gindx:      Indices of objects in each group.  The i-th group will have
;               its object indices stored in gindx(gstart:gstart+gcount-1).
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT PROCEDURES:
;   djs_search_around
;
; REVISION HISTORY:
;   29-May-1997  Written by D. Schlegel, Durham
;   24-Feb-1999  Converted to IDL 5 (DJS).
;-
;------------------------------------------------------------------------------
; Find all neighbors to object #IOBJ and set their flag in QFLAG.
; This procedure then recursively looks for neighbors of those new neighbors.

pro djs_search_around, iobj, nobj, dtheta

   ; This common block is used for variables passed from the main procedure.
   common djs_group_com1, raSort, decSort, qflag

   ; This common block is used to store variables within this procedure
   ; between procedure calls.
   common djs_group_com2, jStart, jEnd, convRA, convDEC, allRA, highDEC

   ; Keep track of index numbers JSTART and JEND that span the
   ; possible declination range for neighbors of object IOBJ.
   ; Note that the assumption is made that JSTART and JEND are always
   ; increasing, and never decreasing.  This is valid since the points
   ; are first sorted in declination.  However, JSTART is set in the
   ; main procedure based upon the first object in the group.
   while ( decSort[jEnd] LT decSort[iobj] + dtheta $
    AND jEnd LT nobj-1 ) do jEnd = jEnd + 1

   ; Keep track of the RA range, DALPHA, for possible neighbors of object IOBJ.
   maxdec = abs( decSort[iobj] ) + dtheta
   if ( maxdec GE highDEC ) then dalpha = allRA + 1 $
    else dalpha = (convRA/convDEC) * dtheta / cos(maxdec*convDEC)

   for jobj=jStart, jEnd do begin
      if (qflag[jobj] EQ 0) then begin  ; This object not yet grouped

         ; See if object numbers IOBJ and JOBJ are neighbors
         if ( abs( raSort[iobj] - raSort[jobj] ) LT dalpha $
          OR  abs( raSort[iobj] - raSort[jobj] + allRA ) LT dalpha $
          OR  abs( raSort[iobj] - raSort[jobj] - allRA ) LT dalpha ) $
          then begin

            adist = djs_diff_angle( raSort[iobj], decSort[iobj], $
             raSort[jobj], decSort[jobj], units=units )
            if ( adist LT dtheta ) then begin
               qflag[jobj] = 1
               djs_search_around, jobj, nobj, dtheta
            endif
         endif

      endif
   endfor

   return
end
;------------------------------------------------------------------------------
function djs_angle_group, ra, dec, dtheta, $
 gstart=gstart, gcount=gcount, gindx=gindx, units=units

   common djs_group_com1, raSort, decSort, qflag
   common djs_group_com2, jStart, jEnd, convRA, convDEC, allRA, highDEC

   ; Need 3 parameters
   if N_params() LT 3 then begin
      print, 'Syntax - iuniq = djs_angle_group( ra, dec, dtheta, [units=units]'
      print, 'Syntax - ngroup = djs_angle_group( ra, dec, dtheta, $'
      print, ' [gstart=gstart, gcount=gcount, gindx=gindx, units=units] )'
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

   nobj = N_elements(ra)
   qflag = intarr(nobj)
   ngroup = 0L
   glast = 0L
   gstart = lonarr(nobj)
   gcount = lonarr(nobj)
   gindx = lonarr(nobj)
   jStart = 0L
   jEnd = 0L

   ; Sort points by DEC
   sindx = sort(dec)
   raSort = ra[sindx]
   decSort = dec[sindx]

   for iobj=0L, nobj-1 do begin
      if (qflag[iobj] EQ 0) then begin  ; Object #iobj not grouped yet

         while ( decSort[jStart] LT decSort[iobj] - dtheta $
          AND jStart LT nobj-1 ) do jStart = jStart + 1

         tflag = qflag
         djs_search_around, iobj, nobj, dtheta
         tflag = qflag - tflag  ; = 1 for objects in this group

         ; listgrp = (unsorted) indices for objects in this group
         listgrp = sindx[ where(tflag EQ 1) ]

         ; Save this group
         gstart[ngroup] = glast
         gcount[ngroup] = N_elements(listgrp)
         gindx[glast:glast+gcount[ngroup]-1] = listgrp
         glast = glast + gcount[ngroup]
         ngroup = ngroup + 1
      endif
   endfor

   ; Trim the group lists to have only as many entries as there are groups
   gstart = gstart[0:ngroup-1]
   gcount = gcount[0:ngroup-1]

   return, ngroup
end
;------------------------------------------------------------------------------
