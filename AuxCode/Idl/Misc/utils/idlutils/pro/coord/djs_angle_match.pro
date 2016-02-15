;------------------------------------------------------------------------------
;+
; NAME:
;   djs_angle_match
;
; PURPOSE:
;   Given two lists of coordinates on a sphere, find matches within an
;   angular distance.  For each entry in list A, find all the entries
;   in list B that lie within an angular distance dtheta.
;   Optionally output up to mmax of these matches per entry, giving
;   the index number of the matches in mindx, and the angular distance
;   in mdist.
;
;   If the lists A and B are different, then the total number of pairs
;   is given by total(mcount).
;   If the lists A and B are the same, then the total number of unique
;   pairs is given by (total(mcount) - N_elements(raA)) / 2.
;
;   This function loops over the objects in each list (sort of), so it's
;   not very fast.
;
; CALLING SEQUENCE:
;   ntot = djs_angle_match( raA, decA, [raB, decB,] dtheta=dtheta, $
;    [ mcount=mcount, mindx=mindx, mdist=mdist, mmax=mmax, units=units ]
;
; INPUTS:
;   raA:        RA of first point(s) in radians/degrees/hours
;   decA:       DEC of first point(s) in radians/degrees
;   dtheta:     Maximum angular distance for points to be considered matches
;
; OPTIONAL INPUTS:
;   raB:        RA of second point(s) in radians/degrees/hours
;   decB:       DEC of second point(s) in radians/degrees
;   mmax:       Maximum number of matches per point.  Default to 1.
;   units:      Set to
;                  degrees - All angles in degrees
;                  hrdeg - RA angles in hours, DEC angles and output in degrees
;                  radians - All angles in radians
;               Default to "degrees".
;
; OUTPUTS:
;   ntot:       Total number of points A with one or more matches in B
;
; OPTIONAL OUTPUTS:
;   mcount:     For each A, number of matches in B.  Vector of length A.
;   mindx:      For each A, indices of matches in B, sorted by their distance.
;               If mmax > 1, this array is of dimensions (mmax, A).
;               For each A, only the values (0:mcount-1,A) are relevant.
;               If mmax = 1, then the return value is a vector.
;               Any unused array elements are set to -1.
;   mdist:      For each A, distance to matches in B, sorted by their distance.
;               If mmax > 1, this array is of dimensions (mmax, A).
;               For each A, only the values (0:mcount-1,A) are relevant.
;               If mmax = 1, then the return value is a vector.
;               Any unused array elements are set to -1.
;
; COMMENTS:
;   By specifying only one set of coordinates (raA, decA), matches are found
;   within that list, but avoiding duplicate matches (i.e., matching 1 to 2
;   and then 2 to 1) and avoiding matching an object with itself (i.e.,
;   matching 1 to 1).  If you wish to include self-matches and duplicates,
;   then call with raB=raA and decB=decA.
;
; PROCEDURES CALLED:
;   djs_diff_angle()
;
; INTERNAL PROCEDURES:
;   djs_angle_1match()
;   djs_angle_2match()
;
; REVISION HISTORY:
;   27-May-1997  Written by David Schlegel, Durham
;   24-Feb-1999  Converted to IDL 5 (DJS)
;   05-Mar-1999  Made the internal routines for more efficient matching
;                within the same coordinate list without duplicates, e.g.
;                by only specifying raA, decA and not raB, decB.
;-
;------------------------------------------------------------------------------
function djs_angle_1match, raA, decA, dtheta=dtheta, $
 mcount=mcount, mindx=mindx, mdist=mdist, mmax=mmax, units=units

   if (NOT keyword_set(units)) then units="degrees"
   if (NOT keyword_set(mmax)) then mmax=1

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
   mcount = lonarr(numA)
   mindx = lonarr(mmax, numA) - 1
   mdist = dblarr(mmax, numA) - 1
   tempindx = lonarr(numA)
   tempdist = dblarr(numA)

   ; Sort points by DEC
;print, 'Sorting list A'
   indxA = sort(decA)

   iStart = 0L
   iEnd = 0L
;print, 'Looking for duplicates'
   for iA=0L, numA-2 do begin ; Don't loop over the last point

      ; Limit search to declination range within "dtheta"
      while ( decA[indxA[iStart]] LT decA[indxA[iA]] - dtheta $
       AND iStart LT numA-1 ) do iStart = iStart + 1
      while ( decA[indxA[iEnd]] LT decA[indxA[iA]] + dtheta $
       AND iEnd LT numA-1 ) do iEnd = iEnd + 1

      ; Avoid double-counting by always forcing iStart > iA
      iStart = iStart > iA + 1

if (iEnd GE iStart) then begin
      nmatch = 0

      maxdec = abs( decA[indxA[iA]] ) + dtheta < 90.
      if (maxdec GE highDEC) then dalpha = allRA + 1 $
       else dalpha = (convRA/convDEC) * dtheta / cos(maxdec*convDEC)

      iBvec = iStart + lindgen(iEnd-iStart+1)

      ; Select objects whose RA falls in the range of "dtheta" about point A
      ii = where( abs( raA[indxA[iA]] - raA[indxA[iBvec]] ) LT dalpha $
              OR  abs( raA[indxA[iA]] - raA[indxA[iBvec]] + allRA ) LT dalpha $
              OR  abs( raA[indxA[iA]] - raA[indxA[iBvec]] - allRA ) LT dalpha, $
              cti )

      if (cti GT 0) then begin
         adist = djs_diff_angle( raA[indxA[iA]], decA[indxA[iA]], $
          raA[indxA[iBvec[ii]]], decA[indxA[iBvec[ii]]], units=units )
         jj = where(adist LT dtheta, ctj)
         ; The following are matches in distances computed by djs_diff_angle.
         if (ctj GT 0) then begin
            tempindx[nmatch:nmatch+ctj-1] = iBvec[ii[jj]]
            tempdist[nmatch:nmatch+ctj-1] = adist[jj]
            nmatch = nmatch + ctj
         endif
      endif

      mcount[indxA[iA]] = min ( [mmax, nmatch] )
      if (nmatch GT 0) then begin
         ; Sort the matches, and keep only the mmax closest ones
         tempsort = sort ( tempdist[0:nmatch-1] )
         mindx[0:mcount[indxA[iA]]-1,indxA[iA]] = $
          indxA[ tempindx[ tempsort[0:mcount[indxA[iA]]-1] ] ]
         mdist[0:mcount[indxA[iA]]-1,indxA[iA]] = $
          tempdist[ tempsort[0:mcount[indxA[iA]]-1] ]
      endif
endif

   endfor

   if (mmax EQ 1) then begin
      mindx = transpose(mindx)
      mdist = transpose(mdist)
   endif

   junk = where(mcount GT 0, ntot)

   return, ntot
end
;------------------------------------------------------------------------------
function djs_angle_2match, raA, decA, raB, decB, dtheta=dtheta, $
 mcount=mcount, mindx=mindx, mdist=mdist, mmax=mmax, units=units

   if (NOT keyword_set(units)) then units="degrees"
   if (NOT keyword_set(mmax)) then mmax=1

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
   mcount = lonarr(numA)
   mindx = lonarr(mmax, numA) - 1 ; Set equal to -1 for no matches
   mdist = dblarr(mmax, numA) - 1 ; Set equal to -1 for no matches
   tempindx = lonarr(numB)
   tempdist = dblarr(numB)

   ; Sort points by DEC
;print, 'Sorting list A'
   indxA = sort(decA)
;print, 'Sorting list B'
   indxB = sort(decB)

   iStart = 0L
   iEnd = 0L
;print, 'Looking for duplicates'
   for iA=0L, numA-1 do begin

      ; Limit search to declination range within "dtheta"
      while ( decB[indxB[iStart]] LT decA[indxA[iA]] - dtheta $
       AND iStart LT numB-1 ) do iStart = iStart + 1
      while ( decB[indxB[iEnd]] LT decA[indxA[iA]] + dtheta $
       AND iEnd LT numB-1 ) do iEnd = iEnd + 1

      nmatch = 0

      maxdec = abs( decA[indxA[iA]] ) + dtheta < 90.
      if (maxdec GE highDEC) then dalpha = allRA + 1 $
       else dalpha = (convRA/convDEC) * dtheta / cos(maxdec*convDEC)

      iBvec = iStart + lindgen(iEnd-iStart+1)

      ; Select objects whose RA falls in the range of "dtheta" about point A
      ii = where( abs( raA[indxA[iA]] - raB[indxB[iBvec]] ) LT dalpha $
              OR  abs( raA[indxA[iA]] - raB[indxB[iBvec]] + allRA ) LT dalpha $
              OR  abs( raA[indxA[iA]] - raB[indxB[iBvec]] - allRA ) LT dalpha, $
              cti )

      if (cti GT 0) then begin
         adist = djs_diff_angle( raA[indxA[iA]], decA[indxA[iA]], $
          raB[indxB[iBvec[ii]]], decB[indxB[iBvec[ii]]], units=units )
         jj = where(adist LT dtheta, ctj)
         ; The following are matches in distances computed by djs_diff_angle.
         if (ctj GT 0) then begin
            tempindx[nmatch:nmatch+ctj-1] = iBvec[ii[jj]]
            tempdist[nmatch:nmatch+ctj-1] = adist[jj]
            nmatch = nmatch + ctj
         endif
      endif

      mcount[indxA[iA]] = min ( [mmax, nmatch] )
      if (nmatch GT 0) then begin
         ; Sort the matches, and keep only the mmax closest ones
         tempsort = sort ( tempdist[0:nmatch-1] )
         mindx[0:mcount[indxA[iA]]-1,indxA[iA]] = $
          indxB[ tempindx[ tempsort[0:mcount[indxA[iA]]-1] ] ]
         mdist[0:mcount[indxA[iA]]-1,indxA[iA]] = $
          tempdist[ tempsort[0:mcount[indxA[iA]]-1] ]
      endif

   endfor

   if (mmax EQ 1) then begin
      mindx = transpose(mindx)
      mdist = transpose(mdist)
   endif

   junk = where(mcount GT 0, ntot)

   return, ntot
end
;------------------------------------------------------------------------------
function djs_angle_match, raA, decA, raB, decB, dtheta=dtheta, $
 mcount=mcount, mindx=mindx, mdist=mdist, mmax=mmax, units=units

   ; Need 5 parameters
   if (N_params() EQ 2) then begin
      ; Call with same RA,DEC
      ntot = djs_angle_1match( raA, decA, dtheta=dtheta, $
       mcount=mcount, mindx=mindx, mdist=mdist, mmax=mmax, units=units)
   endif else if (N_params() EQ 4) then begin
      ; Call with different RA,DEC
      ntot = djs_angle_2match( raA, decA, raB, decB, dtheta=dtheta, $
       mcount=mcount, mindx=mindx, mdist=mdist, mmax=mmax, units=units)
   endif else begin
      print, 'Syntax - ntot = djs_angle_match( raA, decA, [raB, decB,] dtheta=dtheta, $'
      print, ' [ mcount=mcount, mindx=mindx, mdist=mdist, mmax=mmax, units=units ]'
      return, -1
   endelse

   return, ntot
end
;------------------------------------------------------------------------------
