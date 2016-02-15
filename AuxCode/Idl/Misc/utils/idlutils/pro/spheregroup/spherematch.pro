;+
; NAME:
;   spherematch
;
; PURPOSE:
;   Take two sets of ra/dec coords and efficiently match them. It 
;   returns all matches between the sets of objects, and the distances
;   between the objects. The matches are returned sorted by increasing
;   distance. A parameter "maxmatch" can be set to the number of matches 
;   returned for each object in either list. Thus, maxmatch=1 (the default)
;   returns the closest possible set of matches. maxmatch=0 means to
;   return all matches
;
; CALLING SEQUENCE:
;   spherematch, ra1, dec1, ra2, dec2, matchlength, match1, match2, $
;     		       distance12, [maxmatch=maxmatch]
;
; INPUTS:
;   ra1         - ra coordinates in degrees (N-dimensional array)
;   dec1        - dec coordinates in degrees (N-dimensional array)
;   ra2         - ra coordinates in degrees (N-dimensional array)
;   dec2        - dec coordinates in degrees (N-dimensional array)
;   matchlength - distance which defines a match (degrees)
;
; OPTIONAL INPUTS:
;   maxmatch    - Return only maxmatch matches for each object, at
;                 most. Defaults to maxmatch=1 (only the closest
;                 match for each object). maxmatch=0 returns all
;                 matches.
;   estnmatch   - Estimate of the TOTAL number of matches.  If this is 
;                 absent or wrong, the C code is called twice,
;                 doubling execution time!
;
; OUTPUTS:
;   match1     - List of indices of matches in list 1; -1 if no matches
;   match2     - List of indices of matches in list 2; -1 if no matches
;   distance12 - Distance of matches; 0 if no matches
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The code breaks the survey region into chunks of size
;   4*matchlength. Matches are then performed by considering only 
;   objects in a given chunk and neighboring chunks. This makes the
;   code fast.
;
;   The matches are returned sorted by distance.
;
;   If you have a big list and a small list, call with the 
;   BIG LIST FIRST!!!
;   i.e.
;   
;   spherematch, BIGra, BIGdec, SMALLra, SMALLdec, matchlength, $
;                      matchBIG, matchSMALL, distance12
;
;   This method is inherently asymmetric.  Calling in this order will 
;   exploit the asymmetry to reduce memory usage and execution time. 
;
; EXAMPLES:
;
; BUGS:
;   Behavior at poles not well tested.
;
; PROCEDURES CALLED:
;   Dynamic link to spherematch.c
;
; REVISION HISTORY:
;   20-Jul-2001  Written by Mike Blanton, Fermiland
;   01-Mar-2006  estnmatch keyword added - D. Finkbeiner, Princeton
;          estnmatch allows the caller to estimate the number of 
;          matches, so the wrapper can allocate memory for results before
;          calling the C code.  If the estimate is absent or wrong,
;          the code is called a second time (as before). 
;-
;------------------------------------------------------------------------------
pro spherematch, ra1, dec1, ra2, dec2, matchlength, match1, match2, $
                 distance12, maxmatch=maxmatch, chunksize=chunksize, $
                 estnmatch=estnmatch

   ; Need at least 3 parameters
   if (N_params() LT 7) then begin
      print, 'Syntax - spherematch, ra1, dec1, ra2, dec2, matchlength, match1, match2, $'
      print, ' distance12, [maxmatch=]'
      return
  endif

   if (n_elements(maxmatch) eq 0) then begin
       maxmatch=1l
   end else begin
       if (maxmatch lt 0l) then begin
           print,'illegal maxmatch value: '+maxmatch 
           return
       endif
   endelse

   ; Set default return values
   match1 = -1
   match2 = -1
   distance12 = 0

   if(NOT keyword_set(chunksize)) then chunksize=max([4.*matchlength,0.1])

   npoints1 = N_elements(ra1)
   if (npoints1 le 0l) then begin
       print, 'Need array with > 0 elements'
       return
   endif
   if (npoints1 ne N_elements(dec1)) then begin
       print, 'ra1 and dec1 must have same length'
       return
   endif
   if (N_elements(ra2) le 0l) then begin
       print, 'Need array with > 0 elements'
       return
   endif
   npoints2=N_elements(ra2)
   if (npoints2 ne N_elements(dec2)) then begin
       print, 'ra2 and dec2 must have same length'
       return
   endif

   if (matchlength le 0l) then begin
       print, 'Need matchlength > 0'
       return
   endif

   soname = filepath('libspheregroup.'+idlutils_so_ext(), $
    root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')
   onmatch = keyword_set(estnmatch) ? long(estnmatch) : 0L
   onmatch_save = onmatch

; -------- First pass on matching C code. 
   omatch1     = lonarr(onmatch > 1)
   omatch2     = lonarr(onmatch > 1)
   odistance12 = dblarr(onmatch > 1)
   retval = call_external(soname, 'spherematch', $
                          long(npoints1), double(ra1), double(dec1), $
                          long(npoints2), double(ra2), double(dec2), $
                          double(matchlength), double(chunksize), $
                          long(omatch1),long(omatch2), $
                          double(odistance12), long(onmatch))

   if onmatch EQ 0 then begin 
      return
   endif

; -------- Second pass, if we did not allocate enough space before
   if onmatch GT onmatch_save then begin

      omatch1=lonarr(onmatch)
      omatch2=lonarr(onmatch)
      odistance12=dblarr(onmatch)
      retval = call_external(soname, 'spherematch', $
                             long(npoints1), double(ra1), double(dec1), $
                             long(npoints2), double(ra2), double(dec2), $
                             double(matchlength), double(chunksize), $
                             long(omatch1),long(omatch2), $
                             double(odistance12), long(onmatch))
   endif

; -------- trim padding in output arrays
   if onmatch lt n_elements(omatch1) then begin 
      omatch1 = omatch1[0:onmatch-1]
      omatch2 = omatch2[0:onmatch-1]
      odistance12 = odistance12[0:onmatch-1]
   endif 

   ; Retain only desired matches
   sorted=sort(odistance12)
   if (maxmatch gt 0l) then begin
       gotten1=lonarr(npoints1)
       gotten2=lonarr(npoints2)
       nmatch=0l
       for i = 0l, onmatch-1l do begin
           if ((gotten1[omatch1[sorted[i]]] lt maxmatch) and $
               (gotten2[omatch2[sorted[i]]] lt maxmatch)) then begin
               gotten1[omatch1[sorted[i]]]=gotten1[omatch1[sorted[i]]]+1l
               gotten2[omatch2[sorted[i]]]=gotten2[omatch2[sorted[i]]]+1l
               nmatch=nmatch+1l
           end
       end
       gotten1=lonarr(npoints1)
       gotten2=lonarr(npoints2)
       match1=lonarr(nmatch)
       match2=lonarr(nmatch)
       distance12=dblarr(nmatch)
       nmatch=0l
       for i = 0l, onmatch-1l do begin
           if ((gotten1[omatch1[sorted[i]]] lt maxmatch) and $
               (gotten2[omatch2[sorted[i]]] lt maxmatch)) then begin
               gotten1[omatch1[sorted[i]]]=gotten1[omatch1[sorted[i]]]+1l
               gotten2[omatch2[sorted[i]]]=gotten2[omatch2[sorted[i]]]+1l
               match1[nmatch]=omatch1[sorted[i]]
               match2[nmatch]=omatch2[sorted[i]]
               distance12[nmatch]=odistance12[sorted[i]]
               nmatch=nmatch+1l
           end
       end
   endif else begin
       nmatch=onmatch[sorted]
       onmatch=0
       match1=omatch1[sorted]
       omatch1=0
       match2=omatch2[sorted]
       omatch2=0
       if(arg_present(distance12)) then distance12=odistance12[sorted]
       odistance12=0
   endelse

   return
end
;------------------------------------------------------------------------------
