;+
; NAME:
;   spheregroup
;
; PURPOSE:
;   Perform friends-of-friends grouping given ra/dec type coords 
;
; CALLING SEQUENCE:
;   ingroup = spheregroup( ra, dec, linklength, [chunksize=], $
;     [multgroup=], [firstgroup=], [nextgroup=] )
;
; INPUTS:
;   ra         - ra coordinates in degrees (N-dimensional array)
;   dec        - dec coordinates in degrees (N-dimensional array)
;   linklength - linking length for groups (degrees)
;
; OPTIONAL INPUTS:
;   chunksize  - the algorithm breaks the sphere up into a bunch of
;                regions with a characteristic size chunksize
;                (degrees). By default this is max(1.,4*linklength)
;
; OUTPUTS:
;   ingroup    - group number of each object (N-dimensional array);
;                -1 if no groups
;
; OPTIONAL INPUT/OUTPUTS:
;   multgroup  - multiplicity of each group 
;   firstgroup - first member of each group 
;   nextgroup  - index of next member of group for each object
;
; COMMENTS:
;   The code breaks the survey region into chunks which overlap by
;   about linklength. Friends-of-friends is run on each chunk
;   separately. Finally, the bookkeeping is done to combine the
;   results (i.e. joining groups across chunk boundaries). This should
;   scale as area*density^2, 
;
;   It is important that chunksize is >=4.*linklength, and this is
;   enforced.
;
;   firstgroup and nextgroup form a primitive "linked list", which 
;   can be used to step through a particular group, as in the example 
;   below.
;
; EXAMPLES:
;   Group a set of points on a scale of 55'', then step through
;   members of the third group:
;
;   > ingroup=spheregroup(ra,dec,.0152778,multgroup=mult, $
;   > firstgroup=first, nextgroup=next)
;   > indx=firstgroup[2]
;   > for i = 0, multgroup[2] do begin & $
;   > print,ra[indx],dec[indx] & $
;   > indx=nextgroup[indx]  & $
;   > end
;
;   Of course, you could just "print,ra[where(ingroup eq 2)]", but I
;   wanted to demostrate how the linked list worked.
;
; BUGS:
;   Behavior at poles not well tested.
;
; PROCEDURES CALLED:
;   Dynamic link to spheregroup.c
;
; REVISION HISTORY:
;   19-Jul-2001  Written by Mike Blanton, Fermiland
;-
;------------------------------------------------------------------------------
function spheregroup, ra, dec, linklength, chunksize=chunksize, multgroup=multgroup, firstgroup=firstgroup, nextgroup=nextgroup

   ; Need at least 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - ingroup = spheregroup( ra, dec, linklength, [chunksize=], $'
      print, ' [multgroup=], [firstgroup=], [nextgroup=] )' 
      return, -1 
  endif

   if (NOT keyword_set(chunksize)) then begin
       chunksize=max([4.*linklength,1.])
   end else begin
       if (chunksize lt 4.*linklength) then begin
           chunksize=4.*linklength
           print,'chunksize changed to ',chunksize
       endif
   endelse 

   npoints = N_elements(ra)
   if (npoints le 0) then begin
       print, 'Need array with > 0 elements'
       return, -1
   endif

   if (linklength le 0) then begin
       print, 'Need linklength > 0'
       return, -1
   endif
   
   ; Allocate memory for the ingroups array
   ingroup=lonarr(npoints)

   ; Call grouping software
   soname = filepath('libspheregroup.'+idlutils_so_ext(), $
    root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')
   retval = call_external(soname, 'spheregroup', long(npoints), double(ra), $
                          double(dec), double(linklength), double(chunksize), $
                          ingroup)
   
   ; Make multiplicity, etc.
   multgroup=lonarr(npoints)
   nextgroup=lonarr(npoints)-1L
   firstgroup=lonarr(npoints)-1L
   for i = 0l, npoints-1l do multgroup[ingroup[i]]=multgroup[ingroup[i]]+1l
   for i = npoints-1l, 0l, -1l do begin 
       nextgroup[i]=firstgroup[ingroup[i]]
       firstgroup[ingroup[i]]=i
   end
   
   return, ingroup
end
;------------------------------------------------------------------------------
