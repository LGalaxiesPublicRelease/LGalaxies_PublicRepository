;------------------------------------------------------------------------------
;+
; NAME:
;   sortinplace
;
; PURPOSE:
;   Sort an array in-place, without making a second copy of the data
;
; CALLING SEQUENCE:
;   sortinplace, a, isort
;
; INPUTS:
;   a          - Array of any data type
;   isort      - Indices for sorting the A array, such as an index
;                list returned by the SORT() function
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   a          - (Modified)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine is particularly useful when sorting in-place a large
;   array of structures.  More simply, one would execute
;      a = a[isort]
;   but the above command makes a second copy of the array "a".
;
; EXAMPLES:
;
; BUGS:
;   The index list ISORT must be a unique list of indices into
;   the array A, with no duplicate entries and no out-of-bounds
;   values (such as -1).
;
; DATA FILES:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   29-May-2004  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro sortinplace, objs, isort

   nobj = n_elements(objs)
   if (n_elements(isort) NE nobj) then $
    message, 'Number of objects in OBJS and ISORT must agree'
   revsort = lonarr(nobj)
   revsort[isort] = lindgen(nobj)
   qdone = bytarr(nobj)

   for i=0L, nobj-1 do begin
      j = i
      tmp = objs[isort[j]]
      while (qdone[j] EQ 0) do begin
         tmp2 = objs[j]
         objs[j] = tmp
         qdone[j] = 1B
         j = revsort[j]
         tmp = tmp2
      endwhile
   endfor

   return
end
;------------------------------------------------------------------------------
