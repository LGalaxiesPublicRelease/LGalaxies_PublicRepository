;+
; NAME:
;   struct_append
;
; PURPOSE:
;   Append more array elements to a structure.
;
; CALLING SEQUENCE:
;   outstruct = struct_append( struct1, struct2, [ /force ] )
;
; INPUTS:
;   struct1    - First structure; the output structure will match the tags
;                in this, and match the name if it's a named structure.
;   struct2    - Second structure to append to the first.
;
; OPTIONAL INPUTS:
;   force      - If set, then append these two structures, even if one or
;                both are arrays of zeros rather than actual structures.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   If either structure is undefined are set to zero, then return the other
;   one only.  If one of the structures is defined (and non-zero), but not
;   a structure, then a blank structure is put in its place.
;
; EXAMPLES:
;   > a={one:1,two:2}
;   > b={one:11,three:33}
;   > print,struct_append(a,b)
;     {       1       2}{      11       0}
;
; BUGS:
;
; PROCEDURES CALLED:
;   headfits()
;
; REVISION HISTORY:
;   26-Jun-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function struct_append, struct1, struct2, force=force

   num1 = n_elements(struct1)
   num2 = n_elements(struct2)
   if (num1 EQ 0 AND num2 EQ 0) then return, 0
   if (num1 EQ 0) then return, struct2
   if (num2 EQ 0) then return, struct1

   if (NOT keyword_set(force)) then begin
      if (keyword_set(struct1) EQ 0 AND keyword_set(struct2) EQ 0) then $
       return, 0
      if (keyword_set(struct1) EQ 0) then return, struct2
      if (keyword_set(struct2) EQ 0) then return, struct1
   endif

   ; In the event that STRUCT1 is not a structure, but STRUCT2 is,
   ; then use STRUCT2 as the template for the output structure.
   obj1 = struct1[0]
   if (size(struct1, /tname) NE 'STRUCT') then begin
      if (size(struct2, /tname) EQ 'STRUCT') then begin
         obj1 = struct2[0]
      endif else begin
         outstruct = [struct1, struct2]
         return, outstruct
      endelse
   endif
   struct_assign, {junk:0}, obj1 ; Zero-out this new structure
   outstruct = replicate(obj1, num1+num2)

   if (size(struct1, /tname) EQ 'STRUCT') then begin
      outstruct[0:num1-1] = struct1[*]
   endif

   for irow=0L, num2-1 do begin
      if (size(struct2, /tname) EQ 'STRUCT') then begin
         struct_assign, struct2[irow], obj1
         outstruct[num1+irow] = obj1
      endif
   endfor

   return, outstruct
end
;------------------------------------------------------------------------------
