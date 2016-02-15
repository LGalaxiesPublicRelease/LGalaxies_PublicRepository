;-----------------------------------------------------------------------
;+
; NAME:
;   hdr2struct
;
; PURPOSE:
;   Convert a string array with the format of a FITS header 
;       into a single structure
;
; CALLING SEQUENCE:
;   struct = hdr2struct(hdr)
;
; INPUTS:
;   hdr:        FITS-(like) header (array of strings)
;
; OUTPUTS:
;   struct:     Single structure with keyword(s) and comment(s)
;
; OPTIONAL INPUTS:
;
; COMMENTS:
;   Only passes back anonymous structures as it appends as it works
;    through all the header keywords.  Use struct_append to make an
;    array from an array of headers.
;
; PROCEDURES CALLED:
;   headfits()
;
; REVISION HISTORY:
;   17-Jul-2001 Written by Burles
;-
;-----------------------------------------------------------------------
function hdr2struct, hdr

   s = size(hdr)
   nkeys = s[1] 
   tags = strmid(hdr,0,8)
   tags = strtrim(tags,2)

   duplicate = strarr(nkeys)
   order = sort(tags)
   for j=1,nkeys - 1 do begin
     i   = order[j]
     im1 = order[j-1]
     if tags[i] EQ tags[im1] then duplicate[i] = strtrim(duplicate[im1] + 1,2)
   endfor

   ; Now let's fill structure

   i = 0

   value = sxpar(hdr, tags[i], comment = comments)
   this_struct = create_struct(tags[i]+duplicate[i], value[0])
   if comments NE '' then $
       this_struct = create_struct(this_struct, $
         tags[i]+duplicate[i]+'COMMENT', comments[0]) 

   for i=1,nkeys - 1 do begin
     value = sxpar(hdr, tags[i], comment = comments, count=count)
     this_struct = create_struct(this_struct, tags[i]+duplicate[i], value[0])
     if comments NE '' then this_struct = create_struct(this_struct, $
       tags[i]+duplicate[i]+' COMMENT', comments[0])
   endfor

   return, this_struct
end

