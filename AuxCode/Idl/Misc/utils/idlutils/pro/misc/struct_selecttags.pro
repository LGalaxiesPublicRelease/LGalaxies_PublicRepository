;+
; NAME:
;   struct_selecttags
;
; PURPOSE:
;   Trim a structure to a specified list of tags using SELECT and EXCEPT
;
; CALLING SEQUENCE:
;   outstruct = struct_selecttags(tdat, [ select_tags=, except_tags= ])
;
; INPUTS:
;   tdat       - Input structure, which can be an array
;
; OPTIONAL INPUTS:
;   except_tags- Option string array of column names to not return.
;                This can include wildcards.
;                Set to ' ' to not exclude any tags.
;   select_tags- Option string array of column names to return, which
;                takes priority over EXCEPT_TAGS.
;                This can include wildcards.
;
; OUTPUTS:
;   outstruct  - Ouput structure array
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   11-Nov-2003  Written by D. Schlegel, Princeton, broken out of SDSS_READOBJ
;------------------------------------------------------------------------------
function struct_selecttags, tdat1, select_tags=select_tags, $
 except_tags=except_tags

   if (keyword_set(except_tags) OR keyword_set(select_tags)) then begin
      tags = tag_names(tdat1)
      if (keyword_set(except_tags)) then begin
         qkeep = bytarr(n_elements(tags)) + 1b
         for i=0L, n_elements(except_tags)-1 do $
          qkeep = qkeep AND (strmatch(tags,strupcase(except_tags[i])) EQ 0)
      endif
      if (keyword_set(select_tags)) then begin
         qkeep = bytarr(n_elements(tags))
         for i=0L, n_elements(select_tags)-1 do $
          qkeep = qkeep OR (strmatch(tags,strupcase(select_tags[i])) EQ 1)
      endif
      ikeep = where(qkeep, nkeep)
      if (nkeep EQ 0) then begin
         splog, 'No tags selected'
         if keyword_set(outfile) then mwrfits, 0, outfile, /create
         return, 0
      endif
      tdat2 = create_struct(tags[ikeep[0]], tdat1[0].(ikeep[0]))
      for i=1L, nkeep-1 do $
       tdat2 = create_struct(tdat2, tags[ikeep[i]], tdat1[0].(ikeep[i]))

      ; If the input structure is an array, then fill it...
      if (n_elements(tdat1) GT 1) then begin
         tdat2 = replicate(tdat2, n_elements(tdat1))
         struct_assign, tdat1, tdat2
      endif
   endif else begin
      tdat2 = tdat1
   endelse

   return, tdat2
end
;------------------------------------------------------------------------------
