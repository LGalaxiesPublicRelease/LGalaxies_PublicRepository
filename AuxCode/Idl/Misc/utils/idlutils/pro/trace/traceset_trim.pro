;+
; NAME:
;   traceset_trim
;
; PURPOSE:
;   Trim a trace set to a selected number of traces
;
; CALLING SEQUENCE:
;   newset = traceset_trim(tset, [ indx ] )
;
; INPUTS:
;   tset       - Structure containing trace set
;
; OPTIONAL KEYWORDS:
;   indx       - Indices (0-indexed) of trace numbers to select; if not set,
;                then simply return the input structure
;
; OUTPUTS:
;   newset     - Trimmed trace set
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This function returns the input structure (which cannot be an array)
;   with exactly the same structure except for the tag name COEFF.
;   That particular value is trimmed to the values COEFF[*,INDX].
;   In our "trace set" notation, this returns a sub-set of the traces.
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   14-Jan-2004  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function traceset_trim, tset, indx

   ; Need 1 parameter
   if (N_params() LT 1) then begin
      print, 'Syntax - traceset_trim, tset, [ indx ]'
      return, 0
   endif

   if (n_elements(indx) EQ 0) then return, tset

   tags = tag_names(tset)
   for itag=0L, n_tags(tset)-1 do begin
      thisval = tset.(itag)
      if (tags[itag] EQ 'COEFF') then thisval = thisval[*,indx]
      if (itag EQ 0) then newset = create_struct(tags[0], thisval) $
       else newset = create_struct(newset, tags[itag], thisval)
   endfor

   return, newset
end
;------------------------------------------------------------------------------
