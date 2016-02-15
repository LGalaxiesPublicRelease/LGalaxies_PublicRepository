;+
; NAME:
;   findbkpt
;
; PURPOSE:
;   Choose bkpts for b-spline given different constraints
;
; CALLING SEQUENCE:
;   
;   fullbkpt = findbkpt(x, good, bkpt, nord, bkspace=bkspace, nbkpts=nbkpts, $
;                   everyn=everyn, silent=silent)
;
; INPUTS:
;   bkpt       - Breakpoint vector returned by efc
;
; RETURNS:
;   fullbkpt   - The fullbkpt vector required by evaluations with bvalu
;
; OPTIONAL KEYWORDS:
;   bkspace    - Spacing of breakpoints in units of x
;   everyn     - Spacing of breakpoints in good pixels
;   nbkpts     - Number of breakpoints to span x range
;                 minimum is 2 (the endpoints)
;   silent     - Do not produce non-critical messages
;
; OPTIONAL OUTPUTS:
;   bkpt       - breakpoints without padding
;
; COMMENTS:
;   If both bkspace and nbkpts are passed, bkspace is used.
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   none
;
; REVISION HISTORY:
;   10-Mar-2000  Written by Scott Burles, FNAL
;-
;------------------------------------------------------------------------------
function findbkpt, x, good, bkpt, nord, bkspace=bkspace, nbkpts=nbkpts, $
                   everyn=everyn, silent=silent, bkptspread=bkptspread

      ngood = n_elements(good)

      if (NOT keyword_set(bkpt)) then begin
 
         range = (max(x) - min(x))
         startx = min(x)
         if (keyword_set(bkspace)) then begin
            nbkpts = long(range/float(bkspace)) + 1
            if (nbkpts LT 2) then nbkpts = 2
            tempbkspace = double(range/(float(nbkpts-1)))
            bkpt = (findgen(nbkpts))*tempbkspace + startx
         endif else if keyword_set(nbkpts) then begin
            nbkpts = long(nbkpts)
            if (nbkpts LT 2) then nbkpts = 2
            tempbkspace = double(range/(float(nbkpts-1)))
            bkpt = (findgen(nbkpts))*tempbkspace + startx
         endif else if keyword_set(everyn) then begin
            nbkpts = ngood / everyn
            xspot = lindgen(nbkpts)*ngood/(nbkpts-1)
            bkpt = x[good[xspot]]
         endif else message, 'No information for bkpts'
      endif

      bkpt = float(bkpt)

      if (min(x) LT min(bkpt,spot)) then begin
         if (NOT keyword_set(silent)) then $
          print, 'Lowest breakpoint does not cover lowest x value: changing'
         bkpt[spot] = min(x)
      endif

      if (max(x) GT max(bkpt,spot)) then begin
         if (NOT keyword_set(silent)) then $
          print, 'highest breakpoint does not cover highest x value, changing'
         bkpt[spot] = max(x)
      endif

      nshortbkpt = n_elements(bkpt)
      fullbkpt = bkpt

      if (NOT keyword_set(bkptspread)) then bkptspread = 0.1
      bkptspace = (bkpt[1] - bkpt[0]) * bkptspread

      for i=1, nord-1 do $
       fullbkpt = [bkpt[0]-bkptspace*i, fullbkpt, $
        bkpt[nshortbkpt - 1] + bkptspace*i]
 
   return, fullbkpt
end 
