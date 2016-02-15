;+
; NAME:
;   sxcombinepar
;
; PURPOSE:
;   Combine values of specified header cards from many FITS headers.
;
; CALLING SEQUENCE:
;   sxcombinepar, hdrarr, cardname, outhdr, [ func=, weights=, /zeros, $
;    outcard=, _EXTRA=KeywordsForSxaddpar ]
;
; INPUTS:
;   hdrarr     - Array of pointers to FITS headers
;   cardname   - Name(s) of header cards to average
;   outhdr     - Output header
;
; OPTIONAL KEYWORDS:
;   func       - Function to apply:
;                  'average'
;                  'median'
;                  'min'
;                  'max'
;                  'total'
;                Default to 'average'
;   weights    - If set, then weight each of the input headers by these weights;
;                only applicable when the function type is 'average'.
;   zeros      - If set, then include zero values when determining the
;                average or other function.  But never use the zeros
;                returned by SXPAR() if a header is missing that card
;                altogether.
;   outcard    - Card name(s) in output header; if not specified, then use
;                the same name as in CARDNAME.
;   _EXTRA     - Optional keywords for SXADDPAR (such as BEFORE,AFTER,FORMAT).
;
; OUTPUTS:
;   outhdr     - (Modified.)
;
; COMMENTS:
;
; BUGS:
;
; PROCEDURES CALLED:
;   sxaddpar
;   sxpar()
;
; REVISION HISTORY:
;   31-Jan-2002  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro sxcombinepar, hdrarr, cardname, outhdr, func=func, weights=weights, $
 zeros=zeros, outcard=outcard, _EXTRA=KeywordsForSxaddpar

   if (n_params() LT 3) then begin
      print, 'Syntax - sxcombinepar, hdrarr, cardname, outhdr, [ func=, /zeros, outcard= ]'
      return
   endif

   if (NOT keyword_set(func)) then func = 'average'
   if (n_elements(zeros) EQ 0) then zeros = 0
   if (NOT keyword_set(outcard)) then outcard = cardname
   if (NOT keyword_set(weights)) then weights = lonarr(n_elements(hdrarr)) + 1

   if (n_elements(outcard) NE n_elements(cardname)) then $
    message, 'Number of elements in OUTCARD and CARDNAME must agree'
   if (n_elements(hdrarr) NE n_elements(weights)) then $
    message, 'Number of elements in HDRARR and WEIGHTS must agree'

   ;----------
   ; Call this routine recursively if CARDNAME has multiple elements

   ncard = n_elements(cardname)
   if (ncard EQ 0) then return
   if (ncard GT 1) then begin
      for icard=0, ncard-1 do begin
         sxcombinepar, hdrarr, cardname[icard], outhdr, func=func, $
          weights=weights, zeros=zeros, outcard=outcard[icard]
      endfor
      return
   endif

   for ihdr=0, n_elements(hdrarr)-1 do begin
      thisval = sxpar(*hdrarr[ihdr], cardname[0], count=thisct)
      if (thisct GT 0 AND (keyword_set(thisval) or keyword_set(zeros))) $
       then begin $
         if (NOT keyword_set(allval)) then allval = thisval $
          else allval = [allval, thisval]
      endif
   endfor

   nval = n_elements(allval)
   if (nval GT 0) then begin
      case strlowcase(func) of
         'average': outval = total(allval * weights) / total(weights)
         'median' : outval = median(allval)
         'min'    : outval = min(allval)
         'max'    : outval = max(allval)
         'total'  : outval = total(allval)
         else     : message, 'Invalid FUNCTION'
      endcase
      sxaddpar, outhdr, outcard[0], outval, _EXTRA=KeywordsForSxaddpar
   endif

   return
end
;------------------------------------------------------------------------------
