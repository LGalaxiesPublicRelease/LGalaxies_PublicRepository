;+
; NAME:
;   HOGG_MRDFITS
; PURPOSE:
;   Wrapper on MRDFITS() to read in a FITS file one chunk at a time
; CALLING SEQUENCE:
; INPUTS:
;   see MRDFITS
; KEYWORDS:
;   see MRDFITS
;   nrowchunk  - Number of rows to read at a time; default to 2880
;   nchunk     - Number of chunks to read; this keyword takes precedence
;                over NROWCHUNk
; COMMENTS:
;   Useful when "columns" is set, so you can get a couple
;   of columns without holding the whole file in memory
; OUTPUTS:
;   see MRDFITS
; REVISION HISTORY:
;   2002-02-08  written by Hogg
;-
function hogg_mrdfits, file, extension, header, silent=silent, $
 range=range1, structyp=structyp, nrowchunk=nrowchunk1, nchunk=nchunk, $
 _EXTRA=inputs_for_mrdfits

   if (keyword_set(range1) EQ 0 OR keyword_set(nchunk)) then begin
      hdr = headfits(file,exten=extension,silent=silent)
      if (size(hdr,/tname) NE 'STRING') then return, 0
      naxis2= long64(sxpar(hdr,'NAXIS2'))
   endif

   if (keyword_set(range1)) then begin
      range = range1
   endif else begin
      if (NOT keyword_set(silent)) then splog, naxis2
      range= [0LL,naxis2-1LL]
   endelse

   if (keyword_set(nchunk)) then begin
      nrowchunk = ceil(naxis2/double(nchunk))
   endif else begin
      if (keyword_set(nrowchunk1)) then nrowchunk = long(nrowchunk1) $
       else nrowchunk = 2880L
   endelse

   chunkrange = [range[0],((range[0]+nrowchunk-1LL) < range[1])]
   while (chunkrange[1] GE chunkrange[0]) do begin

      if (NOT keyword_set(silent)) then splog, chunkrange
      chunkresult = mrdfits(file, extension, header, $
       range=chunkrange, structyp=structyp, $
       silent=silent, _EXTRA=inputs_for_mrdfits)
      if (NOT keyword_set(result)) then $
       result = replicate(chunkresult[0],range[1]-range[0]+1LL)
      result[chunkrange[0]-range[0]:chunkrange[1]-range[0]] = chunkresult

      chunkrange[0] = chunkrange[0]+nrowchunk
      chunkrange[1] = (chunkrange[1]+nrowchunk) < range[1]
   endwhile

   return, result
end
