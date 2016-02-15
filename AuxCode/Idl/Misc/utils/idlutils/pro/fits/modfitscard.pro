;-----------------------------------------------------------------------
;+
; NAME:
;   modfitscard
;
; PURPOSE:
;   Modify FITS card(s) in a file without changing the data.
;
; CALLING SEQUENCE:
;   modfitscard, filename, card, value, [ comment, /delete, $
;    _EXTRA=KeywordsForSxaddpar ]
;
; INPUTS:
;   filename  - File name(s) to modify; this can be an array of file names,
;               and it can include wildcards
;   card      - Name of FITS card(s) to add or modify
;
; OPTIONAL INPUTS:
;   value     - New value(s) for FITS card; mandatory card if DELETE not set;
;               must have the same number of elements as CARD.
;   comment   - Comment to appear in the card after its value; passed to
;               the routine SXADDPAR.  If specified, it must have the same
;               number of elements as CARD.
;   delete    - If set, then delete all cards CARD from the header;
;               VALUE is ignored if set.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   Modify the value of the DATE keyword in the primary header of all FITS
;   files with '666' or '777' in the file name:
;   IDL> modfitscard, ['*666*.fits','*777*.fits'], 'DATE', '1994-03-23'
;
; BUGS:
;   This routine calls DJS_MODFITS, which allows the size of the header
;   to be changed.
;
;   Wildcards are not supported for CARD, so you cannot say something like
;   IDL> modfitscard, 'test.fits', 'DATE*', '1994-03-23' ; Will not work!
;
; PROCEDURES CALLED:
;   djs_modfits
;   headfits()
;   sxaddpar
;
; REVISION HISTORY:
;   19-Apr-2000  Written by David Schlegel, Princeton.
;-
;-----------------------------------------------------------------------
pro modfitscard, filename, card, value, comment, delete=delete, $
 _EXTRA=KeywordsForSxaddpar

   ; Need at least 3 parameters
   if (N_params() LT 3 AND (N_params() EQ 2 AND NOT keyword_set(delete))) then begin
      print, 'Syntax - modfitscard, filename, card, value, _EXTRA=KeywordsForSxaddpar'
      return
   endif

   ncard = N_elements(card)
   if (ncard NE N_elements(value) AND NOT keyword_set(delete)) then begin
      print, 'Number of elements in CARD and VALUE do not agree'
      return
   endif

   if (keyword_set(comment)) then begin
      if (ncard NE N_elements(comment)) then begin
         print, 'Number of elements in CARD and COMMENT do not agree'
         return
      endif
   endif

   ;----------
   ; Parse the file name list to be all the specified files that are on disk.

   for ii=0, N_elements(filename)-1 do begin
      tempname = findfile(filename[ii], count=ct)
      if (ct GT 0) then begin
         if (keyword_set(fullname)) then fullname = [fullname, tempname] $
          else fullname = tempname
      endif
   endfor
   nfile = N_elements(fullname)

   ;----------
   ; Loop through 1 file at a time, modifying the header.

   for ifile=0, nfile-1 do begin

      hdr = headfits(fullname[ifile])
      if (size(hdr, /tname) NE 'STRING') then begin
         print, 'File does not exist or FITS header is invalid: ', $
          fullname[ifile]
      endif else begin

         if (keyword_set(delete)) then begin
;            sxdelpar, hdr, card
            keyword = strmid(hdr, 0, 8)
            for icard=0, ncard-1 do begin
               cardname = string(strupcase(card[icard])+'        ',format='(a8)')
               ifound = where(keyword EQ cardname)
               if (ifound[0] NE -1) then $
                hdr[ifound] = string(' ', format='(" ",a79)')
            endfor
         endif else begin
            for icard=0, ncard-1 do begin
               if (keyword_set(comment)) then $
                sxaddpar, hdr, card[icard], value[icard], comment[icard], $
                 _EXTRA=KeywordsForSxaddpar $
               else $
                sxaddpar, hdr, card[icard], value[icard], $
                 _EXTRA=KeywordsForSxaddpar
            endfor
         endelse

         djs_modfits, fullname[ifile], 0, hdr
      endelse
   endfor

   return
end 
;-----------------------------------------------------------------------
