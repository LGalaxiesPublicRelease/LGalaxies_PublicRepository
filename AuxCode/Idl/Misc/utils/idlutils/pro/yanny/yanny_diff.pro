;+
; NAME:
;   yanny_diff
;
; PURPOSE:
;   Compare two Yanny files, optionally copying values from one to another.
;
; CALLING SEQUENCE:
;   yanny_diff, filename1, filename2, [ outfile, $
;    cardnames=, stnames=, /verbose, count=, errcode= ]
;
; INPUTS:
;   filename1  - Yanny parameter file from which copy cards
;   filename2  - Yanny parameter file to which copy cards
;
; OPTIONAL INPUTS:
;   outfile    - If specified, then write the modified 2nd file to this file.
;   stnames    - The Yanny structure names to compare between the two
;                input files; if not set, then compare all structures that
;                exist in the 2nd file.  These names are case-insensitive.
;   cardnames  - The tag names in any structure to compare; if not set,
;                then compare all tags that exist in the 2nd file.
;                Note that if this is set to 'FOO', then that tag is
;                compared (if it exists) between all specified structures.
;                These names are case-insensitive.
;   verbose    - If set, then explicitly print each change that is made;
;                if not set, then simply print the number of changes made.
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;   count      - The total number of elements that differ.
;   errcode    - Return 0 upon successful completion, non-zero otherwise.
;                A code of -1 means that there were errors reading the
;                input files.  A code of +1 means that the number of
;                rows or elements in the input structures did not agree.
;                In the latter case, changes are still made to any other
;                structure elements where the numbers do agree.
;
; COMMENTS:
;
; EXAMPLES:
;   Simply compare two opECalib files, printing any differences:
;     IDL> yanny_diff,'opECalib-50000.par','opECalib-51430.par',/verbose
;
;   Compare the same two files, but copy values for the structure
;   ecalib.fullwelldn2 from the 1st file to the 2nd, creating a 3rd file:
;     IDL> yanny_diff,'opECalib-50000.par','opECalib-51430.par', $
;     IDL>  'opECalib-new.par', stnames='ecalib', cardnames='fullwelldn2'
;
; BUGS:
;
; PROCEDURES CALLED:
;   splog
;   yanny_free
;   yanny_read
;   yanny_write
;
; REVISION HISTORY:
;   26-Feb-2002  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro yanny_diff, filename1, filename2, outfile, $
 cardnames=cardnames, stnames=stnames, verbose=verbose, $
 count=count, errcode=errcode

   count = 0L
   errcode = 0L

   if (N_params() LT 2) then begin
      print, 'Syntax - yanny_diff, filename1, filename2, [ outfile, $'
      print, '  cardnames=, stnames=, /verbose ]'
      errcode = -1L
      return
   endif
   if (NOT keyword_set(stnames)) then stnames = ''
   if (NOT keyword_set(cardnames)) then cardnames = ''

   ;----------
   ; Read the two Yanny files

   yanny_read, filename1, pdata1, hdr=hdr1, enums=enums1, structs=structs1, $
    stnames=stnames1, /anonymous, errcode=errcode1
   if (errcode1 NE 0) then begin
      splog, 'Warning: Error reading input file ' + filename1
      errcode = errcode1
      return
   endif

   yanny_read, filename2, pdata2, hdr=hdr2, enums=enums2, structs=structs2, $
    stnames=stnames2, /anonymous, errcode=errcode2
   if (errcode1 NE 0) then begin
      splog, 'Warning: Error reading input file ' + filename2
      errcode = errcode2
      return
   endif

   ;----------
   ; Loop over each structure looking for differences

   for is2=0, n_elements(stnames2)-1 do begin
      is1 = (where(stnames1 EQ stnames2[is2]))[0]
      if ((NOT keyword_set(stnames) $
        OR total(strupcase(stnames) EQ stnames2[is2]) EQ 1) $
       AND is1 NE -1) then begin
         splog, 'Looking at structure ' + stnames2[is2]
         tags1 = tag_names(*pdata1[is1])
         tags2 = tag_names(*pdata2[is2])

         ;----------
         ; Loop over each tag name in this structure

         for ic2=0, n_elements(tags2)-1 do begin
            ic1 = (where(tags1 EQ tags2[ic2]))[0]
            if ((NOT keyword_set(cardnames) $
              OR total(strupcase(cardnames) EQ tags2[ic2]) EQ 1) $
             AND ic1 NE -1) then begin
               if (n_elements((*pdata1[is1]).(ic1)) $
                EQ n_elements((*pdata2[is2]).(ic2))) then begin

                  indx = where((*pdata1[is1]).(ic1) NE (*pdata2[is2]).(ic2), ct)
                  count = count + ct
                  if (ct GT 0) then begin
                     if (keyword_set(verbose)) then begin
                        narr = n_elements( (*pdata1[is1])[0].(ic1) )
                        for ii=0, ct-1 do begin
                           indexnum = indx[ii] MOD narr
                           linenum = (indx[ii] - indexnum) / narr
                           thisname = stnames2[is2] $
                            + '[' + strtrim(string(linenum),2) + '].' $
                            + tags2[ic2]
                           if (narr GT 1) then thisname = thisname $
                            + '[' + strtrim(string(indexnum),2) + ']'

                           splog, 'Changing ' + thisname + ' = ', $
                            ((*pdata1[is1]).(ic1))[indx[ii]], $
                            "->", ((*pdata2[is2]).(ic2))[indx[ii]]
                        endfor
                     endif else begin
                        splog, ct, ' different values for card=' $
                         + tags2[ic2]
                     endelse
                     (*pdata2[is2]).(ic2) = (*pdata1[is1]).(ic1)
                  endif

               endif else begin
                  splog, 'Warning: Number of elements do not agree for card=' $
                   + tags2[ic2]
                  errcode = 1L
               endelse
            endif
         endfor

      endif
   endfor

   ;----------
   ; Write the output file

   if (keyword_set(outfile)) then $
    yanny_write, outfile, pdata2, hdr=hdr2, enums=enums2, structs=structs2, $
     stnames=stnames2

   yanny_free, pdata1
   yanny_free, pdata2

   return
end
;------------------------------------------------------------------------------
