;+
; NAME:
;   yanny_par
;
; PURPOSE:
;   Obtain the value of a parameter in the header of a Yanny file.
;
; CALLING SEQUENCE:
;   result = yanny_par( hdr, keyname, [count=, indx= ] )
;
; INPUTS:
;   hdr        - Header lines in Yanny file, which are usually keyword pairs.
;   keyname    - Keyword name of which to find its corresponding value.
;                If this is an array, then this routine is called recursively
;                and the results appended together.
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;   result     - Value of parameter in header as either a string or an
;                array of strings; return '' if parameter not found
;
; OPTIONAL OUTPUTS:
;   count      - Return the number of parameters found.  There may be more
;                returned values than INDX, if there are several values on
;                the same line.
;   indx       - Index of the lines that match (0-indexed); -1 for no match.
;
; COMMENTS:
;   This routine is meant to be analogous to the Goddard function SXPAR()
;   for reading from FITS headers.
;
; EXAMPLES:
;
; BUGS:
;   Wildcards are not supported for KEYNAME.
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   02-Nov-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function yanny_par, hdr, keyname, count=count, indx=indx

   count = 0L

   if (N_params() LT 2) then begin
      print, 'Syntax - result = yanny_par(hdr, keyname, [count=, indx= ] )'
      return, ''
   endif

   nhead = N_elements(hdr)
   if (nhead EQ 0) then begin
      print, 'Header contains no elements'
      return, ''
   endif

   if (size(keyname, /tname) EQ 'UNDEFINED') then begin
      print, 'KEYNAME is undefined'
      return, ''
   endif

   ; Call this routine recursively if KEYNAME is an array
   nkey = n_elements(keyname)
   if (nkey GT 1) then begin
      for ikey=0, nkey-1 do begin
         res1 = yanny_par(hdr, keyname[ikey], count=count1, indx=indx1)
         if (count1 GT 0) then begin
            count = count + 1
            if (NOT keyword_set(result)) then begin
               result = res1
               indx = indx1
            endif else begin
               result = [result, res1]
               indx = [indx, indx1]
            endelse
         endif
      endfor
      if (NOT keyword_set(result)) then begin
         result = ''
         indx = -1L
      endif
      return, result
   endif

   keylist = strarr(nhead)
   keystring = strarr(nhead)
   for i=0, nhead-1 do begin
      keylist[i] = (str_sep( strtrim(hdr[i],2), ' '))[0]
   endfor

   ; Locate the first keyword that matches
   indx = where(keyname[0] EQ keylist, ct)
   if (ct EQ 0) then return, ''
   if (ct EQ 1) then begin
      j = indx[0]
   endif else begin
      for j=0, ct-1 do begin
         result1 = yanny_par(hdr[indx[j]], keyname[0])
         if (NOT keyword_set(result)) then result = result1 $
          else result = [result, result1]
      endfor
      ; If the result has only 1 element, then return a scalar.
      count = n_elements(result)
      if (count EQ 1) then result = result[0]
      return, result
   endelse

   ; Find the string after the keyword
   ipos = strpos(hdr[j], keylist[j]) + strlen(keylist[j])
   keystring = strtrim( strmid(hdr[j], ipos+1), 2 )

   ; Remove any comments from this string
   ipos = strpos(keystring, '#')
   if (ipos EQ 0) then keystring = '' $
    else if (ipos GT 0) then keystring = strtrim(strmid(keystring, 0, ipos),2)

   ; If any single quote exists, then split the string by looking for
   ; everything within single quotes.
   ; Otherwise, split the string using spaces.
   if (strpos(keystring, "'") NE -1) then begin
      ; The following is a kludge to take strings between successive
      ; pairs of single quotes.

; Below is the 5.3 version ???
;      result = strsplit(keystring, "'", /extract)
;      result = result[ 2*lindgen((n_elements(result)/2) > 1) ]

; Below is the 5.2 version ???
      result = str_sep(keystring, "'")
      result = result[ 2*lindgen((n_elements(result)/2) > 1) + 1 ]

   endif else begin
; Below is the 5.3 version ???
;      result = strsplit(keystring, " ", /extract)

; Below is the 5.2 version ???
      result = str_sep(strcompress(keystring), " ")
   endelse

   ; If the result has only 1 element, then return a scalar.
   count = n_elements(result)
   if (count EQ 1) then result = result[0]

   return, result
end        
;------------------------------------------------------------------------------
