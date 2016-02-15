;+
; NAME:
;   hogg_strsplit
; PURPOSE:
;   split strings on whitespace, except inside double quotes, plus other
;     stuff specialized for yanny_read
; BUGS:
;   demolishes the string
; REVISION HISTORY:
;   2002-10-10  written - Hogg
;-
pro hogg_strsplit, string, output, count, recurse=recurse, verbose=verbose

   if keyword_set(verbose) then splog, 'splitting string >'+string+'<'

   ; Initialize unset variables, putting in first-element dummy value
   ; and splitting out comments
   if (NOT keyword_set(recurse)) then begin
      output= 'NULL'
      count= 0
      string= (strsplit(' '+strcompress(string),'#',/extract))[0]
   endif

   ; Do the dumbest thing, if possible
   if (strcompress(string,/remove_all) EQ '') then return

   ; Do the sedond-dumbest thing, if possible
   if stregex(string,'[\"]') LT 0 then begin

      pos = stregex(string,'\{ *\{ *\} *\}',length=len)
      if (pos GE 0) then begin
         string= strmid(string,0,pos)+' "" '+strmid(string,pos+len)
         hogg_strsplit, string, output, count, /recurse, $
                        verbose=verbose
      endif else begin

         ; just split on spaces
         word = strsplit(string, '[ ,;]+', /regex,/extract)
         output= [output, word]
         count = count + n_elements(word)
      endelse
   endif else begin

      ; split on quotation marks and operate recursively
      ; Find the position and length of the first double-quoted string.
      pos = stregex(string,'\"([^"]*)\"',length=len)
      if (pos GE 0) then begin

         ; Split everyting prior to that quote, appending to OUTPUT
         hogg_strsplit, strmid(string,0,pos), output, count, /recurse, $
                        verbose=verbose

         ; Now add to that the quoted string, but excluding the quotation
         ; marks themselves.
         word = strmid(string,pos+1,len-2)
         output= [output, word]
         count = count + 1

         ; Finally, split everything after the quoted part,
         ; which might contain more quoted strings.
         hogg_strsplit, strmid(string,pos+len), output, count, /recurse, $
                        verbose=verbose
      endif
   endelse

   ; Remove first-element dummy value
   if (NOT keyword_set(recurse)) AND (count GT 0) then begin
      output= output[1:count]
      if keyword_set(verbose) then for i=0,count-1 do print, i,'>'+output[i]+'<'
   endif
   return
end

