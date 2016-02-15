;+
; NAME:
;   lookforgzip
;
; PURPOSE:
;   Look for a gzip-ed file, 
;
; CALLING SEQUENCE:
;   thisfile = lookforgzip( filename, count= )
;
; INPUTS:
;   filename   - Input file name w/out any ".gz" or ".Z" extension
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   thisfile   - Returns input file name with no extension if it exists,
;                otherwise a ".gz" extension if that exists,
;                otherwise a ".Z" extension if that exists,
;                otherwise '' if none of the above exist.
;
; OPTIONAL OUTPUTS:
;   count      - Number of files that matched
;
; COMMENTS:
;   This routine uses FINDFILE().
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   20-Oct-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function lookforgzip, filename, count=ct

   thisfile = findfile(filename, count=ct)
   if (ct GT 0) then return, thisfile

   thisfile = findfile(filename+'.gz', count=ct)
   if (ct GT 0) then return, thisfile

   thisfile = findfile(filename+'.Z', count=ct)
   if (ct GT 0) then return, thisfile

   return, ''
end
;------------------------------------------------------------------------------
