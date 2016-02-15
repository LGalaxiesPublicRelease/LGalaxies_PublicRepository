;-----------------------------------------------------------------------
;+
; NAME:
;   djs_modfits
;
; PURPOSE:
;   Wrapper for MODFITS that allows the header to increase in size,
;   and that works with g-zipped files.
;
; CALLING SEQUENCE:
;   djs_modfits, filename, data, [hdr, exten_no=, /delete_data]
;
; INPUTS:
;   filename  - FITS file name; if the name ends in the suffix ".gz",
;               then the file is g-unzipped first, modified, then re-g-zipped.
;   data      - New data array or structure for extension EXTEN_NO;
;               if this is undefined or zero, then don't modify the data.
;
; OPTIONAL INPUTS:
;   hdr       - New FITS header for extension EXTEN_NO
;   exten_no  - FITS extension number to modify; default to 0.
;   delete_data - If set, then delete this data.  Note that this cannot
;               be accomplished by setting DATA=0, since that simply says
;               to not change the data array/structure (to be consistent
;               with the functionality of MODFITS).
;
; OUTPUTS:
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
;   modfits
;   mrdfits()
;   mwrfits
;   readfits()
;   writefits
;
; INTERNAL PROCEDURES:
;   bitsperpixel()
;
; REVISION HISTORY:
;   17-May-2000  Written by David Schlegel, Princeton.
;-
;-----------------------------------------------------------------------
function bitsperpixel, type

   tvec = [0,1,2,4,4,8,8,0,0,16,0,0,2,4,8,8]
   return, tvec[type]
end
;-----------------------------------------------------------------------
pro djs_modfits, filename, thisdata, hdr, exten_no=exten_no, $
 delete_data=delete_data

   ; Need at least 2 parameters
   if (N_params() LT 2) then begin
      print, 'Syntax - djs_modfits, filename, data, [hdr, exten_no=, /delete_data]'
      return
   endif

   ; If the file appears to be g-zipped, then unzip the file first
   if (strmatch(filename, '*.gz')) then begin
      spawn, 'gunzip '+filename
      thisfile = strmid(filename, 0, strlen(filename)-3)
   endif else begin
      thisfile = filename
   endelse

   if (NOT keyword_set(exten_no)) then exten_no = 0

   ;----------
   ; Make sure the header cards are all exactly 80 characters

   if (keyword_set(hdr)) then begin
      hdr = strmid(hdr+string(' ',format='(a80)'),0,80)
   endif

   ;----------
   ; If we don't think the header size or data size will require a new block,
   ; call MODFITS to modify the header

   qbigger = 0

   if (keyword_set(hdr)) then begin
      hdr1 = headfits(thisfile, exten=exten_no)
      if (NOT keyword_set(hdr1)) then $
       message, 'EXTEN_NO does not exist in file '+thisfile
;      if ((n_elements(hdr)+79)/80 GT (n_elements(hdr1)+79)/80) then qbigger = 1
      if (n_elements(hdr) GT n_elements(hdr1)) then qbigger = 1
   endif

   qstruct = 0
   if (keyword_set(thisdata) OR keyword_set(delete_data)) then begin
      data1 = mrdfits(thisfile, exten_no)
      qstruct = size(thisdata, /tname) EQ 'STRUCT'
      if (qstruct) then begin
         nbytes1 = n_tags(data1, /length)
         nbytes = n_tags(thisdata, /length)
      endif else begin
         nbytes1 = n_elements(data1) * bitsperpixel(size(data1,/type))
         if (keyword_set(delete_data)) then $
          nbytes = 0 $
         else $
          nbytes = n_elements(thisdata) * bitsperpixel(size(thisdata,/type))
      endelse
      if (NOT keyword_set(data1) OR $
       (nbytes+2879)/2880 GT (nbytes1+2879)/2880) then qbigger = 1
   endif

   ;----------
   ; For now, the Goddard routine MODFITS does not work with structures.
   ; So don't use it in that case.

   if ( (qbigger EQ 0) AND (qstruct EQ 0) $
    AND (NOT keyword_set(delete_data)) ) then begin

      modfits, thisfile, thisdata, hdr, exten_no=exten_no

   endif else begin

      ;----------
      ; Read in all the headers and data arrays for FILENAME

      data1 = readfits(thisfile, hdr1)
      pdata = ptr_new(data1)
      phdr = ptr_new(hdr1)
      nhdu = 1

      while (keyword_set(hdr1)) do begin
         hdr1 = 0
         data1 = mrdfits(thisfile, nhdu, hdr1)
         if (keyword_set(hdr1)) then begin
            phdr = [phdr, ptr_new(hdr1)]
            pdata = [pdata, ptr_new(data1)]
            nhdu = nhdu + 1
         endif
      endwhile

      if (exten_no GE n_elements(phdr)) then $
       message, 'EXTEN_NO does not exist in file '+thisfile

      ;----------
      ; Re-write all the headers and data arrays to FILENAME, modifying the
      ; one specified by EXTEN_NO

      if (keyword_set(hdr)) then begin
         ptr_free, phdr[exten_no]
         phdr[exten_no] = ptr_new(hdr)
      endif

      if (keyword_set(delete_data)) then begin
         ptr_free, pdata[exten_no]
         pdata[exten_no] = ptr_new(0)
      endif else if (keyword_set(thisdata)) then begin
         ptr_free, pdata[exten_no]
         pdata[exten_no] = ptr_new(thisdata)
      endif

      writefits, thisfile, *(pdata[0]), *(phdr[0])
      for ihdu=1, nhdu-1 do begin
         mwrfits, *(pdata[ihdu]), thisfile, *(phdr[ihdu])
      endfor

      ;----------
      ; Free memory

      for ihdu=0, nhdu-1 do begin
         ptr_free, phdr[ihdu]
         ptr_free, pdata[ihdu]
      endfor

   endelse

   ;----------
   ; If the file appears to be g-zipped, then re-g-zip the file at the end

   if (strmatch(filename, '*.gz')) then begin
      spawn, 'gzip '+thisfile
   endif

   return
end 
;-----------------------------------------------------------------------
