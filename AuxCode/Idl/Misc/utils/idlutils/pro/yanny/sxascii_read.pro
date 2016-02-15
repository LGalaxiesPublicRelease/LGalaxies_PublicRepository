;+
; NAME:
;   sxascii_read 
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   sxascii_read, filename, [ pstr, hdr=hdr]
;
; INPUTS:
;   filename   - Input file name for Yanny parameter file
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;   pstr       - Only one structure can be returned
;   hdr        - hdr will include the query
;
; COMMENTS:
;   Return 0's if the file does not exist.
;
;   Read and write variables that are denoted INT in the Yanny file
;   as IDL-type LONG, and LONG as IDL-type LONG64.  This is because
;   Yanny files assume type INT is a 4-byte integer, whereas in IDL
;   that type is only 2-byte.
;
; EXAMPLES:
;
; BUGS:
;   Not set up yet to deal with multi-dimensional arrays.
;
;   Not yet tested with char?
;
; PROCEDURES CALLED:
;   needs to compile yanny_read.pro to run
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   16-May-2001  Written by S Burles, FNAL
;   18-Jun-2002  Revisions (better or worse?)
;-
;------------------------------------------------------------------------------
pro sxascii_read, filename, pstr, hdr=hdr

   if (N_params() LT 1) then begin
      print, 'Syntax - sxascii_read, filename, [ pstr, hdr=hdr]'
      return
   endif

   tname = ['char', 'short', 'int', 'long', 'float', 'double']
   ; Read and write variables that are denoted INT in the Yanny file
   ; as IDL-type LONG, and LONG as IDL-type LONG64.  This is because
   ; Yanny files assume type INT is a 4-byte integer, whereas in IDL
   ; that type is only 2-byte.
;   tvals = ['""'  , '0'    , '0'  , '0L'  , '0.0'  , '0.0D'  ]
;   tarrs = ['strarr(' , $
;            'intarr(' , $
;            'intarr(' , $
;            'lonarr(' , $
;            'fltarr(' , $
;            'dblarr(' ]
   tvals = ['""'  , '0'    , '0L'  , '0LL'  , '0.0'  , '0.0D'  ]
   tarrs = ['strarr(' , $
            'intarr(' , $
            'lonarr(' , $
            'lon64arr(' , $
            'fltarr(' , $
            'dblarr(' ]

   hdr = 0
   hdrpresent = ARG_PRESENT(hdr)
   enums = 0
   structs = 0

   ; List of all possible structures
   pcount = 0       ; Count of structure types defined
   pname = 0        ; Name of each structure
   pdata = 0        ; Pointer to each structure
   pnumel = 0       ; Number of elements in each structure

   get_lun, ilun
   openr, ilun, filename, error=err
   if (err NE 0) then begin
      close, ilun
      free_lun, ilun
      return
   endif

;
; Quickly check length
;
   filelength = 0L
   rawline = ''
   while (NOT eof(ilun)) do begin
      readf, ilun, rawline 
      filelength = filelength + 1
   endwhile

   print, filelength
   point_lun, ilun, 0L
   
   qdone = 0
   while (NOT eof(ilun)) do begin

      ; Read the next line
      readf, ilun, rawline 
      words = yanny_getwords(rawline) ; Divide into words and strings
      nword = N_elements(words)
      amp = strpos(rawline,'#')

      if (strmid(rawline,0,7) EQ '#XSTART') then begin
         ; LOOK FOR STRUCTURES TO BUILD 

        fulltag = 0L 
        junk = ''
        readf, ilun, junk, fulltag, format='(a1,d)'
      
          
        readf, ilun, rawline
        rawtags = rawline
        readf, ilun, rawline
        rawtypes = rawline
        readf, ilun, rawline
        rawformats = rawline
        readf, ilun, rawline
        ends = rawline

        if (strmid(ends,0,5) NE '#XEND') then message, 'Confused'

        tags = STR_SEP(rawtags,'"')
        good = where(strmid(tags,0,1) NE ' ' AND tags NE '',ngood)
        if ngood LT 2 then message, 'NO tags?'

        tags = tags[good[1:*]]
        ntags = n_elements(tags)
        print, ntags, fulltag

        ;look for arrays 
       
          
        stags = sort(tags)
        cleantags = strarr(ntags)
        for i=0,ntags - 1 do begin 
          junk = tags[i] 
          brack = strpos(junk,'(') 
          if brack GT 0 then begin
             junk=strmid(junk,0,brack) 
          endif
          cleantags[i] = junk 
        endfor


        place = uniq(cleantags[stags])
;        arrlen = place - [-1,place] 
; 
;        scalar = where(arrlen EQ 1)
;        if scalar[0] NE -1 then $
;             cleantags[stags[place[scalar]]] = tags[stags[place[scalar]]]

        ;Now dump periods and anything proceeding periods
        cleantags2 = strarr(ntags)

        for i=0,ntags - 1 do begin & $
          junk = cleantags[i] & $
          brack = strpos(junk,'.') & $
          while brack NE -1 do begin & $
             junk=strmid(junk,brack+1) & $
             brack = strpos(junk,'.') & $
          endwhile & $
          brack = strpos(junk,'[') & $
          brack2 = strpos(junk,']',/reverse_search) - brack - 1 & $
          if brack NE -1 AND brack2 GT 0 then $
             junk = strmid(junk,0,brack) + strmid(junk,brack+1,brack2) & $
          cleantags2[i] = junk & $
        endfor

        place2 = uniq(cleantags2[sort(cleantags2)])

        nplace = n_elements(place2)
        if n_elements(place) NE nplace then $
           message,'You have conflicting tags'

        place3 = uniq(cleantags2)
        arrlen = place3 - [-1,place3] 
    
   
        ; setup types
  
        types = (str_sep(rawtypes," " ))[1:*] 
        cleantypes = types[place3]
        
        for itag=0, nplace-1 do begin
          i = (where(cleantypes[itag] EQ tname, ct))[0]
          if i EQ -1 then i = 0
          addname = cleantags2[place3[itag]]

          if arrlen[itag] GT 1 then $
             addval = tarrs[i] + strtrim(string(arrlen[itag]),2) + ')' $
          else addval = tvals[i]
   
          if NOT keyword_set(names) then begin
            names = addname
            values = addval
          endif else begin
            names = [names, addname]
            values = [values, addval]
          endelse
 
       endfor
       pstr = mrd_struct(names, values, filelength)
       temppstr = pstr[0]

       qdone = 1
     ; LOOK FOR A STRUCTURE ELEMENT
     ; Only look if some structures already defined
     ; Note that the structure names should be forced to uppercase
     ; such that they are not case-sensitive.
     endif else if (amp EQ 0) then begin
          if (hdrpresent NE 0) then yanny_add_comment, rawline, hdr 
     endif else if (qdone EQ 1) then begin

            ; Split this text line into words
            sline = strcompress( yanny_strip_brackets(rawline) )
            ww = yanny_getwords(sline)

            i = 0 ; Counter for which word we're currently reading
                     ; Skip the 0-th word since it's the structure name.

            ntag = N_tags( temppstr )
            ; Now fill in this structure from the line of text
            for itag=0, ntag-1 do begin
               ; This tag could be an array - see how big it is
               sz = N_elements( temppstr.(itag) )
                  for j=0, sz-1 do begin
                     temppstr.(itag)[j] = ww[i]
                     i = i + 1
                  endfor
            endfor

            pstr[pnumel] = temppstr
            pnumel = pnumel + 1
            if ((pnumel mod 100) EQ 0) then print, pnumel
      endif

   endwhile

   close, ilun
   free_lun, ilun
   pstr = pstr[0:pnumel-1]

   return
end
;------------------------------------------------------------------------------
