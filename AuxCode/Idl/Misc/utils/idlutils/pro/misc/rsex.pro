;+
; NAME:
;       RSEX
;
; PURPOSE:
;       Read in arbitrary SExtractor format catalogs.
;
; INPUTS:
;       A SExtractor-format catalog
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       Returns a structure with all catalog entries, using field
;       names for tagnames. 
;
; COMMON BLOCKS:
;       None.
;
; RESTRICTIONS:
;       None.
;
; PROCEDURE:
;       Use syntax
;       cs = rsex('catalog.cat')
;
; COMMENTS:
;       Uses native header information & data themselves.  
;       Correctly reads longs, strings, and doubles.  
;
; PROCEDURES USED:
;       FILE_SEARCH
;       CREATE_STRUCT
;       VALID_NUM
;       VALID_NUM_ARR
;
; MODIFICATION HISTORY:
; LAM = L. Moustakas
;       LAM '04may02 - fixed problem case of there being an array of
;                      values based on the last header tag position.
;                      will now work with that case after special
;                      check. 
;       LAM '04feb04 - converted from the old lrsex.pro; adapted to
;                      detect and read longs, strings, and doubles.  
;       LAM '04may02 - correctly identifies and reads arrays at the
;                      end of the header
;       LAM '04dec07 - changed findfile use to file_search 
;-


function valid_num_arr, arr, integer=integer
     narr=n_elements(arr)
     nval=intarr(narr)
     for i=0l,narr-1 do nval[i]=valid_num(arr[i],integer=integer)
     return,nval
end


FUNCTION rsex,catalog
; Check that an argument has been passed     
     IF n_params() LE 0 THEN BEGIN 
         print,'cat=rsex(catalog)'
         return,-1
     ENDIF 

; See if the catalog file exists at all
     jnk = file_search(catalog,count=catexist)
     if catexist eq -1 then begin
         message,'file '+catalog+' does not seem to exist',/inf
         return,-1
     endif

; Count catalog lines
     spawn,'wc -l < '+catalog, nlines, /sh
     nlines = long(nlines[0]);-1

; Open the catalog to read
     u0=22
     openr,u0,catalog,error=err
     if err ne 0 then begin
         message,'error opening catalog - ',!error_state.msg,/inf
         close,u0,/force
         return,-1
     endif

; Read in the header
     head = ''
     numb = 0l
     cstr=''
     tag = 1
     while tag do begin         ; while '#' tag is true
         readf,u0,cstr
         if strpos(cstr,'#') ne -1 then begin
             head = [head, strupcase((str_sep(strcompress(cstr),' '))[2])]
             numb = [numb, (str_sep(strcompress(cstr),' '))[1]]
         endif else tag=0
     endwhile
     nhead = n_elements(head)   ; number of header entries
     head = head[1:(nhead-1)]
     numb = numb[1:(nhead-1)]
     nhead = n_elements(head)
     close,u0,/force

; Check header for implicit array entries, and expand relevant tagnames
     tothead=''
     for i=0l,nhead-2 do begin
         tothead = [tothead,head[i]]
         mlen = numb[i+1]-numb[i]-1
         if mlen ne 0 then $
           for j=1,mlen do $
           tothead=[tothead,head[i]+strcompress(j,/remove_all)]
     endfor

; The header info so far
     ntothead = n_elements(tothead)
     tothead=[tothead[1:(ntothead-1)],head[(nhead-1)]]

; Need to check one more thing -- whether the last header field is
; actually the first one of an array...

; Read the first data line.  If it has more entries than the total
; number of header fields so far, we've missed an array at the end. 
     junk=strarr(nhead)
     openr,u0,catalog
     readf,u0,junk
     readf,u0,cstr
     close,u0,/force
     nstr = n_elements(str_sep(strcompress(strtrim(cstr,2)),' '))

     if nstr gt ntothead then begin
         mlen = nstr - ntothead
         for j=1,mlen do $
           tothead=[tothead,head[(nhead-1)]+strcompress(j,/remove_all)]
     endif
     ntothead=n_elements(tothead)

; That should do it!  Onwards, to read the data...
     nbody=nlines-nhead
     junk=strarr(nhead)
     openr,u0,catalog
     readf,u0,junk
     readf,u0,cstr
     close,u0,/force
     body=strarr(nbody)
     openr,u0,catalog
     readf,u0,junk
     readf,u0,body
     close,u0,/force
     bodyslam=strarr(ntothead,nbody)
     for i=0l,nbody-1 do $
       bodyslam[*,i]=str_sep(strcompress(strtrim(body[i],2)),' ')

     carr=str_sep(strcompress(strtrim(cstr,1)),' ')
     tind    = valid_num_arr(carr) ; 0=string 1=number
     tindint = valid_num_arr(carr,/int) ; 0=non-integer 1=integer
     if tind[0] eq 0 then $
       type='' else $
       if tindint[0] eq 1 then $
       type=0l else $
       type=0.d0
     cs=create_struct(tothead[0],type)

     for i=1l,ntothead-1 do begin
         if tind[i] eq 0 then $
           type='' else $
           if tindint[i] eq 1 then $
           type=0l else $
           type=0.d0
         cs=create_struct(cs,tothead[i],type)
     endfor
     cs=replicate(cs,nbody)

     for i=0l,ntothead-1 do for j=0l,nbody-1 do cs[j].(i)=bodyslam[i,j]
     
     return,cs
END 
