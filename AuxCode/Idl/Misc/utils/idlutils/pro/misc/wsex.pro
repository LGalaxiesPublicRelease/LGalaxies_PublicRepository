pro wsex,cat,outfile=outfile
;+
; NAME:
;       WSEX
;
; PURPOSE:
;       Write out arbitrary SExtractor format catalogs.
;
; INPUTS:
;       A SExtractor-format catalog
;
; OUTPUTS:
;       Prints out a catalog file for the given structure
;
; COMMON BLOCKS:
;       None.
;
; RESTRICTIONS:
;       None.
;
; PROCEDURE:
;       Use syntax
;       wset,catalog,outfile='catalog.cat'
;
; COMMENTS:
;       To do -- 
;       . better error checking, e.g. for existing output filename,
;       etc. 
;
; PROCEDURES USED:
;       VALID_NUM
;
; MODIFICATION HISTORY:
;       L. Moustakas '04feb11 - created
;-

; Check that an argument has been passed     
     if n_params() le 0 then begin 
         print,'wsex,cat,outfile=''catalog.cat'''
         return
     endif 

     if n_elements(outfile) eq 0 then begin
; herein, assign an output catalog name
         outfile='wsexout.cat'
     endif

     nbody = n_elements(cat)
     nhead = n_tags(cat)
     tags  = tag_names(cat)

     u0=25
     openw,u0,outfile
     for i=0l,nhead-1 do $
       printf,u0,'# '+string(i+1,format='(i4)')+' '+strupcase(tags[i])

     tind    = intarr(nhead)
     tindint = intarr(nhead)

     for i=0l,nhead-1 do tind[i]    = valid_num((cat[0].(i))[0])
     for i=0l,nhead-1 do tindint[i] = valid_num((cat[0].(i))[0],/int)

     formstr='(a20)'
     formlon='(i10)'
     formdbl='(f14.7)'
     formarr=strarr(nhead)
     for i=0l,nhead-1 do $
       if tind[i] eq 0 then $
       formarr[i]=formstr else $
       if tindint[i] eq 1 then $
       formarr[i]=formlon else $
       formarr[i]=formdbl

     for k=0l,nbody-1 do begin 
         tstr=''
         for i=0l,nhead-1 do tstr=tstr+string(cat[k].(i),format=formarr[i])
         printf,u0,tstr
     endfor
     close,u0,/force
 end  
