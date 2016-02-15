;+
; NAME:
;	catfits
; PURPOSE:
;       Read and concatenate a bunch of binary FITS tables
;
; CALLING SEQUENCE:
;	struc = catfits(flist)
;
; INPUTS:
;       filst  - string array of FITS filenames to read and concatenate
;		for each column of data to be read.  Allowed letters are 
; OUTPUTS:
;	struc  - IDL structure
;
; REVISION HISTORY:
;	Written  2001-Nov-28 D. Finkbeiner
;-

function catfits, flist
  
  n = n_elements(flist)
  for i=0, n-1 do begin 
     a = mrdfits(flist[i], 1, /silent)
     
     if n_elements(aa) EQ 0 then aa = a $
     else begin 
; test whether structures are compatible:
        if n_tags(a) GT n_tags(aa) then begin 
           temp = replicate(a[0], n_elements(aa))
           struct_assign, aa, temp
           aa = temporary(temp)
        endif 
        if n_tags(a) EQ n_tags(aa) AND $
          total(tag_names(a) NE tag_names(aa)) EQ 0 then begin 
           aa = [aa, a]
        endif else begin 
           temp = replicate(aa[0], n_elements(a))
           struct_assign, a, temp
           aa = [aa, temp]
        endelse
     endelse
;     aa = n_elements(aa) EQ 0 ? a : [aa, a]
  endfor

  return, aa
end
