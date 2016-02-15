;+
; NAME:
;	  lle_sm_mult
; PURPOSE: (one line)
;   multiply two sparse matrices kind of efficiently for LLE routines
; DESCRIPTION:
;   rather crappy sparse matrix format for NxM matrix:
;      .N - number of rows
;      .M - number of columns
;      .VALS - each nonzero value
;      .NINDX - row of each nonzero value 
;      .MINDX - column of each nonzero value 
;   but it handles nonsquare matrices
; CATEGORY:
;       Numerical
; CALLING SEQUENCE:
;   cc= lle_sm_mult(aa,bb)
; INPUTS:
;   aa - [N,M] sparse matrix input
;   bb - [M,P] sparse matrix input
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;   cc - [N,P] sparse matrix output (should be equal to idl's aa#bb) 
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
;   always at double precision
; PROCEDURE:
; MODIFICATION HISTORY:
;	  Blanton 2003-05-26 
;-
pro lle_sm_mult_compress,numlist,vallist

isort=sort(numlist)
numlist=numlist[isort]
vallist=vallist[isort]
iuniq=[0L,uniq(numlist)+1]
sumval=dblarr(n_elements(iuniq)-1)
for k=0L, n_elements(iuniq)-2 do begin
    sumval[k]=total(vallist[iuniq[k]:iuniq[k+1]-1L],/double)
endfor
numlist=numlist[iuniq[0:n_elements(iuniq)-2]]
vallist=sumval

end
;
function lle_sm_mult, a_in, b_in

a=(a_in)
b=lle_sm_transpose(b_in)

iinnera=[0,uniq(a.mindx)+1]
ninnera=n_elements(iinnera)-1L
mindx_innera=a.mindx[iinnera[0:ninnera-1L]]

jinnerb=[0,uniq(b.mindx)+1]
ninnerb=n_elements(jinnerb)-1L
mindx_innerb=b.mindx[jinnerb[0:ninnerb-1L]]

j=0L
count=0L
compresscount=30L
for i=0L,ninnera-1L do begin
    contloop=1
    while(j lt ninnerb and contloop) do begin
        contloop=mindx_innerb[j] lt mindx_innera[i]
        if(contloop) then j=j+1L
    endwhile 
    if(j lt ninnerb) then begin
        if(mindx_innerb[j] eq mindx_innera[i]) then begin
            tmp_outn=a.nindx[iinnera[i]:iinnera[i+1]-1L]
            tmp_outm=b.nindx[jinnerb[j]:jinnerb[j+1]-1L]
            tmp_num=replicate(1L,n_elements(tmp_outn))#tmp_outm*a.n+ $
              tmp_outn#replicate(1L,n_elements(tmp_outm))
            tmp_val=(replicate(1.,n_elements(tmp_outn))# $
                     b.vals[jinnerb[j]:jinnerb[j+1]-1L])* $
              (a.vals[iinnera[i]:iinnera[i+1]-1L]# $
               replicate(1.,n_elements(tmp_outm)))
            tmp_num=reform(tmp_num,n_elements(tmp_num))
            tmp_val=reform(tmp_val,n_elements(tmp_val))
            if(n_elements(numlist) eq 0) then begin
                numlist=tmp_num
                vallist=tmp_val
            endif else begin
                numlist=[numlist,tmp_num]
                vallist=[vallist,tmp_val]
            endelse
            if((count mod compresscount) eq 0) then $
              lle_sm_mult_compress,numlist,vallist
            count=count+1L
        endif
    endif
endfor
if(n_elements(numlist) eq 0) then stop
lle_sm_mult_compress,numlist,vallist

nindx=numlist mod a.n
mindx=numlist/a.n

c={n:a.n, $
   m:b.n, $
   nindx:nindx, $
   mindx:mindx, $
   vals:vallist}

return,c

end
