;+
; NAME:
;   shuffle_indx
; PURPOSE:
;   yield an index to randomly rearrange an array 
; CALLING SEQUENCES:
;   indx= shuffle_indx(num [, seed=, num_sub=])
; INPUTS:
;   num - number of elements in an array
; OPTIONAL INPUTS:
;   seed - seed to pass to randomu
;   num_sub - only get the first num_sub elements of the rearranged array
; OUTPUTS:
;   indx - [num] or [num_sub] index of the shuffled rearrangement
; BUGS:
; Written MRB 2003-03-02
;-
function shuffle_indx, num, num_sub=num_sub, seed=seed

if(n_params() ne 1) then begin
    print, 'Usage: indx= shuffle_indx(num [, seed=, num_sub= ])'
    return, -1
endif

if(n_elements(num_sub) eq 0) then num_sub=num

indx=lindgen(num)
ran=randomu(seed,num_sub)
for i=0L, num_sub-1L do begin
    iswap=i+floor(ran[i]*double(num-i))
    tmp=indx[i]
    indx[i]=indx[iswap]
    indx[iswap]=tmp
endfor

return,indx[0:num_sub-1]

end
