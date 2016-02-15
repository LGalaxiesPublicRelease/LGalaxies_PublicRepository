
;
;  for now assume ideriv = 0
;
;
function bvalu, x, fullbkpt, coeff, ideriv=ideriv

     if NOT keyword_set(ideriv) then ideriv = 0L

     out = x*0.
     xs = sort(x)
     n = n_elements(coeff)
     k = n_elements(fullbkpt) - n
     nx = n_elements(x)
     work = fltarr(nx,k)
     left = fltarr(nx,k)
     right = fltarr(nx,k)
     t = fullbkpt

     i = intrv(x[xs], t, k)
     km1 = k - 1
     imk = i-k+1
     kmider = k - ideriv

     for j=0,k-1 do work[*,j] = coeff[imk+j]

     if ideriv GT 0 then begin
       for j=1,ideriv do begin
         kmj = k - j
         fkmj = kmj
         for jj=0,kmj-1 do begin 
            IHI = I + JJ
            IHMKMJ = IHI - KMJ
            WORK[*,JJ] = (WORK[*,JJ+1]-WORK[*,JJ])/(T[IHI]-T[IHMKMJ])*FKMJ
         endfor
       endfor
     endif

     if ideriv EQ km1 then begin
       out[xs] = work[*,0]
       return, out
     endif

     for j=0,k-1 do begin
       left[*,j] = t[i+j+1] - x[xs]
       right[*,j] = x[xs] - t[i-j] 
     endfor

     for j=ideriv+1,k - 1 do begin
       ilo = k-j-1
       for jj=0,k-j-1 do begin
         work[*,jj] = (work[*,jj+1] * right[*,ilo] + work[*,jj] * left[*,jj]) $
                        / (right[*,ilo] + left[*,jj])
         ilo = ilo - 1
       endfor
     endfor    

     out[xs] = work[*,0]
     return, out
end

