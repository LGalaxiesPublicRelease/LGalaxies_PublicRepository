
;
;	This is a slow version of slatec_bvalu, needs to be sped
;       up before using as full replacement.
;
;
function bvalu_slow, x, fullbkpt, coeff

      n = n_elements(coeff)
      k = n_elements(fullbkpt) - n
      nx = n_elements(x)
      indx = intrv(x, fullbkpt, k)

      bf1 = bsplvn(fullbkpt, k, x, indx[p])
      answer = x*0.0

      spot = lindgen(k) - k + 1
      for i=k-1L,n-1 do begin
        inside = where(indx EQ i)
        if inside[0] NE -1 then answer[inside] = bf1[inside,*] # coeff[spot+i]
      endfor

      return, answer
end

