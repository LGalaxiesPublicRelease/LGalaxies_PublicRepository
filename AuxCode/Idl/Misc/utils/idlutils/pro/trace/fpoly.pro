function fpoly, x, ncoeff

  if N_params() LT 2 then begin
        print,'Syntax - result = FPOLY( x, ncoeff)
        return,0
  endif

  if ncoeff LT 1 then message, $
        'ERROR - Order of polynomial must be at least 1'
  N = N_elements(x)

  rr = (x[*]*0 + 1) # replicate(1,ncoeff)
  for i=1, ncoeff-1 do rr[*,i] = rr[*,i-1] * x

  return, rr
end

   
