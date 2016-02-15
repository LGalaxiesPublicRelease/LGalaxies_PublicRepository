function gauss_kernel, sigma, nsub=nsub, exp=exp, hpix=hpix, xl=xl

    if NOT keyword_set(nsub) then nsub = 5


    if NOT keyword_set(hpix) then begin
      if keyword_set(exp) then hpix = 10.0*sigma + 1 $
      else hpix = 5.0*sigma + 1
    endif

    x1 = 1.0+findgen(hpix)
    xl = [reverse(-x1),0.0,x1]
    sneaky = (findgen(nsub) - (nsub-1)/2.0)/nsub

    x = sneaky # replicate(1,2*hpix+1) + xl ## replicate(1,nsub)

    if keyword_set(exp) then expl = exp(-abs(x)/sigma) $
    else expl = exp(-x^2/(2.0*sigma^2))

    kernel = total(expl,1)/ total(expl)

    return, kernel
end



