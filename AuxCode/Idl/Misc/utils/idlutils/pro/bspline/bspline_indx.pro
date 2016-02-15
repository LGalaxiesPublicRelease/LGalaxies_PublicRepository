pro bspline_indx, x, lowerbkpt, upperbkpt, lower, upper

      ; first find new lower bkpt

      nx = n_elements(x)

      lower = upper + 1
      if lower GE nx then return
      if x[lower] GT upperbkpt then return

      while x[lower] LT lowerbkpt do lower = lower + 1

      upper = lower 
      if (upper GE nx-1) then return

      while (x[upper+1] LE upperbkpt) do begin
           upper = upper + 1
           if (upper GE nx-1) then return
      endwhile
 
return
end
      
