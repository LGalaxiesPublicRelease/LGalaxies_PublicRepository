; Make r_V tables for SFD dustmap web site (for Peregrine)
; 2004-Mar-24  D. Finkbeiner

pro sdss_rv_table, source=source, zsource=zsource, fname=fname

  if NOT keyword_set(zsource) then zsource = 0
  if NOT keyword_set(source) then source = 'Galaxy'
  if NOT keyword_set(fname) then ilun = -1 else begin 
     openw, ilun, fname, /get_lun
  endelse
  
  rv = dindgen(41)/10.+2

  printf, ilun, 'Source: ', source, ' at redshfit z = ', zsource
  printf, ilun, '  R_V       u      g      r      i      z'
  for i=0, n_elements(rv)-1 do begin 
     aval = dust_sdssfilter(1.d ,source=source, rv=rv[i])
     printf, ilun, rv[i], aval/.0184, format='(F6.3,F9.3,4F7.3)'
  endfor 

  if ilun NE -1 then free_lun, ilun


  return
end
