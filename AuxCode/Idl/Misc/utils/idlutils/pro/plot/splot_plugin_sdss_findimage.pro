; plugin to call sdss_findimage from splot
; Finkbeiner

function splot_plugin_sdss_findimage, coord0, coord1

  info = sdss_findimage(coord0, coord1)
  if keyword_set(info) then begin 
     if n_elements(info) GT 20 then begin 
        loc_string = strjoin(strcompress(string(info[0:19].run)))+'...'
     endif else begin 
        loc_string = strjoin(strcompress(string(info.run)))
     endelse 
  endif else loc_string = 'no run'
  
  return, loc_string

end
