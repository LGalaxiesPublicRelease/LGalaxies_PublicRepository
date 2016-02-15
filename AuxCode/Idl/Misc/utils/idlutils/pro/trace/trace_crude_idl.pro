;
;  Don't call this directly, used as a utility routine for trace_crude
;

pro trace_crude_idl, image, invvar, radius, xstart, ypass, $
                   xset, xerr, maxerr, maxshift, maxshift0

  ntrace = n_elements(ypass)
  ny = (size(image))[2]

  FOR itrace = 0,ntrace-1 DO BEGIN

;
;  Recenter INITIAL Row for all traces simultaneously
;
    xinit = xstart[itrace]
    iy = ypass[itrace]
    xfit = trace_fweight(image, xinit, iy, invvar=invvar, $
                            radius=radius, xerr=xfiterr, /idl)
  
    xshift = (((xfit - xinit) < maxshift0) > (-maxshift)) * $
              (xfiterr LT maxerr)
 
    xset[iy,itrace] = xinit + xshift
    xerr[iy,itrace] = xfiterr * (xfiterr LT maxerr)  + $
                        999.0 * (xfiterr GE maxerr)
 
;    /* LOOP FROM INITIAL (COL,ROW) NUMBER TO LARGER ROW NUMBERS */
    for iy=ypass[itrace]+1, ny-1 do begin
      xinit = xset[iy-1, itrace]
      xfit = trace_fweight(image, xinit, iy, invvar=invvar, $
                              radius=radius, xerr=xfiterr, /idl)
      
      xshift = (((xfit - xinit) < maxshift0) > (-maxshift)) * $
                (xfiterr LT maxerr)

      xset[iy,itrace] = xinit + xshift
      xerr[iy,itrace] = xfiterr * (xfiterr LT maxerr)  + $
                        999.0 * (xfiterr GE maxerr)
    endfor

;      /* LOOP FROM INITIAL (COL,ROW) NUMBER TO SMALLER ROW NUMBERS */
    for iy=ypass[itrace]-1, 0, -1 do begin
      xinit = xset[iy+1, itrace]
      xfit = trace_fweight(image, xinit, iy, invvar=invvar, $
                              radius=radius, xerr=xfiterr, /idl)
      
      xshift = (((xfit - xinit) < maxshift0) > (-maxshift)) * $
                (xfiterr LT maxerr)

      xset[iy,itrace] = xinit + xshift
      xerr[iy,itrace] = xfiterr * (xfiterr LT maxerr)  + $
                        999.0 * (xfiterr GE maxerr)
    endfor

  ENDFOR
        
  return
end

