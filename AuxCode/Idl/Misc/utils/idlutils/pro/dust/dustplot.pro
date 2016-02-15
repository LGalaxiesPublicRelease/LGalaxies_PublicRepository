;------------------------------------------------------------------------------
;+
; NAME:
;   dustplot
;
; PURPOSE:
;   Make a PostScript plot of the dust maps in a rectalinear projection.
;
; CALLING SEQUENCE:
;   dustplot
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   latrange:   Latitude range; default to [-30,30]
;   lonrange:   Longitude range; default to [-30,30]
;   rangestr:   3-element string with plot limits and units;
;               default to ['0.0', '1.0', 'A(B)']
;               Valid entries for RANGESTR[2] are:
;                  A(B)
;                  A(I)
;                  E(B-V)
;                  I100
;                  T
;   csys:       Coordinate system:
;               'gal': Galactic coordinates (default)
;               'equ1950': Equatorial coordinates, epoch 1995
;               'equ2000': Equatorial coordinates, epoch 2000
;               'ecl': Ecliptic coordinates
;   tspace:     Tick spacing for grid overlay; default to 5 deg
;   pixpdeg:    Pixels per degree; default to 500 pixels across the image
;   filename:   Output file name; default to "test.ps"
;   encap:      If set, then produce encapsulated PostScript
;   ctnum:      Color table number; default to 23 (Purple-Red + Stripes)
;   invert:     Invert color map
;   colorbar:   Plot color bar on bottom of page
;   nonames:    If set, then disable the title string on the top,
;               and our names on the bottom
;
; OUTPUTS:
;
; PROCEDURES CALLED:
;   adstring()
;   djs_laxisgen()
;   dust_getval()
;   euler
;   jprecess
;
; REVISION HISTORY:
;   Written D. Schlegel, 18 June 1999, Princeton
;-
;------------------------------------------------------------------------------
pro dustplot, lrange, brange, rangestr=rangestr, csys=csys, $
 tspace=tspace, pixpdeg=pixpdeg, $
 filename=filename, encap=encap, ctnum=ctnum, invert=invert, $
 logplot=logplot, colorbar=colorbar, nonames=nonames

   if (NOT keyword_set(filename)) then filename = 'test.ps'
   if (size(encap,/tname) EQ 'UNDEFINED') then encap = 0
   if (size(ctnum,/tname) EQ 'UNDEFINED') then ctnum = 23
   if (size(invert,/tname) EQ 'UNDEFINED') then invert = 0
   if (size(logplot,/tname) EQ 'UNDEFINED') then logplot = 0
   color = 1
   if (size(colorbar,/tname) EQ 'UNDEFINED') then colorbar = 1

   if (NOT keyword_set(lrange)) then lrange = [-10, 10]
   if (NOT keyword_set(brange)) then brange = [-10, 10]
   if (NOT keyword_set(csys)) then csys = 'gal'
   if (NOT keyword_set(tspace)) then tspace = 5
   if (NOT keyword_set(pixpdeg)) then begin
      pixpdeg = 500.0 / min([ abs(lrange[1]-lrange[0]), $
       abs(brange[1]-brange[0]) ])
   endif

;lrange = [-5, 5]
;brange = [-14, -4]
;tspace = 1.0 ; tick spacing in degrees
;   rangestr = ['0.0', '2.0', 'A(B)']
;   rangestr = ['0.0', '3.0', 'A(I)']
   if (NOT keyword_set(rangestr)) then $
    rangestr = ['0.0', '1.0', 'A(B)']

   ; Create maps of latitude, longitude in the coordinate system specified
   ; by CSYS.
   nx = abs(lrange[1] - lrange[0]) * pixpdeg
   ny = abs(brange[1] - brange[0]) * pixpdeg
   lonmap = djs_laxisgen([nx+1,ny+1], iaxis=0) * float(lrange[1]-lrange[0])/nx $
    + lrange[0]
   latmap = djs_laxisgen([nx+1,ny+1], iaxis=1) * float(brange[1]-brange[0])/ny $
    + brange[0]

   ; Coordinate transformation
   case csys of
   'equ1950': $
      begin
         cstring = 'B1950'
         jprecess, lonmap, latmap, ltemp, btemp  ; B1950 -> J2000
lonmap = 0
latmap = 0
         euler, ltemp, btemp, lmap, bmap, 1  ; J2000 -> Galactic
ltemp = 0
btemp = 0
      end
   'equ2000': $
      begin
         cstring = 'J2000'
         euler, lonmap, latmap, lmap, bmap, 1  ; J2000 -> Galactic
lonmap = 0
latmap = 0
      end
   'ecl': $
       begin
         cstring = 'ecliptic'
         euler, lonmap, latmap, lmap, bmap, 5  ; Ecliptic -> Galactic
lonmap = 0
latmap = 0
      end
   'gal': $
      begin
         cstring = 'Galactic'
         lmap = lonmap
         bmap = latmap
      end
   else: message, 'Unsupported CSYS'
   endcase

   ; Read the dust maps and convert to A(B), A(I), E(B-V), I100, or T.

   case rangestr[2] of
      'A(B)': image = 4.325 * dust_getval(lmap, bmap, /interp, /noloop)
      'A(I)': image = 1.940 * dust_getval(lmap, bmap, /interp, /noloop)
      'E(B-V)': image = dust_getval(lmap, bmap, /interp, /noloop)
      'I100': image = dust_getval(lmap, bmap, /interp, /noloop, map='I100')
      'T': image = dust_getval(lmap, bmap, /interp, /noloop, map='T')
      else: message, 'Unknown RANGESTR'
   endcase

   ; Set up postscript device
   set_plot, 'ps'
   xoffs = 0.5
   yoffs = 1.25
   scale = 1
   xsize=7.5
   ysize=8.5
   device, file=filename, bits=8, xsize=xsize, ysize=ysize, $
    xoffs=xoffs, yoffs=yoffs, /inch, encapsul=encap, color=color, scale=scale
   loadct, 0

   ; Set up "data" coordinate system
   plot, [0,xsize], [0,ysize], /nodata, xmargin=[0,0], ymargin=[0,0], $
    xstyle=5, ystyle=5

   ; Print title (caption)
   !p.font = 0
   cs = 1.5
   caption = 'SFD Dust Map - ' + cstring
   if (NOT keyword_set(nonames)) then begin
      xyouts, 0.25,8.10, caption, charsize=2.25, /data
      xyouts, 0.25,-0.5, 'Schlegel, Finkbeiner, & Davis (1998)', charsize=cs
   endif

   ; Load color table and reset stretch
   loadct, ctnum

  if (colorbar) then begin
      grayline = bindgen(256,8)
      if (keyword_set(invert)) then grayline = 255b-grayline
      tv, grayline, .25, .25, xsize=7, ysize=0.25, /inch
      oplot, [.25, 7.25, 7.25, .25, .25], [.25, .25, .5, .5, .25]
      for i=1, 9 do begin
         xtmp = 0.25 + 7.00*i/10
         oplot, [xtmp,xtmp], [.25, .5]
      endfor

      xyouts, 0.25, 0., rangestr[0], /data, charsize=cs, align=0.0
      xyouts, 7.25, 0., rangestr[1], /data, charsize=cs, align=1.0

      unitstr = rangestr[2]
      if (keyword_set(logplot)) then unitstr = 'log '+unitstr
      xyouts, 3.75, 0., unitstr, /data, charsize=cs, align=0.5
   endif

   ; Make byte array
   minval = float(rangestr[0])
   maxval = float(rangestr[1])
   if (keyword_set(logplot)) then begin
      bimage = bytscl(alog(image), min=alog(minval), max=alog(maxval))
   endif else begin
      bimage = bytscl(image, min=minval, max=maxval)
   endelse
   if (keyword_set(invert)) then bimage = 255b-bimage

   ; Display byte array
   tv, bimage, 0.75, 0.75, xsize=6.5, ysize=7, /inch

   ; Label only every ILABELX, ILABELY grid lines in X and Y; at most 11 labels
   nlabelx = fix(abs(lrange[1]-lrange[0])/tspace) + 1
   nlabely = fix(abs(brange[1]-brange[0])/tspace) + 1
   ilabelx = fix( (nlabelx+10) / 11 )
   ilabely = fix( (nlabely+10) / 11 )

   ; Plot grid
   lpad1 = (0.75/6.5)*(lrange[1]-lrange[0]) ; Can be negative if RA descends
   lpad2 = (0.25/6.5)*(lrange[1]-lrange[0]) ; Can be negative if RA descends
   bpad = 0.5 * (1.5/7.0)*abs(brange[1]-brange[0])
   plot, [lrange[0]-lpad1,lrange[1]+lpad2], [brange[0]-bpad,brange[1]+bpad], $
    /nodata, xstyle=5, ystyle=5, xmargin=[0,0], ymargin=[0,0], /noerase, /data
;contour, bimage, levels=findgen(11)*25.5, /noerase, $
; position=[lrange[0],brange[0],lrange[1],brange[1]]
   tspacex = abs(tspace)
   if (lrange[1] LT lrange[0]) then tspacex = -tspacex
   for ix=0, nlabelx-1 do begin
      ltemp = lrange[0] + ix * tspacex
      ; The following hack is to fix an IDL bug!!??
      if (lrange[0] GT lrange[0]) then lthis = ltemp $
       else lthis = lrange[1] + lpad2-lpad1 + ix * abs(tspacex)
      oplot, [lthis,lthis], brange
      if (ix MOD ilabelx EQ 0) then begin
         if (csys EQ 'equ1950' OR csys EQ 'equ2000') then begin
            stemp = strmid( adstring(ltemp, 0.0), 1, 8)
         endif else begin
            stemp = strtrim(string(ltemp),2)
         endelse
         xyouts, lthis, brange[0]-0.02*(brange[1]-brange[0]), stemp, $
          /data, charsize=1.0, align=1.0
      endif
   endfor
   for iy=0, nlabely-1 do begin
      btemp = min(brange) + iy * abs(tspace)
      ; The following hack is to fix an IDL bug!!??
      if (lrange[0] GT lrange[0]) then lthis = lrange $
       else lthis = lrange+lpad2-lpad1
      oplot, lthis, [btemp,btemp]
      if (iy MOD ilabely EQ 0) then begin
         if (csys EQ 'equ1950' OR csys EQ 'equ2000') then begin
            stemp = strmid( adstring(0.0, btemp), 13, 9)
         endif else begin
            stemp = strtrim(string(btemp),2)
         endelse
         ; The following hack is to fix an IDL bug!!??
         if (lrange[0] GT lrange[0]) then $
          lthis = lrange[0]-0.005*(lrange[1]-lrange[0]) $
         else $
          lthis = lrange[1] + lpad2-lpad1 -0.005*(lrange[1]-lrange[0])
         xyouts, lthis, btemp, stemp, $
          /data, charsize=1.0, align=1.0
      endif
   endfor

   device, /close
   set_plot, 'X'

  return
end
;------------------------------------------------------------------------------
