;+
; NAME:
;   djs_mosaic_rgb
; PURPOSE:
;   Make color (RGB) images from 3 FITS files, and overplot text and points
; CALLING SEQUENCE:
;   djs_mosaic_rgb, prefix
; INPUTS:
;   prefix        - start of filenames
; OPTIONAL KEYWORDS:
;   resizefactor  - factor by which to scale the JPG relative to the
;                   fits files; default 0.5
;   stretch       - factor by which to multiply the default RGB
;                   scales; default 1.0
;   rotation      - integer to pass to IDL rotate command; default 0
;   xtext         - X position for text [NTEXT]
;   ytext         - Y position for text [NTEXT]
;   text          - Text [NTEXT]
;   colortext     - Colors for text; dimensioned either [3] or [3,NTEXT],
;                   where [255,0,0] is red, [0,255,0] is green, etc.
;   xplot         - X position for points [NPOINTS]
;   yplot         - Y position for points [NPOINTS]
;   colorplot     - Colors for points; dimensioned either [3] or [3,NPOINTS],
;                   where [255,0,0] is red, [0,255,0] is green, etc.
;   _EXTRA        - Keywords to pass to XYOUTS and PLOT, such as
;                   CHARSIZE, PSYM, etc
; COMMENTS:
;   Stretch is constant in f_lambda
; EXAMPLES:
;   djs_mosaic_rgb, 'marla-001'
; BUGS:
;   Memory issues with asinh etc.
;   The current implementation is very slow if COLORTEXT or COLORPLOT
;     are 2-dimensional arrays.
; REVISION HISTORY:
;   2003-11-24  written - Hogg
;   2004-01-03  Modified (generalized) by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro djs_mosaic_rgb, prefix, resizefactor=resizefactor, rotation=rotation, $
 stretch=stretch, $
 xtext=xtext, ytext=ytext, text=text, colortext=colortext, $
 xplot=xplot, yplot=yplot, colorplot=colorplot, _EXTRA=extra

   if NOT keyword_set(resizefactor) then resizefactor = 0.50
   if NOT keyword_set(stretch) then stretch = 1.0
   if NOT keyword_set(prefix) then begin
      print, 'PREFIX must be specified'
      return
   endif
   if NOT keyword_set(rotation) then rotation = 0
   scales = [20.0,7.8,5.7,4.9,4.0]*stretch
   nonlinearity = 3.0
   origin = fltarr(3)+0.05
   filterlist = [2]
   if (keyword_set(colortext)) then begin
      if (size(colortext, /n_dimen) EQ 1) then ncolortext = 1 $
       else ncolortext = (size(colortext, /dimens))[1]
   endif
   nplot = n_elements(xplot)
   if (keyword_set(colorplot)) then begin
      if (size(colorplot, /n_dimen) EQ 1) then ncolorplot = 1 $
       else ncolorplot = (size(colorplot, /dimens))[1]
   endif

   for ifilt=0, n_elements(filterlist)-1 do begin
      filternum = filterlist[ifilt]

      Rfilter= filtername(filternum+1)
      Gfilter= filtername(filternum)
      Bfilter= filtername(filternum-1)

      Rfile= prefix+'-'+Rfilter+'.fits*'
      Gfile= prefix+'-'+Gfilter+'.fits*'
      Bfile= prefix+'-'+Bfilter+'.fits*'

      im= rotate(mrdfits(Rfile),rotation)
      tmp= size(im,/dimensions)
      nxorig= tmp[0]
      nyorig= tmp[1]
      nx= floor(nxorig*resizefactor)
      ny= floor(nyorig*resizefactor)

      RGBim= fltarr(nx,ny,3)
      RGBim[*,*,0]= congrid(temporary(im),nx,ny,/interp)
      RGBim[*,*,1]= congrid(rotate(mrdfits(Gfile),rotation),nx,ny,/interp)
      RGBim[*,*,2]= congrid(rotate(mrdfits(Bfile),rotation),nx,ny,/interp)

      scales1= scales[[filternum+1,filternum,filternum-1]]
      RGBim = nw_scale_rgb(temporary(RGBim),scales=scales1)
      RGBim = nw_arcsinh(temporary(RGBim),nonlinearity=nonlinearity)
      RGBim = nw_cut_to_box(temporary(RGBim),origin=origin)

      RGBim = nw_float_to_byte(temporary(RGBim))

      ;----------
      ; Create overlay images of text + points

      if (ifilt EQ 0 AND (keyword_set(text) OR nplot GT 0)) then begin
         set_plot, 'z'
         device, z_buffering=0, set_resolution=[nx,ny]
         charimg = bytarr(nx,ny,3)

         ; Add any text...
         if (keyword_set(text)) then begin
            ntext = n_elements(text)
            erase
            for i=0L, ntext-1 do begin
print,'Text ',i,ntext
               xyouts, xtext[i]*resizefactor, ytext[i]*resizefactor, $
                text[i], /device, _EXTRA=extra, color=1
               if (ncolortext GT 1 OR i EQ ntext-1) then begin
                  tmpimg = tvrd(0, 0, nx, ny)
                  if (keyword_set(colortext)) then $
                   color = colortext[*,i<(ncolortext-1)] $
                  else $
                   color = [255,255,255] ; default to white
                  for j=0, 2 do $
                   charimg[*,*,j] = charimg[*,*,j] + tmpimg * color[j]
                  erase
               endif
            endfor
         endif

         ; Add any points...
         if (nplot GT 0) then begin
            erase
            for i=0L, nplot-1 do begin
print,'Point ',i,nplot
	       plot, [xplot[i]*resizefactor], [yplot[i]*resizefactor], $
                _EXTRA=extra, color=1, $
                xrange=[0,nx-1], yrange=[0,ny-1], $
                xmargin=[0,0], ymargin=[0,0], xstyle=5, ystyle=5, /noerase
               if (ncolorplot GT 1 OR i EQ nplot-1) then begin
                  tmpimg = tvrd(0, 0, nx, ny)
                  if (keyword_set(colorplot)) then $
                   color = colorplot[*,i<(ncolorplot-1)] $
                  else $
                   color = [255,255,255] ; default to white
                  for j=0, 2 do $
                   charimg[*,*,j] = charimg[*,*,j] + tmpimg * color[j]
                  erase
               endif
            endfor
         endif

         ; We only replace pixels that are changed...
         set_plot, 'X'
         qmask = total(charimg NE 0, 3)
         ipix = where(qmask NE 0, npix)
      endif

      if (keyword_set(npix)) then begin
         for j=0, 2 do begin
            tmpimg1 = RGBim[*,*,j]
            tmpimg2 = charimg[*,*,j]
            tmpimg1[ipix] = (tmpimg1[ipix] + float(tmpimg2[ipix])) < 255B
            RGBim[*,*,j] = tmpimg1
         endfor
      endif

      write_jpeg, prefix+'-'+Rfilter+Gfilter+Bfilter+'.jpg', $
       RGBim,TRUE=3,QUALITY=75
   endfor

   return
end
;------------------------------------------------------------------------------
