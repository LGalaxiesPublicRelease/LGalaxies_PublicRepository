pro imcontour, im, hdr, TYPE=type, PUTINFO=putinfo, XTITLE=xtitle,  $
      YTITLE=ytitle, SUBTITLE = subtitle, XDELTA = xdelta, YDELTA = ydelta, $
      ANONYMOUS_ = dummy_,_EXTRA = extra, XMID = xmid, YMID = ymid
;+
; NAME:
;       IMCONTOUR
; PURPOSE:
;       Make a contour plot labeled with astronomical coordinates.
; EXPLANATION:
;       The type of coordinate display is controlled by the keyword TYPE
;       Set TYPE=0 (default) to measure distances from the center of the image
;       (IMCONTOUR will decide whether the plotting units will be in
;       arc seconds, arc minutes, or degrees depending on image size.)
;       Set /TYPE for standard RA and Dec labeling
;
;       By using the /NODATA keyword, IMCONTOUR can also be used to simply
;       provide astronomical labeling of a previously displayed image.
; CALLING SEQUENCE
;       IMCONTOUR, im, hdr,[ /TYPE, /PUTINFO, XDELTA = , YDELTA =, _EXTRA = 
;                            XMID=, YMID= ]
;
; INPUTS:
;       IM - 2-dimensional image array
;       HDR - FITS header associated with IM, string array, must include
;               astrometry keywords.   IMCONTOUR will also look for the
;               OBJECT and IMAGE keywords, and print these if found and the 
;               PUTINFO keyword is set.
;
; OPTIONAL PLOTTING KEYWORDS:
;       /TYPE - the type of astronomical labeling to be displayed.   Either set
;               TYPE = 0 (default), distance to center of the image is
;               marked in units of Arc seconds, arc minutes, or degrees
;
;               TYPE = 1 astronomical labeling with Right ascension and 
;               declination.
;
;       /PUTINFO - If set, then IMCONTOUR will add information about the image
;               to the right of the contour plot.  Information includes image
;               name, object, image center, image center, contour levels, and
;               date plot was made
;
;       XDELTA, YDELTA - Integer scalars giving spacing of labels for TYPE=1.  
;               Default is to label every major tick (XDELTA=1) but if 
;               crowding occurs, then the user might wish to label every other
;               tick (XDELTA=2) or every third tick (XDELTA=3)
;
;       XMID, YMID - Scalars giving the X,Y position from which offset distances
;               will be measured when TYPE=0.   By default, offset distances 
;               are measured from the center of the image.
;
;       Any keyword accepted by CONTOUR may also be passed through IMCONTOUR
;       since IMCONTOUR uses the _EXTRA facility.     IMCONTOUR uses its own
;       defaults for the XTITLE, YTITLE XMINOR, YMINOR, and SUBTITLE keywords
;       but these may be overridden.
;
; NOTES:
;       (1) The contour plot will have the same dimensional ratio as the input
;           image array
;       (2) To contour a subimage, use HEXTRACT before calling IMCONTOUR
;       (3) Use the /NODATA keyword to simply provide astronomical labeling
;           of a previously displayed image.
;       (4) The IMCONTOUR display currently does not indicate the image 
;           rotation in any way, but only specifies coordinates along the 
;           edges of the image 
;
; EXAMPLE:
;       Overlay the contour of an image, im2, with FITS header, h2, on top
;       of the display of a different image, im1.   Use RA, Dec labeling, and
;       seven equally spaced contour levels.    The use of a program like
;       David Fanning's TVIMAGE  http://www.dfanning.com/programs/tvimage.pro
;       is suggested to properly overlay plotting and image coordinates.  The
;       /Keep_aspect_ratio keyword must be used.
;
;       IDL> tvimage,im1,/keep_aspect, position = pos
;       IDL> imcontour,im2,h2,nlevels=7,/Noerase,/TYPE,position = pos
;
; PROCEDURES USED:
;       CHECK_FITS, EXTAST, GETROT, TICPOS, TICLABEL, TIC_ONE, TICS, XYAD
;       CONS_RA(), CONS_DEC(), ADSTRING()
;
; REVISION HISTORY:
;       Written   W. Landsman   STX                    May, 1989
;       Fixed RA,Dec labeling  W. Landsman             November, 1991
;       Fix plottting keywords  W.Landsman             July, 1992
;       Recognize GSSS headers  W. Landsman            July, 1994
;       Removed Channel keyword for V4.0 compatibility June, 1995
;       Add _EXTRA CONTOUR plotting keywords  W. Landsman  August, 1995
;       Add XDELTA, YDELTA keywords  W. Landsman   November, 1995
;       Use SYSTIME() instead of !STIME                August, 1997
;       Remove obsolete !ERR system variable W. Landsman   May 2000 
;       Added XMID, YMID keywords to specify central position (default is still
;          center of image)  W. Landsman               March 2002     
;       Recognize Galactic coordinates, fix Levels display when /PUTINFO set
;           W. Landsman                May 2003
;       Correct conversion from seconds of RA to arcmin is 4 not 15.
;       	M. Perrin					July 2003
;       Fix integer truncation which appears with tiny images WL  July 2004
;       
;-
  On_error,2                                 ;Return to caller

  if N_params() LT 2 then begin             ;Sufficient parameters?
      print,'Syntax - imcontour, im, hdr, [ /TYPE, /PUTINFO, XDELTA=, YDELT= '
      print,'                               XMID=, YMID = ]'
      print,'         Any CONTOUR keyword is also accepted by IMCONTOUR'  
     return
  endif

  ;Make sure header appropriate to image
  check_fits, im, hdr, dimen, /NOTYPE, ERRMSG = errmsg    
  if errmsg NE '' then message,errmsg

; Set defaults if keywords not set

  if not keyword_set( TYPE ) then type = 0
  if not keyword_set( XDELTA ) then xdelta = 1
  if not keyword_set( YDELTA ) then ydelta = 1
 
  if not keyword_set(XMINOR) then $
       if !X.MINOR EQ 0 then xminor = 5 else xminor = !X.MINOR

  if not keyword_set(YMINOR) then $
       if !Y.MINOR EQ 0 then yminor = 5 else yminor = !Y.MINOR

  EXTAST, hdr, astr, noparams      ;Extract astrometry from header
  if noparams LT 0 then $                       ;Does astrometry exist?
      message,'FITS header does not contain astrometry'
  if strmid( astr.ctype[0], 5, 3) EQ 'GSS' then begin
        hdr1 = hdr
        gsss_STDAST, hdr1
        extast, hdr1, astr, noparams
  endif
  sexig = strmid(astr.ctype[0],0,4) EQ 'RA--'
 
; Adjust plotting window so that contour plot will have same dimensional 
; ratio as the image

  xlength = !D.X_VSIZE &  ylength = !D.Y_VSIZE
  xsize = fix( dimen[0] )  &   ysize = fix( dimen[1] )
  xsize1 = xsize-1 & ysize1 = ysize-1
  xratio = xsize / float(ysize)
  yratio = ysize / float(xsize)
  if N_elements(XMID) EQ 0 then xmid = xsize1/2.
  if N_elements(YMID) EQ 0 then ymid = ysize1/2.

  if ( ylength*xratio LT xlength ) then begin

    xmax = 0.15 + 0.8*ylength*xratio/xlength
    pos = [ 0.15, 0.15, xmax, 0.95 ]

  endif else begin

     xmax = 0.95
     pos = [ 0.15, 0.15, xmax, 0.15+ 0.8*xlength*yratio/ylength ]

  endelse

  if !X.TICKS GT 0 then xtics = abs(!X.TICKS) else xtics = 8
  if !Y.TICKS GT 0 then ytics = abs(!Y.TICKS) else ytics = 8

  pixx = float(xsize)/xtics            ;Number of X pixels between tic marks
  pixy = float(ysize)/ytics            ;Number of Y pixels between tic marks

  getrot,hdr,rot,cdelt               ;Get the rotation and plate scale

  xyad,hdr,xmid,ymid,ra_cen,dec_cen         ;Get coordinates of image center
  if sexig then ra_dec = adstring(ra_cen,dec_cen,1)       ;Make a nice string

; Determine tic positions and labels for the different type of contour plots

  if type NE 0 then begin                  ;RA and Dec labeling

     xedge = [ 0, xsize1, 0]          ;X pixel values of the four corners
     yedge = [ 0, 0, ysize1]          ;Y pixel values of the four corners

     xy2ad, xedge, yedge, astr, a, d
 
     pixx = float(xsize)/xtics          ;Number of X pixels between tic marks
     pixy = float(ysize)/ytics          ;Number of Y pixels between tic marks

; Find an even increment on each axis
     tics, a[0], a[1], xsize, pixx, raincr, RA=sexig  ;Find an even increment for RA
     tics, d[0], d[2], ysize, pixy, decincr    ;Find an even increment for Dec

; Find position of first tic on each axis
     tic_one, a[0], pixx, raincr, botmin, xtic1, RA= sexig  ;Position of first RA tic
     tic_one, d[0], pixy, decincr,leftmin,ytic1       ;Position of first Dec tic

     nx = fix( (xsize1-xtic1)/pixx )             ;Number of X tic marks
     ny = fix( (ysize1-ytic1)/pixy )             ;Number of Y tic marks

     if sexig then ra_grid = (botmin + findgen(nx+1)*raincr/4.) else $ 
                   ra_grid = (botmin + findgen(nx+1)*raincr/60.)
     dec_grid = (leftmin + findgen(ny+1)*decincr/60.)

     ticlabels, botmin, nx+1, raincr, xlab, RA=sexig, DELTA=xdelta
     ticlabels, leftmin, ny+1, decincr, ylab,DELTA=ydelta

     xpos = cons_ra( ra_grid,0,astr )     ;Line of constant RA
     ypos = cons_dec( dec_grid,0,astr)   ;Line of constant Dec

     if sexig then begin 
        xunits = 'Right Ascension'
        yunits = 'Declination'
     endif else begin
        xunits = 'Longitude'
        yunits = 'Latitude'
     endelse                          

  endif else begin ; label with distance from center.
     ticpos, xsize1*cdelt[0], xsize, pixx, incrx, xunits     
     numx = fix(xmid/pixx)              ;Number of ticks from left edge
     ticpos, ysize1*cdelt[1], ysize, pixy, incry, yunits
     numy = fix(ymid/pixy)             ;Number of ticks from bottom to center
     nx = numx + fix((xsize1-xmid)/pixx)    ;Total number of X ticks 
     ny = numy + fix((ysize1-ymid)/pixy)    ;Total number of Y ticks  
     xpos = xmid + (findgen(nx+1)-numx)*pixx
     ypos = ymid + (findgen(ny+1)-numy)*pixy
     xlab = string(indgen(nx+1)*incrx - incrx*numx,'(I3)')
     ylab = string(indgen(ny+1)*incry - incry*numy,'(I3)')
  
  endelse

; Get default values of XTITLE, YTITLE, TITLE and SUBTITLE

  if not keyword_set(PUTINFO) then putinfo = 0

  if N_elements(xtitle) EQ 0 then $
  if !X.TITLE eq '' then xtitle = xunits else xtitle = !X.TITLE

  if N_elements(ytitle) EQ 0 then $
      if !Y.TITLE eq '' then ytitle = yunits else ytitle = !Y.TITLE

  if (not keyword_set( SUBTITLE) ) and (putinfo LT 1) then $
      if sexig then $
      subtitle = 'Center:  R.A. '+ strmid(ra_dec,1,13)+'  Dec ' + $
               strmid(ra_dec,13,13) else $
     subtitle = 'Center:  Longitude '+ strtrim(string(ra_cen,'(f6.2)'),2) + $
                          ' Latitude ' + strtrim(string(dec_cen,'(f6.2)'),2)

  if (not keyword_set( SUBTITLE) ) then subtitle = !P.SUBTITLE
   
  contour,im, $
         XTICKS = nx, YTICKS = ny, POSITION=pos, XSTYLE=1, YSTYLE=1,$
         XTICKV = xpos, YTICKV = ypos, XTITLE=xtitle, YTITLE=ytitle, $
         XTICKNAME = xlab, YTICKNAME = ylab, SUBTITLE = subtitle, $
         XMINOR = xminor, YMINOR = yminor, _EXTRA = extra

;  Write info about the contour plot if desired

  if putinfo GE 1 then begin

     xmax = xmax + 0.01

     ypos = 0.92
     object = sxpar( hdr, 'OBJECT', Count = N_object )
     if N_object GT 0  then begin 
           xyouts, xmax, ypos, object, /NORM
           ypos = ypos-0.05
     endif

     name = sxpar( hdr, 'IMAGE', Count = N_image )
     if N_image GT 0 then begin 
           xyouts,xmax,ypos,name, /NORM
           ypos = ypos - 0.05
     endif

     xyouts, xmax, ypos,'Center:',/NORM
     ypos = ypos - 0.05
     if sexig then begin
     xyouts, xmax, ypos, 'R.A. '+ strmid(ra_dec,1,13),/NORM
     xyouts, xmax, ypos-0.05, 'Dec '+  strmid(ra_dec,13,13),/NORM
     endif else begin
     xyouts, xmax, ypos, 'Longitude: '+ strtrim(string(ra_cen,'(f6.2)'),2),/NORM
     xyouts, xmax, ypos-0.05,  $
             'Latitude: '+  strtrim(string(dec_cen,'(f6.2)'),2),/NORM
     endelse
     ypos = ypos - 0.1
     xyouts, xmax, ypos, 'Image Size', /NORM
     xyouts, xmax, ypos-0.05, 'X: ' + strtrim(xsize,2), /NORM
     xyouts, xmax, ypos-0.1, 'Y: ' + strtrim(ysize,2), /NORM
     xyouts, xmax, ypos- 0.15, strmid(systime(),4,20),/NORM
     xyouts, xmax, ypos - 0.2, 'Contour Levels:',/NORM

    sv = !D.NAME
    set_plot,'null'
    contour,im, _EXTRA = extra, PATH_INFO = info
    set_plot,sv

    ypos = ypos - 0.25
    val = info.value
    val = val[uniq(val,sort(val))]
     nlevels = N_elements(val)
     for i = 0,(nlevels < 7)-1 do $
          xyouts,xmax,ypos-0.05*i,string(i,'(i2)') + ':' + $
                              string(val[i]), /NORM

  endif
  
  return                                          
  end                                         
