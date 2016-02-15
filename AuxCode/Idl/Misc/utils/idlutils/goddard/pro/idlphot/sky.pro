pro sky,image,skymode,skysig, SILENT=silent, CIRCLERAD = circlerad, $
           ReadNoise = ReadNoise, Highbad = highbad, NAN = nan
;+
; NAME:
;       SKY
; PURPOSE:
;       Determine the sky level in an image using the the procedure MMM
; EXPLANATION:
;       Approximately 10000 uniformly spaced pixels are selected for the
;       computation.  Adapted from the DAOPHOT routine of the same name.
;
; CALLING SEQUENCE:
;       SKY, image, [ skymode, skysig, HIGHBAD= , READNOISE= ,/SILENT 
;                                      /NAN, CIRCLERAD= ]
; INPUTS:
;       IMAGE - One or two dimensional array
;
; OPTIONAL OUTPUT ARRAYS:
;       SKYMODE - Scalar, giving the mode of the sky pixel values of the 
;               array IMAGE, as determined by the procedure MMM.
;       SKYSIG -  Scalar, giving standard deviation of sky brightness
;
; INPUT KEYWORD PARAMETERS:
;	CIRCLERAD - Use this keyword to have SKY only select pixels within
;		specified pixel radius of the center of the image.  If 
;		CIRCLERAD =1, then the radius is set equal to half the image
;		width.   Can only be used with square images.
;       HIGHBAD - scalar value of the (lowest) "bad" pixel level (e.g. cosmic 
;                rays or saturated pixels) If not supplied, then there is 
;                assumed to be no high bad pixels.
;       /NAN - This keyword must be set to  ignore NaN values when computing 
;              the sky.
;              Note that the CIRCLERAD, HIGHBAD and /NAN are not exclusive, e.g.
;               one can set both /NAN and CIRCLERAD
;       /SILENT - If this keyword is supplied and non-zero, then SKY will not
;               display the sky value and sigma at the terminal
;
;       READNOISE - Scalar giving the read noise (or minimum noise for any 
;                pixel).     Normally, MMM determines the (robust) median by 
;                averaging the central 20% of the sky values.     In some cases
;                where the noise is low, and pixel values are quantized a
;                larger fraction may be needed.    By supplying the optional
;                read noise parameter, MMM is better able to adjust the
;                fraction of pixels used to determine the median. 
;
; PROCEDURE:
;       A grid of points, not exceeding 10000 in number, is extracted
;       from the image array.  The mode of these pixel values is determined
;       by the procedure MMM.   In a 2-d array the grid is staggered in each
;       row to avoid emphasizing possible bad columns
;
; PROCEDURE CALLS:
;       MMM, DIST_CIRCLE
; REVISION HISTORY:
;       Written, W. Landsman   STX Co.            September, 1987     
;       Changed INDGEN to LINDGEN                 January, 1994
;       Fixed display of # of points used         March, 1994
;       Stagger beginning pixel in each row, added NSKY, READNOISE, HIGHBAD
;          W. Landsman        June 2004
;      Adjustments for unbiased sampling  W. Landsman June 2004
;      Added /NAN keyword, put back CIRCLERAD keyword W. Landsman July 2004
;-
; On_error,2              ;Return to caller
 maxsky = 10000          ;Maximum # of pixels to be used in sky calculation

 if N_params() eq 0 then begin
        print,'Syntax - sky, image, [ skymode, skysig , HIGHBAD= '
        print, '                    READNOISE = , /NAN, CIRCLERAD = , /SILENT ]'
        return
 endif

 checkbad = (N_elements(highbad) GT 0) or keyword_set(circlerad) or $
              keyword_set(nan)                          
 s = size(image)      
 nrow = s[1]
 if s[0] EQ 1 then ncol = 1 else begin                      
    if s[0] NE 2 then message, $
          'ERROR - Input array (first parameter) must be 1 or 2 dimensional'
    ncol = s[2]
 endelse
 if keyword_set(circlerad) then if ncol ne nrow then message, $
       'ERROR - The CIRCLERAD keyword only applies to a 2-d square array'
        
 if checkbad then begin 
          mask = replicate(1b, nrow, ncol)
          if N_elements(highbad) GT 0 then mask = mask and (image LT highbad)
          if keyword_set(nan) then mask = mask and finite(image)
          if keyword_set(circlerad) then begin
                  if circlerad EQ 1 then rad = nrow/2 else rad = long(circlerad)
                  dist_circle,drad, nrow
                  mask = mask and (temporary(drad) LT rad)
           endif
          npts = total(mask)
 endif else  npts = N_elements(image)

 skyvec = fltarr(maxsky)
     istep = npts/maxsky +1
    nstep = (nrow/istep)
 
    jj = 0
    index0 = istep*lindgen(nstep) 
    if nstep GT 1 then begin 
          i0 = (nrow-1 - max(index0)  - istep)/2 > 0  ;Adjust margin for symmetry
          index0  = index0 + i0
    endif

; The beginning index in each row is staggered to avoid emphasizing possible
; bad columns
    for i=0, Ncol-1 do begin
        index  = index0 + (i mod istep)  
        row = image[*,i]
        if checkbad then begin         
            g = where(mask[*,i],ng)
            case ng of 
            0: goto, Done
            Nrow: 
            else: row = row[g]
            endcase
          endif else ng = nrow
          imax = value_locate( index, ng-1) > 0
          ix = index[0:imax] < (ng-1)
          skyvec[jj] = row[ix]
          jj = jj + imax + 1
 DONE:
  endfor    
  skyvec = skyvec[0:jj-1] 
 

 MMM, skyvec, skymode, skysig, readnoise = readnoise, highbad = highbad, $
      nsky = nsky

 skymode = float(skymode)  &  skysig = float(skysig)
 if not keyword_set(SILENT) then begin
        print,'Number of points used to find sky = ',nsky
        print,'Approximate sky value for this frame = ',skymode
        print,'Standard deviation of sky brightness = ',skysig
 endif

 return
 end
