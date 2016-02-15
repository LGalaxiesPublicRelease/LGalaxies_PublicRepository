pro find, image, x, y, flux, sharp, roundness, hmin, fwhm, roundlim, sharplim,$
                      PRINT = print, SILENT=silent
;+
; NAME:
;	FIND
; PURPOSE:
;	Find positive brightness perturbations (i.e stars) in an image 
; EXPLANATION:
;	Also returns centroids and shape parameters (roundness & sharpness).
;	Adapted from 1986 STSDAS version of DAOPHOT.
;
; CALLING SEQUENCE:
;	FIND, image, [ x, y, flux, sharp, round, hmin, fwhm, roundlim, sharplim 
;		PRINT= , /SILENT ]
;
; INPUTS:
;	image - 2 dimensional image array (integer or real) for which one
;		wishes to identify the stars present
;
; OPTIONAL INPUTS:
;	FIND will prompt for these parameters if not supplied
;
;	hmin -  Threshold intensity for a point source - should generally 
;		be 3 or 4 sigma above background
;	fwhm  - FWHM to be used in the convolve filter
;	sharplim - 2 element vector giving low and high cutoff for the
;		sharpness statistic (Default: [0.2,1.0] ).   Change this
;		default only if the stars have significantly larger or 
;		or smaller concentration than a Gaussian
;	roundlim - 2 element vector giving low and high cutoff for the
;		roundness statistic (Default: [-1.0,1.0] ).   Change this 
;		default only if the stars are significantly elongated.
;
; OPTIONAL INPUT KEYWORDS:
;	/SILENT - Normally, FIND will write out each star that meets all
;		selection criteria.   If the SILENT keyword is set and 
;		non-zero, then this printout is suppressed.
;	PRINT - if set and non-zero then FIND will also write its results to
;		a file find.prt.   Also one can specify a different output file 
;		name by setting PRINT = 'filename'.
;
; OPTIONAL OUTPUTS:
;	x - vector containing x position of all stars identified by FIND
;	y-  vector containing y position of all stars identified by FIND
;	flux - vector containing flux of identified stars as determined
;		by a Gaussian fit.  Fluxes are NOT converted to magnitudes.
;	sharp - vector containing sharpness statistic for identified stars
;	round - vector containing roundness statistic for identified stars
;
; NOTES:
;	(1) The sharpness statistic compares the central pixel to the mean of 
;       the surrounding pixels.   If this difference is greater than the 
;       originally estimated height of the Gaussian or less than 0.2 the height of the
;	Gaussian (for the default values of SHARPLIM) then the star will be
;	rejected. 
;
;       (2) More recent versions of FIND in DAOPHOT allow the possibility of
;       ignoring bad pixels.    Unfortunately, to implement this in IDL
;       would preclude the vectorization made possible with the CONVOL function
;       and would run extremely slowly.
; PROCEDURE CALLS:
;	GETOPT()
; REVISION HISTORY:
;	Written W. Landsman, STX  February, 1987
;	ROUND now an internal function in V3.1   W. Landsman July 1993
;	Change variable name DERIV to DERIVAT    W. Landsman Feb. 1996
;	Use /PRINT keyword instead of TEXTOUT    W. Landsman May  1996
;	Changed loop indices to type LONG       W. Landsman Aug. 1997
;	Converted to IDL V5.0   W. Landsman   September 1997
;       Replace DATATYPE() with size(/TNAME)   W. Landsman Nov. 2001
;       Fix problem when PRINT= filename   W. Landsman   October 2002
;       Fix problems with >32767 stars   D. Schlegel/W. Landsman Sep. 2004
;-
;
 On_error,2                         ;Return to caller
 compile_opt idl2

 npar   = N_params()
 if npar EQ 0 then begin
    print,'Syntax - FIND, image,' + $
          '[ x, y, flux, sharp, round, hmin, fwhm, roundlim, sharplim'
    print,'                      PRINT= , /SILENT ]'
    return
 endif

 maxbox = 13 	;Maximum size of convolution box in pixels 

; Get information about the input image 

 type = size(image)
 if ( type[0] NE 2 ) then message, $
     'ERROR - Image array (first parameter) must be 2 dimensional'
 n_x  = type[1] & n_y = type[2]
 message,  $
    'Input Image Size is '+strtrim(n_x,2) + ' by '+ strtrim(n_y,2),/INF

;Determine if hardcopy output is desired
 doprint = keyword_set( PRINT)
 if not keyword_set( SILENT ) then silent = 0
 if ( N_elements(fwhm) NE 1 ) then $
           read, 'Enter approximate FWHM: ', fwhm

 radius = 0.637*FWHM > 2.001             ;Radius is 1.5 sigma
 radsq = radius^2
 nhalf = fix(radius) < (maxbox-1)/2   	;
 nbox = 2*nhalf + 1	;# of pixels in side of convolution box 
 middle = nhalf          ;Index of central pixel

 lastro = n_x - nhalf
 lastcl = n_y - nhalf
 sigsq = ( fwhm/2.35482 )^2
 mask = bytarr( nbox, nbox )   ;Mask identifies valid pixels in convolution box 
 c = fltarr( nbox, nbox )      ;c will contain Gaussian convolution kernel

 dd = indgen(nbox-1) + 0.5 - middle	;Constants need to compute ROUND
 dd2 = dd^2
 w = 1. - 0.5*(abs(dd)-0.5) / (middle-.5)   
 ir = (nhalf-1) > 1

 row2 = (findgen(Nbox)-nhalf)^2

 for i = 0, nhalf do begin
	temp = row2 + i^2
	c[0,nhalf-i] = temp         
        c[0,nhalf+i] = temp                           
 endfor

 mask = fix(c LE radsq)     ;MASK is complementary to SKIP in Stetson's Fortran
 good = where( mask, pixels)  ;Value of c are now equal to distance to center

 c = c*mask               
 c[good] = exp(-0.5*c[good]/sigsq)	;Make c into a Gaussian kernel
 sumc = total(c)
 sumcsq = total(c^2) - sumc^2/pixels
 sumc = sumc/pixels
 c[good] = (c[good] - sumc)/sumcsq
 c1 = exp(-.5*row2/sigsq)
 sumc1 = total(c1)/nbox
 sumc1sq = total(c1^2) - sumc1
 c1 = (c1-sumc1)/sumc1sq
 sumc = total(w)                         ;Needed for centroid computation

 print,'RELATIVE ERROR computed from FWHM',sqrt(total(c[good]^2))
 if N_elements(hmin) NE 1 then read, $
    'Enter minimum value above background for threshold detection: ',hmin

 if N_elements(sharplim) NE 2 then begin
      print,'Enter low and high cutoffs, press [RETURN] for defaults:'
GETSHARP:   
      ans = ''
      read, 'Image Sharpness Statistic (DEFAULT = 0.2,1.0): ', ans   
      if ans EQ '' then sharplim = [0.2,1.0] else begin
         sharplim = getopt(ans,'F')
          if N_elements(sharplim) NE 2 then begin  
              message, 'ERROR - Expecting 2 scalar values',/CON
              goto, GETSHARP     
          endif
      endelse                                                      

GETROUND: 
  ans = ''
  read, 'Image Roundness Statistic [DEFAULT = -1.0,1.0]: ',ans
  if ans EQ '' then roundlim = [-1.,1.] else begin
      roundlim = getopt( ans, 'F' )
      if N_elements( roundlim ) NE 2 then begin
           message,'ERROR - Expecting 2 scalar values',/CON
           goto, GETROUND   
      endif
 endelse
 endif 

 message,'Beginning convolution of image', /INF

 h = convol(float(image),c)    ;Convolve image with kernel "c"

    h[0:nhalf-1,*] = 0 & h[n_x-nhalf:n_x-1,*] = 0
    h[*,0:nhalf-1] = 0 & h[*,n_y-nhalf:n_y-1] = 0

 message,'Finished convolution of image', /INF

 mask[middle,middle] = 0	;From now on we exclude the central pixel
 pixels = pixels -1      ;so the number of valid pixels is reduced by 1
 good = where(mask)      ;"good" identifies position of valid pixels
 xx= (good mod nbox) - middle	;x and y coordinate of valid pixels 
 yy = fix(good/nbox) - middle    ;relative to the center
 offset = yy*n_x + xx
SEARCH: 			    ;Threshold dependent search begins here

 index = where( h GE hmin, nfound)  ;Valid image pixels are greater than hmin
 if nfound EQ 0 then begin          ;Any maxima found?

    message,'ERROR - No maxima exceed input threshold of ' + $
             string(hmin,'(F9.1)'),/CON
    goto,FINISH    

 endif

 for i= 0L, pixels-1 do begin                             

	stars = where (h[index] GE h[index+offset[i]], nfound)
        if nfound LT 0 then begin  ;Do valid local maxima exist?
             message,'ERROR - No maxima exceed input threshold of ' + $
                     string(hmin,'(F9.1)'),/CON
             goto,FINISH  
        endif
	index = index[stars]

 endfor 
 
 ix = index mod n_x              ;X index of local maxima
 iy = index/n_x                  ;Y index of local maxima
 ngood = N_elements(index)       
 message,strtrim(ngood,2)+' local maxima located above threshold',/INF

 nstar = 0L       	;NSTAR counts all stars meeting selection criteria
 badround = 0L & badsharp=0L  &  badcntrd=0L
 if (npar GE 2) or (doprint) then begin 	;Create output X and Y arrays? 
  	x = fltarr(ngood) & y = x
 endif

 if (npar GE 4) or (doprint) then begin   ;Create output flux,sharpness arrays?
 	flux = x & sharp = x & roundness = x
 endif

 if doprint then begin	;Create output file?

         if ( size(print,/TNAME) NE 'STRING' ) then file = 'find.prt' $
                                         else file = print
         message,'Results will be written to a file ' + file,/INF
         openw,lun,file,/GET_LUN
	printf,lun,' Program: FIND '+ systime()
	printf,lun,format='(/A,F7.1)',' Threshold above background:',hmin
	printf,lun,' Approximate FWHM:',fwhm
	printf,lun,format='(2(A,F6.2))',' Sharpness Limits: Low', $
                sharplim[0], '  High',sharplim[1]
	printf,lun,format='(2(A,F6.2))',' Roundness Limits: Low', $
                roundlim[0],'  High',roundlim[1]
	printf,lun,format='(/A,i6)',' No of sources above threshold',ngood

 endif                      

 if not SILENT then $
  print,format='(/8x,a)','     STAR      X      Y     FLUX     SHARP    ROUND'

;  Loop over star positions; compute statistics

 for i = 0L,ngood-1 do begin   
     temp = float(image[ix[i]-nhalf:ix[i]+nhalf,iy[i]-nhalf:iy[i]+nhalf])
     d = h[ix[i],iy[i]]                  ;"d" is actual pixel intensity        

;  Compute Sharpness statistic

     sharp1 = (temp[middle,middle] - (total(mask*temp))/pixels)/d
     if ( sharp1 LT sharplim[0] ) or ( sharp1 GT sharplim[1] ) then begin
	badsharp = badsharp + 1
	goto, REJECT             ;Does not meet sharpness criteria
     endif

;   Compute Roundness statistic

     dx = total( total(temp,2)*c1)   
     dy = total( total(temp,1)*c1)
     if (dx LE 0) or (dy LE 0) then begin
         badround = badround + 1
	 goto, REJECT           ;Cannot compute roundness
     endif

     around = 2*(dx-dy) / ( dx + dy )    ;Roundness statistic
     if ( around LT roundlim[0] ) or ( around GT roundlim[1] ) then begin
	badround = badround + 1
	goto,REJECT           ;Does not meet roundness criteria
     endif

; Find X centroid

     derivat = shift(temp,-1,0) - temp
     derivat = total( derivat[0:nbox-2,middle-ir:middle+ir],2)
     sumd = total(w*derivat)
     sumxd = total(w*dd*derivat)
     sumxsq = total(w*dd2) 

     if ( sumxd GE 0. ) then begin
	badcntrd = badcntrd + 1
	goto,REJECT           ;Cannot compute X centroid
     endif

     dx =sumxsq*sumd/(sumc*sumxd)
     if abs(dx) GT nhalf then begin
      	 badcntrd = badcntrd + 1
	 goto,REJECT           ;X centroid too far from local X maxima
     endif

     xcen = ix[i]-dx               ;Convert back to big image coordinates

; Find Y centroid                 

     derivat = shift(temp,0,-1) - temp 
     derivat = total( derivat[middle-ir:middle+ir,0:nbox-2], 1 )
     sumd = total( w*derivat )
     sumxd = total( w*dd*derivat )
     sumxsq = total( w*dd2 )
     if (sumxd GE 0) then begin
	  badcntrd = badcntrd + 1
	  goto, REJECT  
     endif

     dy = sumxsq*sumd/(sumc*sumxd)
     if ( abs(dy) GT nhalf ) then begin
	badcntrd = badcntrd + 1
	goto,REJECT 
     endif
     
     ycen = iy[i] - dy

;  This star has met all selection criteria.  Print out and save results

   if not SILENT then $
      print,FORM = '(12x,i5,2f7.1,f9.1,2f9.2)', $ 
            nstar, xcen, ycen, d, sharp1, around

   if (npar GE 2) or (doprint) then begin
              x[nstar] = xcen & y[nstar] = ycen
   endif

   if ( npar GE 4 ) or (doprint) then begin
	flux[nstar] = d & sharp[nstar] = sharp1 & roundness[nstar] = around
   endif
   
   nstar = nstar+1

REJECT: 
  
 endfor

 nstar = nstar-1		;NSTAR is now the index of last star found

 if doprint then begin
  printf,lun,' No. of sources rejected by SHARPNESS criteria',badsharp
  printf,lun,' No. of sources rejected by ROUNDNESS criteria',badround
  printf,lun,' No. of sources rejected by CENTROID  criteria',badcntrd
 endif
 
  print,' No. of sources rejected by SHARPNESS criteria',badsharp
  print,' No. of sources rejected by ROUNDNESS criteria',badround
  print,' No. of sources rejected by CENTROID  criteria',badcntrd

  if nstar LT 0 then return               ;Any stars found?

  if (npar GE 2) or (doprint) then begin
	x=x[0:nstar]  & y = y[0:nstar]
  endif

  if (npar GE 4) or (doprint) then begin
	flux= flux[0:nstar] & sharp=sharp[0:nstar]  
        roundness = roundness[0:nstar]
  endif

 if doprint then begin                
   printf,lun, $
      format = '(/8x,a)','     STAR       X       Y     FLUX     SHARP    ROUND'
	for i = 0L, nstar do $
	   printf,lun,format='(12x,i5,2f8.2,f9.1,2f9.2)', $
	              i+1, x[i], y[i], flux[i], sharp[i], roundness[i]
        free_lun, lun
 endif

FINISH:

 if SILENT then return

 print,form='(A,F8.1)',' Threshold above background for this pass was',hmin
 ans = ''
 read,'Enter new threshold or [RETURN] to exit: ',ans
 ans = getopt(ans,'F')              
 if ans GT 0. then begin
       hmin = ans
       goto, SEARCH   
 endif

 return                                      
 end
