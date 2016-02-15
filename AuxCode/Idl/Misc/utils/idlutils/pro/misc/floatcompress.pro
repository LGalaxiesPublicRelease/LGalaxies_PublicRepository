;----------------------------------------------------------------------------
;+
; NAME:
;   floatcompress
;
; PURPOSE:
;   Make floating-point data more compressible by trimming binary digits. 
;   The routine keeps the ndig most signifcant binary digits. 
;   If keyword nsig is passed, the algorithm rounds to the nearest
;   power of two less than nsig*sigma, where sigma is evaluated from
;   the passed array (with 5 sigma outlier rejection). 
;
; CALLING SEQUENCE:
;   out = floatcompress(data, ndig=ndig, nsig=nsig)
;
; INPUTS:
;   data       - input data (type float or double)
;                  WARNING: input data array is nuked to save memory
;                  on large arrays. 
;
; OPTIONAL KEYWORDS:
;   ndig       - number of binary significant digits to keep
;   nsig       - number of sigma at which to round data
;
; OUTPUTS:
;   out        - output data array with ndig significant binary digits
;                kept and the rest zeroed. 
;
; COMMENTS:
;   This function does not compress the data in an array, but fills
;   unnecessary digits of the IEEE floating point representation with
;   zeros.  This makes the data more compressible by standard
;   compression routines such as compress or gzip. 
; 
;   The default is to retain 10 binary digits instead of the usual 23
;   bits (or 52 bits for double precision), introducing a fractional
;   error strictly less than 1/1024).  This is adequate for most
;   astronomical images, and results in images that compress a factor
;   of 2-4 with gzip. 
;
; EXAMPLES:
;   image = readfits('map.fits')              ; read in FITS image
;   outimage = floatcompress(image,ndig=8)    ; keep 8 binary digits
;   writefits,'mapsmall.fits',outimage        ; write image
;   
;   Then from the UNIX shell
;   > gzip -8 mapsmall.fits
;
;   (level 8 gzip is slower but more effective than average this is
;       good for files that will be zipped once and unzipped many times)
; 
; PERFORMANCE:
;   On the typical maps of the ISM, gzip -8 compression factor is
;   ~2.1.  Mileage may vary.  For some images, a factor of 4-5 is possible.
;
; BUGS:
;   None known, but it is possible that there are floating point
;   values that are corrupted due to round off errors.  Results should
;   be double-checked.   
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   05-Jul-2000 Written by Doug Finkbeiner (UC Berkeley)
;   16-Sep-2000 Put in current format and commented -DPF
;   22-Sep-2000 Added nsig keyword
;   22-Jun-2002 Deals with Infs and NaNs - DPF
;-
;----------------------------------------------------------------------------;
FUNCTION floatcompress, data, ndig=ndig, nsig=nsig

; check data type
  IF (size(data, /type) NE 4) AND (size(data, /type) NE 5) THEN BEGIN 
     print, 'This routine is for use with FLOAT or DOUBLE types ONLY!'
     return, data
  ENDIF 
  
; set default ndig (number of binary significant digits)
  IF NOT keyword_set(ndig) THEN ndig = 10

; determine sigma with iterated sigma clipping
  IF keyword_set(nsig) THEN BEGIN  
     IF size(data, /n_dim) EQ 2 THEN BEGIN 
        sx = (size(data))[1]
        sy = (size(data))[2]

        IF ((sx/8)*8 EQ sx) AND ((sy/8)*8 EQ sy) THEN $
          sig = djsig(rebin(data, sx/8, sy/8), sigrej=3.0, maxiter= 100) ELSE $
          sig = djsig(data, sigrej=3.0, maxiter= 100)
     ENDIF ELSE BEGIN 
        sig = djsig(data, sigrej=3.0, maxiter= 100)
     ENDELSE 
    

     step = 2.^ceil(alog(sig*nsig)/alog(2))
     print, 'FLOATCOMPRESS: max roundoff error, sig: ', step/2, sig
     data = round(temporary(data)/step)*step
  ENDIF 

; replace zeros and nonfinite values with ones temporarily. 
  wzer = where((data EQ 0) OR finite(data) EQ 0, zcount)
  IF zcount NE 0 THEN BEGIN
     temp = data[wzer]
     data[wzer] = 1.
  ENDIF

; compute log base 2
  log2 = ceil(alog(abs(data))/alog(2.))  ; exponent part

  mant = round(temporary(data)/2.0^(log2-ndig))/(2.0^ndig) ; mantissa, truncated
  out = temporary(mant)*2.0^log2     ; multiply 2^exponent back in 

  IF zcount NE 0 THEN BEGIN
     out[wzer] = temp
  ENDIF
     

  return, out
END

