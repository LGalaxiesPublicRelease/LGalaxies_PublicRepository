PRO medarr, inarr, outarr, mask, output_mask
;+
; NAME:
;       MEDARR
; PURPOSE:
;       Compute the median at each pixel across a set of 2-d images
; EXPLANATION:
;       Each pixel in the output array contains  the median of the 
;       corresponding pixels in the input arrays.   Useful, for example to 
;       combine a stack of CCD images, while removing cosmic ray hits.
;
;       This routine became partially obsolete in V5.6 with the introduction
;       of the DIMENSION keyword to the intrinsic MEDIAN() function.   However,
;       it is  still useful if a input mask is needed (though it is much 
;       faster to set invalid pixels to NaN values.)
; CALLING SEQUENCE:
;       MEDARR, inarr, outarr, [ mask, output_mask ]
; INPUTS:
;       inarr  -- A three dimensional array containing the input arrays to 
;                 combine together.  Each of the input arrays must be two 
;                 dimensional and must have the same dimensions.  These arrays
;                 should then be stacked together into a single 3-D array,
;                 creating INARR.
;
; OPTIONAL INPUT:
;       mask   -- Same structure as inarr, byte array with 1b where
;                 pixels are to be included, 0b where they are to be
;                 excluded.    For floating point images, it is much faster to 
;                 set masked pixels in inarr equal to !VALUES.F_NAN (see below),
;                 rather than use the mask parameter.
;                
; OUTPUTS:
;       outarr -- The output array.  It will have dimensions equal to the
;                 first two dimensions of the input array.
;
; OPTIONAL OUTPUT:
;       output_mask -- Same structure as outarr, byte array with 1b
;                      pixels are valid, 0b where all the input pixels
;                      have been masked out.
; RESTRICTIONS:
;        Prior to V5.6, this procedure was *SLOW* because it had to loop over 
;        each pixel of the image.   See notes below about an alternative with 
;        CALL_EXTERNAL.
;
; EXAMPLE:
;       Suppose one wants to combine three floating point 1024 x 1024 bias 
;       frames which have been read into the IDL variables im1,im2,im3
;
;       IDL> bigim = fltarr(1024,1024,3)        ;Create big array to hold images
;       IDL> bigim(0,0,0) = im1 & bigim(0,0,1) = im2 & bigim(0,0,2) = im2  
;       IDL> medarr, bigim, avgbias
;
;       The variable avgbias will be the desired 1024x 1024 float image.
; PROCEDURE:
;       A scalar median function over the third dimension is looped over 
;       each pixel of the first two dimensions.   The /EVEN keyword is used
;       with MEDIAN (which averages the two middle values), since this avoids 
;       biasing the output for an even number of images.
;
;       Any values set to NAN (not a number) are ignored when computing the
;       median.    If all values for a pixel location are NAN, then the median
;       is also returned as NAN.
;
;       MEDARR is also available as a C procedure linked to IDL via
;       CALL_EXTERNAL (but without the mask parameter).   The callable C 
;       version is 2-3 times faster for large  (~ 500 x 500 x 7) images.   
;       Contact W. Landsman (landsman@mpb.gsfc.nasa.gov) for the C program
; MODIFICATION HISTORY:
;       Written by Michael R. Greason, STX, 12 June 1990.
;       Don't use MEDIAN function for even number of images.
;          W. Landsman Sep 1996
;       Mask added.  RS Hill, HSTX, 13 Mar. 1997
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Use /EVEN keyword to MEDIAN    W. Landsman  September 1997
;       Rearranged code for faster execution   W. Landsman January 1998
;       Faster execution for odd number of images   W. Landsman July 2000
;       V5.4 fix for change in SIZE() definition of undefined variable 
;                W. Landsman/E. Young   May 2001
;       Use MEDIAN(/DIMEN) for V5.6 or later   W. Landsman   November 2002
;       Use keyword_set() instead of ARG_present() to test for presence of mask
;           parameter  D. Hanish/W. Landsman   June 2003
;       Use MEDIAN(/EVEN) when mask not set V5.6 or later W. Landsman Feb 2004
;-
 On_error,2
;                       Check parameters.

 if N_params() LT 2 then begin                  ; # parameters.
        print, "Syntax -  MEDARR, inputarr, outputarr [, maskarr, output_mask]"
        return
 endif
 
 s = size(inarr)
 if s[0] NE 3 then $                    ; Input array size.
        message, "Input array must have 3 dimensions"
 if !VERSION.RELEASE GE '5.6' and (N_elements(mask) EQ 0) then begin
        outarr = median(inarr,dimension=3,/even)
        return
 endif

;                       Create the output array.
 ncol = s[1]
 nrow = s[2]
 narr = s[3]
 type = s[s[0] + 1]
 outarr = make_array( dimen = [ncol,nrow], /NOZERO, TYPE = type )
 if N_params() GT 2 then $
        output_mask = make_array (dimen = [ncol,nrow], VALUE = 1b)
 even = (narr mod 2) EQ 0

;                       Combine the input arrays into the output array.

 mask_given = 0b
 if keyword_set(mask) then begin
    sm = size(mask)
    if N_elements(mask) LT 4 then $ 
           message,'Input mask not valid... must have 3 dimensions'
    w = where(sm[0:3] eq s[0:3], cw)
    if cw eq 4 then begin
       mask_given = 1b 
    endif else begin
       message,'Mask not valid... must be same shape as input cube.'
    endelse
 endif

 if not mask_given then begin
; If the /EVEN keyword is not needed, then it is faster not to use it

     if even then begin     
       for j = 0l, nrow-1 do $
            for i = 0l, ncol-1 do outarr[i,j] = median(inarr[i,j,*],/EVEN)
      endif else begin 
        for j = 0l, nrow-1 do $
            for i = 0l, ncol-1 do outarr[i,j] = median(inarr[i,j,*])
     endelse


 endif else begin

 for j = 0l, (nrow-1) do begin    
        for i = 0l, (ncol-1) do begin
                good_pixels = 1b   
                       wmask = where(mask[i,j,*],cwm)
                       if cwm gt 0 then begin
                          marr = inarr[i,j,wmask] 
                       endif else begin
                          good_pixels = 0b
                          output_mask[i,j] = 0b
                       endelse
  
                if good_pixels then outarr[i,j] = median(marr,/EVEN)
          
        endfor
 endfor
 endelse 

 return
 end
