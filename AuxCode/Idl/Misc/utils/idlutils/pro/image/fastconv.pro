;+
; NAME:
;    FASTCONV
; PURPOSE:
;    Perform a convolution faster by binning both the image and
;    kernel down by a factor of BINFACTOR. 
;
; CALLING SEQUENCE:
;     Fastcon, image, kernel, binfactor
; INPUTS:
;     image - input array to be convolved
;     kernel - array to convolve image with - e.g. a Gaussian
;     binfactor - factor to bin down by (must divide both image and 
;                 kernel dimensions)
; OUTPUTS:
;
;
; REVISION HISTORY:
;   Written D. Finkbeiner, 3 Sept 96
; 2  May 1997     Add nodisplay keyword
; 3  May 1997     Add edge_wrap keyword
; 30 March 1998   Allow non-square arrays (introduced bug)
; 24 April 1998   Bug found - failed to divide by binfactor before
;                  rebinning.  Bug fixed. 
; 29 June 1998    Add disc keyword to allow disc smoothing (DPF)
;-

FUNCTION fastconv,big_image,fwhm,binfactor, nodisplay=nodisplay, $
                  edge_wrap=edge_wrap, silent=silent, disc=disc

  verbose = NOT keyword_set(silent)
  
;image_size = round(sqrt(n_elements(big_image))/binfactor)
  isizex = (size(big_image))[1]/binfactor
  isizey = (size(big_image))[2]/binfactor
  
  IF verbose THEN print,'Binning image to ',isizex,' by ',isizey
  image=rebin(big_image,isizex,isizey)
  
  IF keyword_set(disc) THEN BEGIN 
      rad = (float(fwhm)/2.)/binfactor
      ksize = round(2*rad)+3
      xx=shift(dist(ksize),ksize/2,ksize/2)
      kernel = float(xx LT (rad+1))
      kernel=kernel/total(kernel)
  ENDIF ELSE BEGIN 
      sigma=fwhm/(2.355*binfactor)
      ksize=round(4*sigma+1)*2
      xx=shift(dist(ksize),ksize/2,ksize/2)
      kernel=exp(-xx^2/(2*sigma^2))
      kernel=kernel/total(kernel)
  ENDELSE 

  IF verbose THEN print,'kernel is ',ksize,' by ',ksize
  IF keyword_set(nodisplay) EQ 0 THEN $
    surface,kernel
  IF verbose THEN BEGIN 
      print,'Time estimate: ',2.22e-09*isizex*isizey*ksize^2,' minutes.'
      print,'Performing convolution... begun at ',systime()
  ENDIF 
  IF keyword_set(edge_wrap) THEN $
    smooth_image=convol(image,kernel, /edge_wrap) $
  ELSE $ 
    smooth_image=convol(image,kernel) 
  
  IF verbose THEN print,'Done at ',systime()
  
  return,smooth_image
END


