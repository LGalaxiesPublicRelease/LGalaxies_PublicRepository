; D. Finkbeiner
; 28 Feb 1999
; Returns Galactic coordinates (l,b) for healpix projection

; always get HEALPIX (l,b) via this routine for uniformity
; of sign conventions. 
; 17 Sep 2002 added nest keyword - DPF
PRO healgen_lb, nside, l, b, nest=nest, double=double

  healgen, nside, theta, phi, nest=nest, double=double
  l = phi*!radeg
  b = 90.-theta*!radeg
  
  return
END
