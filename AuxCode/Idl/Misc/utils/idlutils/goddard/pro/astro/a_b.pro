function A_b,l,b
;+
; NAME:
;     A_b
; PURPOSE:
;     Compute B band interstellar extinction according to the RC2.
; EXPLANATION:
;     The predicted B band extinction is computed as a function of  
;     Galactic position  using the 21 parameter function given by
;     deVaucouleurs in the 2nd Reference Catalog of Galaxies (RC2).   Note 
;     that this formula was not used for the RC3 and that reddenings
;     were instead obtained from the Burstein-Heiles 21 cm maps.
;
; CALLING SEQUENCE:
;     result = A_b( l2, b2)
;
; INPUT PARAMETERS
;     l2 = Galactic longitude (degrees), scalar or vector
;     b2 = Galactic latitude  (degrees), scalar or vector
;
; OUTPUT PARAMETERS
;     RESULT - Interstellar extinction Ab in magnitudes, same number of 
;             elements as input l2 and b2 parameters
;
; NOTES:
;     The controversial aspect of the deVaucouleurs reddening curve
;     is that it predicts an extinction of about 0.2 at the poles 
;
;     The parameters used here differ from the ones printed in the RC2
;     but are the ones actually used for entries in the catalog
;     (see Rowan-Robinson 1985) 
;
;     This procedure is mainly of historical interest only, and reddening
;     is now better determined using dust maps, such as those available at
;     http://astro.berkeley.edu/davis/dust/index.html
; REVISION HISTORY
;     Written by R. Cornett and W. Landsman, STX October 1987
;     Vectorized code      W. Landsman   STX    December 1992
;     Converted to IDL V5.0   W. Landsman   September 1997
;-
  On_error,2
  if N_params() LT 2 then begin
        print,'Syntax -- result = A_b(gal_long, gal_lat) '
        print,'Galactic longitude and latitude in degrees (scalar or vector)'
        return, -1
  endif

  lr = l/!radeg  & br = b/!radeg
; compute the RC2 'C'
   c = 1./sin(br+(0.25/!radeg)-(1.7/!radeg)*sin(lr)-(1./!radeg)*cos(3.*lr))
                                    
   npts = N_elements( lr)
   ab = fltarr( npts)

   pos = where( b GE 0, Npos )
   if Npos GT 0 then begin

      lrpos = lr[pos] 

       sn =      0.1948*cos(lrpos)    + 0.0725*sin(lrpos)     $
              + 0.0953*cos(2.*lrpos) - 0.0751*sin(2.*lrpos)   $
              + 0.0936*cos(3.*lrpos) + 0.0639*sin(3.*lrpos)   $
              + 0.0391*cos(4.*lrpos) + 0.0691*sin(4.*lrpos)
   
      ab[pos] = 0.19*( 1.+ sn * cos(br[pos] )) * abs( c[pos] )
    
  endif 

  neg = where( b LT 0, NNeg)
  if Nneg GT 0 then begin

     lrneg = lr[neg]

     ss =     0.1749*cos(lrneg) - 0.01112*sin(lrneg)         $
           + 0.1438*cos(2.*lrneg) - 0.0180*sin(2.*lrneg)     $
           - 0.0897*cos(3.*lrneg) - 0.0013*sin(3.*lrneg)     $
           + 0.0568*cos(4.*lrneg) + 0.0433*sin(4.*lrneg)
  
     ab[neg] = 0.21 * ( 1.+ ss * cos(br[neg]) ) * abs(c[neg]) 

  endif

  return,ab
  end
