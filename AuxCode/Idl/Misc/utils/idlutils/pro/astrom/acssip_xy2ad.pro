;+                                                                  
;  NAME:
;    acssip_xy2ad
;
;  PURPOSE:
;    Read new ACS _flt.fits header to convert pixel locations into
;    RA/Dec pairs relative to CRPIX1,2.  Will only work with headers
;    that have A_p_q, B_p_q keywords (after Aug 2004), and are only found
;    in extensions 1 (CCD #2) and 4 (CCD #1)
;
;  CALLING SEQUENCE:
;    acssip_xy2ad,  pix_col, pix_row, hdr, ra, dec
;
;  INPUT:
;    pix_col     : any dimension array containing pixel column positions
;    pix_row     : any dimension array containing pixel row positions
;    hdr         : ACS header, designed for WFC, from HDU #1 or #4 in _flt files
;    
;  OUTPUT:
;    ra          : right ascension in arcseconds in the tangent plane from 
;                  the pixel desginated by (CRPIX1-1, CRPIX2-1)
;    dec         : declination in arcseconds
;
;  SUBROUTINES CALLED:
;    sxpar()
;
;  COMMENTS: 
;    Based on the reference of "The SIP Convention for Representing 
;      Distortion in FITS Image Headers"
;    by Shupe et al. 2005.
;
;  BUGS:
;    Doesn't do the tangent-plane projection correctly; i.e., this
;      assumes that tan(theta)=theta.
;
;  REVISION HISTORY
;    Implemented by S. Burles during the lost year of 05. 
;
;-

pro acssip_xy2ad, x, y, hdr, rab, decb

      xd = x - (sxpar(hdr, 'CRPIX1') - 1.0d)
      yd = y - (sxpar(hdr, 'CRPIX2') - 1.0d)
     
      f = xd*0.0d
      g = yd*0.0d

      a_order = sxpar(hdr, 'A_ORDER')
      b_order = sxpar(hdr, 'B_ORDER')

      norder = a_order > b_order

      for p=0,norder do begin
        for q = 0,norder - p do begin  
           a = sxpar(hdr, string('A_',p,'_',q,format='(a,i1,a,i1)')) 
           f = f + a*(xd^p) * (yd^q) 
           b = sxpar(hdr, string('B_',p,'_',q,format='(a,i1,a,i1)')) 
           g = g + b*(xd^p) * (yd^q) 
        endfor
      endfor

      xd = xd + f
      yd = yd + g

      rab  = 3600.0*(xd*sxpar(hdr, 'CD1_1') + yd*sxpar(hdr,'CD1_2'))  
      decb = 3600.0*(xd*sxpar(hdr, 'CD2_1') + yd*sxpar(hdr,'CD2_2'))  

    return
    end

