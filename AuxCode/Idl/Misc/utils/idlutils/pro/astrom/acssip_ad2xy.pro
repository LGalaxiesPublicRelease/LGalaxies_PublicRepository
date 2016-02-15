;+                                                                  
;  NAME:
;    acssip_ad2xy
;
;  PURPOSE:
;    Read new ACS _flt.fits header to invert distortion matrix and
;    find distorted x & y positions given tangent RA/Dec in arcseconds
;    relative to the center pixel (CRPIX1,2). Will only work with headers 
;    that have A_p_q, B_p_q keywords (after Aug 2004), 
;    and are only found in extensions 1 (CCD #2) and 4 (CCD #1).
;
;  CALLING SEQUENCE:
;    acssip_xy2ad,  ra, dec, hdr, xd, yd
;
;  INPUT:
;    ra          : right ascension in arcseconds in the tangent plane from 
;                  the pixel desginated by (CRPIX1-1, CRPIX2-1)
;    dec         : declination in arcseconds
;    hdr         : ACS header, designed for WFC, from HDU #1 or #4 in _flt files
;    
;  OUTPUT:
;    xd          : distorted pixel column in native CCD frame
;    yd          : distorted pixel row in native CCD frame
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

pro acssip_ad2xy, ra, dec, hdr, x, y, niter=niter

      if NOT keyword_set(niter) then niter=5L

      cd = 3600.0* [[sxpar(hdr, 'CD1_*')],[sxpar(hdr, 'CD2_*')]]
      icd = invert(cd)

      xdu = icd[0,0] * ra + icd[1,0] * dec
      ydu = icd[0,1] * ra + icd[1,1] * dec

      xd = xdu
      yd = ydu

      for ii = 1,niter do begin

        f = xdu*0.0d
        g = ydu*0.0d

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

        diff_x = (xdu - xd)
        diff_y = (ydu - yd)

        xd = xd - (f - diff_x) 
        yd = yd - (g - diff_y)

      endfor

      x = xd + (sxpar(hdr, 'CRPIX1') - 1.)
      y = yd + (sxpar(hdr, 'CRPIX2') - 1.)

      return
    end

