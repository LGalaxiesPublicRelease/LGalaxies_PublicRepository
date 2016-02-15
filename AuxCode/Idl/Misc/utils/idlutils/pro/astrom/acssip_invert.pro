;+                                                                  
;  NAME:
;    acssip_invert
;
;  COMMENTS: 
;    Based on the reference of "The SIP Convention for Representing 
;      Distortion in FITS Image Headers", by Shupe et al. 2005.
;    Iterates the distortion correction to invert the distortion pixel map
;
;  EXAMPLE:
;
;    adxy, hdr, ra, dec, xl, yl
;    acssip_invert, xl, yl, hdr, x, y
;
;  REVISION HISTORY:
;    Implemented by S. Burles during the lost year of 05. 
;-

pro acssip_invert, xl, yl, hdr, x, y, niter=niter

      if NOT keyword_set(niter) then niter=5L

      xdu = xl - (sxpar(hdr, 'CRPIX1') - 1.)
      ydu = yl - (sxpar(hdr, 'CRPIX2') - 1.)

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

