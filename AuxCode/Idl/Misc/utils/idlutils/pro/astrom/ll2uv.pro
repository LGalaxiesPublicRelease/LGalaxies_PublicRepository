FUNCTION ll2uv, lon_lat, double=double
;+                                                                  
;  NAME:
;    ll2uv
;
;  PURPOSE:
;    To convert from longitude/latitude to unit vectors
;
;  CALLING SEQUENCE:
;    uv = ll2uv(lon_lat)
;
;  INPUT:
;    (n,2) longitude/latitude array
;
;  OUTPUT:
;    (n,3) unit vector array
;
;  SUBROUTINES CALLED:
;    None
;
;  REVISION HISTORY
;
; SPR 9923 26-AUG-1992  Change unit vector output array to float type
; J.M. Gales
; 29-March-2001  Added double keyword - Doug Finkbeiner, Princeton
;
;-
;
	n_ll = n_elements(lon_lat) / 2
        if keyword_set(double) then $
          vector = dblarr(n_ll,3) $
        else $
          vector = fltarr(n_ll,3)
        
	d2r = atan(1.0d0) / 45

	lon = lon_lat(*,0) * d2r
	lat = lon_lat(*,1) * d2r

	cos_lat = cos(lat)
	sin_lat = sin(lat)
	cos_lon = cos(lon)
	sin_lon = sin(lon)

	vector(*,0) = cos_lat * cos_lon
	vector(*,1) = cos_lat * sin_lon
	vector(*,2) = sin_lat

RETURN, vector
END
;DISCLAIMER:
;
;This software was written at the Cosmology Data Analysis Center in
;support of the Cosmic Background Explorer (COBE) Project under NASA
;contract number NAS5-30750.
;
;This software may be used, copied, modified or redistributed so long
;as it is not sold and this disclaimer is distributed along with the
;software.  If you modify the software please indicate your
;modifications in a prominent place in the source code.  
;
;All routines are provided "as is" without any express or implied
;warranties whatsoever.  All routines are distributed without guarantee
;of support.  If errors are found in this code it is requested that you
;contact us by sending email to the address below to report the errors
;but we make no claims regarding timely fixes.  This software has been 
;used for analysis of COBE data but has not been validated and has not 
;been used to create validated data sets of any type.
;
;Please send bug reports to CGIS@COBECL.DNET.NASA.GOV.


