FUNCTION uv2ll, vector
;+                                                                  
;  NAME:
;    uv2ll
;
;  PURPOSE: To convert from unit vectors to longitude/latitude
;
;  CALLING SEQUENCE:
;    lon_lat = uv2ll(uv)
;
;  INPUT:
;    (n,3) unit vector array
;
;  OUTPUT:
;    (n,2) longitude/latitude array
;
;  SUBROUTINES CALLED:
;    None
;
;  REVISION HISTORY
;
;  SPR 10476  Add Documentation  J.M. Gales  01/21/93
;-

	n_vec = n_elements(vector) / 3
	lon_lat = dblarr(n_vec,2)
	; get number of input vectors
	; generate lon/lat output array

	r2d = 45 / atan(1.0d0)
	; calculate radian to degree conversion

	lon_lat(*,0) = atan(double(vector(*,1)),double(vector(*,0)))
	lon_lat(*,1) = asin(double(vector(*,2)))
	; longitude is arc tangent of (y/x)
	; longitude is arc sine of z

	lon_lat = r2d * lon_lat
	; convert to degrees

	lon_lat(*,0) =  lon_lat(*,0) * (lon_lat(*,0) GE 0) + $
			(360 + lon_lat(*,0)) * (lon_lat(*,0) LT 0)
	; if longitude is less than 0 add 360 to it


RETURN, lon_lat
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


