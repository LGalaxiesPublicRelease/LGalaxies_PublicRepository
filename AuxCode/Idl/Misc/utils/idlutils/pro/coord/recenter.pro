;-----------------------------------------------------------------------
;  Rotate to new coordinates centered on (clng,clat) such that locally
;  the new lines of latitude are parallel to the old lines of latitude.
;  This is done by first rotating the coordinate sphere about the
;  z-axis by -clng and then rotating the great circle nlng=0 by -nlat.
;  Input:  lng=  longitude
;          lat=  latitude
;          clng= central longitude
;          clat= central latitude
;  Return: nlng= new longitude
;          nlat= new latitude
;  All coordinates are in degrees.
pro recenter, lng, lat, clng, clat, nlng, nlat

;  Set the following variables:
;  rag  = lng (old RA) for the new pole
;  decg = lat (old DEC) for the new pole
;  pgl  = 90 degrees + nlng (new RA) where lat=0 crosses nlat=0
;                                    (where old DEC=0 crosses new DEC=0)
;  Reference Lang's Astrophysical Formulae pg 505 to visualize.
   rag = clng
   decg = 90.d0 + clat
   pgl = 180.d0
   ; The "0*lng" is to convert the expression from a scalar to an array
   sptr1, lng-rag, 0*lng+90.d0-decg, 90.d0-lat, pgq, gqp, gq
   nlng = djs_angpos(pgl-pgq)
   nlat = 90.d0 - gq
   return
end
;-----------------------------------------------------------------------
