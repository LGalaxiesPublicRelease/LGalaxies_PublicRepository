;-----------------------------------------------------------------------
;  Solve a spherical triangle given two sides and the included angle.
;  This routine is called by the coordinate rotation routines.
pro sptr1,a,ab,ac,b,c,bc

   DRADEG = 180.d0/!dpi
   z = 0.d0

   sa = sin(a/DRADEG)
   ca = cos(a/DRADEG)
   sab = sin(ab/DRADEG)
   cab = cos(ab/DRADEG)
   sac = sin(ac/DRADEG)
   cac = cos(ac/DRADEG)
   cbc = cab*cac + sab*sac*ca

   cs = cos(0.5d0*a/DRADEG) * sin((ab-ac)/DRADEG)
   ss = sin(0.5d0*a/DRADEG) * sin((ab+ac)/DRADEG)
   sbcsq = cs*cs + ss*ss + sab*sab*sac*sac*sa*sa
   sbc = sqrt(sbcsq)

   sb = sac*sa
   cb = cac*sab - cab*sac*ca
   sc = sab*sa
   cc = cab*sac - cac*sab*ca

   ; Special cases
   indx1 = where ( (sb EQ z AND cb EQ z) OR (sc EQ z AND cc EQ z) )
   indx2 = where ( (sb EQ z AND cb EQ z) OR (sc EQ z AND cc EQ z) $
    OR (abs(sab)+abs(sac) LE abs(sa)) )
   if (indx1[0] NE -1) then begin
      sb[indx1] = sa[indx1]
      cb[indx1] = ca[indx1]
      sc[indx1] = 0*sc[indx1]
      cc[indx1] = 0*cc[indx1] - 1.d0
   endif
; DJS 07-OCT-1998, the below appears to be a mistake when used with RECENTER
; and was removed
;   if (indx2[0] NE -1) then begin
;      cb[indx2] = -cab[indx2]*ca[indx2]
;      cc[indx2] = cab[indx2]
;   endif

   b = atan(sb,cb) * DRADEG
   c = atan(sc,cc) * DRADEG
   bc = atan(sbc,cbc) * DRADEG
   return
end
;-----------------------------------------------------------------------
