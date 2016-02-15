PRO ang2pix_ring, nside, theta, phi, ipring
;*********************************************************************************
;+
; ANG2PIX_RING, Nside, Theta, Phi, Ipring
;
;        renders the RING scheme pixel number Ipring for a pixel which, given the
;        map resolution parameter Nside, contains the point on the sphere
;        at angular coordinates Theta and Phi
;
; INPUT
;    Nside     : determines the resolution (Npix = 12* Nside^2)
;	SCALAR
;    Theta : angle (along meridian), in [0,Pi], theta=0 : north pole,
;	can be an ARRAY
;    Phi   : angle (along parallel), in [0,2*Pi]
;	can be an ARRAY of same size as theta
;
; OUTPUT
;    Ipring  : pixel number in the RING scheme of HEALPIX pixelisation in [0,Npix-1]
;	can be an ARRAY of same size as Theta and Phi
;    pixels are numbered along parallels (ascending phi), 
;    and parallels are numbered from north pole to south pole (ascending theta)
;
;
; HISTORY
;    June-October 1997,  Eric Hivon & Kris Gorski, TAC, 
;            original ang_pix
;    Feb 1999,           Eric Hivon,               Caltech
;            name changed to ang2pix_ring
;
;-
;*********************************************************************************

ns_max = 8192L

if N_params() ne 4 then begin
    print,' syntax = ang2pix_ring, nside, theta, phi, ipix'
    stop
endif
if (N_ELEMENTS(nside) GT 1) then begin
	print,'nside should be a scalar in ang2pix_ring'
	stop
endif
if (nside lt 1) or (nside gt ns_max) then stop, 'nside out of range'

np = N_ELEMENTS(theta)
np1 = N_ELEMENTS(phi) 
if (np NE np1) then begin
	print,'inconsistent theta and phi in ang2pix_ring'
	stop
endif
ipring = LONARR(np)
;------------------------------------------------------------
nside  = LONG(nside)
pion2 = !DPI * 0.5d0
twopi = !DPI * 2.d0
nl2   = 2*nside
nl4   = 4*nside
npix  = (3L*nside)*(4L*nside)
ncap  = nl2*(nside-1L)

cth0 = 2.d0/3.d0

cth_in = COS(DOUBLE(theta))
phi_in = phi MOD twopi
phi_in = phi + (phi LE 0.d0)*twopi

pix_eqt = WHERE(cth_in LE  cth0 AND cth_in GT -cth0, n_eqt) ; equatorial strip
pix_np  = WHERE(cth_in GT  cth0, n_np)                      ; north caps
pix_sp  = WHERE(cth_in LE -cth0, n_sp)                      ; south pole

IF (n_eqt GT 0) THEN BEGIN ; equatorial strip ----------------
      tt = phi_in(pix_eqt) / pion2

      jp = LONG(nside*(0.5d0 + tt - cth_in(pix_eqt)*0.75d0)) ; increasing edge line index
      jm = LONG(nside*(0.5d0 + tt + cth_in(pix_eqt)*0.75d0)) ; decreasing edge line index

      ir = (nside + 1) + jp - jm ; in {1,2n+1} (ring number counted from z=2/3)
      k =  ( (ir MOD 2) EQ 0)   ; k=1 if ir even, and 0 otherwise

      ip = LONG( ( jp+jm+k + (1-nside) ) / 2 ) + 1 ; in {1,4n}
      ip = ip - nl4*(ip GT nl4)

      ipring(pix_eqt) = ncap + nl4*(ir-1) + ip - 1
      tt = 0 & jp = 0 & jm = 0 & ir = 0 & k = 0 & ip = 0
ENDIF

IF (n_np GT 0) THEN BEGIN ; north polar caps ------------------------

      tt = phi_in(pix_np) / pion2
      tp = tt MOD 1.d0
      tmp = SQRT( 3.d0*(1.d0 - ABS(cth_in(pix_np))) )

      jp = LONG( nside * tp          * tmp ) ; increasing edge line index
      jm = LONG( nside * (1.d0 - tp) * tmp ) ; decreasing edge line index

      ir = jp + jm + 1         ; ring number counted from the closest pole
      ip = LONG( tt * ir ) + 1 ; in {1,4*ir}
      ir4 = 4*ir
      ip = ip - ir4*(ip GT ir4)

      ipring(pix_np) =        2*ir*(ir-1) + ip - 1
ENDIF ; -------------------------------------------------------

IF (n_sp GT 0) THEN BEGIN ; south polar caps ------------------------

      tt = phi_in(pix_sp) / pion2
      tp = tt MOD 1.d0
      tmp = SQRT( 3.d0*(1.d0 - ABS(cth_in(pix_sp))) )

      jp = LONG( nside * tp          * tmp ) ; increasing edge line index
      jm = LONG( nside * (1.d0 - tp) * tmp ) ; decreasing edge line index

      ir = jp + jm + 1         ; ring number counted from the closest pole
      ip = LONG( tt * ir ) + 1 ; in {1,4*ir}
      ir4 = 4*ir
      ip = ip - ir4*(ip GT ir4)

      ipring(pix_sp) = npix - 2*ir*(ir+1) + ip - 1
ENDIF ; -------------------------------------------------------

return
end ; ang2pix_ring

;=======================================================================
; The permission to use and copy this software and its documentation, 
; without fee or royalty is limited to non-commercial purposes related to 
; Boomerang, Microwave Anisotropy Probe (MAP) and
; PLANCK Surveyor projects and provided that you agree to comply with
; the following copyright notice and statements,
; and that the same appear on ALL copies of the software and documentation.
;
; An appropriate acknowledgement has to be included in any
; publications based on work where the package has been used
; and a reference to the homepage http://www.tac.dk/~healpix
; should be included
;
; Copyright 1997 by Eric Hivon and Kris Gorski.
;  All rights reserved.
;=======================================================================
