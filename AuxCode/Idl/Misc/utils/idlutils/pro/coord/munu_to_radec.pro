;+
; NAME:
;   munu_to_radec
;
; PURPOSE:
;   Convert from SDSS great circle coordinates to equatorial coordinates.
;
; CALLING SEQUENCE:
;   munu_to_radec, mu, nu, [ ra, dec, stripe=, node=, incl=, phi= ]
;
; INPUTS:
;   mu         - Mu coordinate, scalar or array (degrees)
;   nu         - Nu coordinate, scalar or array (degrees)
;
; OPTIONAL INPUTS:
;   stripe     - Stripe number for SDSS coordinate system.  If specified,
;                the NODE,INCL are ignored; scalar or array with same
;                dimensions as MU.
;   node       - Node of great circle on the J2000 celestial equator (degrees),
;                scalar or array with same dimensions as MU.
;   incl       - Inclination of great circle relative to the J2000
;                celestial equator (degrees); scalar or array with same
;                dimensions as MU.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   ra         - Right ascension (J2000 degrees)
;   dec        - Declination (J2000 degrees)
;   phi        - Counter-clockwise position angle w.r.t. north for an arc
;                in the +nu direction.
;
; COMMENTS:
;   Either STRIPE or NODE,INCL must be specified.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   cirrange
;   stripe_to_incl()
;
; REVISION HISTORY:
;   20-Feb-2002  Written by M. Blanton, NYU
;   03-Oct-2002  Modified by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro munu_to_radec, mu, nu, ra, dec, stripe=stripe, node=node, incl=incl, $
 phi=phi

   if (n_params() LT 2) then begin
       print, 'Syntax - munu_to_radec, mu, nu, ra, dec [, stripe=, node=, incl=]'
       return
   endif

   if (keyword_set(stripe)) then begin
      node = 95.d
      incl = stripe_to_incl(stripe)
   endif else begin
      if (n_elements(node) EQ 0 AND n_elements(incl) EQ 0) then $
       message, 'Must specify either STRIPE or NODE,INCL'
   endelse
   if (n_elements(mu) NE n_elements(nu)) then $
    message, 'Number of elements in MU and NU must agree'

   r2d = 180.d / !dpi
   d2r = !dpi / 180.d

   sinnu = sin(nu*d2r)
   cosnu = cos(nu*d2r)
   sini = sin(incl*d2r)
   cosi = cos(incl*d2r)
   sinmu = sin((mu-node)*d2r)
   cosmu = cos((mu-node)*d2r)

   xx = cosmu * cosnu
   yy = sinmu * cosnu * cosi - sinnu * sini
   zz = sinmu * cosnu * sini + sinnu * cosi

   ra = r2d * atan(yy,xx) + node
   dec = r2d * asin(zz)

   cirrange, ra

   if (arg_present(phi)) then begin
      ; See my notebook on 14 Oct 2002 for this derivation.
      phi = r2d $
       * atan(cosmu * sini, (-sinmu * sinnu * sini + cosnu * cosi)*cosnu)

      ; Compute the rotation angle numerically, since the returned value
      ; from the above seems to be incorrect???
;      dradeg = 180.d0 / !dpi
;      radec_to_munu, ra, dec, mu1, nu1, node=node, incl=incl
;      radec_to_munu, ra, dec+0.001d0, mu2, nu2, node=node, incl=incl
;      phi = atan(mu2-mu1,nu2-nu1) * dradeg
   endif

   return
end
;------------------------------------------------------------------------------
