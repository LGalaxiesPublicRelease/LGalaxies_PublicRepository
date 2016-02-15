;+
; NAME:
;   radec_to_munu
;
; PURPOSE:
;   Convert from equatorial coordinates to SDSS great circle coordinates.
;
; CALLING SEQUENCE:
;   radec_to_munu, ra, dec, [ mu, nu, stripe=, node=, incl=, phi= ]
;
; INPUTS:
;   ra         - Right ascension (J2000 degrees)
;   dec        - Declination (J2000 degrees)
;
; OPTIONAL INPUTS:
;   stripe     - Stripe number for SDSS coordinate system.  If specified,
;                the NODE,INCL are ignored; scalar or array with same
;                dimensions as MU.
;   node       - Node of great circle on the J2000 celestial equator (degrees),
;                scalar or array with same dimensions as MU.
;   incl       - Inclination of great circle relative to the J2000
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   mu         - Mu coordinate, scalar or array (degrees)
;   nu         - Nu coordinate, scalar or array (degrees)
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
pro radec_to_munu, ra, dec, mu, nu, stripe=stripe, node=node, incl=incl, $
 phi=phi

   if (n_params() LT 2) then begin
       print, 'Syntax - radec_to_munu, ra, dec, [ mu, nu, stripe=, node=, incl=]'
       return
   endif

   if (keyword_set(stripe)) then begin
      node = 95.d
      incl = stripe_to_incl(stripe)
   endif else begin
      if (n_elements(node) EQ 0 OR n_elements(incl) EQ 0) then $
       message, 'Must specify either STRIPE or NODE,INCL'
   endelse
   if (n_elements(ra) NE n_elements(dec)) then $
    message, 'Number of elements in RA and DEC must agree'

   r2d = 180.d / !dpi

   if (n_elements(node) LT 1 OR n_elements(incl) LT 1) then begin
      node = 95.d
      incl = stripe_to_incl(stripe)
   endif

   sinra = sin((ra-node)/r2d)
   cosra = cos((ra-node)/r2d)
   sindec = sin(dec/r2d)
   cosdec = cos(dec/r2d)
   sini = sin(incl/r2d)
   cosi = cos(incl/r2d)

   x1 = cosdec * cosra
   y1 = cosdec * sinra
   z1 = sindec
   x2 = x1
   y2 = y1 * cosi + z1 * sini
   z2 = -y1 * sini + z1 * cosi
   mu = r2d * atan(y2,x2) + node
   nu = r2d * asin(z2)
   cirrange, mu

   if (arg_present(phi)) then begin
      ; See my notebook on 14 Oct 2002 for this derivation.
      sinmu = sin((mu-node)/r2d)
      cosmu = cos((mu-node)/r2d)
      sinnu = sin(nu/r2d)
      cosnu = cos(nu/r2d)
      phi = r2d $
       * atan(cosmu * sini, (-sinmu * sinnu * sini + cosnu * cosi)*cosnu)

   endif

   return
end
;------------------------------------------------------------------------------
