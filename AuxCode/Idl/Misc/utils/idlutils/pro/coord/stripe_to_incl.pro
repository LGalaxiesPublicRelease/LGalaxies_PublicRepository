;+
; NAME:
;   stripe_to_incl
;
; PURPOSE:
;   Convert from SDSS stripe number to an inclination relative to the equator.
;
; CALLING SEQUENCE:
;   inc = stripe_to_incl(stripe)
;
; INPUTS:
;   stripe     - Stripe number for SDSS coordinate system.  If specified,
;                the NODE,INCL are ignored; scalar or array.
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   incl       - Inclination of great circle relative to the J2000
;                celestial equator (degrees); scalar or array with same
;                dimensions as STRIPE.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   stripe_to_eta()
;
; REVISION HISTORY:
;   20-Feb-2002  Written by M. Blanton, NYU
;   03-Oct-2002  Modified by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function stripe_to_incl, stripe

   if (n_params() NE 1) then $
    message, 'Wrong number of arguments'

   dec_center = 32.5d
   etacenter = stripe_to_eta(stripe)
   incl = etacenter + dec_center

   return, incl
end
;------------------------------------------------------------------------------
