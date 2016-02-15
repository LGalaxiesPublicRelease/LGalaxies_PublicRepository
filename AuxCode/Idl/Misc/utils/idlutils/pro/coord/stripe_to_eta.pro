;+
; NAME:
;   stripe_to_eta
;
; PURPOSE:
;   Convert from SDSS great circle coordinates to equatorial coordinates.
;
; CALLING SEQUENCE:
;   eta = stripe_to_eta(stripe)
;
; INPUTS:
;   stripe     - Stripe number for SDSS coordinate system.  If specified,
;                the NODE,INCL are ignored; scalar or array.
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   eta        - Eta in SDSS (lambda,eta) coordinate system (degrees);
;                scalar or array with same dimensions as STRIPE.
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
;
; REVISION HISTORY:
;   20-Feb-2002  Written by M. Blanton, NYU
;   03-Oct-2002  Modified by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function stripe_to_eta, stripe

   if (n_params() NE 1) then $
    message, 'Wrong number of arguments'

   stripe_sep = 2.5d
   eta = stripe * stripe_sep - 57.5d - 180.d * (stripe GT 46)

   return, eta
end
;------------------------------------------------------------------------------
