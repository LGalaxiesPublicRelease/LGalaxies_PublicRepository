;+
; NAME:
;   eta_to_stripe
; PURPOSE:
;   find the stripe which an eta value is in; hardwired to what astrotools 
;   v5_6 does; except it uses lambda to deal with southern stripes
; CALLING SEQUENCE:
;   etalambda_to_stripe, eta, stripe
; INPUTS:
;   eta      value of eta (survey lat) in deg
;   lambda   value of lambda (survey long) in deg
; OUTPUTS:
;   stripe   Survey Stripe #
; BUGS:
;   Location of the survey center is hard-wired, not read from astrotools.
; REVISION HISTORY:
;   2002-Feb-20  written by Blanton (NYU)
;-
pro eta_to_stripe,eta,lambda,stripe

stripe_separation=2.5D
if(n_elements(eta) eq 1) then begin
  if(abs(lambda) lt 90.D) then begin
		; north
    stripe=long((eta+58.75D)/stripe_separation)
  endif else begin
		; south
    stripe=long((eta+(58.75D)+180.D)/stripe_separation)
  endelse
endif else begin
  nindx=where(abs(lambda) lt 90.D,ncount)
  sindx=where(abs(lambda) ge 90.D,scount)
	stripe=lonarr(n_elements(eta))
	if(ncount gt 0) then $
  	stripe[nindx]=floor((eta[nindx]+58.75D)/stripe_separation)
	if(scount gt 0) then $
    stripe[sindx]=floor((eta[sindx]+(58.75D)+180.)/stripe_separation)
endelse

end

