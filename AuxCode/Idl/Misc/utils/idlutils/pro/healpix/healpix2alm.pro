;+
; NAME:
;   healpix2alm
;
; PURPOSE:
;   Compute spherical harmonic transform (Alm) of a healpix map
;
; CALLING SEQUENCE:
;   Alm = healpix2alm(data, lmax=lmax)
;
; INPUTS:
;   data   - healpix array [npix, nmap] (must be real)
;   
; KEYWORDS:
;   lmax   - maximum l to compute (Default???)
;
; OUTPUTS:
;   Alm    - dcomplex array of Alm[l, m, nmap] (lmax by lmax)
;
; OPTIONAL OUTPUTS:
;   
; RESTRICTIONS:
;   
; EXAMPLES:
;   
; COMMENTS:
;  Note - the Ylm are not even perfectly orthonormal on the healpix sphere
;   print,total(spher_harm(theta,phi,5,3,/doub)*conj(spher_harm(theta,phi,9,3,/doub)))*!dpi*4,format='(F20.10)'
;   
; REVISION HISTORY:
;   2003-Feb-19  Written by Douglas Finkbeiner, Princeton
;   2003-Nov-12  Can call with multiple maps at a time - DPF & NP
;
;----------------------------------------------------------------------
;
; SUBROUTINE: 
;   healpix_fqm
;
; PURPOSE:
;   Compute the product f(phi)*exp(-i*m*phi)
;
; CALLING SEQUENCE:
;   Fqm = healpix_fqm(lmax, nside, qind, phi, data, dumb=dumb)
;
; INPUTS:
;   lmax   - max l to compute
;   nside  - healpix Nside
;   qind   - Index of last pixel in each ring, i.e. uniq(theta)
;   phi    - phi coordinate of healpix pixel centers
;   data   - healpix map
;
; KEYWORDS:
;   dumb   - do it the slow way (for testing)
;
; OUTPUTS:
;   Fqm    - dcomplex array of f(phi)*exp(-i*m*phi) 
;             for each q (ring number) and m. 
;
; COMMENTS:
;   The "dumb" way does it brute force.  Only use this as a check
;   The "smart" way uses an FFT and is 20 times faster for Nside=64
;   
;----------------------------------------------------------------------

function healpix_fqm, lmax, nside, qind, phi, data, dumb=dumb

; -------- Preliminaries
  i   = sqrt(dcomplex(-1, 0))
  nq  = 4*nside-1
  Fqm = dcomplexarr(nq, lmax+1)
  
; -------- Set up phase array for FFT method - healpix ring offsets in phi
  phase = dblarr(nq)
  phase[0:nside-2] = !dpi/(dindgen(nside-1)+1)/4.
  phase[nq-nside+1: nq-1] = reverse(!dpi/(dindgen(nside-1)+1)/4.)
  q = lindgen(nside+1)*2+nside-1
  phase[q] = !dpi/nside/4

; -------- Loop over rings (indexed on q)
  for q=0, nq-1 do begin 
     q0    = (q eq 0) ? 0 : qind[q-1]+1
     q1    = qind[q]
     phiq  = phi[q0:q1]
     dataq = data[q0:q1]
     ll = ((q1-q0) < lmax)
     if keyword_set(dumb) then begin 
        for m=0, ll do begin 
           Fqm[q, m] = total(exp(-i*m*phiq)*dataq)
        endfor 
     endif else begin 
        Fqm[q, 0:ll] = (fft(dataq))[0:ll]*(q1-q0+1)
        if phase[q] ne 0 then for m=0, ll do Fqm[q, m] = Fqm[q, m]*exp(-i*phase[q]*m)
     endelse
  endfor

  return, Fqm
end



function healpix2alm, data, lmax=lmax, mmax=mmax

  t1 = systime(1)
  if NOT keyword_set(data) then begin  ; test setup
     nside = 64L
     npix = 12L*nside*nside
     data = randomn(555., npix, /double)
  endif else begin 
     nmap = size(data, /n_dim) eq 1 ? 1 : (size(data, /dimens))[1]
     nside = round(sqrt(n_elements(data)/nmap/12))
     npix = 12L*nside*nside
  endelse

  ringwt = healpix_ring_weight(nside, iring=iring)
  wt = ringwt[iring]

; -------- default lmax (Nyquist on equatorial ring)
  if NOT keyword_set(lmax) then lmax = nside*2
  if NOT keyword_set(mmax) then mmax = lmax 

  healgen, nside, theta, phi, /double
  qind = uniq(theta)

  Alm = dcomplexarr(lmax+1, lmax+1, nmap)

; -------- Generate exp(i*m*phi) part
  Fqm = dcomplexarr(4*nside-1, lmax+1, nmap)
  for imap=0L, nmap-1 do begin 
     Fqm[*, *, imap] = healpix_fqm(lmax, nside, qind, phi, data[*, imap]*wt)
  endfor

  splog, 'got Fqm'+(nmap gt 1 ? 's' : '')
  t2 = systime(1)
  xq = cos(theta[qind])
  xqhalf = xq[0:nside*2-1]

  lfac = sqrt((2*dindgen(lmax+1)+1)/(4*!dpi))
  revind = reverse(lindgen(nside*2-1))
  for m=0, mmax do begin 
     tt = systime(1)
     mqladvance, xqhalf, m, Mql_1, Mql, lmax=lmax
     
     Fqn = Fqm[*, m, *]
     for l=m, lmax do begin 
        sign = 1-((l+m) mod 2)*2
        Mq_thisl = Mql[*, l]
        Mqlfull = [Mq_thisl, sign*Mq_thisl[revind]]

        for imap=0L, nmap-1 do $
          Alm[l, m, imap] = lfac[l]*total(Mqlfull*Fqn[*, 0, imap])
     endfor 

     Mql_1 = temporary(Mql)
     if (m mod 64) eq 0 then splog, m, systime(1)-tt
  endfor 

  splog, systime(1)-t1, ' seconds Total'
  splog, 'Fqm: ', t2-t1, 'Alm ', systime(1)-t2

  Alm = Alm*(4*!dpi/npix)

  return, Alm
end
