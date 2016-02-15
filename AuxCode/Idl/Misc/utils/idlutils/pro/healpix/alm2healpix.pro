;+
; NAME:
;   alm2healpix
;
; PURPOSE:
;   Compute healpix map from spherical harmonic transform (Alm)
;
; CALLING SEQUENCE:
;   map = alm2healpix(nside, alm, lmax=lmax)
;
; INPUTS:
;   nside   - healpix nside (number of pixels is 12*nside^2)
;   alm     - dcomplex array of Alm coefficients (Alm[l, m])
;                Of course this is zero for l<m
;   
; KEYWORDS:
;   lmax   - maximum l to compute (Default determined by alm array size)
;
; OUTPUTS:
;   map    - healpix map
;
; OPTIONAL OUTPUTS:
;   
; RESTRICTIONS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   Something is a little funny (part in 1 million) at the poles
;
; REVISION HISTORY:
;   2003-Mar-14  Written by Douglas Finkbeiner, Princeton
;   2003-Nov-12  Can do many maps at once, sped up - DPF & NP
;
;----------------------------------------------------------------------

function fqm2healpix, nside, qind, phi, Fqm

  i = sqrt(dcomplex(-1, 0))

  nq  = 4*nside-1
  nm  = (size(Fqm, /dimens))[1]
; maybe make this double??
  map = dcomplexarr(n_elements(phi))

; -------- Set up phase array for FFT method
  phase = dblarr(nq)
  phase[0:nside-2] = !dpi/(dindgen(nside-1)+1)/4.
  phase[nq-nside+1: nq-1] = reverse(!dpi/(dindgen(nside-1)+1)/4.)
  q = lindgen(nside+1)*2+nside-1
  phase[q] = !dpi/nside/4

; -------- Loop over rings
  for q=0, nq-1 do begin 
     q0    = (q eq 0) ? 0 : qind[q-1]+1
     q1    = qind[q]
     phiq  = phi[q0:q1]
     lq    = (q1-q0)+1
     Fm = dcomplexarr(lq) ; correct length
     lq2 = (lq/2) < nm
     m  = dindgen(lq2)
     Fm[0:lq2-1] = reform(Fqm[q, 0:lq2-1], lq2)*exp(i*phase[q]*m)
     Fm[0] = Fm[0]/2   ; hack
     map[q0:q1] = fft(Fm, /inv)
     
  endfor

  return, map
end



function alm2healpix, nside, alm, lmax=lmax

; -------- find number of maps to generate
  nmap = size(alm, /n_dim) eq 2 ? 1 : (size(alm, /dimens))[2]

; -------- check keywords
  if NOT keyword_set(lmax) then lmax = (size(alm, /dimens))[0]-1
  if NOT keyword_set(mmax) then mmax = lmax 

  t1 = systime(1)
  npix = 12L*nside*nside
  healgen, nside, theta, phi, /double
  qind = uniq(theta)
  xq   = cos(theta[qind])

; -------- Initialize Fqm array (q is ring index; Fqm is FFT of each ring)
  nq  = 4*nside-1
  Fqm = dcomplexarr(nq, lmax+1, nmap)

  xqhalf = xq[0:nside*2-1] ; only compute half the indices
  lfac = sqrt((2*dindgen(lmax+1)+1)/(4*!dpi))
  revind = reverse(lindgen(nside*2-1))
  mlqind = [lindgen(nside*2), revind]
  for m=0, mmax do begin 
     tt = systime(1)
     mqladvance, xqhalf, m, Mql_1, Mql, lmax=lmax

     Mlq = transpose(Mql) ; 5 sec just for transpose!! Can we do this better?

     sign = 1-((indgen(lmax+1)+m) mod 2)*2

     ; suppress m index
     Fqn = dcomplexarr(nq, 1, nmap, /nozero)
     aln = alm[m:lmax, m, *]

     ; loop over rings and sum on l index
     for iq=0L, nq-1 do begin 

        if check_math() NE 0 then stop
        Ml = Mlq[m:lmax, mlqind[iq]]*lfac[m:lmax]
        if iq GE 2*nside then Ml = Ml*sign[m:lmax]

        for imap=0L, nmap-1 do begin 
           Fqn[iq, 0, imap] = total(Ml*aln[*,0,imap])
           cm = check_math()
           if (cm NE 32) and (cm NE 0) then stop
        endfor 
     endfor
     ; fill Fqm array
     Fqm[*, m, *] = Fqn
     Mql_1 = temporary(Mql)
     if (m mod 64) eq 0 then splog, m, systime(1)-tt
  endfor 

  map = dblarr(n_elements(phi), nmap)
  for imap=0L, nmap-1 do $
    map[*, imap] = double(fqm2healpix(nside, qind, phi, Fqm[*, *, imap]))

  print, 'total time', systime(1)-t1

; we only sum over half the alms, so we need a factor of 2
  map = 2*map

  return, map
end
