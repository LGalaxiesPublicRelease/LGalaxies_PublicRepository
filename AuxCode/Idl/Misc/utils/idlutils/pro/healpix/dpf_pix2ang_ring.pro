;+
; NAME:
;   dpf_pix2ang_ring
;
; PURPOSE:
;   Compute coordinates for HEALPix pixel numbers
;
; CALLING SEQUENCE:
;   dpf_pix2ang_ring, nside, ipix, theta, phi, double=double
;
; INPUTS:
;   nside   - healpix nside
;   ipix    - pixel numbers
; 
; OUTPUTS:
;   theta   - angle from north pole [radians]
;   phi     - longitude angle [radians]
;
; EXAMPLES:
;   
; COMMENTS:
;   Calling syntax is same as pix2ang_ring and agrees to machine 
;     precision. 
;   This routine has somewhat better performance when called for the
;     full sky than pix2ang_ring. 
;   
; REVISION HISTORY:
;   2003-Dec-06  Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
pro dpf_pix2ang_ring, nside, ipix, theta, phi, double=double

; -------- number in each ring

  t1 = systime(1)
  nside = long(nside)
  npix = 12L*nside*nside
  tind = lindgen(npix)
  nq = [(lindgen(nside)+1)*4, lonarr(2*nside-1)+(4*nside), $
        reverse((lindgen(nside)+1)*4)]
  thetaq = dblarr(4L*nside-1)
  bound = 2*nside*(1+nside)

  qind  = lindgen(nside)
  qindrev = nside-qind
  xoffs = 2L*nside*(nside+1)
; -------- ring0 is the lowest pixel number in each ring
  ring0 = [2*qind*(qind+1), xoffs+(4*nside)*lindgen(2*nside-1), $
           npix-2*qindrev*(qindrev+1)]

; -------- initialize arrays for output
  if keyword_set(double) then begin 
     theta = dblarr(npix, /nozero)
     phi   = dblarr(npix, /nozero)
  endif else begin 
     theta = fltarr(npix, /nozero)
     phi   = fltarr(npix, /nozero)
  endelse 


  np = tind LT bound
  w_np = where(np, n_np)
  if n_np gt 0 then begin 
     qind = lindgen(nside)
     thetaq[qind] = acos(1-((qind+1)/double(nside))^2/3.d)

     q = long((sqrt(0.25d + 0.5d*tind[w_np])-0.5d))
     theta[w_np] = thetaq[q]
     phi[w_np] = (tind[w_np]-ring0[q]+0.5d)/(nq/(2*!dpi))[q]
  endif 

  w_eq = where((np eq 0) and (tind LT (npix-bound)), n_eq)
  if n_eq gt 0 then begin 
     qind = lindgen(2*nside-1)+nside
     thetaq[qind] = acos((2-(qind+1)/double(nside))*2.d/3.d)

     q = (tind[w_eq]-bound)/(4*nside)+nside
     theta[w_eq] = thetaq[q]

     phi[w_eq] = (tind[w_eq]-ring0[q]+(q AND 1)*0.5d)/(4*nside/(2*!dpi))[q]
  endif 

  w_sp = where(tind GE (npix-bound), n_sp)
  if n_sp gt 0 then begin 
     qind = lindgen(nside)+(3*nside-1)
     thetaq[qind] =acos(-1+((qind+1-4*nside)/double(nside))^2/3.d)

     q = (4*nside-2)-long(sqrt(0.5d*((0.5d + npix-1)-tind[w_sp]))-0.5d)
     theta[w_sp] = thetaq[q]
     phi[w_sp] = (tind[w_sp]-ring0[q]+0.5d)/(nq/(2*!dpi))[q]
  endif 

; -------- Now lookup descired pixels from table 

  theta = theta[ipix]
  phi = phi[ipix]

;  print, 'time', systime(1)-t1

  return
end
