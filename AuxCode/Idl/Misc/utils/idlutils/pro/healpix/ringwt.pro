; take a stab at computing ring weights for healpix integrations
; Finkbeiner 5 Dec

pro ringwt, theta, wt

  nside = 32L
  lmax = 2*nside
  npix = 12L*nside*nside
  healgen, nside, theta, phi, /double
  
  gorski = healpix_ring_weight(nside)

  uind = uniq(theta)
  thetaq = theta[uind]
  nq = [(lindgen(nside)+1)*4, lonarr(2*nside-1)+(4*nside), $
        reverse((lindgen(nside)+1)*4)]
  cq = cos(thetaq)


  lfac = sqrt(2*lindgen(lmax+1)+1)

  ct = cos(theta)
  wt = dblarr(n_elements(theta))+1
  fn0 = legendre(ct,0,/doub)*lfac[0]
  fn1 = legendre(ct,1,/doub)*lfac[1]


  for k=1, 15 do begin 
     for l=2L, lmax do begin 
        leg = legendre(ct,l,/doub)*lfac[l]

        if ((l mod 2) eq 0) then begin 
           coeff = total(leg*fn0*wt)/npix
           wt = wt-coeff*leg/fn0
        endif 
        
        print, 'even', l, coeff
     endfor 
     for l=2L, lmax do begin 
        leg = legendre(ct,l,/doub)*lfac[l]

        if ((l mod 2) eq 0) then begin 
           coeff = total(leg*fn1*wt)/npix
           wt = wt-coeff*leg/fn1
        endif 
 
;        fn = ((l mod 2) eq 0) ? fn0 : fn1
;        coeff = total(leg*fn*wt)/npix
;        wt = wt-coeff*leg/fn
        
        print, l, coeff
     endfor 
     splot, wt[uind], /yno, xr=[0, nside]
     soplot, gorski, color='red'
  endfor 

  for i=1, 10 do begin 

     l1 = round(randomu(iseed)*lmax)
     l2 = round(randomu(iseed)*lmax)
     coeff = total(legendre(ct, l1, /doub)*legendre(ct, l2, /doub))/npix $
       *lfac[l1]*lfac[l2]
     print, l1, l2, coeff

  endfor 

  return
end

  
function polymatrix, nside, wt=wt


  lmax = 4*nside-2
  npix = 12L*nside*nside
  healgen, nside, theta, phi, /double
  
  uind = uniq(theta)
  thetaq = theta[uind]
  nq = [(lindgen(nside)+1)*4, lonarr(2*nside-1)+(4*nside), $
        reverse((lindgen(nside)+1)*4)]

  if keyword_set(wt) then nq = nq*wt

  cq = cos(thetaq)

  lfac = sqrt(2*lindgen(lmax+1)+1)
  M = dblarr(n_elements(thetaq), lmax+1)

  for l=0L, lmax do begin 
     M[*, l]= legendre(cq, l, /doub)*lfac[l] * sqrt(nq)
  endfor 

  return, M
end


pro tryit

  nside = 16L
  wt=healpix_ring_weight(nside)
  lmax = 4*nside
  npix = 12L*nside*nside
  m=polymatrix(nside, wt=wt)

  p1 = m#transpose(m)/npix                ; ntheta by ntheta
  p2 = transpose(m)#m/npix           ; npoly by npoly

  ii = lindgen(4L*nside-1)


  return
end
