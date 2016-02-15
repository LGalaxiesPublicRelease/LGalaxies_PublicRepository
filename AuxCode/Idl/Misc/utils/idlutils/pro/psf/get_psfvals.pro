
; pass in psfs
; return components, coefficients
pro stamps2pca, psfs, comp=comp, coeff=coeff, recon=recon

  sz = size(psfs, /dimen)

  data0 = transpose(reform(psfs, sz[0]*sz[1], sz[2]))
  means = total(data0, 2)/(sz[0]*sz[1])

  data = data0 ; make a copy
  allcomp = transpose(pcomp(data, /double, /standard))
  
  ncomp = 9   ; keep ncomp components
  comp = allcomp[*, 0:ncomp-1]
  norm = 1./sqrt(total(comp^2, 1))
  for i=0L, ncomp-1 do comp[*, i] = comp[*, i]*norm[i]

  coeff = comp##data0

  if arg_present(recon) then begin 
     recon = reform(transpose(coeff)##comp, sz[0], sz[1], sz[2])
     for i=0L, sz[2]-1 do recon[*, *, i] = recon[*, *, i]+means[i]
  endif

  return
end


function spread_stack, stack, npad

  sz = size(stack, /dim)
  box = sz[0]
  if sz[1] NE sz[0] then stop
  np = sz[2]

  if NOT keyword_set(npad) then npad = 0
  mask = bytarr(box, box)
  mask[npad:box-npad-1, npad:box-npad-1] = 1B

  nx   = ceil(sqrt(np))
  npix = nx*box
  arr  = fltarr(npix, npix)

  for i=0, np-1 do begin 
     stamp = stack[*, *, i]

     ix = i mod nx
     iy = i / nx
     arr[ix*box:(ix+1)*box-1, iy*box:(iy+1)*box-1] = stamp*mask

  endfor

  return, arr
end


function psf_par

  par = {boxrad:  12, $  ; box radius of PSF cutout (2*boxrad+1, 2*boxrad+1)
         fitrad:  9, $   ; radius of region used in PSF fit
         cenrad:  1}   ; region used to center PSF
         
; nsigma ?
; nfaint?
; use ivar in psf_polyfit
; check condition of matrix in psf_polyfit


  return, par
end



pro stamp_renorm, stamp, sivar, par

  nstamp = (size(stamp, /dimen))[2]
  for istamp=0L, nstamp-1 do begin 
     norm = psf_norm(stamp[*, *, istamp], par.cenrad)
     stamp[*, *, istamp] = stamp[*, *, istamp]/norm
     sivar[*, *, istamp] = sivar[*, *, istamp]*norm^2
  endfor 

  return
end



pro psfplot, sub, x, coeff=coeff, bin=bin

  binfac = keyword_set(bin) ? 2:1
  x0 = 6*binfac
  x1 = 12*binfac

  w=where(abs(sub[9*binfac,9*binfac,*]) lt 1)

  yr = [-1, 1]*0.005

  !p.multi = [0, x1-x0+1, x1-x0+1]
  xx = findgen(max(x))
  
  for j=x1, x0, -1 do for i=x0, x1 do begin 
     plot, x[w], sub[i, j, w], ps=3, yr=yr, chars=.1
     if keyword_set(coeff) then oplot, xx, poly(xx, coeff[i, j, *])
  endfor 

  return
end



function mockimage, psf

  pad = 20
  nx = 2048
  ny = 1489
  image = fltarr(nx, ny)
  npix = (size(psf, /dim))[0]
  brad = (npix-1)/2


  for i=0, 200 do begin
     stamp = psf*randomu(iseed)*100
     dcen = randomu(iseed, 2)-0.5
     stamp = sshift2d(stamp, -dcen)
     stamp = stamp+sqrt(stamp > 0)*randomn(iseed, npix, npix)
     cen = randomu(iseed, 2)*([nx, ny]-2*pad)+pad
     image[cen[0]-brad:cen[0]+brad, cen[1]-brad:cen[1]+brad] = $
       image[cen[0]-brad:cen[0]+brad, cen[1]-brad:cen[1]+brad] + stamp
  endfor

  image = image+randomn(iseed, nx, ny)*0.1

  return, image
end



; subtract PSF stars from image and see how well we did!
function image_psf_sub, image_in, psfs, px, py, dx, dy

  image = image_in
  rad = 1
  sz = size(image, /dim)
  boxrad = 9
  box    = boxrad*2+1
  nstar  = n_elements(px) 

  x0 = (px-boxrad) > 0
  x1 = (px+boxrad) < (sz[0]-1)
  y0 = (py-boxrad) > 0
  y1 = (py+boxrad) < (sz[1]-1)
  
  sx0 = x0-px+boxrad
  sx1 = x1-px+boxrad
  sy0 = y0-py+boxrad
  sy1 = y1-py+boxrad

  for i=0L, nstar-1 do begin 
     stamp    = psfs[*, *, i] ; assume this is already sinc shifted
     model    = stamp[sx0[i]:sx1[i], sy0[i]:sy1[i]]
     subimage = image[x0[i]:x1[i], y0[i]:y1[i]] 
     poly_iter, model, subimage, 1, 3, fit

     image[x0[i]:x1[i], y0[i]:y1[i]] = image[x0[i]:x1[i], y0[i]:y1[i]] -fit


  endfor

  return, image
end



function psf_chisq, stamps, stampivar, par, dx=dx, dy=dy

  nstamp = (size(stamps, /dimen))[2]
  if ~(keyword_set(dx) && keyword_set(dy)) then begin
     dx = fltarr(nstamp)
     dy = fltarr(nstamp)
  endif 

  if n_elements(dx) NE nstamp then message, 'index bug!'

  npix = (size(stamps, /dimen))[0]
  npad = par.boxrad-par.fitrad
  mask = bytarr(npix, npix)
  mask[npad:npix-npad-1, npad:npix-npad-1] = 1B
  chisq = fltarr(nstamp)

  for i=0L, nstamp-1 do begin 
     w = where(mask AND (stampivar[*, *, i] NE 0), nmask)
     shiftstamp =  sshift2d(stamps[*, *, i], [dx[i], dy[i]])
     chisq[i] = total(shiftstamp[w]^2 * (stampivar[*, *, i])[w])/nmask
  endfor

  return, chisq
end


; doesn't use psf to get radius -- maybe should!
; changes stamps and sub!!!
pro psf_zero, stamps, stampivar, sub, psf, par, faint, chisq

  nstamp = (size(stamps, /dimen))[2]
  npix = (size(stamps, /dimen))[0]
  rad = 3

;  chisq = fltarr(nstamp)
  xbox = djs_laxisnum([npix, npix], iaxis = 0)
  ybox = djs_laxisnum([npix, npix], iaxis = 1)

  mask = bytarr(npix, npix)
  npad = par.boxrad-par.fitrad
  mask[npad:npix-npad-1, npad:npix-npad-1] = 1B

  xcen = (npix-1)/2
  ycen = (npix-1)/2
  mask = mask AND (((xbox-xcen)^2+(ybox-ycen)^2) GT rad^2)

  for i=0L, nstamp-1 do begin 
     fmask = mask
     for j=0L, faint[i].npsf-1 do begin 
        fmask = fmask AND ((xbox-faint[i].px[j])^2+ $
                           (ybox-faint[i].py[j])^2) GT rad^2
     endfor
     wmask = where(fmask, ngoodpix)
     if ngoodpix LT 10 then message, 'we have a problem'
     djs_iterstat, (sub[*, *, i])[wmask], mean=mean, sigma=sigma

     stamps[*, *, i] = stamps[*, *, i]-mean
     sub[*, *, i] = sub[*, *, i]-mean

  endfor

  return
end


pro chisq_cut, chisq, x, y, dx, dy, stamps, stampivar, psfs, arr
  
; -------- keep stars with chisq LT 1 but keep at least half of them.
  chisqval = median(chisq) > 1
  w = where(chisq LT chisqval, nstar)

  splog, 'Cutting to', nstar, ' stars.'

  chisq     = chisq[w]
  x         = x[w]
  y         = y[w]        
  dx        = dx[w]       
  dy        = dy[w]       
  stamps    = stamps[*, *, w]
  stampivar = stampivar[*, *, w]
  psfs      = psfs[*, *, w]     
  arr       = arr[*, *, w]      

  return
end

; input image information, fit paramters
; output PSf fit coeeficients and structure array of star positions used.

function psf_fit_coeffs, image, ivar, satmask, par, status=status

; -------- highest order fit to attempt
  ndeg = 3

; -------- timer
  t1 = systime(1)
  status = 0L

; -------- look for negative pixels
  wneg = where(image LE 0, nneg)
  if nneg GT 0 then begin 
     status = status OR 1
     ivar[wneg] = 0
  endif

; -------- determine badpixels mask (pixels not to use for stamps)
  badpixels = smooth(float(ivar EQ 0), par.fitrad*2+1, /edge) GT $
    ((1./par.fitrad)^2/10)

; -------- find some stars (not a complete list)
  psf_findstars, image, ivar, par.boxrad, clean, x, y, $
    satmask=satmask, badpixels=badpixels

; -------- cut out postage stamps around them
  stamps = psf_stamps(clean, ivar, x, y, par, /shift, dx=dx, dy=dy, $
                      stampivar=stampivar)
  nstamp = n_elements(dx) 

; -------- get median psf
  if nstamp GT 1 then begin
     median_psf = djs_median(stamps, 3)
     

     psfs0 = stamps
     for i=0L, nstamp-1 do psfs0[*, *, i] = median_psf

     psf_multi_fit, stamps, stampivar, psfs0, par, sub, faint, nfaint=3
     psf_zero, stamps, stampivar, sub, median_psf, par, faint, chisq
     stamp_renorm, stamps, stampivar, par

; -------- first chisq cut.  Some garbage (diffraction spikes, etc.)
;          may have leaked in - let's get it before we go on. 
; or not     
;     chisq = psf_chisq(stamps-psfs0, stampivar, par, dx=dx, dy=dy)
     
     psf_multi_fit, stamps, stampivar, psfs0, par, sub1, faint, nfaint=3
; maybe recenter here
     recon = sub1+psfs0
     
; -------- do spatial PSF fit
     cf = psf_polyfit(recon, stampivar, x, y, par, ndeg=ndeg, /reject, cond=cond)
     
; -------- if condition number is bad try again
     while max(cond) GT 1000 do begin 
        ndeg--
        splog, 'Condition Number:', max(cond), '  Falling back to ndeg =', ndeg
        if ndeg EQ -1 then stop
        cf = psf_polyfit(recon, stampivar, x, y, par, ndeg=ndeg, /reject, cond=cond)
     endwhile
     
     
     psfs1 = psf_eval(x, y, cf, par.cenrad)
     
     psf_multi_fit, stamps, stampivar, psfs1, par, sub2, faint, nfaint=3
     
;  stamps2pca, psfs1, comp=comp, coeff=coeff, recon=recon
     
     chisq = psf_chisq(sub2, stampivar, par, dx=dx, dy=dy)
     
     sub3 = sub2
     
     chisq_cut, chisq, x, y, dx, dy, stamps, stampivar, psfs1, sub3
     
     recon = sub3+psfs1
     cf = psf_polyfit(recon, stampivar, x, y, par, ndeg=ndeg, /reject, cond=cond)
     psfs2 = psf_eval(x, y, cf, par.cenrad)
  endif else begin 
     cond = 1.0
     cf = stamps
     ndeg = 0
     chisq = 0.0
     status = status OR 2
  endelse

print, '--------> CONDITION NUMBER ', max(cond)

; -------- maybe make room for some QA stuff in here...
  delvarx, result
  result = {nstar: n_elements(x), $
            x: x+dx, $
            y: y+dy, $
            chisq: chisq, $
            ndeg: ndeg, $
            coeff: cf, $
            boxrad: par.boxrad, $
            fitrad: par.fitrad, $
            cenrad: par.cenrad, $
            condition: max(cond), $
            status: status}

  splog, 'Time: ', systime(1)-t1
  
  return, result
end



pro callit

  run = 273
  camcol = 3
  field = 500

  run = 4828
  camcol = 3
  field = 258


; -------- get parameters
  par = psf_par()

; -------- read image
  fname = sdss_name('idR', run, camcol, field, filter='i')
  sdss_readimage, fname, image, ivar, satmask=satmask, /silent

; -------- do the fit
  pstr = psf_fit_coeffs(image, ivar, satmask, par)

; -------- see how we did
  chisq = psf_chisq(sub2, stampivar, par, dx=dx, dy=dy)

  sind = sort(chisq)
  atv, spread_stack(sub2[*, *, sind], 3)


; center up
; compute psf0 (median)
; subtract 3 faint stars
; get new zero
; recenter
; fit psf1
; subtract 5 faint stars
; get new zero  (with tilt?)
; fit psf2
; toss anything that doesn't fit well within fitrad. 

  return
end



pro validate_all
  
  flist = file_search('psC*fit')
  for i=0L, n_elements(flist)-1 do begin
     print, i, '   ', flist[i]
     pstr = mrdfits(flist[i], 1, /silent)
     psf_validate, pstr
  endfor

  return
end




; loop over some SDSS fields, see if we crash...
pro tryit

  run = 273
  camcol = 1
  fstart = sdss_fieldrange(run, fend=fend)

; -------- get parameters
  par = psf_par()

  for ifilt=0, 4 do begin 
     filtname = filtername(ifilt)
     
     for field=fstart, fend do begin 
        print, run, camcol, field

; -------- read image
        fname = sdss_name('idR', run, camcol, field, filter=filtname)
        sdss_readimage, fname, image, ivar, satmask=satmask, /silent

; -------- do the fit
        pstr = psf_fit_coeffs(image, ivar, satmask, par)

; -------- write outputs
        outname = string('psCoeff-', run, '-'+filtname, camcol, '-', $
                         field, '.fit', $
                         format='(A,I6.6,A,I1,A,I4.4,A)')
        splog, 'Writing ', outname
        mwrfits, pstr, outname, /create

     endfor
  endfor


  return
end



pro foobar

  raw = readfits('/scr/apache20/dfink/ctio/oct98/fits/n7obj/obj7100.fits')
  bias = readfits('bias.fits')
  flat = readfits('final_flat.fits')

  image = mask*(raw-bias)/(flat + (mask eq 0))



  return
end
