;+
; NAME:
;   extract_profmean
; PURPOSE:
;   Extract a photoesque radial profile
; CALLING SEQUENCE
;   extract_profmean, image, center, profmean, proferr $
;     [,profradius=profradius] [,invvar=invvar] 
; INPUTS:
;   image       - [nx,ny] image
;   center      - [2] center of extraction (LONG)
; OPTIONAL INPUTS:
;   invvar      - [nx,ny] inverse variance image; default to unity
;   profradius  - [nrad] defining profile, in pixels; default to PHOTO aps
; KEYWORDS:
; OUTPUTS:
;   nprof         - number of "good" profmean values (based on image
;                   size alone)
;   profmean      - [nrad] annular fluxes
;   profmean_ivar - [nrad] uncertainties
;   qstokes       - if present, calculates Stokes Q parameter in each circle
;                   (not within annuli)
;   ustokes       - if present, calculates Stokes U parameter in each circle
;                   (not within annuli)
; OPTIONAL INPUT/OUTPUTS:
;   cache         - cache storing photfracs for re-use
; COMMENTS:
;   Image must be centered on the center pixel. 
; DEPENDENCIES:
;   idlutils
; BUGS:
; REVISION HISTORY:
;   2002-09-04  Written - Blanton
;   2002-09-12  Modified to use djsphot - Hogg
;-
pro extract_profmean, image, center, profmean, profmean_ivar, $
                      qstokes=qstokes, ustokes=ustokes, $
                      profradius=profradius, invvar=invvar, $
                      nprof=nprof, cache=cache, area=area

; set defaults
if NOT keyword_set(profradius) then $
  profradius= [  0.564190,   1.692569,   2.585442,   4.406462, $
                 7.506054,  11.576202,  18.584032,  28.551561, $ 
                 45.503910,  70.510155, 110.530769, 172.493530, $
                 269.519104, 420.510529, 652.500061]
if NOT keyword_set(invvar) then invvar= 0.0*image+1.0

; choose precise and imprecise regions
cutoff= 50.0
iexact= where(profradius LT cutoff,nexact)
iquick= where(profradius GE cutoff,nquick)

; create output arrays
nrad= n_elements(profradius)
center= reform(center,2)

; check if we have cached the profradius stuff, if not remake  
if(n_tags(cache) eq 0) then begin
    for i=0, nexact-1 do begin
        dimen=(2L*(long(profradius[iexact[i]]*2.))+1L)>11L
        exact_photfrac, dimen/2L, dimen/2L, profradius[iexact[i]], $
          xdimen=dimen, ydimen=dimen, fracs=fracs, xpixnum=xpixnum, $
          ypixnum=ypixnum, pixnum=pixnum
        cache1={dimen:dimen, $
                radius:profradius[iexact[i]], $
                fracs:ptr_new(fracs), $
                mxx:ptr_new(fltarr(n_elements(fracs))), $
                myy:ptr_new(fltarr(n_elements(fracs))), $
                mxy:ptr_new(fltarr(n_elements(fracs))), $
                xpixnum:ptr_new(xpixnum), $
                ypixnum:ptr_new(ypixnum), $
                pixnum:ptr_new(pixnum) $
               }
        xx2=(float(xpixnum-(dimen/2L)))^2
        yy2=(float(ypixnum-(dimen/2L)))^2
        xy2=(float(ypixnum-(dimen/2L)))*(float(xpixnum-(dimen/2L)))
        rpix2=xx2+yy2
        *cache1.mxx=xx2/rpix2
        *cache1.myy=yy2/rpix2
        *cache1.mxy=xy2/rpix2
        ii=where(rpix2 eq 0.,iicount)
        if(iicount gt 0) then  begin
            (*cache1.mxx)[ii]=1.
            (*cache1.myy)[ii]=1.
            (*cache1.mxy)[ii]=1.
        endif
        if(n_tags(cache) eq 0) then $
          cache=cache1 $
        else $
          cache=[cache,cache1]
    endfor
    for i=0, nquick-1 do begin
        dimen=(2L*(long(profradius[iquick[i]]*2.))+1L)>11L
        quick_photfrac, dimen/2L, dimen/2L, profradius[iquick[i]], $
          xdimen=dimen, ydimen=dimen, fracs=fracs, xpixnum=xpixnum, $
          ypixnum=ypixnum, pixnum=pixnum
        cache1={dimen:dimen, $
                radius:profradius[iquick[i]], $
                fracs:ptr_new(fracs), $
                mxx:ptr_new(fltarr(n_elements(fracs))), $
                myy:ptr_new(fltarr(n_elements(fracs))), $
                mxy:ptr_new(fltarr(n_elements(fracs))), $
                xpixnum:ptr_new(xpixnum), $
                ypixnum:ptr_new(ypixnum), $
                pixnum:ptr_new(pixnum) $
               }
        xx2=(float(xpixnum-(dimen/2L)))^2
        yy2=(float(ypixnum-(dimen/2L)))^2
        xy2=(float(ypixnum-(dimen/2L)))*(float(xpixnum-(dimen/2L)))
        rpix2=xx2+yy2
        *cache1.mxx=xx2/rpix2
        *cache1.myy=yy2/rpix2
        *cache1.mxy=xy2/rpix2
        ii=where(rpix2 eq 0.,iicount)
        if(iicount gt 0) then  begin
            (*cache1.mxx)[ii]=1.
            (*cache1.myy)[ii]=1.
            (*cache1.mxy)[ii]=1.
        endif
        if(n_tags(cache) eq 0) then $
          cache=cache1 $
        else $
          cache=[cache,cache1]
    endfor
endif

; now for each circular region:
nx=(size(image,/dim))[0]
ny=(size(image,/dim))[1]
area=fltarr(n_elements(cache)+1L)
profflux=fltarr(n_elements(cache)+1L)
if(arg_present(profmean_ivar)) then $
  profflux_ivar=fltarr(n_elements(cache)+1L)
if(arg_present(ustokes) OR arg_present(qstokes)) then begin
    ustokes=fltarr(n_elements(cache))
    qstokes=fltarr(n_elements(cache))
endif
for i=0L, n_elements(cache)-1L do begin
    xpixnum=(*cache[i].xpixnum)-cache[i].dimen/2L+center[0]
    ypixnum=(*cache[i].ypixnum)-cache[i].dimen/2L+center[1]
    ikeep=where(xpixnum ge 0 and $
                xpixnum lt nx and $
                ypixnum ge 0 and $
                ypixnum lt ny, nkeep)
    if(nkeep eq 0) then begin
        area[i+1L]=0.
        profflux[i+1L]=0.
    endif else begin
        pixnum= ypixnum[ikeep]*nx+xpixnum[ikeep]
        fracs= (*cache[i].fracs)[ikeep]
        mxx= (*cache[i].mxx)[ikeep]
        myy= (*cache[i].myy)[ikeep]
        mxy= (*cache[i].mxy)[ikeep]
        iuse=where(invvar[pixnum] gt 0., nuse)
        if(nuse gt 0) then begin
            area[i+1L]=total(fracs[iuse],/double)
            profflux[i+1L]=total(fracs*image[pixnum[iuse]],/double)
            if(n_elements(profflux_ivar) gt 0) then $
              profflux_ivar[i+1L]=1./total(fracs/invvar[pixnum[iuse]],/double)
            if(arg_present(ustokes) OR arg_present(qstokes)) then begin
                mxx=total(fracs*image[pixnum[iuse]]*mxx[iuse], /double)/ $
                  profflux[i+1L]
                myy=total(fracs*image[pixnum[iuse]]*myy[iuse], /double)/ $
                  profflux[i+1L]
                mxy=total(fracs*image[pixnum[iuse]]*mxy[iuse], /double)/ $
                  profflux[i+1L]
                ustokes[i]=2.*mxy
                qstokes[i]=mxx-myy
            endif
        endif
    endelse
endfor

; deal with annuli
for i=nrad,1L,-1L do begin
    area[i,*]= area[i,*]-area[i-1,*]
    profflux[i,*]= profflux[i,*]-profflux[i-1,*]
    if(n_elements(profflux_ivar) gt 0) then $
      profflux_ivar[i,*]= 1./(1./profflux_ivar[i,*]-1./profflux_ivar[i-1,*])
endfor

radii=[0,cache.radius]
fullarea=!DPI*(radii[1:n_elements(cache)]^2-radii[0:n_elements(cache)-1L]^2)
inotenough=where(area[1:n_elements(cache)]/fullarea le 0.95, nnotenough)
ienough=where(area[1:n_elements(cache)]/fullarea gt 0.95, nenough)


profmean=fltarr(n_elements(cache))
if(arg_present(profmean_ivar)) then $
  profmean_ivar=fltarr(n_elements(cache))

nprof=n_elements(cache)
if(nnotenough gt 0) then begin 
    nprof=inotenough[0]
    if(n_elements(ustokes) gt 0) then ustokes[inotenough]=0.
    if(n_elements(qstokes) gt 0) then qstokes[inotenough]=0.
endif
for i=0, nenough-1 do begin
    profmean[ienough[i]]=profflux[ienough[i]+1]/area[ienough[i]+1]
    if(arg_present(profmean_ivar)) then $
      profmean_ivar[ienough[i]]= $
      profflux_ivar[ienough[i]+1]*area[ienough[i]+1]^2
endfor
area=area[1:n_elements(cache)]

; check for very negative annuli
if(n_elements(profmean_ivar) eq 0) then $
  profmean_ivar=fltarr(n_elements(cache))+1.
ineg=where(profmean[ienough]*sqrt(profmean_ivar[ienough]) lt -100., nneg)
if(nneg gt 0) then nprof=ineg[0]

if(NOT arg_present(cache)) then begin
    for i=0, n_elements(cache)-1 do begin
        ptr_free,cache[i].fracs
        ptr_free,cache[i].xpixnum
        ptr_free,cache[i].ypixnum
        ptr_free,cache[i].pixnum
        ptr_free,cache[i].mxx
        ptr_free,cache[i].myy
        ptr_free,cache[i].mxy
    endfor
    cache=0
endif

end
