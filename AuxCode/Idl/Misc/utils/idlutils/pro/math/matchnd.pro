;+
; NAME:
;   matchnd
; PURPOSE:
;   match two sets of points in N dimensions
; CALLING SEQUENCE:
;   matchnd, x1, x2, distance [, m1=, m2=, d12=, nmatch= ]
; INPUTS:
;   x1 - [M,N1] positions in M-dimensions
;   x2 - [M,N2] positions in M-dimensions
;   distance - match distance
; OPTIONAL INPUTS:
;   maxmatch - MRB: please explain!
;   /silent  - don't splog anything
; OUTPUTS:
;   m1 - [nmatch] matches to x1
;   m2 - [nmatch] matches to x2
;   d12 - [nmatch] distance between matches
;   nmatch - number of matches
; COMMENTS:
;   This code is BETA! Use at your own risk.
; REVISION HISTORY:
;   12-Oct-2005  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro matchnd, x1, x2, distance, m1=m1, m2=m2, d12=d12, nmatch=nmatch, $
             maxmatch=maxmatch, nd=nd, silent=silent

if(n_elements(maxmatch) eq 0) then maxmatch=1

if((size(x1))[0] eq 1) then begin
    if(NOT keyword_set(nd)) then begin
        mm=1
        nn=n_elements(x1)
        x1=reform(x1, 1, nn)
    endif else begin
        mm=nd
        nn=n_elements(x1)/mm
        x1=reform(x1, mm, nn)
    endelse 
endif else begin
    mm=(size(x1,/dim))[0]
    nn=(size(x1,/dim))[1]
endelse

maxdiff=-1.
for i=0L, mm-1L do begin
    if(max(x1[i,*])-min(x1[i,*]) gt maxdiff) then $
      maxdiff=max(x1[i,*])-min(x1[i,*])
endfor
maxn=long((10.*n_elements(x1)+n_elements(x2))^(1./float(mm))) > $
  (long(1000000.^(1./float(mm))))
binsize=(maxdiff/float(maxn)) > (distance)
if (not keyword_set(silent)) then splog, binsize

distance2=distance^2

gridnd, x1, ix=ix1, binsize=binsize, grid=grid1, ngrid=ngrid1, nx=nx, $
  igrid=igrid1, xminmax=xminmax, nd=nd
gridnd, x2, ix=ix2, binsize=binsize, grid=grid2, ngrid=ngrid2, nx=nx, $
  igrid=igrid2, xminmax=xminmax, nd=nd

ifilled=where(ngrid1 gt 0, nfilled)
nmatch=0L
nmax=nn*100L
m1=lonarr(nmax)
m2=lonarr(nmax)
d12=dblarr(nmax)
for i=0L, nfilled-1L do begin
    i1tmp=*(grid1[ifilled[i]])
    ii1=lonarr(mm)
    itmp=ifilled[i]
    for j=0L, mm-1L do begin
        ii1[j]=itmp mod nx[j]
        itmp=itmp/nx[j]
    endfor
    iist=(ii1-1L) > 0L
    iind=(ii1+1L) < (nx-1L)
    ii=lonarr(mm)
    nii=1L
    for j=0L, mm-1L do $
      nii=(iind[j]-iist[j]+1L)*nii
    ii2=lonarr(mm)
    tmpnmatch=0L
    tmpm1=-1L
    tmpm2=-1L
    tmpd12=-1.
    for j=0L, nii-1L do begin
        igrid2=0L
        ibase=1L
        for k=0L, mm-1L do begin
            igrid2=igrid2+ibase*(ii2[k]+iist[k])
            ibase=ibase*nx[k]
        endfor
        
        if(ngrid2[igrid2] gt 0) then begin 
            i2tmp=*(grid2[igrid2])
            for k=0L, ngrid2[igrid2]-1L do begin
                dtmp=fltarr(mm, n_elements(i1tmp))
                for l=0L, mm-1L do $
                  dtmp[l,*]=(x1[l,i1tmp]-x2[l,i2tmp[k]])
                dtmp2=total(dtmp^2,1)
                imtmp=where(dtmp2 lt distance2, nmtmp)
                if(nmtmp gt 0) then begin
                    if(tmpnmatch gt 0) then begin
                        tmpm1=[tmpm1, i1tmp[imtmp]]
                        tmpm2=[tmpm2, replicate(i2tmp[k], nmtmp)]
                        tmpd12=[tmpd12, sqrt(dtmp2[imtmp])]
                    endif else begin
                        tmpm1=i1tmp[imtmp]
                        tmpm2=replicate(i2tmp[k], nmtmp)
                        tmpd12=sqrt(dtmp2[imtmp])
                    endelse
                    tmpnmatch=tmpnmatch+nmtmp
                endif
            endfor
        endif

        iidim=0L
        nnx=iind-iist+1L
        ii2[iidim]=(ii2[iidim]+1L) mod nnx[iidim]
        while(ii2[iidim] eq 0 and iidim lt mm-1L) do begin
            iidim=iidim+1L
            ii2[iidim]=(ii2[iidim]+1L) mod nnx[iidim]
        endwhile
    endfor

    if(tmpnmatch gt 0) then begin
        isort=sort(tmpm1)
        tmpd12=tmpd12[isort]
        tmpm1=tmpm1[isort]
        tmpm2=tmpm2[isort]
        iuniq=uniq(tmpm1)
        istart=0L
        for j=0L, n_elements(iuniq)-1L do begin
            iend=iuniq[j]
            icurr=istart+lindgen(iend-istart+1L)
            isort=sort(tmpd12[icurr])
            icurr=icurr[isort]
            if(n_elements(icurr) gt maxmatch AND $
               maxmatch gt 0) then $
                icurr=icurr[0:maxmatch-1L]
            if(nmatch+n_elements(icurr) gt nmax) then begin
                m1new=lonarr(nmax*2L)
                m2new=lonarr(nmax*2L)
                d12new=lonarr(nmax*2L)
                m1new[0:nmatch-1]=m1[0:nmatch-1]
                m2new[0:nmatch-1]=m2[0:nmatch-1]
                d12new[0:nmatch-1]=d12[0:nmatch-1]
                m1=m1new
                m2=m2new
                d12=d12new
                m1new=0
                m2new=0
                d12new=0
                nmax=nmax*2L
            endif
            m1[nmatch:nmatch+n_elements(icurr)-1L]=tmpm1[icurr]
            m2[nmatch:nmatch+n_elements(icurr)-1L]=tmpm2[icurr]
            d12[nmatch:nmatch+n_elements(icurr)-1L]=tmpd12[icurr]
            nmatch=nmatch+n_elements(icurr)
            istart=iend+1L
        endfor
    endif
endfor

if(nmatch gt 0) then begin
    m1=m1[0:nmatch-1]
    m2=m2[0:nmatch-1]
    d12=d12[0:nmatch-1]
    
    isort=sort(m1*(max(m2)+1L)+m2)
    m1=m1[isort]
    m2=m2[isort]
    d12=d12[isort]
endif

end
