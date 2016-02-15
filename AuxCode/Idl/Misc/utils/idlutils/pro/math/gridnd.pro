;+
; NAME:
;   gridnd
; PURPOSE:
;   determine positions on a grid, for easy access
; CALLING SEQUENCE:
;   gridnd, x [, ix=, grid=, ngrid=, binsize=, igrid= ]
; INPUTS:
;   x - [M,N] positions in M-dimensions
; OPTIONAL INPUTS:
;   binsize - grid size (default 1.)
; OUTPUTS:
;   grid - [nx[0]*nx[1]*...*nx[M-1]] set of pointers to particles in
;          each cell
;   ngrid - [nx[0]*nx[1]*nx[2]*...*nx[M-1]] number of particles in
;           each cell
;   igrid - [N] position in grid
;   xminmax - [2,M] limits of grid
;   nx - [M] dimensions of grid
; COMMENTS:
;   This code is BETA! Use at your own risk.
; REVISION HISTORY:
;   12-Oct-2005  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro gridnd, x, ix=ix, binsize=binsize, grid=grid, ngrid=ngrid, $
            xminmax=xminmax, nx=nx, igrid=igrid, nd=nd

if(NOT keyword_set(binsize)) then binsize=1.

if((size(x))[0] eq 1) then begin
    if(NOT keyword_set(nd)) then begin
        mm=1
        nn=n_elements(x)
        x1=reform(x, 1, nn)
    endif else begin
        mm=nd
        nn=n_elements(x)/mm
        x=reform(x, mm, nn)
    endelse 
endif else begin
    mm=(size(x,/dim))[0]
    nn=(size(x,/dim))[1]
endelse

if(NOT keyword_set(xminmax)) then begin
    xminmax=fltarr(2,mm)	
    nx=lonarr(mm)
    for i=0L, mm-1L do begin
        xminmax[*,i]=minmax(x[i,*])
        xminmax[0,i]=float(long(xminmax[0,i]/binsize)-1.)*binsize
        xminmax[1,i]=float(long(xminmax[1,i]/binsize+1L)+1.)*binsize
        nx[i]=long((xminmax[1,i]-xminmax[0,i])/binsize)
    endfor
endif

ix=long(x)*0L
for i=0L, mm-1L do $
  ix[i,*]=long(float(nx[i])*(x[i,*]-xminmax[0,i])/ $
               (xminmax[1,i]-xminmax[0,i]))

ntot=1L
for i=0L, mm-1L do ntot=ntot*nx[i]

igrid=lonarr(nn)
ngrid=lonarr(ntot)
grid=ptrarr(ntot)
ibase=1L
for i=0L, mm-1L do begin
  igrid=igrid+ibase*ix[i,*]
  ibase=ibase*nx[i]
endfor
isort=sort(igrid)
iuniq=uniq(igrid[isort])
istart=0L
for i=0L, n_elements(iuniq)-1L do begin
    iend=iuniq[i]
    icurr=isort[istart:iend]
    if(igrid[icurr[0]] ge 0L and $
       igrid[icurr[0]] le ntot-1L) then begin
        ngrid[igrid[icurr[0]]]=n_elements(icurr)
        grid[igrid[icurr[0]]]=ptr_new(icurr[sort(icurr)])
    endif
    istart=iend+1L
endfor

end
