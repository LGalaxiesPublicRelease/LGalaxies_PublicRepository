;-----------------------------------------------------------------------
;+
; NAME:
;   hogg_oplot_covar
; PURPOSE:
;   Wrapper for oplot that plots covariance ellipses.
; CALLING SEQUENCE:
;   hogg_oplot_covar, x,y,covar,...
; INPUT:
;   x,y     - [N] arrays of points
;   covar   - [2,2,N] array of symmetric covariance matrices
; KEYWORDS:
;   fill    - make filled rather than outlined ellipses
; OPTIONAL INPUTS:
;   nsigma  - plot n-sigma ellipse (default to 1)
;   color   - color or vector of [N] colors
;   thick   - line thickness or vector of [N] line thicknesses
;   [keywords for IDL plot]
; COMMENTS:
;   Allows color to be a [N] array.
; REVISION HISTORY:
;   2002-04-11  written by Hogg
;-
;-----------------------------------------------------------------------
pro hogg_oplot_covar,x,y,incovar,color=color,nsigma=nsigma,fill=fill, $
                     thick=thick,_EXTRA=KeywordsForPlot

; set defaults, etc
if NOT keyword_set(color) then color= !p.color
if NOT keyword_set(thick) then thick= !p.thick
ndata= n_elements(x)
x= reform([x],ndata)
y= reform([y],ndata)
if (n_elements(incovar) EQ 2L*2L*ndata) then begin
    covar= reform(incovar,2,2,ndata)
endif else begin
    covar= dblarr(2,2,ndata)
    for ii=0,ndata-1 do covar[*,*,ii]= incovar
endelse
ncolor= N_elements(color)
if ncolor EQ 1 then icolor= intarr(ndata)+djs_icolor(color) $
else icolor= djs_icolor(color)
nthick= N_elements(thick)
if nthick EQ 1 then ithick= fltarr(ndata)+thick $
else ithick= thick
if NOT keyword_set(nsigma) then nsigma= 1D0

; set up useful vectors
npoint= 100
theta= 2D0*!PI*dindgen(npoint+2)/double(npoint)
xx= cos(theta)
yy= sin(theta)

; loop over points
for ii=0L,ndata-1L do begin

; symmetrize and eigensolve covariance matrix
    tcovar= double(covar[*,*,ii])
    infinite= where(finite(tcovar) EQ 0,ninf)
    if (ninf GT 0) then tcovar[*,*]= 0D0
    tcovar= 5d-1*(tcovar+transpose(tcovar))
    eval= sqrt(eigenql(tcovar,eigenvectors=evec,/double))

; make ellipses
    xxx= x[ii]+xx*nsigma*eval[0]*evec[0,0]+yy*nsigma*eval[1]*evec[0,1]
    yyy= y[ii]+xx*nsigma*eval[0]*evec[1,0]+yy*nsigma*eval[1]*evec[1,1]
    if keyword_set(fill) then $
      polyfill, xxx,yyy,color=icolor[ii],noclip=0 $
    else $
      oplot, xxx,yyy,color=icolor[ii],thick=ithick[ii],_EXTRA=KeywordsForPlot
endfor
return
end
