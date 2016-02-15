;+
; NAME:
;   hogg_plothist
; PURPOSE:
;   plot histogram of weighted points
; INPUTS:
;   x           - data values
; OPTIONAL INPUTS:
;   weight      - weighting for data points; default unity
;   npix        - number of bins in range
;   xrange      - x range; default minmax(x)
;   yrange      - y range; default chosen dumbly!
;   [etc]       - extras passed to "plot" command
; KEYWORDS:
;   overplot    - overplot, don't plot anew
;   ploterr     - plot Poisson error bars too
;   log         - take log_10 before plotting
;   meanweight  - plot the mean of the weights in each bin rather than
;                 the total of the weights; don't divide by binwidth
;   totalweight - plot the total weight; don't divide by binwidth
;   dontplot    - just calculate, do not plot anything
; OPTIONAL OUTPUTS:
;   xvec        - [npix] vector of x values of grid pixel centers
;   hist        - the histogram itself (ie, the total weight in each
;                 bin divided by the binwidth).
;   err         - the Poisson uncertainties on each point (ie, the
;                 sqrt of the sum of the squares of the weights,
;                 divided by the binwidth).
; COMMENTS:
;   Divides total weight in each bin by binwidth, unless /totalweight
;   or /meanweight is set.
; BUGS:
;   Doesn't check inputs.
;   Super-slow!
; REVISION HISTORY:
;   2002-12-14  written -- Hogg
;-
pro hogg_plothist, x,weight=weight, $
                   xrange=xrange,yrange=yrange,npix=npix, $
                   linestyle=linestyle, $
                   xvec=xvec,hist=hist,err=err, $
                   overplot=overplot,ploterr=ploterr,log=log, $
                   meanweight=meanweight,totalweight=totalweight, $
                   dontplot=dontplot, _EXTRA=KeywordsForPlot

; set defaults
ndata= n_elements(x)
if not keyword_set(weight) then weight= dblarr(ndata)+1.0
if not keyword_set(npix) then npix= ceil(0.3*sqrt(ndata)) > 10
if not keyword_set(xrange) then begin
    if keyword_set(overplot) then xrange= !X.CRANGE else xrange= minmax(x)
endif

; check inputs
; [tbd]

; snap points to grid
xvec= xrange[0]+(xrange[1]-xrange[0])*(dindgen(npix)+0.5)/double(npix)
xgrid= reform(floor(npix*(x-xrange[0])/(xrange[1]-xrange[0])),ndata)

; make and fill histogram
num= dblarr(npix)
hist= dblarr(npix)
err= dblarr(npix)
inxgrid= where(xgrid GE 0 AND xgrid LT npix,ninxgrid)
for ii=0L,ninxgrid-1 do begin
    num[xgrid[inxgrid[ii]]]= num[xgrid[inxgrid[ii]]] $
      +1.0
    hist[xgrid[inxgrid[ii]]]= hist[xgrid[inxgrid[ii]]] $
      +weight[inxgrid[ii]]
    err[xgrid[inxgrid[ii]]]= err[xgrid[inxgrid[ii]]] $
      +(weight[inxgrid[ii]])^2
endfor
if keyword_set(meanweight) then begin
    hist= hist/(num+(num eq 0.0))
    err= err/(num+(num eq 0.0))^2
endif
if ((NOT keyword_set(meanweight)) AND $
    (NOT keyword_set(totalweight))) then begin
    hist= hist*double(npix)/abs(xrange[1]-xrange[0])
    err= sqrt(err)*double(npix)/abs(xrange[1]-xrange[0])
endif

; take log?
if keyword_set(log) then begin
    err= (err/hist)/alog(10.0)
    hist= alog10(hist)

; set y range
    if not keyword_set(yrange) then yrange= [-1.5,0.2]*max(hist)
endif else begin
    if not keyword_set(yrange) then yrange= [-0.1,1.1]*max(hist)
endelse

; plot
if (NOT keyword_set(dontplot)) then begin
    if NOT keyword_set(overplot) then begin
        plot, xrange,0.0*xrange,psym=0, $
          xrange=xrange,yrange=yrange,/xstyle,/ystyle, $
          _EXTRA=KeywordsForPlot,thick=1
    endif
    oplot, xvec,hist,psym=10,linestyle=linestyle
    if keyword_set(ploterr) then $
      djs_oploterr, xvec,hist,yerr=err,psym=0
endif

end
