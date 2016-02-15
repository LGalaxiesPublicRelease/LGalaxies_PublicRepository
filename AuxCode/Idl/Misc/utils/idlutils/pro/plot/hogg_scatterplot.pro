;+
; NAME:
;   hogg_scatterplot
; PURPOSE:
;   plot greyscale scatterplot with contours
; COMMENTS:
;   Doesn't overplot -- only plots.  This is because the IDL tvscl blots
;     out any other information on the plot.
;   Compares cumulated grid to the *total* weight -- ie, including points
;     outside the range (which is what you want; trust me).
; INPUTS:
;   x,y         - data values
; OPTIONAL INPUTS:
;   weight      - weighting for data points; default unity
;   xnpix       - width of greyscale grid in pixels; default 0.3*sqrt(N)
;   ynpix       - height of greyscale grid in pixels; same default
;   xrange      - x range; default minmax(x)
;   yrange      - y range; default minmax(y)
;   levels      - contour levels; default in source code
;   quantiles   - quantiles to plot on conditional plot; default
;                 [0.25,0.5,0.75]
;   cthick      - thickness for contours
;   exponent    - stretch greyscale at exponent power; default 1.0
;   satfrac     - fraction of pixels to saturate in greyscale; default 0
;   darkest     - darkest shade at saturation; default 127; lower darker
;   outpsym    - NEEDS DOCUMENTATION
;   outcolor   - NEEDS DOCUMENTATION
;   outsymsize - NEEDS DOCUMENTATION
;   grid        - input grid for running with "/usegrid"
;   [etc]       - extras passed to "plot" command
; KEYWORDS:
;   conditional     - normalize each column separately
;   labelcont       - label contours with numbers
;   internal_weight - use only the points in the image to determine
;                     contours
;   nocontours      - don't plot the contours
;   nogreyscale     - don't plot the greyscale
;   adqgreyscale    - do greyscale in the "ADQ" style (only works when
;                     conditional is set)
;   outliers        - NEEDS DOCUMENTATION
;   meanweight      - plot the mean of the weight values that land in
;                     that pixel, rather than the total; don't use
;                     with /conditional!
;   usegrid         - use input grid rather than compute from x,y;
;                     over-rules x,y inputs; needs matched
;                     xnpix,ynpix; brittle
; OPTIONAL OUTPUTS:
;   xvec        - [xnpix] vector of x values of grid pixel centers
;   yvec        - [ynpix] vector of y values of grid pixel centers
;   grid        - the greyscale grid [xnpix,ynpix] that was plotted
;   cumimage    - the cumulated grid [xnpix,ynpix] that was contoured
;   outquantiles - the plotted quantiles (when /conditional is set)
;   ioutliers   - indices of outlier points
; COMMENTS:
;   When output, the grid is in units of unit_weight, not in 
;   unit_weight per unit_x per unit_y (as you would want to do if 
;   you wanted to directly compare two results using different
;   resolution grids); the user will have to convert to that themselves.
; BUGS:
;   Doesn't check inputs.
;   Ought to specify saturation not as a fraction of pixels, but as a fraction
;     of the total weight (ie, saturate inside a particular, specifiable
;     confidence region).  This mod is trivial.
;   Ought to specify min and max grey levels, and contour colors.
;   Contour thicknesses hard-coded to unity.
; DEPENDENCIES:
;   hogg_histogram
;   plus much, much more
; REVISION HISTORY:
;   2002-12-04  written --- Hogg
;-
pro hogg_scatterplot, xxx,yyy,weight=weight, $
                      xnpix=xnpix,ynpix=ynpix, $
                      xrange=xrange,yrange=yrange, $
                      levels=levels,quantiles=quantiles, $
                      cthick=cthick, $
                      exponent=exponent, $
                      satfrac=satfrac, $
                      darkest=darkest, $
                      internal_weight=internal_weight, $
                      conditional=conditional, $
                      labelcont=labelcont, $
                      nocontours=nocontours,nogreyscale=nogreyscale, $
                      adqgreyscale=adqgreyscale, $
                      xvec=xvec,yvec=yvec,grid=grid, $
                      cumimage=cumimage,outquantiles=outquantiles, $
                      outliers=outliers, outpsym=outliers_psym, $
                      outcolor=outliers_color, outsymsize=outliers_symsize, $
                      ioutliers=ioutliers, $
                      meanweight=meanweight, $
                      usegrid=usegrid, $
  overplot=overplot,$
                      _EXTRA=KeywordsForPlot

if(n_params() lt 2) then begin
    print,'Syntax - hogg_scatterplot, x, y [, weight=, xnpix=, ynpix=, '
    print,'          xrange=, yrange=, levels=, quantiles=, cthick=, '
    print,'          exponent=, satfrac=, darkest=, /internal_weight, '
    print,'          /conditional, labelcont=, xvec=, yvec=, grid=, '
    print,'          cumimage= ]'
    print,'(also takes keywords associated with plot)'
    return
endif

; set defaults
ndata= n_elements(xxx)
if not keyword_set(weight) then weight= dblarr(ndata)+1.0
if not keyword_set(xnpix) then xnpix= ceil(0.3*sqrt(ndata)) > 10
if not keyword_set(ynpix) then ynpix= ceil(0.3*sqrt(ndata)) > 10
if not keyword_set(xrange) then xrange= minmax(xxx)
if not keyword_set(yrange) then yrange= minmax(yyy)
if not keyword_set(levels) then levels= errorf(0.5*(dindgen(3)+1))
nlevels= n_elements(levels)
if not keyword_set(quantiles) then quantiles= [0.25,0.5,0.75]
nquantiles= n_elements(quantiles)
if not keyword_set(satfrac) then satfrac= 0.0
if not keyword_set(exponent) then exponent= 1.0
if not keyword_set(darkest) then darkest= 127.0
if not keyword_set(cthick) then cthick= !P.THICK
if (n_elements(cthick) EQ 1) then begin
    if keyword_set(conditional) then cthick= replicate(cthick,nquantiles) $
    else cthick= replicate(cthick,nlevels)
endif
if keyword_set(adqgreyscale) then nogreyscale=1

; check inputs
; [tbd]

; cram inputs into correct form
x= reform(xxx,ndata)
y= reform(yyy,ndata)

; make axes
  plot, [0],[0],xrange=xrange,yrange=yrange,/xstyle,/ystyle, $
    _EXTRA=KeywordsForPlot,/nodata, noerase=keyword_set(overplot)
;if (not keyword_set(overplot)) then $ ; jm07mar
;  plot, [0],[0],xrange=xrange,yrange=yrange,/xstyle,/ystyle, $
;  _EXTRA=KeywordsForPlot,/nodata

; snap points to grid
deltax= (xrange[1]-xrange[0])/double(xnpix)
deltay= (yrange[1]-yrange[0])/double(xnpix)
xvec= xrange[0]+deltax*(dindgen(xnpix)+0.5)
yvec= yrange[0]+deltay*(dindgen(ynpix)+0.5)

; make and fill 1-d grid first, if necessary
if keyword_set(conditional) then begin
    colnorm= hogg_histogram(x,xrange,xnpix,weight=weight)
endif

; make and fill 2-d grid
; (this puts the grid in units of the weights, not per unitx per
; unity)
if (not keyword_set(usegrid)) then begin
    grid= hogg_histogram(transpose([[x],[y]]),[[xrange],[yrange]], $
                         [xnpix,ynpix],weight=weight)
    if keyword_set(meanweight) then begin
        dgrid= hogg_histogram(transpose([[x],[y]]), $
                              [[xrange],[yrange]],[xnpix,ynpix])
        grid= grid/(dgrid+(dgrid EQ 0.0))
    endif
endif

; renormalize columns, if necessary
if keyword_set(conditional) then begin
    zeroindx= where(grid EQ 0.0,nzeroindx)
    grid= grid/(colnorm#(dblarr(ynpix)+1))
    if nzeroindx GT 0 then grid[zeroindx]= 0.0
endif

; compute quantiles, if necessary
xgrid= floor(xnpix*(x-xrange[0])/(xrange[1]-xrange[0]))
if keyword_set(conditional) then begin
    outquantiles= dblarr(xnpix,nquantiles)
    for ii=0L,xnpix-1 do begin
        inii= where(xgrid EQ ii,ninii)
        if ninii GT 0 then begin
            outquantiles[ii,*]= weighted_quantile(y[inii],weight[inii], $
                                                  quant=quantiles)
        endif
    endfor

; otherwise cumulate image
endif else begin
    cumindex= reverse(sort(grid))
    cumimage= dblarr(xnpix,ynpix)
    cumimage[cumindex]= total(grid[cumindex],/cumulative) 
; renormalize the cumulated image so it really represents fractions of the
; *total* weight
    if(NOT keyword_set(internal_weight)) then $
      cumimage= cumimage/total(weight) $
    else $
      cumimage= cumimage/total(grid)
endelse

if keyword_set(outliers) then begin
    ilow=where(cumimage gt max(levels), nlow)
    if(nlow gt 0) then $
      grid[ilow]=0.
endif

; scale greyscale
if NOT keyword_set(nogreyscale) then begin
    mingrey= 255.0
    maxgrey= darkest
    maxgrid= grid[(reverse(sort(grid)))[ceil(satfrac*xnpix*ynpix)]]
    mingrid= 0.0
    tvgrid= mingrey+(maxgrey-mingrey)*((grid-mingrid)/(maxgrid-mingrid))^exponent
    tvgrid= (tvgrid < mingrey) > maxgrey

; plot greyscale; re-draw the axes
    tvimage, tvgrid, overplot=overplot, _EXTRA=KeywordsForPlot
    if (n_elements(KeywordsForPlot) ne 0L) then begin
       KeywordsForPlot1 = KeywordsForPlot
       if tag_exist(KeywordsForPlot1,'XTITLE') then KeywordsForPlot1.xtitle = ''
       if tag_exist(KeywordsForPlot1,'YTITLE') then KeywordsForPlot1.ytitle = ''
       if tag_exist(KeywordsForPlot1,'XTICKNAME') then KeywordsForPlot1.xtickname = replicate(' ',10)
       if tag_exist(KeywordsForPlot1,'YTICKNAME') then KeywordsForPlot1.ytickname = replicate(' ',10)
    endif
    plot, [0],[0],xrange=xrange,yrange=yrange,/xstyle,/ystyle, $
      _EXTRA=KeywordsForPlot1,/nodata,/noerase,color=djs_icolor('default')
;   tv, tvgrid,xrange[0],yrange[0],/data, $
;     xsize=(xrange[1]-xrange[0]),ysize=(yrange[1]-yrange[0])
endif

; plot quantiles, if necessary
if keyword_set(conditional) then begin
    if keyword_set(adqgreyscale) then begin
        for ii=1L,nquantiles-1L do begin
            shade= floor(256.999*(0.5+0.5*(abs(quantiles[ii-1]-0.5)+ $
                                           abs(quantiles[ii]-0.5))))
            for jj=0L,xnpix-1L do begin
                xf= xvec[jj]+0.5*[-deltax,-deltax,deltax,deltax]
                yf= [outquantiles[jj,ii-1],outquantiles[jj,ii], $
                     outquantiles[jj,ii],outquantiles[jj,ii-1]]
                polyfill, xf,yf,color=shade,noclip=0
            endfor
        endfor
    endif
    for ii=0L,nquantiles-1 do begin
        oplot, xvec,outquantiles[*,ii],psym=10,thick=cthick[ii]
    endfor

; otherwise overplot contours
endif else begin
    if NOT keyword_set(labelcont) then labelcont=0
    if (NOT keyword_set(nocontours)) then begin
        contour, cumimage,xvec,yvec,levels=levels,/overplot, $
          c_labels=lonarr(n_elements(levels))+labelcont,c_thick=cthick
    endif
endelse

if keyword_set(outliers) then begin
    if(n_elements(outliers_psym) eq 0) then $
      outliers_psym=6
    if(n_elements(outliers_color) eq 0) then $
      outliers_color='red'
    if(n_elements(outliers_symsize) eq 0) then $
      outliers_symsize=0.13
    iin=where((x-xrange[0])/(xrange[1]-xrange[0]) gt 0 and $
              (x-xrange[1])/(xrange[1]-xrange[0]) lt 0 and $
              (y-yrange[0])/(yrange[1]-yrange[0]) gt 0 and $
              (y-yrange[1])/(yrange[1]-yrange[0]) lt 0, nin)
    if(nin gt 0) then begin
        pval=interpolate(cumimage, $
                         (x[iin]-xvec[0])/(xvec[xnpix-1L]-xvec[0])* $
                         float(xnpix-1L), $
                         (y[iin]-yvec[0])/(yvec[ynpix-1L]-yvec[0])* $
                         float(ynpix-1L))
        ioutliers=where(pval gt max(levels), noutliers)
        if(noutliers gt 0) then begin
            ioutliers= iin[ioutliers]
            djs_oplot, x[ioutliers], y[ioutliers], $
              psym=outliers_psym, color=outliers_color, $
              symsize=outliers_symsize
        endif
    endif
endif

;; re-plot axes (yes, this is a HACK)
;if (not keyword_set(overplot)) then begin
;;  !P.MULTI[0]= !P.MULTI[0]+1
;   plot, [0],[0],xrange=xrange,yrange=yrange,/xstyle,/ystyle, $
;     _EXTRA=KeywordsForPlot,/nodata,/noerase
;endif

end
