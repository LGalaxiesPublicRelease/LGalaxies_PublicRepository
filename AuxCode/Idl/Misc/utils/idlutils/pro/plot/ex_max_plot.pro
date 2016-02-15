;+
; NAME:
;   ex_max_plot
; PURPOSE:
;   plot ex_max outputs or other N-dimensional data sets
; INPUTS:
;   weight       [N] array of data-point weights
;   point        [d,N] array of data points - N vectors of dimension d
;   amp          [M] array of gaussian amplitudes
;   mean         [d,M] array of gaussian mean vectors
;   var          [d,d,M] array of gaussian variance matrices
;   psfilename   name for PostScript file; if no filename is given, then the
;                  plots will simply be sent to the currently active device
; OPTIONAL INPUTS:
;   nsig         number of sigma for half-width of each plot; default 5
;   label        [d] array of axis labels; default 'x_i'
;   contlevel    confidence levels for contouring; defaults in source code
;   range        [2,d] array of plotting ranges
;   textlabel    [q] vector of text labels
;   textpos      [d,q] array of text label positions
;   box          [2,d] array of coordinates for d-dimensional box
;   xdims,ydims  indices of data dimensions to use on each x and y axis
;   npix_x,npix_y  number of pixels in x and y dimensions of each panel
;   sqrt         sqrt stretch on image
;   log          logarithmic stretch on image
;   axis_char_scale  size of characters on labels
;   overpoints   [d,P] set of points to overplot a red box on each panel
;   model_npix_factor  for gaussians, use this factor to scale down pixel size
;   panellabels  label string for each panel (size of xdims, ydims arrays)
;   quantfrac    vector of fractions at which to plot quantiles on conditional
;                   panels
;   default_font  font command to send to set font for plotting
;   pthick       thickness of plot lines
;   yrangefudge  fudge range on histograms (default 1.)
; KEYWORDS:
;   nomodel      don't show model fits as greyscales or histograms
;   nodata       don't show data as greyscales or histograms
;   noellipse    don't show ellipses on 2-d plots
;   bw           show ellipses in b/w
;   conditional  plot the conditional distribution of y on x 
;   nogreyscale  don't plot the greyscale
; OUTPUTS:
; OPTIONAL OUTPUTS:
;   quantile     output array of quantile positions; read the source code
;   quantneff    the effective number of data points contributing to the medians
;   twodimages   set of all of the 2-dimensional projections 
;                (doesn't save the gaussian model unless
;                model_npix_factor==1)
; BUGS:
;   Greyscales need work?
; DEPENDENCIES:
; REVISION HISTORY:
;   2001-10-22  written - Hogg
;   2002-03-22  many added features - Blanton
;-
pro ex_max_plot, weight,point,amp,mean,var,psfilename,nsig=nsig, $
                 label=label,contlevel=contlevel,nomodel=nomodel, $
                 nodata=nodata, noellipse=noellipse,range=range,bw=bw, $
                 textlabel=textlabel,textpos=textpos,box=box,sqrt=sqrt, $
                 log=log,weight2=weight2,point2=point2, $
                 conditional=conditional, xdims=xdims, ydims=ydims, $
                 axis_char_scale=axis_char_scale, npix_x=npix_x, $
                 npix_y=npix_y,overpoints=overpoints, $
                 model_npix_factor=model_npix_factor, $
                 panellabels=panellabels, panellabelpos=panellabelpos, $
                 sigrejimage=sigrejimage, paneluse=paneluse, $
                 quantfrac=quantfrac, quantile=quantile, quantneff=quantneff, $
                 default_font=default_font,twodimages=twodimages, $
                 psym_overpoints=psym_overpoints, $
                 sizepanellabel=sizepanellabel, $
                 pthick=pthick,nogreyscale=nogreyscale, $
                 yrangefudge=yrangefudge

; set defaults
if(NOT keyword_set(model_npix_factor)) then model_npix_factor= 4.0
if(NOT keyword_set(npix_x)) then npix_x= 32L
if(NOT keyword_set(npix_y)) then npix_y= 32L
if(NOT keyword_set(pthick)) then pthick= 2.
if(NOT keyword_set(yrangefudge)) then yrangefudge= 1.

; recast range
range = double(range)

; check dimensions
ndata= n_elements(weight)       ; N
ngauss= n_elements(amp)         ; M
dimen= n_elements(point)/n_elements(weight) ; d
splog, ndata,' data points,',dimen,' dimensions,',ngauss,' gaussians'

; which dimensions should we look at?
if(NOT keyword_set(xdims)) then xdims=lindgen(dimen)
if(NOT keyword_set(ydims)) then ydims=lindgen(dimen)
if(NOT keyword_set(default_font)) then default_font='!3'
xdimen=n_elements(xdims)
ydimen=n_elements(ydims)

; cram inputs into correct format
point= reform(double(point),dimen,ndata)
weight= reform(double(weight),ndata)
amp= reform(double([[amp]]),ngauss)
mean= reform(double(mean),dimen,ngauss)
var= reform(double(var),dimen,dimen,ngauss)

; by default show quartiles in conditional case
if (keyword_set(conditional) and (not keyword_set(quantfrac))) then $
  quantfrac=[0.25,0.5,0.75]

; setup arrays for quantiles in conditional case
if keyword_set(quantfrac) then begin
  nquant= n_elements(quantfrac)
  quantile= dblarr(xdimen,ydimen,nquant,npix_x,2)
  quantneff= dblarr(xdimen,ydimen,npix_x)
endif

; if there is a second set of data, deal with it
if(keyword_set(weight2)) then begin
    ndata2= n_elements(weight2) 
    weight2= reform(double(weight2),ndata2)
    if(NOT keyword_set(point2)) then point2=point
    point2=reform(double(point2),dimen,ndata2)
endif

; set defaults
if NOT keyword_set(label) then $
  label= 'x!d'+strcompress(string(lindgen(dimen)),/remove_all)
if NOT keyword_set(nsig) then nsig= 5d
if NOT keyword_set(contlevel) then contlevel= [0.01,0.05,0.32,2.0/3.0]

; invert all matrices
invvar= reform(dblarr(dimen*dimen*ngauss),dimen,dimen,ngauss)
for j=0L,ngauss-1 do invvar[*,*,j]= invert(var[*,*,j],/double)

; compute mean and variance of the whole sample for plot ranges
amp1= total(weight)
mean1= total(weight##(dblarr(dimen)+1D)*point,2)/amp1
var1= 0d
for i=0L,ndata-1 do begin
    delta= point[*,i]-mean1
    var1= var1+weight[i]*delta#delta
endfor
var1= var1/amp1
if NOT keyword_set(range) then begin
  range= dblarr(2,dimen)
  for d1= 0,dimen-1 do range[*,d1]= mean1[d1]+[-nsig,nsig]*sqrt(var1[d1,d1])
endif

; set up 2-d images to save
twodimages=dblarr(npix_x,npix_y,xdimen,ydimen)

; save system plotting parameters for later restoration
bangP= !P
bangX= !X
bangY= !Y

; setup postscript file
!P.FONT= -1
!P.BACKGROUND= djs_icolor('white')
!P.COLOR= djs_icolor('black')
xsize= 7.5 & ysize= 7.5
if keyword_set(psfilename) then begin
  set_plot, "PS"
  device, file=psfilename,/inches,xsize=xsize,ysize=ysize, $
    xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
endif
!P.THICK= pthick
!P.CHARTHICK= !P.THICK
!P.CHARSIZE= 1.0
if(NOT keyword_set(axis_char_scale)) then axis_char_scale= 1.75
tiny= 1.d-4
!P.PSYM= 0
!P.LINESTYLE= 0
!P.TITLE= ''
!X.STYLE= 1
!X.THICK= 0.5*pthick
!X.CHARSIZE= tiny
!X.MARGIN= [1,1]*0.0
!X.OMARGIN= [6,6]*axis_char_scale-!X.MARGIN
!X.RANGE= 0
!X.TICKS= 0
!Y.STYLE= 1
!Y.THICK= !X.THICK
!Y.CHARSIZE= !X.CHARSIZE
!Y.MARGIN= 0.6*!X.MARGIN
!Y.OMARGIN= 0.6*!X.OMARGIN
!Y.RANGE= 0
!Y.TICKS= !X.TICKS
!P.MULTI= [xdimen*ydimen,xdimen,ydimen]
xyouts, 0,0,default_font
;  djs_xyouts, 0,0,'!BBlanton & Hogg (NYU)',/device,color='grey'

; make useful vectors for plotting
colorname= ['red','green','blue','magenta','cyan','dark yellow', $
            'purple','light green','orange','navy','light magenta','yellow green']
ncolor= n_elements(colorname)
theta= 2.0D *double(!PI)*dindgen(31)/30.0D
x= cos(theta)
y= sin(theta)

; loop over all pairs of dimensions
for id2=ydimen-1,0L,-1 do begin
    for id1=0L,xdimen-1 do begin
        d1=xdims[id1]
        d2=ydims[id2]

        if(d1 lt 0 or d2 lt 0) then begin
            !P.MULTI[0]=!P.MULTI[0]-1L
        endif else begin 

            if(keyword_set(paneluse)) then begin
                paneluse=reform(paneluse,ndata,xdimen,ydimen)
                panelindx=where(paneluse[*,id1,id2] gt 0)
                panelweight=weight[panelindx]
                panelpoint=point[*,panelindx]
            endif else begin
                panelweight=weight
                panelpoint=point
            endelse
            pointamp1= total(panelweight)
            
            if(keyword_set(paneluse2)) then begin
                if(keyword_set(weight2) and keyword_set(point2)) then begin
                    paneluse2=reform(paneluse2,ndata,xdimen,ydimen)
                    panelindx2=where(paneluse2[*,id1,id2] gt 0)
                    panelweight2=weight2[panelindx]
                    panelpoint2=point2[*,panelindx]
                endif
            endif else begin
                if(keyword_set(weight2) and keyword_set(point2)) then begin
                    panelweight2=weight2
                    panelpoint2=point2
                endif
            endelse
            
; are we on one of the plot edges?
            xprevblank=0
            yprevblank=0
            xnextblank=0
            ynextblank=0
            if(id1 gt 0) then $
              if(xdims[id1-1] eq -1) then $
              xprevblank=1
            if(id1 lt xdimen-1) then $
              if(xdims[id1+1] eq -1) then $
              xnextblank=1
            if(id2 gt 0) then $
              if(ydims[id2-1] eq -1) then $
              ynextblank=1
            if(id2 lt xdimen-1) then $
              if(ydims[id2+1] eq -1) then $
              yprevblank=1
            leftside= 0B
            if (!P.MULTI[0] EQ 0) OR $
              (((!P.MULTI[0]-1) MOD xdimen) EQ (xdimen-1) OR $
               xprevblank eq 1) then leftside= 1B
            topside= 0B
            if (!P.MULTI[0] EQ 0) OR $
              (floor(float(!P.MULTI[0]-1)/xdimen) EQ (ydimen-1) OR $
               yprevblank eq 1) then topside= 1B
            rightside= 0B
            if (((!P.MULTI[0]-1) MOD xdimen) EQ 0 OR $
                xnextblank eq 1) then rightside= 1B
            bottomside= 0B
            if (floor(float(!P.MULTI[0]-1)/xdimen) EQ 0 OR $
                ynextblank eq 1) then bottomside= 1B

; set axis label properties
            !X.CHARSIZE= tiny
;      if bottomside then !X.CHARSIZE= 1.0 ; axis_char_scale
            !Y.CHARSIZE= tiny
;      if leftside AND (d1 NE d2) then !Y.CHARSIZE= 1.0 ; axis_char_scale
            !X.TITLE= ''
            !Y.TITLE= ''

; set plot range and make axes
            !X.RANGE= range[*,d1]
            !Y.RANGE= range[*,d2]
            nticks=6/axis_char_scale
            xinterval= hogg_interval(!X.RANGE,nticks=nticks)
            yinterval= hogg_interval(!Y.RANGE,nticks=nticks)
            if d1 EQ d2 then begin
                xfrac=(panelpoint[d1,*]-!X.RANGE[0])/(!X.RANGE[1]-!X.RANGE[0])
                indxra=where(xfrac gt 0. and xfrac le 1.,countra)
                if(countra eq 0) then begin
                    !Y.RANGE=[-0.1,1.1]
                endif else begin
                    totweight=total(panelweight[indxra])
                    totmean=total(panelweight[indxra]*panelpoint[d1,indxra])/ $
                      totweight
                    totdelta=panelpoint[d1,indxra]-totmean
                    totsig=sqrt(total(totdelta^2*panelweight[indxra]) $
                                /totweight)
                    !Y.RANGE= yrangefudge*[-0.1,1.1]*totweight/(totsig)
                endelse
                yinterval= 1000*(!Y.RANGE[1]-!Y.RANGE[0])
            endif
            djs_plot,[0],[1],/nodata, $
              xtickinterval=xinterval,ytickinterval=yinterval

; make extra axis labels where necessary
            if bottomside then begin
                axis,!X.RANGE[0],!Y.RANGE[0],xaxis=0,xtickinterval=xinterval, $
                  xtitle=label[d1],xcharsize=axis_char_scale
            endif
            if topside then begin
                axis,!X.RANGE[0],!Y.RANGE[1],xaxis=1,xtickinterval=xinterval, $
                  xtitle=label[d1],xcharsize=axis_char_scale
            endif
            if leftside AND (d1 NE d2) then begin
                axis,!X.RANGE[0],!Y.RANGE[0],yaxis=0,ytickinterval=yinterval, $
                  ytitle=label[d2],ycharsize=axis_char_scale
            endif
            if rightside AND (d1 NE d2) then begin
                axis,!X.RANGE[1],!Y.RANGE[0],yaxis=1,ytickinterval=yinterval, $
                  ytitle=label[d2],ycharsize=axis_char_scale
            endif

; reset image and set x and y pixel centers
            usenpix_x=npix_x
            usenpix_y=npix_y
            if (((d2 LT d1) OR keyword_set(nodata)) AND $
                (NOT keyword_set(nomodel))) then begin
                usenpix_x= usenpix_x*long(model_npix_factor)
                usenpix_y= usenpix_y*long(model_npix_factor)
            endif
            delta_x= (range[1,d1]-range[0,d1])/usenpix_x
            delta_y= (range[1,d2]-range[0,d2])/usenpix_y
            image= dblarr(usenpix_x,usenpix_y)
            ximg= range[0,d1]+delta_x*(dindgen(usenpix_x)+0.5)
            yimg= range[0,d2]+delta_y*(dindgen(usenpix_y)+0.5)

; increment greyscale image with data
            if NOT keyword_set(nodata) then begin
                if ((d2 GE d1) OR $
                    (keyword_set(nomodel))) then begin
                    xx= floor(double(usenpix_x)* $
                              (panelpoint[d1,*]-range[0,d1])/ $
                              (range[1,d1]-range[0,d1]))
                    yy= floor(double(usenpix_y)* $
                              (panelpoint[d2,*]-range[0,d2])/ $
                              (range[1,d2]-range[0,d2]))
;                    if d1 eq d2 then yy= xx
                    for i=0L,n_elements(panelweight)-1 do begin
                        if xx[i] GE 0 AND xx[i] LT usenpix_x AND $
                          yy[i] GE 0 AND yy[i] LT usenpix_y then $
                          image[xx[i],yy[i]]= image[xx[i],yy[i]]+panelweight[i]
                    endfor

; if we are in the lower quadrant and
; there is a second set of data, use it
                    if((d2 le d1) AND keyword_set(panelweight2)) then begin
                        xx= floor(double(usenpix_x)* $
                                  (panelpoint2[d1,*]-RANGE[0,d1])/ $
                                  (RANGE[1,d1]-RANGE[0,d1]))
                        yy= floor(double(usenpix_y)* $
                                  (panelpoint2[d2,*]-RANGE[0,d2])/ $
                                  (RANGE[1,d2]-RANGE[0,d2]))
                        image2= dblarr(usenpix_x,usenpix_y)
;                        if d1 EQ d2 then yy= xx
                        for i=0L,n_elements(panelweight2)-1 do begin
                            if xx[i] GE 0 AND xx[i] LT usenpix_x AND $
                              yy[i] GE 0 AND yy[i] LT usenpix_y then $
                              image2[xx[i],yy[i]]= image2[xx[i],yy[i]]+panelweight2[i]
                        endfor
                    endif 
                endif
                image= image/abs(delta_x*delta_y)
                
                if(keyword_set(image2)) then begin
                    image2= image2/abs(delta_x*delta_y)
                    if(d2 lt d1) then image=image2
                endif
            endif

; begin loop over gaussians
            for j=0L,ngauss-1 do begin
; get eigenvalues and eivenvectors of this 2x2
                var2d= [[var[d1,d1,j],var[d1,d2,j]],[var[d2,d1,j],var[d2,d2,j]]]
                tr= trace(var2d)
                det= determ(var2d,/double)
                eval1= tr/2.0+sqrt(tr^2/4.0-det)
                eval2= tr/2.0-sqrt(tr^2/4.0-det)
                evec1= [var2d[1,0],eval1-var2d[0,0]]
                evec1= evec1/(sqrt(transpose(evec1)#evec1))[0]
                evec2= [evec1[1],-evec1[0]]
                evec1= evec1*2.0*sqrt(eval1)
                evec2= evec2*2.0*sqrt(eval2)
; increment greyscale image with gaussians
                if(NOT keyword_set(nogreyscale)) then begin
                    if (((d2 LT d1) OR keyword_set(nodata)) AND $ 
                        (NOT keyword_set(nomodel))) then begin
                        invvar2d= invert(var2d,/double)
                        for xxi=0L,usenpix_x-1 do for yyi=0L,usenpix_y-1 do $
                          image[xxi,yyi]= image[xxi,yyi]+ $
                          amp[j]/sqrt(det)/2.0/!PI* $
                          exp(-0.5*([mean[d1,j],mean[d2,j]]- $
                                    [ximg[xxi],yimg[yyi]])# $
                              invvar2d#([mean[d1,j],mean[d2,j]]- $
                                        [ximg[xxi],yimg[yyi]]))
                    endif
                endif
            endfor

; save greyscale image
            if(NOT keyword_set(nogreyscale)) then begin
                if(usenpix_x eq npix_x AND usenpix_y eq npix_y) then $
                  twodimages[*,*,id1,id2]=image
            endif

; plot greyscale image
            if (d1 NE d2) then begin
                if (NOT keyword_set(nodata) OR NOT keyword_set(nomodel)) then begin
                    loadct,0,/silent
                    if(NOT keyword_set(nogreyscale)) then begin
                    if(keyword_set(conditional)) then begin
                        xx= floor(double(usenpix_x)* $
                                  (panelpoint[d1,*]-range[0,d1])/ $
                                  (range[1,d1]-range[0,d1]))
                        yy= floor(double(usenpix_y)* $
                                  (panelpoint[d2,*]-range[0,d2])/ $
                                  (range[1,d2]-range[0,d2]))
                        amp1col=dblarr(usenpix_x)
                        avgcol=total(image,1)
                        avgcol=avgcol/total(avgcol)
                        if(keyword_set(quantfrac)) then begin
                            quantuse=lonarr(nquant,usenpix_x) 
                        endif
                        for ii=0L, usenpix_x-1L do begin
                            indx=where(xx eq ii, countin)
                            if(countin gt 0 and total(image[ii,*]) gt 0) $
                              then begin
                                amp1col[ii]=total(image[ii,*]^2)/ $
                                  total(image[ii,*]*avgcol)
                                image[ii,*]=image[ii,*]/amp1col[ii]
                                if(keyword_set(quantfrac)) then begin
                                    quantile[id1,id2,*,ii,1]= $
                                      weighted_quantile(panelpoint[d2,indx], $
                                                        panelweight[indx], $
                                                        quant=quantfrac)
                                    quantuse[*,ii]=1L
                                    quantneff[id1,id2,ii]= (total(panelweight[indx],/double))^2/total((panelweight[indx])^2,/double)
                                endif
                            endif
                        endfor
                    endif

; mess with "outimage" for plotting
                    outimage=image
                    if(keyword_set(sqrt)) then outimage=sqrt(outimage)
                    if(keyword_set(log)) then begin
                        notzero=where(outimage ne 0.d,nzcount)
                        iszero=where(outimage le 0.d,zcount)
                        if(nzcount gt 0) then begin
                            outimage[notzero]=alog10(outimage[notzero])
                            scale=0.00*(max(outimage[notzero])-min(outimage[notzero]))
                            if(zcount gt 0) then begin
                                outimage[iszero]=min(outimage[notzero])-scale
                            endif
                        endif
                    endif

; rescale outimage
                        if(keyword_set(sigrejimage)) then begin
                            tvrange= [0.0,djsig(outimage,sigrej=sigrejimage*2)]
                        endif else begin
                            tvrange= [0.0,2.0*max(outimage)]
                        endelse
                        tvimage= 255-((floor(256*(outimage-tvrange[0])/ $
                                             (tvrange[1]-tvrange[0])) > 0) $
                                      < 255)
                    endif

; plot greyscale image
                    if(NOT keyword_set(nogreyscale)) then $
                      tv, tvimage,!X.CRANGE[0],!Y.CRANGE[0],/data, $
                      xsize=(!X.CRANGE[1]-!X.CRANGE[0]), $
                      ysize=(!Y.CRANGE[1]-!Y.CRANGE[0]) 

; re-make axes (yes, this is a HACK)
                    !P.MULTI[0]= !P.MULTI[0]+1
                    djs_plot,[0],[1],/nodata, $
                      xtickinterval=xinterval,ytickinterval=yinterval

; cumulate the image
                    cumindex= reverse(sort(image))
                    cumimage= dblarr(usenpix_x,usenpix_y)
                    cumimage[cumindex]= total(image[cumindex],/cumulative) 

; renormalize the cumulated image so it really represents fractions of the
; total
                    cumimage= 1.0-cumimage*abs(delta_x*delta_y)/pointamp1

; contour
                    if(NOT keyword_set(conditional)) then begin
;                        contour, cumimage/max(cumimage),ximg,yimg,levels=contlevel, $
;                          thick=3,/overplot,color=djs_icolor('white')
                        contour, cumimage/max(cumimage),ximg,yimg,levels=contlevel, $
                          thick=1,/overplot
                    endif else begin
                        ylineupper=dblarr(usenpix_x)
                        ylinelower=dblarr(usenpix_x)
                        xline=RANGE[0,d1]+ $
                          ((dindgen(usenpix_x)+0.5)*(RANGE[1,d1]-RANGE[0,d1])/ $
                           double(usenpix_x))
                        if(keyword_set(quantfrac)) then begin
                            for m=0L, nquant-1L do begin
                                quantile[id1,id2,m,*,0]= xline
                                outmedindx=where(quantuse[m,*] ne 0,outcount)
                                if(outcount gt 0) then begin
;                                    djs_oplot,quantile[id1,id2,m,outmedindx,0], $
;                                      quantile[id1,id2,m,outmedindx,1],thick=4,psym=10, $
;                                      color=djs_icolor('white')
                                    djs_oplot,quantile[id1,id2,m,outmedindx,0], $
                                      quantile[id1,id2,m,outmedindx,1],thick=1,psym=10
                                endif
                            endfor
                        endif else begin
                            for c=0, n_elements(contlevel)-1L do begin
                                for ii=0L, usenpix_x-1L do begin
;column=smooth(image[ii,*],3,/edge_truncate)
                                    column=image[ii,*]
                                    if(total(column) gt 0.) then begin
                                        cumindex=reverse(sort(column))
                                        cumimage=dblarr(usenpix_y)
                                        cumimage[cumindex]=total(column[cumindex], $
                                                                 /cumulative) 
                                        cumimage=1.-cumimage*abs(delta_y)
                                        maximage=max(cumimage,imaximage)
                                        ilower=imaximage-lindgen(imaximage+1L)
                                        nlower=n_elements(ilower)
                                        contindx=where(cumimage[ilower] gt  $
                                                       contlevel[c] and $
                                                       cumimage[ilower-1L] lt  $
                                                       contlevel[c], $
                                                       ncont)
                                        if(ncont gt 0) then begin
                                            sp=(contlevel[c]- $
                                                cumimage[ilower[contindx[0]]-1L])/ $
                                              (cumimage[ilower[contindx[0]]]- $
                                               cumimage[ilower[contindx[0]]-1L])
                                            poslower=ilower[contindx[0]]-1L+sp
                                            ylinelower[ii]=RANGE[0,d2]+ $
                                              poslower*(RANGE[1,d2]-RANGE[0,d2])/ $
                                              double(usenpix_y)
                                        endif else begin
                                            if(ii gt 0) then ylinelower[ii]= $
                                              RANGE[0,d2]
                                        endelse
                                        iupper=imaximage+lindgen(usenpix_y- $
                                                                 imaximage)
                                        nupper=n_elements(iupper)
                                        contindx=where(cumimage[iupper] gt  $
                                                       contlevel[c] and $
                                                       cumimage[iupper+1L] lt  $
                                                       contlevel[c], $
                                                       ncont)
                                        if(ncont gt 0) then begin
                                            sp=(contlevel[c]- $
                                                cumimage[iupper[contindx[0]]])/ $
                                              (cumimage[iupper[contindx[0]]+1L]- $
                                               cumimage[iupper[contindx[0]]])
                                            posupper=iupper[contindx[0]]+sp
                                            ylineupper[ii]=RANGE[0,d2]+ $
                                              posupper*(RANGE[1,d2]- $
                                                        RANGE[0,d2])/ $
                                              double(usenpix_y)
                                        endif else begin
                                            if(ii gt 0) then ylineupper[ii]=RANGE[1,d2]
                                        endelse
                                    endif
                                    
                                endfor
                                djs_oplot, xline, ylinelower,thick=4,psym=10, $
                                    color=djs_icolor('white')
                                djs_oplot, xline, ylineupper,thick=4,psym=10, $
                                    color=djs_icolor('white')
                                djs_oplot, xline, ylinelower,thick=1,psym=10
                                djs_oplot, xline, ylineupper,thick=1,psym=10
                            endfor
                        endelse
                    endelse
                endif

                
; put on extra text labels, if asked
                if keyword_set(textlabel) AND keyword_set(textpos) AND d1 NE d2 then begin
                    ilabel= where((textpos[d1,*] GT min(!X.CRANGE)) AND $
                                  (textpos[d1,*] LT max(!X.CRANGE)) AND $
                                  (textpos[d2,*] GT min(!Y.CRANGE)) AND $
                                  (textpos[d2,*] LT max(!Y.CRANGE)),nlabel)
                    if nlabel GT 0 then begin
                        xyouts, textpos[d1,ilabel],textpos[d2,ilabel],$
                          '!D'+textlabel[ilabel],noclip=0,align=0.5,charsize=0.3
                    endif
                endif

                if keyword_set(overpoints) then begin
                    if(d1 ne d2) then begin
                        if(n_elements(psym_overpoints) eq 0) then $
                          psym_overpoints=6
                        djs_oplot,overpoints[d1,*],overpoints[d2,*], $
                          psym=psym_overpoints,color='red', thick=4
                    endif 
                endif

; plot ellipses
                if NOT keyword_set(noellipse) then begin
                    var2d= var[[d1,d2],[d1,d2],*]
                    if keyword_set(bw) then begin
                        hogg_oplot_covar, [mean[d1,*]],[mean[d2,*]],var2d,nsigma=2, $
                          color='white',thick=3.0*!P.THICK
                        hogg_oplot_covar, [mean[d1,*]],[mean[d2,*]],var2d,nsigma=2, $
                          color='black'
                    endif else begin
                        hogg_oplot_covar, [mean[d1,*]],[mean[d2,*]],var2d,nsigma=2, $
                          color=colorname[lindgen(ngauss) MOD ncolor]
                    endelse
                endif

; plot superimposed box, if asked
                if keyword_set(box) then begin
                    box_x= [box[0,d1],box[1,d1],box[1,d1],box[0,d1],box[0,d1]]
                    box_y= [box[0,d2],box[0,d2],box[1,d2],box[1,d2],box[0,d2]]
                    djs_oplot,box_x,box_y,color='white',thick=3.0*!P.THICK
                    djs_oplot,box_x,box_y,color='black'
                endif
            endif

; if we are on the diagonal...
            if (d1 EQ d2) then begin

; plot data and model histograms
                djs_oplot,RANGE[*,d1],[0,0],psym=0,xstyle=5,ystyle=5,color='grey'
                if NOT keyword_set(nodata) then begin
                    scale=abs(delta_y)
                    if(keyword_set(image2)) then begin
                        yhist2= total(image2,2)*scale
                        djs_oplot,ximg,yhist2,psym=10,thick=2*!P.THICK, $
                          color='grey'
                    endif
                    if(keyword_set(image)) then begin
                        yhist= total(image,2)*scale
                        djs_oplot,ximg,yhist,psym=10,thick=2*!P.THICK
                    endif
                endif

                if(keyword_set(overpoints)) then $
                  for io=0L, n_elements(overpoints[d1,*])-1L do $
                  djs_oplot,[overpoints[d1,io],overpoints[d1,io]], $
                  [!Y.RANGE[0],!Y.RANGE[1]], thick=6, color='red'
            endif               ; end if d1 EQ d2  
            
            if(keyword_set(panellabels)) then begin
                panellabels=reform(panellabels,xdimen,ydimen)
                if(NOT keyword_set(panellabelpos)) then begin
                    xlabelpos=!X.CRANGE[0]+0.1*(!X.CRANGE[1]-!X.CRANGE[0])
                    ylabelpos=!Y.CRANGE[1]-0.1*(!Y.CRANGE[1]-!Y.CRANGE[0])
                endif else begin
                    xlabelpos=!X.CRANGE[0]+panellabelpos[0]* $
                      (!X.CRANGE[1]-!X.CRANGE[0])
                    ylabelpos=!Y.CRANGE[0]+panellabelpos[1]* $
                      (!Y.CRANGE[1]-!Y.CRANGE[0])
                endelse
                xyouts, xlabelpos, ylabelpos, panellabels[id1,id2], $
                  charsize=sizepanellabel, charthick=4, color=255
                xyouts, xlabelpos, ylabelpos, panellabels[id1,id2], $
                  charsize=sizepanellabel, charthick=1
            endif
                                ; end loops and close file
        endelse 
    endfor
endfor
if keyword_set(psfilename) then device, /close

                                ; restore system plotting parameters
!P= bangP
!X= bangX
!Y= bangY
    
return
end
