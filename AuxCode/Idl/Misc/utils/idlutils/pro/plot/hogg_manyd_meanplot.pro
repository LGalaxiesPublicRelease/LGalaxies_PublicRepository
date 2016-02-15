;+
; NAME:
;   hogg_manyd_meanplot
; PURPOSE:
;   plot N+1-dimensional data sets, in the mean!
; INPUTS:
;   weight       [N] array of data-point weights
;   point        [d,N] array of data points - N vectors of dimension d
;   zdim         index of quantity to be averaged in the plots
;   psfilename   name for PostScript file; if no filename is given, then the
;                  plots will simply be sent to the currently active device
; OPTIONAL INPUTS:
;   nsig         number of sigma for half-width of each plot; default 5
;   dbin         [d] array of bin widths for averaging
;   label        [d] array of axis labels; default 'x_i'
;   range        [2,d] array of plotting ranges
;   xdims,ydims  indices of data dimensions to use on each x and y axis
;   axis_char_scale size of characters on labels
;   default_font font command to send to set font for plotting
; KEYWORDS:
;   [etc]        [options for hogg_meanplot, see documentation]
; OUTPUTS:
; OPTIONAL OUTPUTS:
; BUGS:
;   WAY too much duplicated code with hogg_manyd_scatterplot.
; DEPENDENCIES:
; REVISION HISTORY:
;   2003-01-12  translated from hogg_manyd_scatterplot
;-
pro hogg_manyd_meanplot, weight,point,zdim,psfilename, $
                         nsig=nsig,dbin=dbin, $
                         label=label,range=range, $
                         xdims=xdims,ydims=ydims, $
                         axis_char_scale=axis_char_scale, $
                         default_font=default_font, $
                         manyd_mean=manyd_mean, $
                         input_mean=input_mean, $
                         onedrange=onedrange, $
                         _EXTRA=KeywordsForHoggMeanplot

; check dimensions
ndata= n_elements(weight)       ; N
dimen= n_elements(point)/n_elements(weight) ; d
splog, ndata,' data points,',dimen,' dimensions'

; set defaults
if NOT keyword_set(label) then $
  label= 'x!d'+strcompress(string(lindgen(dimen)),/remove_all)
if NOT keyword_set(nsig) then nsig= 5d
if NOT keyword_set(axis_char_scale) then axis_char_scale= 1.75

; which dimensions should we look at?
if(NOT keyword_set(xdims)) then xdims=lindgen(dimen)
if(NOT keyword_set(ydims)) then ydims=lindgen(dimen)
if(NOT keyword_set(default_font)) then default_font='!3'
xdimen=n_elements(xdims)
ydimen=n_elements(ydims)

; set manyd_mean
if(NOT keyword_set(input_mean)) then begin
    manyd_mean1={range:dblarr(2,2), $
                 label:strarr(2), $
                 dbin:dblarr(2), $
                 number:ptr_new(), $
                 weight:ptr_new(), $
                 mean:ptr_new(), $
                 scatter:ptr_new()}
    manyd_mean=replicate(manyd_mean1,xdimen,ydimen)
endif

; cram inputs into correct format
point= reform(double(point),dimen,ndata)
weight= reform(double(weight),ndata)
quantity= reform(double(point[zdim,*]),ndata)

; compute mean and variance of the whole sample for plot ranges
if NOT keyword_set(range) then begin
    amp1= total(weight)
    mean1= total(weight##(dblarr(dimen)+1D)*point,2)/amp1
    var1= 0d
    for i=0L,ndata-1 do begin
        delta= point[*,i]-mean1
        var1= var1+weight[i]*delta#delta
    endfor
    var1= var1/amp1
    range= dblarr(2,dimen)
    for d1= 0,dimen-1 do range[*,d1]= mean1[d1]+[-nsig,nsig]*sqrt(var1[d1,d1])
endif

; save system plotting parameters for later restoration
bangP= !P
bangX= !X
bangY= !Y

; setup postscript file
xsize= 7.5 & ysize= 7.5
if keyword_set(psfilename) then begin
    set_plot, "PS"
    device, file=psfilename,/inches,xsize=xsize,ysize=ysize, $
      xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
endif
hogg_plot_defaults, axis_char_scale=axis_char_scale
tiny= 1.d-4
!X.CHARSIZE= tiny
!Y.CHARSIZE= !X.CHARSIZE
!P.MULTI= [xdimen*ydimen,xdimen,ydimen]
xyouts, 0,0,default_font

; loop over all pairs of dimensions
for id2=ydimen-1L,0L,-1 do begin
    for id1=0L,xdimen-1 do begin
        d1=xdims[id1]
        d2=ydims[id2]
        if d1 lt 0 or d2 lt 0 then begin
            !P.MULTI[0]=!P.MULTI[0]-1L
        endif else begin 

; set axis label properties
            !X.CHARSIZE= tiny
            !Y.CHARSIZE= tiny
            !X.TITLE= ''
            !Y.TITLE= ''

; are we on one of the plot edges?
; NB: must run this check before plotting!
            xprevblank=0
            yprevblank=0
            xnextblank=0
            ynextblank=0
            if (id1 gt 0) then if (xdims[id1-1] eq -1) then xprevblank=1
            if (id1 lt xdimen-1) then if (xdims[id1+1] eq -1) then xnextblank=1
            if (id2 gt 0) then if (ydims[id2-1] eq -1) then ynextblank=1
            if (id2 lt xdimen-1) then if (ydims[id2+1] eq -1) then yprevblank=1
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

; set dxbin and dybin
            if keyword_set(dbin) then begin
                dxbin= dbin[d1]
                dybin= dbin[d2]
            endif else begin
                dxbin= 0
                dybin= 0
            endelse

; plot
            if d1 NE d2 then begin
                if(NOT keyword_set(input_mean)) then begin
                    manyd_mean[id1,id2].range[*,0]=range[*,d1]
                    manyd_mean[id1,id2].range[*,1]=range[*,d2]
                    manyd_mean[id1,id2].label[0]=label[d1]
                    manyd_mean[id1,id2].label[1]=label[d2]
                    manyd_mean[id1,id2].dbin[0]=dxbin
                    manyd_mean[id1,id2].dbin[1]=dybin
                endif else begin
                    bin_mean=*manyd_mean[id1,id2].mean
                    bin_number=*manyd_mean[id1,id2].number
                endelse
                hogg_meanplot, point[d1,*],point[d2,*],quantity, $
                  weight=weight, axis_char_scale=axis_char_scale, $
                  dxbin=dxbin,dybin=dybin, bin_mean=bin_mean, $
                  bin_number=bin_number, bin_scatter=bin_scatter, $
                  bin_weight=bin_weight, input_mean=input_mean, $
                  xrange=range[*,d1],yrange=range[*,d2], $
                  _EXTRA=KeywordsForHoggMeanplot
                if(NOT keyword_set(input_mean)) then begin
                    manyd_mean[id1,id2].mean=ptr_new(bin_mean)
                    manyd_mean[id1,id2].scatter=ptr_new(bin_scatter)
                    manyd_mean[id1,id2].number=ptr_new(bin_number)
                    manyd_mean[id1,id2].weight=ptr_new(bin_weight)
                endif
            endif else begin
; HACK: placeholder
                if(NOT keyword_set(input_mean)) then begin
                    manyd_mean[id1,id2].range[*,0]=range[*,d1]
                    manyd_mean[id1,id2].range[*,1]=range[*,d2]
                    manyd_mean[id1,id2].label[0]=label[d1]
                    manyd_mean[id1,id2].label[1]=label[d2]
                    manyd_mean[id1,id2].dbin[0]=dxbin
                    manyd_mean[id1,id2].dbin[1]=dybin
                endif else begin
                    bin_mean=*manyd_mean[id1,id2].mean
                    bin_number=*manyd_mean[id1,id2].number
                endelse
                blanton_oned_meanplot, point[d1,*],quantity, $
                  weight=weight,axis_char_scale=axis_char_scale, dxbin=dxbin, $
                  bin_mean=bin_mean,bin_number=bin_number, $
                  bin_scatter=bin_scatter, bin_weight=bin_weight, $
                  input_mean=input_mean, xrange=range[*,d1], $
                  yrange=onedrange, _EXTRA=KeywordsForHoggMeanplot
                if(NOT keyword_set(input_mean)) then begin
                    manyd_mean[id1,id2].mean=ptr_new(bin_mean)
                    manyd_mean[id1,id2].weight=ptr_new(bin_weight)
                    manyd_mean[id1,id2].number=ptr_new(bin_number)
                    manyd_mean[id1,id2].scatter=ptr_new(bin_scatter)
                endif
            endelse

; make axis labels afterwards
            if bottomside then begin
                axis,!X.CRANGE[0],!Y.CRANGE[0],xaxis=0, $
                  xtitle=label[d1],xcharsize=axis_char_scale
            endif
            if topside then begin
                axis,!X.CRANGE[0],!Y.CRANGE[1],xaxis=1, $
                  xtitle=label[d1],xcharsize=axis_char_scale
            endif
            if leftside AND (d1 NE d2) then begin
                axis,!X.CRANGE[0],!Y.CRANGE[0],yaxis=0, $
                  ytitle=label[d2],ycharsize=axis_char_scale
            endif
            if rightside AND (d1 NE d2) then begin
                axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
                  ytitle=label[d2],ycharsize=axis_char_scale
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
