;+
; NAME:
;   hogg_make_image_label
; PURPOSE:
;   make text label for RGB images
; CALLING SEQUENCE:
;   label= hogg_make_image_label(title,naxis1, [subtitle1=,subtitle2=])
; EXAMPLES:
;   label= hogg_make_image_label('M51',2048, $
;     subtitle1='data from the Sloan Digital Sky Survey', $
;     subtitle2='image by Hogg, Wherry, Blanton, Finkbeiner, Schlegel')
; INPUTS:
;   title      - main title text
;   naxis1     - width in pixels to make the label
; OPTIONAL INPUTS:
;   subtitle1  - first subtitle
;   subtitle2  - second subtitle; defaults to idlutils plug
; OUTPUT:
;   label      - [naxis1,n2] image containing text
; COMMENTS:
; BUGS:
;   - Format and label padding all hard-wired.
;   - Makes insecure temporary LaTeX file etc.
;   - Relies on horrifying UNIX bitmapping code.
;   - Slightly odd behavior if subtitles 1 and 2 are both empty.
; DEPENDENCIES:
;   LaTeX etc
;   pstopnm etc
; REVISION HISTORY:
;   2003-12-05  written - Hogg
;-
function hogg_make_image_label, title,naxis1, $
                                subtitle1=subtitle1,subtitle2=subtitle2
if not keyword_set(subtitle1) then subtitle1='~'
if not keyword_set(subtitle2) then $
  subtitle2='made with \texttt{idlutils} (http:$/\!/$skymaps.info)'
prefix= 'tmp_hogg_make_image_label'
openw, wlun,prefix+'.tex',/get_lun
printf, wlun,'\documentclass[10pt]{article}'
printf, wlun,'\usepackage{color}'
printf, wlun,'\setlength{\headheight}{0in}'
printf, wlun,'\setlength{\headsep}{0in}'
printf, wlun,'\setlength{\topmargin}{0in}'
printf, wlun,'\setlength{\oddsidemargin}{-1.0truein}'
printf, wlun,'\setlength{\textheight}{1.0truein}'
printf, wlun,'\setlength{\textwidth}{8.5truein}'
printf, wlun,'\setlength{\parskip}{0ex}'
printf, wlun,'\setlength{\parsep}{0ex}'
printf, wlun,'\pagestyle{empty}'
printf, wlun,'\pagecolor{black}\color{white}'
printf, wlun,'\begin{document}'
printf, wlun,'\noindent'
printf, wlun,'\textsf{'
printf, wlun,'\textbf{\Large '+title+' }'
printf, wlun,'\hfill'
printf, wlun,'\footnotesize'
printf, wlun,'\begin{tabular}{r}'
printf, wlun,subtitle1+' \\'
printf, wlun,subtitle2
printf, wlun,'\end{tabular}'
printf, wlun,'}'
printf, wlun,'\end{document}'
close, wlun
free_lun, wlun
cmd= '\rm -fv '+prefix+'.ppm'
splog, cmd
spawn, cmd
cmd= '\latex '+prefix
splog, cmd
spawn, cmd
spawn, cmd
cmd= '\dvips -E -Ppdf '+prefix
splog, cmd
spawn, cmd
factor= 2
cmd= '\pstopnm -ppm -stdout -xborder 0 -yborder 0'+ $
  ' -ysize '+strtrim(string(factor*naxis1),2)+ $
  ' '+prefix+' > '+prefix+'.ppm'
splog, cmd
spawn, cmd
label= read_image(prefix+'.ppm')
tmp= size(label,/dimensions)
nx= tmp[1]
ny= tmp[2]
label= rotate(reform(label[0,*,*],nx,ny),6)
label= float(label)/255.0
if (nx NE (nx/2)*2) then begin
    label= [[label],[fltarr(ny)]]
    nx= nx+1
endif
label= rebin(label,ny/factor,nx/factor)
if nx GT (0.1*factor*naxis1) then begin
    splog, 'warning: re-sized label will not be full-width'
    newfactor= (0.03*float(naxis1))/float(nx)
    label= congrid(label,ny*newfactor,nx*newfactor,/interp)
endif
help, label
return, label
end
