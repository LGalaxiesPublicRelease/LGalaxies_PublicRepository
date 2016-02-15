;;
; 
; Copyright (C) 2006 Patricio Rojo
; 
; This program is free software; you can redistribute it and/or
; modify it under the terms of the GNU General Public License
; as published by the Free Software Foundation; either version 2
; of the License, or (at your option) any later version.
; 
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
; 
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
;   MA  02110-1301, USA.
; 
; 
;;

;+
; NAME:
;	GETFROMWIN
;
; PURPOSE:
;	This function reads coordinate array from window.
;
; CATEGORY:
;	Pato's input routine.
;
; CALLING SEQUENCE:
;
;	data = GETFROMWIN(Array[, Xaxis], /YPAIRS, /XPAIRS,
;                   	  /ABSORPSEARCH, /EMISSIONSEARCH, 
;                         /COORDINATE, MODENAME=MODENAME, TITLE=TITLE,
;                         SORTFLAG=SORTFLAG, COLORMARK=COLORMARK,
;                         KEYALTERNATIVES={delete:xx, add:xx,
;                         reset:xx, exit:xx, default:xx, zoom:xx,
;                         zoomx:xx, zoomy:xx, refit:xx, fitdeg:xx,
;                         center:xx, tag:xx, print:xx, stop:xx,
;                         undo:xx}, 
;                         DEFAULT=DEFAULT, /DODEFAULT, /CENTERDEF,
;                         /NONINTERACTIVE, /NOADD, /NOTAG,
;                         /FITDEFAULT, /NOFIT, FITCOLOR=FITCOLOR,
;                         FITDEG=FITDEG)
;
; INPUTS:
;	Array:  Array to show in the window
;
; OPTIONAL INPUTS:
;	Xcoo:   X coordinate array. It will use pixel value if omitted.
;
; KEYWORD PARAMETERS:
;       KEYALTERNATIVES: Structure with alternative keys to accept for
;                 "delete", "add", "reset", "exit", "default", "zoom",
;                 "zoomx", "zoomy", "refit", "fitdeg", "center",
;                 "tag", "print", "break", and "undo" events. Defaults are 'd',
;                 'a', 'r', 'x', 'e', 'z', 'x', 'y', 'f', 'g', 'c',
;                 't', 'p', 'b' and 'u', respectively.
;       MODENAME: This is the method name that will be displayed in
;                 the lower left corner of the plotting window.
;       TITLE:    Plot title
;       SORTFLAG: It is set to 1 if you want an x-axis sort, or to 2
;                 if y-axis
;       COLORMARK: Color of the user marks. Default is green
;                 (255L*256)
;       NONINTERACTIVE: Set this keyword to just produce plot and
;                 optionally process the defaults. No user
;                 interaction.
;       NOADD:    Set this flag to disable adding of new points.
;       NOTAG:    Set this flag to disable add & taging of new
;                 points.
;	DEFAULT:  Default values if DODEFAULT is set or "default"
;                 event key is pressed, has to have format according
;                 to the selected return option.
;       DODEFAULT: Flag that indicates to start with default values.
;       CENTERDEF: Set this flag to center default points given any of
;                 the two peak search methods. Setting this activates
;                 DODEFAULT as well.
;       KEEPZOOM: Do not reset view once you quit interactive.
;       ALLOWBREAK: Set this flag to allow interruption of interactive
;                 window, whenever key for break is pressed. Used for
;                 debugging.
;
;   Fitting keywords for xpair mode:
;       NOFIT:    Forbids polynomial fitting.
;       FITDEG:   (I/O) Degree of polynomial fit. Default is 2
;       FITCOLOR: Color to use for the fitted line. Default is 255
;                 (Red).
;       FITDEFAULT: Set this flag to fit default range if
;                 given. Setting this activates DODEFAULT as well.
;
;   Return format options (only one of these flags can be specified
;   per run, if more than then the last is used):
;	YPAIRS:   Set this flag to return paired y coordinate data
;                 values. 
;       XPAIRS:	  Set this flag to return paired x coordinate data
;                 values. 
;       COORDINATES: Set this flag to return each x,y coordinate
;                 pressed
;       EMISSIONSEARCH: Set this flag to find the center of emission
;                 lines in the array.
;       ABSORPSEARCH: Set this flag to find the center of absorption
;                 lines in the array.
;       DATASELECT: Set this flag to look for X, Y values that match
;                 the input data. A proximity search is performed.
;
; OUTPUTS:
;	This function returns cursor position over window data
;	coordinate according to return option
;
; OPTIONAL OUTPUTS:
;       FITCOEF:  Return coefficients of polynomial fit.
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; TODO:
;
; PROCEDURE:
;	Algorithm description.
;
; EXAMPLE:
;	Take values given some default points
;
;		F = GETFROMWIN(array, def=[[1,2],[4,5]])
;
; MODIFICATION HISTORY:
; 	Written by:	Pato Rojo, Cornell.  2005 Mar
;			pato@astro.cornell.edu
;-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function polyn_getfw, c, x
compile_opt idl2, hidden

;build polynomial
ret = 0
for i=n_elements(c)-1, 0, -1 do ret = x*ret + c[i]

return, ret

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function center1peak_getfw, pos, absorp, x, y, fitmethod
compile_opt idl2, hidden

posi = polyn_getfw(poly_fit(x, dindgen(n_elements(x)), 2), pos[0])
posi = [long(posi), long(posi)+1]

;make it always emission for the following analysis
if absorp then dat = -y $
else dat = y

;is it at the right of the peak?
right = dat[posi[0]] gt dat[posi[1]]
cr = right

;check end of peaks
beg = posi[0]
while beg gt 0 do begin
    if dat[beg-1] gt dat[beg] then begin
        if ~ cr then break
    endif else if cr then cr = 0
    beg--
endwhile
fin = posi[1]
while fin lt n_elements(x)-1 do begin
    if dat[fin+1] gt dat[fin] then begin
        if cr then break
    endif else if ~ cr then cr = 1
    fin++
endwhile

switch fitmethod of
    0:begin
;spline centering
        ii = (x[fin]-x[beg])*dindgen((fin-beg)*100+1) / (100.0*(fin-beg)) + $
          x[beg]
        c = spl_init(x[beg:fin], dat[beg:fin], /double)
        dd = spl_interp(x[beg:fin], dat[beg:fin], c, ii, /double)
        fl = max(dd, cnt)
        cnt = ii[cnt]
        break
    end
    else: message, "Fit method '" + string(fitmethod) + "' is not enabled"
endswitch

if absorp then fl = -fl


return, [cnt, fl]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro centerpeak_getfw, coo, absorp, x, y, fitmethod
compile_opt idl2, hidden

;get dimensions and prepare
if absorp ne 1 and absorp ne 0 then return
n = (size(coo, /dim))[0]
if n eq 0 then return

;transform x value into indices values.
for i=0, n-1 do $
  coo[i, *] = center1peak_getfw(reform(coo[i,*]), absorp, x, y, fitmethod)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro addmarks, xaxis, array, coo, mode, ncolor, bordrange=bordrange
compile_opt idl2, hidden

if ~ keyword_set(bordrange) then bordrange = 20.0

xpos = !x.crange[0] + (!x.crange[1] - !x.crange[0]) / bordrange
ypos = !y.crange[0] + (!y.crange[1] - !y.crange[0]) / bordrange
n = (size(coo, /dim))[0]

if n lt 1 then return

oldcol = !p.color
!p.color = ncolor

switch mode of
    5:
    0: begin                    ;data mode
        oplot, coo[*,0], coo[*,1], psym=7
        break
    end
    1: begin                    ;xpair
        np = (size(coo, /dim))[0]
        if np and 1 then np--
        if np gt 0 then begin
            idx = makeidxrange(reform(coo[0:np-1, 0], 2, $
                                      np/2), xaxis)
            oplot, xaxis[idx], array[idx], color=ncolor, psym=7
        endif
        usersym, [0, 0], [2, -2]
        for i=0, n-2, 2 do $
          oplot, [coo[i, 0], coo[i+1, 0]], [ypos, ypos]
        oplot, coo[*, 0], replicate(ypos, n), psym=8
        break
    end
    2: begin                    ;ypair
        usersym, [2, -2], [0, 0]
        for i=0, n-2, 2 do $
          oplot, [xpos, xpos], [coo[i, 1], coo[i+1, 1]]
        oplot, repliocate(xpos, n), coo[*, 1], psym=8
        break
    end
    3:
    4: begin                    ;peaksearch
        if mode eq 3 then usersym, [0, 0], [1, 5] $
        else usersym, [0, 0], [-1, -5]
        oplot, coo[*, 0], coo[*, 1], psym=8
        break
    end
    else: message, $
      "getfromwin() mode '" + string(mode) + "' is not implemented"
endswitch

!p.color = oldcol

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function defaultset, def, tags, mode, xrange, yrange
compile_opt idl2, hidden

;Check formats
if ~ keyword_set(def) then return, 0
if mode gt 5 then $
  message, "Mode '" + string(mode) + $
  "' is not an implemented data collection mode in getfromwin()"
if size(def, /n_dim) ne 2 then message, $
  "Format of DEFAULT's keyword is not correct in getfromwin()"
datapairs = mode eq 1 or mode eq 2
if datapairs then dim2 = 0 $
else dim2 = 1
if dim2 ge 0 and (size(def, /dim))[dim2] ne 2 then $
  message, "Format of DEFAULT's keyword is not correct in getfromwin()"

;look for insiders
if datapairs then begin
    range = [[xrange], [yrange]]
    cond = def[0, *] ge range[0, mode-1] and $
      def[0, *] le range[1, mode-1] and $
      def[1, *] ge range[0, mode-1] and $
      def[1, *] le range[1, mode-1]
endif else begin
    xcond = (def[*, 0] ge xrange[0] and def[*, 0] le xrange[1])
    ycond = 1
    if mode eq 0 or mode eq 5 then $
      ycond = (def[*, 1] ge yrange[0] and def[*, 1] le yrange[1])
    cond = xcond and ycond
endelse

;choose insiders and returns if there are none
idx = where(cond, n)
if idx[0] eq -1 then begin
    tags = 0
    return, 0
endif

;prepare return
if datapairs then ret = [[reform(def[*, idx], n*2)], [replicate(0, n*2)]] $
else ret = def[idx, *]

;prepare tags
if datapairs then nt = n_elements(def) $
else nt = (size(def, /dim))[0]
case n_elements(tags) of
    0: tags = datapairs? strarr(n*2): strarr(n)
    nt: tags = datapairs? reform(tags[*, idx], n*2): tags[idx]
    else: message, "Format of default TAGS is not valid"
endcase

return, ret

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function closestpoint_getfw, x, y, xaxis, array
compile_opt idl2, hidden

cnt = n_elements(xaxis)
idx = (sort(total(([[xaxis], [array]] - [x, y] ##  $
                   replicate(1, cnt))^2, 2) ))[0]
return, [xaxis[idx[0]], array[idx[0]]]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function in_struct, struct, tag
compile_opt idl2, hidden

if ~ keyword_set(struct) or size(struct, /type) ne 8 then $
  return, 0
return, (where(strpos(tag_names(struct), tag) gt -1))[0]  ne -1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro doplots, xaxis, array, title, psym, coo, printlabel, fineaxis, $
             fitcoef, fitcolor, colormarks, retmode, modename, $
             plotyrng, plotxrng

erase
if keyword_set(plotyrng) and keyword_set(plotxrng) then $
  plot, xaxis, array, title=title, xrange=plotxrng, yrange=plotyrng, $
  psym=psym, /xstyle, /ystyle $
else plot, xaxis, array, title=title, psym=psym

if keyword_set(fitcoef) then $
  oplot, fineaxis, polyn_getfw(fitcoef, fineaxis), color=fitcolor

addmarks, xaxis, array, coo, retmode, colormarks, bordrange=2
if printlabel then xyouts, 0.01, 0.01, "Interactive Window. " + modename, $
  color=255L+256*255L, /normal, charsize=1.5

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function getfromwin, usrarray, usraxis, $
                     default=default, dodefault=dodefault, $
                     keyalternatives=keyalternatives, $
                     xpairs=xpairs, ypairs=ypairs, $
                     coordinates=coordinates, emmisionsearch=emissionsearch, $
                     absorpsearch=absorpsearch, dataselect=dataselect, $
                     modename=modename, twowin=twowin, $
                     noninteractive=noninteractive, noadd=noadd, $
                     title=title, posttitle=posttitle, sortflag=sortflag, $
                     colormark=colormark, verbose=verbose, $
                     fitdeg=fitdeg, fitcoef=fitcoef, nofit=nofit, $
                     fit=fit, fitcolor=fitcolor, fitdefault=fitdefault, $
                     centerdef=centerdef, keepzoom=keepzoom, $
                     tags=tags, notag=notag, tagname=tagname, $
                     convolvewidth=convolvewidth, $
                     allowbreak=allowbreak
compile_opt idl2
noerasestatus = !p.noerase
!p.noerase = 1

;only one of the following should have been specified
retmode = 1 
if keyword_set(coordinates)    then retmode = 0
if keyword_set(xpairs)         then retmode = 1
if keyword_set(ypairs)         then retmode = 2
if keyword_set(emissionsearch) then retmode = 3
if keyword_set(absorpsearch)   then retmode = 4
if keyword_set(dataselect)     then retmode = 5

if ~ keyword_set(modename) then begin
    case retmode of
        0: modename = "Coordinate return"
        1: modename = "Pair of X elements return"
        2: modename = "Pair of Y elements return"
        3: modename = "Emission peaks return"
        4: modename = "Absorption peaks return"
        5: modename = "Selected data points"
        else: message, "Mode " + string(retmode) + $
          ", doesn't have a default name. Aborting"
    endcase
endif

;fit keywords
if ~ keyword_set(fitcolor) then fitcolor = 255 ;red
if ~ keyword_set(fitdeg) then fitdeg = 2
if ~ keyword_set(fit) then fit = 0
if ~ keyword_set(nofit) then nofit = 0
if keyword_set(fitdefault) then dodefault = 1
if retmode ne 1 and retmode ne 5 then nofit = 1  ;only fit on horizontal
                                                 ;ranges or selected points

;simple keywords
if retmode eq 5 then psym = 5 else psym = 0
if ~ keyword_set(tagname) then tagname = "tag"
if ~ keyword_set(notag) then notag = 0
if ~ keyword_set(noadd) then noadd = 0
if keyword_set(centerdef) then dodefault = 1
if ~ keyword_set(title) then title = ""
if ~ keyword_set(posttitle) then posttitle = title
if ~ keyword_set(colormarks) then colormarks = 256L*255 ;green
if keyword_set(verbose) then verbosity = 1 else verbosity = 0

;Axis, arrays and double plot
if ~ keyword_set(twowin) then twowin = 0
if twowin then begin
    npl = 2
    if size(usrarray, /n_dim) gt 1 and $
      (size(usrarray, /dim))[1] gt 1 then array = usrarray[*, 0:npl-1] $
    else array = usrarray[*,0] # [1, 1]
    if n_params() lt 2 then $
      xaxis = dindgen((size(array, /dim))[0]) # [1, 1] $
    else if size(usraxis, /n_dim) ne 2 then $
      xaxis = usraxis[*,0] # [1, 1] $
    else xaxis = usraxis[*, 0:npl-1]
endif else begin
    npl = 1
    array = usrarray[*,0]
    if n_params() lt 2 then $
      xaxis = dindgen((size(usrarray, /dim))[0], npl) $
    else xaxis = usraxis[*,0]
endelse
for i=0, npl-1 do begin
    idx = sort(xaxis[*, i])
    xaxis[*, i] = xaxis[idx, i]
    array[*, i] = array[idx, i]
endfor

;sortflag could be 1 for x sorting, or 2 for y-sorting, otherwise
;there is no sorting
if ~ keyword_set(sortflag) then sortflag = 0

;set alternate keyboard mapping
if ~ in_struct(keyalternatives, "DEFAULT") then defaultkeys = "e" $
else defaultkeys = keyalternatives.default
if ~ in_struct(keyalternatives, "DELETE") then deletekeys = "d" $
else deletekeys = keyalternatives.delete
if ~ in_struct(keyalternatives, "FITDEG") then fitdegkeys = "g" $
else fitdegkeys = keyalternatives.fitdeg
if ~ in_struct(keyalternatives, "CENTER") then centerkeys = "c" $
else centerkeys = keyalternatives.center
if ~ in_struct(keyalternatives, "REFIT") then refitkeys = "f" $
else refitkeys = keyalternatives.refit
if ~ in_struct(keyalternatives, "PRINT") then printkeys = "p" $
else printkeys = keyalternatives.print
if ~ in_struct(keyalternatives, "RESET") then resetkeys = "r" $
else resetkeys = keyalternatives.reset
if ~ in_struct(keyalternatives, "ZOOMX") then zoomxkeys = "x" $
else zoomxkeys = keyalternatives.zoomx
if ~ in_struct(keyalternatives, "ZOOMY") then zoomykeys = "y" $
else zoomykeys = keyalternatives.zoomy
if ~ in_struct(keyalternatives, "BREAK") then breakkeys = "b" $
else breakkeys = keyalternatives.break
if ~ in_struct(keyalternatives, "EXIT") then exitkeys = "q" $
else exitkeys = keyalternatives.exit
if ~ in_struct(keyalternatives, "ZOOM") then zoomkeys = "z" $
else zoomkeys = keyalternatives.zoom
if ~ in_struct(keyalternatives, "UNDO") then undokeys = "u" $
else undokeys = keyalternatives.undo
if ~ in_struct(keyalternatives, "TAG") then tagkeys = "t" $
else tagkeys = keyalternatives.tag
if ~ in_struct(keyalternatives, "ADD") then addkeys = "a" $
else addkeys = keyalternatives.add
if ~ in_struct(keyalternatives, "CONV") then convkeys = "k" $
else convkeys = keyalternatives.conv

erase
plot, xaxis[*,0], array[*,0], title=title, psym=psym
fullxrng = !x.crange
fullyrng = !y.crange

fineaxis = (indgen(200) * (fullxrng[1] - fullxrng[0]) / 200.0 + $
  fullxrng[0]) # replicate(1.0, npl)

if keyword_set(convolvewidth) then begin
    t = convolvewidth
    dw = xaxis[1,0] - xaxis[0,0]
    nk = long(30.0 * t / dw)
    if nk gt (size(array, /dim))[0] then begin
        print, 'IW: Convolution profile is bigger than array. ' + $
          'Ignoring convolution request'
    endif else begin
        xx = (dindgen(nk) - nk/2) * dw / t
        prof = 1.0d / (sqrt(2*!pi) * t / dw) * exp(-0.5*xx^2)
        array = convol(array[*, 0], prof,1)
    endelse
endif

if keyword_set(dodefault) then begin
    coo = defaultset(default, tags, retmode, fullxrng, fullyrng)
    if retmode eq 5 then for i=0, (size(coo, /dim))[0]-1 do $
      coo[i, *] = closestpoint_getfw(coo[i,0], coo[i,1], xarray, data)

    if keyword_set(fitdefault) then begin
        np = (size(coo, /dim))[0]
        if nofit then print, $
          "IW: Polynomial fit is disabled in this window" $
        else if np ge 0 then begin
            if np gt 0 then begin
                if np and 1 then np--
                idxrng = makeidxrange(reform(coo[0:np-1, 0], 2, $
                                             np/2), xaxis)
            endif else idxrng = lindgen(n_elements(xaxis))
            fitcoef = reform(poly_fit(xaxis[idxrng], $
                                      array[idxrng], fitdeg))
        endif
    endif
    if keyword_set(centerdef) then $
      centerpeak_getfw, coo, retmode-3, xaxis, array, 0

endif

if keyword_set(noninteractive) then interactive = 0 else interactive = 1

doplots, xaxis, array, title, psym, coo, interactive,         $
  fineaxis, fitcoef, fitcolor, colormarks, retmode, modename, $
  plotyrng, plotxrng


while interactive do begin
    t = 0
    k = get_kbrd(1)

    cursor, x, y, /nowait, /data
    addition = 0

    if      strpos(defaultkeys, k) ge 0 then presskey = 1  $
    else if strpos(deletekeys, k) ge 0  then presskey = 2  $
    else if strpos(resetkeys, k) ge 0   then presskey = 3  $
    else if strpos(exitkeys, k) ge 0    then presskey = 4  $
    else if strpos(undokeys, k) ge 0    then presskey = 5  $
    else if strpos(addkeys, k) ge 0     then presskey = 6  $
    else if strpos(zoomkeys, k) ge 0    then presskey = 7  $
    else if strpos(zoomxkeys, k) ge 0   then presskey = 8  $
    else if strpos(zoomykeys, k) ge 0   then presskey = 9  $
    else if strpos(refitkeys, k) ge 0   then presskey = 10 $
    else if strpos(fitdegkeys, k) ge 0  then presskey = 11 $
    else if strpos(centerkeys, k) ge 0  then presskey = 12 $
    else if strpos(tagkeys, k) ge 0     then presskey = 13 $
    else if strpos(printkeys, k) ge 0   then presskey = 14 $
    else if strpos(breakkeys, k) ge 0   then presskey = 15 $
    else if strpos(convkeys, k) ge 0    then presskey = 16 $
    else presskey = -1

;prepare undo if appropiate
    if presskey eq 5 and keyword_set(undov) then begin
        x = undoc[0]
        y = undoc[1]
        t = undot
        k = undov
    endif
;exiting
    if presskey eq 4 then break

    switch presskey of
        1: begin                ;defaults
            if keyword_set(default) then $
              coo = defaultset(default, tags, retmode, $
                               fullxrng, fullyrng)
            if retmode eq 5 then for i=0, (size(coo, /dim))[0]-1 do $
              coo[i, *] = closestpoint_getfw(coo[i,0], coo[i,1], xarray, data)

            break
        end
        2: if keyword_set(coo) then begin ;deleting points
            cnt = (size(coo, /dim))[0]
            case retmode of
                0: rem = (sort(total((coo - [x, y] ##  $
                                     replicate(1, cnt))^2, 2) ))[0]
                1: rem = (sort((coo[*, 0] - x)^2))[0]
                2: rem = (sort((coo[*, 1] - y)^2))[0]
                4: rem = (sort((coo[*, 0] - x)^2))[0]
                5: rem = (sort((coo[*, 0] - x)^2))[0]
                else: message, "getfromwin() mode '" + string(retmode) + $
                  "' is not implemented"
            endcase
            undot = tags[rem[0]]
            undoc = reform(coo[rem[0], *])
            undov = strmid(addkeys, 0, 1)
            if cnt eq 1 then begin
                coo = 0
                tags = 0
            endif else begin
                idx = where(lindgen(cnt) ne rem[0])
                coo = coo[idx, *]
                tags = tags[idx, *]
;rotate if this deletion is breaking a range
                if (retmode eq 1 or retmode eq 2) and $
                  rem[0] lt (size(coo, /dim))[0]-1 then begin
                    if rem[0] and 1 then rem[0]--
                    coo[rem[0]:*, *] = shift(coo[rem[0]:*, *], -1, 0)
                    tags[rem[0]:*] = shift(tags[rem[0]:*], -1)
                endif
            endelse
            erase
            break
        endif
        3: begin                ;reset
            coo = 0
            zoomwin = 0
            zoomxrng = 0
            zoomyrng = 0
            if twowin then begin
                npl = 2
                if size(usrarray, /n_dim) gt 1 and $
                  (size(usrarray, /dim))[1] gt 1 then $
                    array = usrarray[*, 0:npl-1] $
                else array = usrarray[*,0] # [1, 1]
                if n_params() lt 2 then $
                  xaxis = dindgen((size(array, /dim))[0]) # [1, 1] $
                else if size(usraxis, /n_dim) ne 2 then $
                  xaxis = usraxis[*,0] # [1, 1] $
                else xaxis = usraxis[*, 0:npl-1]
            endif else begin
                npl = 1
                array = usrarray[*,0]
                if n_params() lt 2 then $
                  xaxis = dindgen((size(usrarray, /dim))[0], npl) $
                else xaxis = usraxis[*,0]
            endelse
            for i=0, npl-1 do begin
                idx = sort(xaxis[*, i])
                xaxis[*, i] = xaxis[idx, i]
                array[*, i] = array[idx, i]
            endfor
            erase
            break
        end
        6: begin                ;add
            if noadd then break
            if x lt !x.crange[0] or y lt !y.crange[0] or $
              x gt !x.crange[1] or y gt !y.crange[1] then $
              break
            t = ""
            if ~ keyword_set(coo) then begin
                coo = [[x], [y]]
                tags = [t]
            endif else begin
                coo = [[coo[*,0], x], [coo[*,1], y]]
                tags = [tags, t]
            endelse
            undot = t
            undoc = [x, y]
            undov = strmid(deletekeys, 0, 1)
            addition = 1
            break
        end
        7:                      ;zoom
        8:                      ;zoomx
        9: begin                ;zoomy
            if keyword_set(zoomwin) then begin
                if x eq zoomwin[0] or presskey eq 9 then plotxrng = fullxrng $
                else plotxrng = [zoomwin[0], x]
                if y eq zoomwin[1] or presskey eq 8 then plotyrng = fullyrng $
                else plotyrng = [zoomwin[1], y]
                zoomwin = 0
            endif else zoomwin = [x, y]
            erase
            break
        end
        11:begin                ;fitdeg
            if nofit then begin 
                print, "IW: Polynomial fit is disabled in this window"
                break
            endif else begin
                print, "IW: Enter new polynomial degree (default: " + $
                  string(fitdeg, format='(i0)') + "): ", format='(a,$)'
                read, fitdeg, format='(i)'
            endelse
        end
        10:begin                ;refit (note that it is connected with
                                ;previous case: don't put anything; in
                                ;between)
            np = (size(coo, /dim))[0]
            if nofit then print, $
              "IW: Polynomial fit is disabled in this window" $
            else begin
                if np gt 0 then begin
                    if np and 1 then np--
                    idxrng = makeidxrange(reform(coo[0:np-1, 0], 2, $
                                                 np/2), xaxis)
                endif else if np eq 0 then $
                  idxrng = lindgen(n_elements(xaxis)) $
                else break
                fitcoef = reform(poly_fit(xaxis[idxrng], $
                                          array[idxrng], fitdeg))
                erase
            endelse
            break
        end
        12: begin               ;Centering in all peaks
            print, "IW: Centering... ", format='(a,$)'
            centerpeak_getfw, coo, retmode-3, xaxis, array, 0
            print, 'done'
            erase
            break
        end
        13: begin               ;adding and tagging a point
            if notag then break
            if x lt !x.crange[0] or y lt !y.crange[0] or $
              x gt !x.crange[1] or y gt !y.crange[1] then $
              break
            if ~ keyword_set(t) then begin
                print, "IW: Enter new " + tagname + " for this point (", $
                  x, ", ", y, ")", format='(a,g9.3,a,g8.2,a,$)'
                t = ""
                read, t, format='(a)'
            endif
            if ~ keyword_set(coo) then begin
                coo = [[x], [y]]
                tags = [t]
            endif else begin
                coo = [[coo[*,0], x], [coo[*,1], y]]
                tags = [tags, t]
            endelse
            undoc = [x, y]
            undot = t
            undov = strmid(deletekeys, 0, 1)
            addition = 1
            break
        end
        14: begin               ;print point and tag information
            n = (size(coo, /dim))[0]
            for i=0, n-1 do $
              print, tags[i], "(", coo[i,0], ",", coo[i,1], ")", $
              format='(a-25,a,g11.6,a,g11.6,a)'
            break
        end
        15: begin               ;breaking
            if keyword_set(allowbreak) then begin 
                print, "IW: Interrupting getfromwin by user request"
                stop
            endif
            break
        end
        16: begin               ;convolving
            t = ""
            print, 'IW: Enter convolution width', format='(a,$)'
            read, t, format='(a)'
            
            dw = xaxis[1,0] - xaxis[0,0]
            nk = long(30.0 * t / dw)
            if nk gt (size(array, /dim))[0] then begin
                print, 'IW: Convolution profile is bigger than array. ' + $
                  'Not posible'
            endif else begin
                xx = (dindgen(nk) - nk/2) * dw / t
                prof = 1.0d / (sqrt(2*!pi) * t / dw) * exp(-0.5*xx^2)
                array = convol(array[*, 0], prof,1)
                erase
            endelse
            break
        end
        else: if verbosity then print, "IW: Key '" + k + $
          "' is not an enabled function of getfromwin()"
    endswitch

;some fixes after point addition
    i = (size(coo, /dim))[0] - 1
    if addition then begin
        case retmode of
;focus in closest point
            5: coo[i,*] = closestpoint_getfw(x, y, xaxis, array)
;center if we are peak searching
            3: coo[i,*] = center1peak_getfw(reform(coo[i,*]), retmode-3, $
                                            xaxis, array, 0)
            4: coo[i,*] = center1peak_getfw(reform(coo[i,*]), retmode-3, $
                                            xaxis, array, 0)
            else: break
        endcase
    endif

;plotting
    doplots, xaxis, array, title, psym, coo, 1, $
      fineaxis, fitcoef, fitcolor, colormarks, retmode, modename, $
      plotyrng, plotxrng
endwhile

!p.noerase = noerasestatus
doplots, xaxis, array, posttitle, psym, coo, 0, $
  fineaxis, fitcoef, fitcolor, colormarks, retmode, modename, $
  keyword_set(keepzoom)?plotyrng:0, keyword_set(keepzoom)?plotxrng:0

case n_elements(coo) of
    0: return, keyword_set(default) ? default: 0
    1: return, 0
    else:
endcase

;discard last point if it is a ranged mode
cnt = (size(coo, /dim))[0]
if (cnt and 1) and ((retmode eq 1) or (retmode eq 2)) then begin
  coo = coo[0:cnt-2, *]
  tags = tags[0:cnt-2, *]
endif

;reordering of output
case retmode of
    0: ret = coo
    1: ret = reform(coo[*, 0], 2, (size(coo, /dim))[0]/2)
    2: ret = reform(coo[*, 1], 2, (size(coo, /dim))[0]/2)
    3: ret = coo
    4: ret = coo
    5: ret = coo
endcase
if retmode eq 1 or retmode eq 2 then $
  tags = reform(tags, 2, (size(coo, /dim))[0]/2)

;sorting
if (sortflag eq 1 or sortflag eq 2) then begin
    switch retmode of
        sortflag: begin
            for i = 0, cnt/2 - 1 do begin
                si = sort(ret[*, i])
                ret[*, i] = ret[si, i]
                tags[*, i] = tags[si, i]
            endfor
            si = sort(ret[0, *])
            ret = ret[*, si]
            tags = tags[*, si]
            break
        end
        3:
        4:
        5:
        0: begin
            si = sort(ret[*, sortflag-1])
            ret = ret[si, *]
            tags = tags[si, *]
            break
        end
    endswitch
endif

return, ret

end

