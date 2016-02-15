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
;	FXCOR
;
; PURPOSE:
;	This function computes the cross-correlation in fourier space
;
; CATEGORY:
;	Pato's general functions.
;
; CALLING SEQUENCE:
;
;	Result = FXCOR(template, target)
;
; INPUTS:
;	template:	array to which to compare the correlation
;       target:         array to correlate with the template
;
; OPTIONAL INPUTS:
;       Showfourier:    Flag to show fourier space transform.
;       Showccorr:      Flag to show cross crorelation.
;       Showorig:       Flag to show original spectra
;       Nopad:          Flag to indicate that padding is not to be
;                       performed. Otherwise, padding is such that
;                       there is at least one full power of two worth
;                       of zeroes; this is to ensure there is no
;                       overlaping.
;       Linearscale:    Flag to perform the best linear scale of the
;                       spectra before correlating them 
;       Maxprec:        Maximum precision of fit. Default is 0.01
;
; OPTIONAL OUTPUTS:
;       Convwidth:      Convolution width of function that makes
;                       target to resolution of template
;       Factor:         Factor that makes target of template
;                       magnitude. 
;                       Target = Factor * Template x 
;                                Gauss(shift, Convwidth)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	This function returns the shift that most closely makes
;	target to match template.
;
; SIDE EFFECTS:
;          Disturbing noises
;
; RESTRICTIONS:
;
; PROCEDURE:
;	This follows algorithm descrobed in paper by Johyn Tonry and
;	Marc Davis: "A Survey of Galaxy Redshifts. 1. Data Reduction
;	Techniques". AJ 84, pags. 1511-1525.
;
; EXAMPLE:
;	
;
;		ret = FXCOR(template, target)
;
; MODIFICATION HISTORY:
; 	Written by:	Pato Rojo, Cornell.  2005 Mar.
;			pato@astro.cornell.edu
;-

function fxcor, ut, ug, $
                showfourier=showfourier, showccorr=showccorr, $
                showorig=showorig, showfilt=showfilt, showarr=showarr, $
                nopad=nopad, maxshift=maxshift, preinterp=preinterp, $
                maxprec=maxprec, linearscale=linearscale, $
                convwidth=convwidth, factor=factor, $
                filtername=filtername, cuton=cuton, cutoff=cutoff
compile_opt idl2

if ~ keyword_set(showarr) or n_elements(showarr) lt 4 then begin
    showarr = bytarr(5)
    if keyword_set(showorig) then showarr[0] = 1
    if keyword_set(showfourier) then showarr[1] = 1
    if keyword_set(showfilt) then showarr[2] = 1
    if keyword_set(showccorr) then showarr[3] = 1
endif

if ~ keyword_set(linearscale) then linearscale = 0
if ~ keyword_set(maxprec) then maxprec = .01
if ~ keyword_set(preinterp) then preinterp = 0

if keyword_set(nopad) then pad = 0 else pad = 1


message, "Please check this function before using it. It was not working always"

un = n_elements(ut)
if un ne n_elements(ug) then $
  message, "Template and target arrays don't have the same number " + $
  "of elements in 'fxcor.pro'"
if ~ keyword_set(maxshift) then maxshift = un

xo = dindgen(un)
;do a preinterpolation if requested
if preinterp gt 1 then begin
    itpn = un * preinterp
    x = dindgen(itpn) / preinterp
    coef = spl_init(xo, ug, /double)
    itpg = spl_interp(xo, ug, coef, x, /double)
    coef = spl_init(xo, ut, /double)
    itpt = spl_interp(xo, ut, coef, x, /double)
    dx = 1.0 / preinterp
endif else begin
    dx = 1
    x = xo
    itpn = un
    itpg = ug
    itpt = ut
endelse

;pad to avoid cycling
if pad then begin
    padstp = 0
    n = 2^long(alog(itpn)/alog(2)+1.9999)
    g = dblarr(n)
    t = dblarr(n)
    g[padstp:padstp+itpn-1] = itpg
    t[padstp:padstp+itpn-1] = itpt
endif else begin
    padstp = 0
    n = itpn
    g = itpg
    t = itpt
endelse
positk = long(n/2.0 + 1)

x = dindgen(n) * dx
x[positk:*] -= n*dx

if keyword_set(linearscale) then begin
    fitlin = linfit(g, t)
    g      = fitlin[0] + g * fitlin[1]
endif

;find standard deviations
st = total(t*t, /double) / n
sg = total(g*g, /double) / n

;find fourier transforms and cross-correlation
ft = n * fft(t, /double)
fg = n * fft(g, /double)

k = x / (double(n) * dx)

if keyword_set(filtername) and strupcase(filtername) ne 'NONE' then begin
    if ~ keyword_set(cuton) then cuton = 0.01
    if ~ keyword_set(cutoff) then cutoff = 0.4
    ;approximate to the closest frequency
    icuton = (sort(abs(k[0:positk-1] - cuton)))[0]
    icutoff = (sort(abs(k[0:positk-1] - cutoff)))[0]

    nj = (icutoff - icuton) + 1
    jj = dindgen(nj)
    case strupcase(filtername) of
        'WELCH': filter = 1.0 - ((jj - 0.5 * (nj - 1)) / (0.5 * (nj + 1)))^2
        'HANNING': filter = 0.5 * (1.0 - cos(2 * !pi * (jj + 1) / (nj + 1)))
        'UNIFORM': filter = dblarr(nj) + 1.0
        'NONE': filter = 0
    endcase
endif
filtmask = dblarr(positk)
if ~ keyword_set(filter) then filtmask +=1 $
else filtmask[icuton:icutoff] = filter

fltt = ft * [filtmask, rotate(filtmask[1:n-positk], 2)]
fltg = fg * [filtmask, rotate(filtmask[1:n-positk], 2)]

;find real crosscorrelation, and interpolate to desired resolution,
;find the peak by just looking for the maximum
fc = fltg * conj(fltt) / (n * sg * st)
c = fft(fc, /inverse, /double) / n

rc = shift(c, n/2)
rx = shift(x, n/2)

xrange = rx[n-1]-rx[0]
nn = long(n/double(maxprec)) + 1
rxx = findgen(nn) * xrange / double(nn) + rx[0]

coef = spl_init(rx, rc, /double)
rcc = spl_interp(rx, rc, coef, rxx, /double)

ims = [(rotate(where(-rxx[0:positk-1] gt maxshift), 2))[0], $
       nn/2+(where(xx[nn/2:*] gt maxshift))[0]]
if ims[0] eq -1 then ims[0] = 0
if ims[1] eq -1 then ims[1] = nn-1

rdelv = max(rcc[ims[0]:ims[1]], rdeli)
rdel = rxx[rdeli+ims[0]]

delv = max(rcc, deli)
del = rxx[deli]

if rdel ne del then print, $
  "fxcor() detected a shift(" + string(del, format='(g0)') + $
  ") outside the given boundaries"

;OUT OF MEMORY! It doesn't works, it centers at zero.
;kk = xx / double(n)
;k = x / double(n)
;shf = shift(ft, -positk)
;coef = spl_init(k, shf, /double)
;fft = spl_interp(k, shf, coef, kk, /double)
;shg = shift(fg, -positk)
;coef = spl_init(k, shg, /double)
;ffg = spl_interp(k, shg, coef, kk, /double)

;ffc = ffg * conj(fft) / (n * sg * st)
;cc = shift(fft(shift(ffc, nfr), /inverse, /double) / n, -nfr)

;this is to obtain convolution width or factor, bt they are not
;working. Note that the order of k, ft, fg, x, c,fc might be
;completely changed. You have been warned:D
if keyword_set(factor) or keyword_set(convwidth) then begin
    sig = dindgen(n)

;expo[k, sig]
    kk = k # replicate(1.0, n)
    ssig = sig ## replicate(1.0, n)
    ffg = fg # replicate(1.0, n)
    fft = ft # replicate(1.0, n)

    expo = (2 * !pi * ssig * kk)^2 / (2 * n * n)
    expoi = (2 * kk * del * !pi / n)
    
    alp = total(exp(- complex(expo, expoi)) *  $
                (conj(fft) * ffg * exp(complex(0, 2*expoi)) + $
                 fft * conj(ffg)), 1) /                            $
      (2.0 * total(abs(fft)^2 * exp(-2*expo), 1))
    fcn = alp * total(kk * kk * abs(fft)^2 * exp(-2*expo), 1) -     $
      total(kk * kk * exp(-complex(expo, expoi)) * $
            (conj(fft) * ffg * exp(complex(0, 2*expoi)) +    $
             fft * conj(ffg)), 1)

    message, "Sorry but dispersion and factor fitting are not enabled"

    convwidth = sig[0]          ;sig that nullifies fcn
    factor    = alp[0]          ;alp for respective sig
endif


;plotting desired graphics
nshow = total(showarr)
if nshow gt 0 then begin
    !p.multi = [nshow*2, 2, nshow]
    erase
    if showarr[0] then begin
        plot, t[padstp:padstp+itpn], title='TEMPLATE array', $
          charsize=2, /xstyle, xrange=[0,itpn-1]
        plot, g[padstp:padstp+itpn], title='TARGET array', $
          charsize=2, /xstyle, xrange=[0,itpn-1]
    endif
    if showarr[1] then begin
        plot, k[0:positk-1], abs(ft[0:positk-1]), $
          title='Fourier transform of TEMPLATE', charsize=2, $
          /xstyle, xrange=[-0.01, 0.51], /ylog
        tr = !y.crange
        plot, k[0:positk-1], abs(fg[0:positk-1]), $
          title='Fourier transform of TARGET', charsize=2, $
          /xstyle, xrange=[-0.01, 0.51], /ylog
        gr = !y.crange
    endif
    if showarr[2] then begin
        if ~ keyword_set(tr) then $
          tr = alog(minmax((abs(fltt))[where(abs(fltt) gt 0)]))
        plot, k[0:positk-1], abs(fltt[0:positk-1]), $
          title='Fourier transform of TEMPLATE after filter', charsize=2, $
          /xstyle, xrange=[-0.01, 0.51], /ylog, yrange=10^tr
        if ~ keyword_set(gr) then $
          gr = alog(minmax((abs(fltg))[where(abs(fltg) gt 0)]))
        plot, k[0:positk-1], abs(fltg[0:positk-1]), $
          title='Fourier transform of TARGET after filter', charsize=2, $
          /xstyle, xrange=[-0.01, 0.51], /ylog, yrange=10^gr
    endif
    if showarr[3] then begin
        plot, rxx, rcc, title='Cross-correlation. Peak:' + string(del), $
          charsize=2, /xstyle, xrange=[-maxshift, maxshift]
        oplot, x, c, psym=6
        oplot, [del, del], [!y.crange[0], delv]
        lx = lindgen(itpn)
        lxx = lindgen(itpn/double(maxprec) + 1)
        gg_coef = spl_init(lx, g[0:itpn-1], /double)
        gg = spl_interp(lx, g[0:itpn-1], gg_coef, lxx+del, /double)
        plot, abs(t-gg), title='Difference after shifting', charsize=2, $
          /xstyle, xrange=[0,itpn-1], yrange=minmax(g)-max(g)/2
    endif
    !p.multi = 0
endif

;stop

return, del

end
