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
;;

;+
; NAME:
;	READMODEL
;
; PURPOSE:
;	This function reads a spectra from a datafile and convolve it
;	down to the specified resolution.
;
; CATEGORY:
;	Pato's data reduction potpourri
;
; CALLING SEQUENCE:
;
;	Data = READMODEL(File, Wav)
;
; INPUTS:
;	File:	File containing the model. Either a FITS or text file.
;       Wav:    Wavelength scale. If not specified the same scale in
;               the model file is used.
;
;	
; KEYWORD PARAMETERS:
;	MODCOLS: A 2-element array containing the columns with the
;	         wavelength scale and the model, respectively.
;       DELTAWAV: Width of the convolution kernel.
;       RESOLUTION: Spectral resolution with reference to the middle
;                wavelength. 
;       VELOCITY: Relative velocity (in units of the speed of light)
;                to shift the spectrum. Negative is towards Earth.
;       SHOW:    Shows a plot.
;       
;
; OUTPUTS:
;	This function returns the spectrum from file File, in the
;	wavelength scale Wav after convolving according to DELTAWAV or
;	RESOLUTION.
;
; MODIFICATION HISTORY:
; 	Written by:	Pato Rojo, Cornell.  2006
;			pato@astro.cornell.edu
;-

function readmodel, modelfile, wav,                                  $
                    modcols=modcols, quiet=quiet,                    $
                    deltawav=deltawav, resolution=resolution,        $
                    velocity=veloverc, show=show

if keyword_set(resolution) then                         $
  if (where(resolution gt 50))[0] eq -1 then            $
  message, "Check use of RESOLUTION keyword in READMODEL(). " + $
  "It was used to identify DELTAWAV before, but now that has" + $
  "its own keyword, and RESOLUTION has its true meaning as mean(wav)/DELTAWAV"

on_error, 2
if n_params() lt 2 then     $
  message, "Usage: model = READMODEL(modelfile, wav, modcols=modcols, " + $
  "resol=resol, veloc=veloc, /show))"

if ~ keyword_set(modcols)    then modcols    = [0,1]
if ~ keyword_set(quiet)      then quiet      = 0
if ~ keyword_set(resolution) and ~ keyword_set(deltawav) then $
  noconv = 1 else noconv = 0

if ~ quiet then $
  print, "Getting theoretical model from file '" + modelfile + "':"

;Read model
rpos = strlen(modelfile)-5
if ((strpos(modelfile,'.fits') eq rpos           $
     or strpos(modelfile,'.FITS') eq rpos        $
     or strpos(modelfile,'.fits.gz') eq rpos-3   $
     or strpos(modelfile,'.FITS.GZ') eq rpos-3)) then begin
    if ~ quiet then $
      print, ' Reading model from FITS. Channels ', modcols[0], $
      '(equisp wav), ', modcols[1], '(data)...', format='(a,i0,a,i0,a,$)'
    model = (readfits(modelfile,h,/silent))[*,[modcols]]
endif else begin
    if ~ quiet then $
      print,' Reading it from column ', modcols[1]+1, $
      ' (assuming equispaced wavenumber in column ', modcols[0]+1, ')...', $
      format='(a,i0,a,i0,a,$)'
    model = (arrayfromfile(modelfile))[*,modcols]
endelse

if n_params() eq 1 then wav = rotate(1e4 / model[*,0], 2)

;Confirm wav and find extremes from it
nx = (size(wav, /dim))[0]
goodspec = bytarr(size(reform(wav[0,*,*]), /dim)) + 1 

ncyc = (size(goodspec, /dim))[0]
npos = n_elements(goodspec) / ncyc
velarray   = dblarr(ncyc, npos)
if keyword_set(veloverc) then velarray += veloverc
if n_elements(wav)/nx ne n_elements(velarray) then $
  message, string("VELOCITY(", n_elements(velarray), " spectra) and WAV(", $
                  n_elements(wav)/nx, " spectra) needs to be consistent" + $
                  "with one another (vel could be fixed as well)", $
                  format='(a,i0,a,i0,a)')
on_error, 0

;Define resolution element
exwn  = rotate(1e4 / minmax((reform(wav, nx, ncyc*npos)) $
                            [*, where(reform(goodspec, ncyc*npos) eq 1)]), 2)
mwav = total(wav,1)/nx
exdwav = keyword_set(resolution)?1e4/(2.35482*min(mwav)*min(resolution)): $
  (keyword_set(deltawav)?max(deltawav):0)

;discard out-of-range
valrng = where(model[*,0] gt exwn[0]-50*exdwav and $
               model[*,0] lt exwn[1]+50*exdwav)
if valrng[0] eq -1 then valrng = lindgen((size(model,/dim))[0])
model = model[valrng, *]

iw = model[1,0] - model[0,0]
valrng = 0
modwav = rotate(1e4 / model[*,0], 2)
if ~ quiet then print, ' done'

if ~ quiet then $
  print," Convolving and interpolating model to data's sampling...", $
  format='(a,$)'

;prepare convolution element
if ~ noconv then begin
    case size(resolution, /n_dim) of
        0: dwav_res = keyword_set(resolution)? $
          1e4/(2.35482*mwav*resolution) : (mwav*0+iw/5.0)
        1: dwav_res = 1e4/(2.35482*mwav*(resolution # replicate(1.0, npos)))
        2: dwav_res = 1e4/(2.35482*mwav*resolution)
    endcase
    case size(deltawav, /n_dim) of
        0: dwav = keyword_set(deltawav)?(mwav*0+deltawav):dwav_res
        1: dwav = deltawav # replicate(1.0, npos)
        2: dwav = deltawav
    endcase

    nk = long(20.0*dwav/iw)
endif
;stop

;do convolution
wavmodel = dblarr(nx, ncyc, npos)
for j=0, npos-1 do for i=0, ncyc-1 do begin

    mdl = model[*,1]
    if noconv then cnvmod = mdl $
    else begin
        if nk[i,j] gt 0 then begin
            xx = (dindgen(nk[i,j]) - nk[i,j]/2) * iw / dwav[i,j]
            telprof = 1.0d / (sqrt(2*!pi) * dwav[i,j] / iw) * exp(-0.5*xx^2)
            
            mdl = convol(mdl, telprof, 1)
        endif
        cnvmod = rotate(mdl, 2)
    endelse

;stop
    c = spl_init(modwav, cnvmod, /double)
    wavmodel[*, i, j] = spl_interp(modwav, cnvmod, c, $
                                   wav[*,i,j] * $
                                   (1.0d - velarray[i,j]), /double)
endfor
if ~ quiet then print,' done'

if keyword_set(show) then $
  plot, modwav, cnvmod, xtitle='Wavelength', ytitle='Convolved Model'

return, wavmodel

end

