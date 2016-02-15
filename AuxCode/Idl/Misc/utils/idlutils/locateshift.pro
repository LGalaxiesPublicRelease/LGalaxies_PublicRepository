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
;	LOCATESHIFT
;
; PURPOSE:
;	This function locates the shift at which two arrays look the
;	most similar.
;
; CATEGORY:
;	Pato's data reduction potpourri
;
; CALLING SEQUENCE:
;
;	shift = LOCATESHIFT(Dat, Oref)
;
; INPUTS:
;       Dat:    Data array
;       Oref:   Reference array for comparison
;
; OUTPUTS:
;	This function returns the lag where the reference and the
;	array are the most similar.
;
; PROCEDURE:
;	Algorithm description.
;
; EXAMPLE:
;	Describe example here
;
;		F = LOCATESHIFT(Data, ref)
;
; MODIFICATION HISTORY:
; 	Written by:	Pato Rojo, Cornell.  2006
;			pato@astro.cornell.edu
;-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function dens_solve, sft, rms, nsofar, nsamp, nsampisodd, lastdelta, $
                     nsft, maxsample, oldx, conv, xextreme
compile_opt idl2, hidden

;get the lowest nsamp points, and get a weighted low
lowi = (sort(rms[0:nsofar]))[0:nsamp-1]
xextreme = total(sft[lowi]/rms[lowi]) / total(1.0/rms[lowi])

;if converged from previous value or maxsampling reached then stop
if abs(xextreme - oldx) le conv*oldx then return, 0
if nsofar + nsamp + 1 ge maxsample then return, 1
refsft = sft[(sort(abs(sft-xextreme)))[0]]

;get newest nsample value around bin containing candidate
nsft = refsft + (dindgen(nsamp)-(nsamp-1)/2.0) * lastdelta
lastdelta /= 2.0
if nsampisodd then begin
    if xextreme lt refsft then nsft -= lastdelta $
    else nsft += lastdelta
endif
 
;store old value
oldx = xextreme

return, -1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function parab_fit, sft, rms, error, noerror, nsft, maxsample, conv, $
                    xextreme
compile_opt idl2, hidden

message, 'Sorry, but parab_fit alternative is not working'
;find fit accordingly
if keyword_set(noerror) then $
  coef = poly_fit(sft[0:nsofar], rms[0:nsofar], 2) $
else $
  coef = poly_fit(sft[0:nsofar], rms[0:nsofar], 2, $
                  meas=error[0:nsofar])

xextreme = - coef[1] / (2.0*coef[2])
yextreme = coef[0] - coef[2]*xextreme*xextreme

;if yextreme is a maximum, complain
if (where(rms gt yextreme))[0] eq -1 then $
  message, "yextreme is a maximum. Aborting"

;if it is the first tiome then reference is the central point,
;otherwise is the last point
refsft=sft[nsofar]

;if found satisfactory minimum then exit
if abs(refsft - xextreme) le conv then return, 0

;throw away the highest point, and the candidate center will be the
;new point
idx = reverse(sort(rms[0:nsofar]))
nsofar--
sofar--
sft[0:nsofar] = sft[idx[1:*]]
error[0:nsofar] = error[idx[1:*]]
nsft=[xextreme]

if sofar ge maxsample-1 then return, 1

return, -1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function rmsshift, dat, ref, sft, datcoefs, xdat, xref, lin
;return rms values for sft not computed.
compile_opt idl2, hidden

rms = double(sft)
mmdat = minmax(xdat)
mmref = minmax(xref)

for i=0, n_elements(sft)-1 do begin
    sf = sft[i]

    ;subset to avoid extrapolation
    maxv = min([mmdat[1]+sf, mmref[1]])
    minv = max([mmdat[0]+sf, mmref[0]])
    if maxv le minv then $
      message, "Requested shift of " + string(sf) + $
      " leaves no common x axis of target and reference to compare"

    refi = where((xref-sf le maxv) and (xref-sf ge minv))
    newref = ref[refi]
    newx = (xref - sf)[refi]

   ;interpolate to the new value
    newdat = spl_interp(xdat, dat, datcoefs, newx, /double)

    ;linear scale to the best possible
    if lin then begin
        fitlin = linfit(newdat, newref)
        newdat = fitlin[0] + newdat * fitlin[1]
    endif

    ;get RMS
    rms[i] = sqrt(mean((newdat-newref)^2))
endfor

return, rms

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function locateshift, dat, oref, $
                      inishift=inishift, conv=conv, maxsample=maxsample, $
                      noerror=noerror, osamp=osamp, parab=parab, $
                      nsamp=nsamp, status=s, numbershift=sofar, $
                      linearscale=linearscale
compile_opt idl2

;number of samples
if not keyword_set(nsamp) then nsamp = 5
if nsamp and 1 then nsampisodd = 1 else nsampisodd = 0
;If not set then assume that ideal shift is within a pixel
if not keyword_set(inishift)  then inishift = 1
;if not set, then only stop on the exact result.
if not keyword_set(conv)      then conv = double(1e-7)
;maximum number of shift samples that are taken
if not keyword_set(maxsample) then maxsample = 90
;If using density increase method
if not keyword_set(parab) then parab = 0

if not keyword_set(linearscale) then linearscale = 0

nr = n_elements(oref)
xref = dindgen(nr)

;do sampling at half point of the reference array
if keyword_set(osamp) then begin
    ox=xref
    nr = 2*nr + 1
    xref = dindgen(nr) / 2.0
    datcoefs = spl_init  (ox, oref, /double)
    ref      = spl_interp(ox, oref, datcoefs, xref)
endif else begin
    ref = oref
endelse

;precompute spline coefficients
xdat = dindgen(n_elements(dat))
datcoefs   = spl_init(xdat, dat, /double)
lastdelta = inishift/2.0

;make nsft according to nsamp
nsft = (dindgen(nsamp)-(nsamp-1)/2.0) * lastdelta
if not nsampisodd then nsft = [nsft-lastdelta/2.0, lastdelta*(nsamp/2)]
n=n_elements(nsft)
nsft[nsamp/2:n-2] = nsft[nsamp/2+1:n-1]
nsft[n-1] = 0

rms = dblarr(maxsample)
sft = dblarr(maxsample)
error = dblarr(maxsample)
sofar = 0
cnt = 0
oldx = 0

;loop until desired result is found
while 1 do begin

;compute rms and add it to sample
    nsofar = sofar+n_elements(nsft)-1
    nrms = rmsshift(dat, ref, nsft, datcoefs, xdat, xref, linearscale)

    exact = where(nrms eq 0)
    if exact[0] ne -1 then begin
        xextreme = nsft[exact[0]]
        s = 0
        break
    endif
;add new values to the sample
    rms[sofar:nsofar] = nrms
    sft[sofar:nsofar] = nsft
    if parab then begin
        error *= 2.0
        error[sofar:nsofar] = sqrt(nrms)
    endif else begin
        idx = sort(sft[0:nsofar])
        rms[0:nsofar] = rms[idx]
        sft[0:nsofar] = sft[idx]
    endelse
    sofar = nsofar+1

;if trying parabolic fit
    if parab then s = parab_fit(sft, rms, error, noerror,       $
                                nsft, maxsample, conv,          $
                                xextreme)                       $
;if trying density increase method
    else s = dens_solve(sft, rms, nsofar, nsamp, nsampisodd,    $
                        lastdelta, nsft, maxsample, oldx, conv, $
                        xextreme)

    if s ge 0 then break

endwhile

return, xextreme

end
