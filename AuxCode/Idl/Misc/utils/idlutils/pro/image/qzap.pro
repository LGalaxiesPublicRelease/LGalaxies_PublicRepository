;+
; NAME:
;   qzap.pro
;
; PURPOSE:
;   Remove cosmic rays from a 2-D image.
;
; CALLING SEQUENCE:
;   qzap, name, outname, [ outmaskname, skyfiltsize=skyfiltsize, $
;    boxsize=boxsize, nsigma=nsigma, /nofluxratio, maxiter=maxiter, $
;    fluxcompare=fluxcompare, nrings=nrings, path=path, nzap=nzap ]
;
; INPUTS:
;   name       - 2-D image array, or name of input FITS file.
;   outname    - Output image array, or name of output FITS file.
;
; OPTIONAL INPUTS:
;   outmaskname- Output mask array, or name of output FITS file.
;   skyfiltsize- Boxsize for computing local sky value; default to 15.
;   boxsize    - Boxsize for computing local median; default to 5.
;   nsigma     - Rejection threshhold in sigma; default to 4.
;   fluxratio  - Comparison value for identifying cosmics; default to 0.15
;   maxiter    - Number of zapping iterations; default to 2.
;   nofluxcompare - Set to disable the flux comparison algorithm, which
;                is the "black magic" heart of this routine.
;   nrings     - Radius of cosmic ray neighbors to also zap; default to 1.
;   path       - Input/output path name
;
; OPTIONAL OUTPUTS:
;   NZAP       - Number of pixels zapped.
;
; COMMENTS:
;   Based on the tried and true IRAF QZAP routine by Mark Dickinson.
;   Results from IDL qzap.pro and IRAF QZAP are found to be virtually
;   identical.
;
; PROCEDURES CALLED:
;   djs_iterstat
;
; REVISION HISTORY:
;   20-Aug-1999  Written by Cullen Blake & David Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro qzap, name, outname, outmask, skyfiltsize=skyfiltsize, $
 boxsize=boxsize, nsigma=nsigma, fluxratio=fluxratio, maxiter=maxiter, $
 nofluxcompare=nofluxcompare, nrings=nrings, path=path, nzap=nzap

   if (NOT keyword_set(skyfiltsize)) then skyfiltsize=15
   if (NOT keyword_set(boxsize)) then boxsize=5
   if (NOT keyword_set(nsigma)) then nsigma=4
   if (NOT keyword_set(fluxratio)) then fluxratio=1
   if (NOT keyword_set(maxiter)) then maxiter=2
   if (NOT keyword_set(nofluxcompare)) then nofluxcompare=0
   if (NOT keyword_set(nrings)) then nrings=1
   if (NOT keyword_set(path)) then path = ''

   if (size(name, /tname) EQ 'STRING') then outimg=readfits(path+name) $
    else outimg=name

   dims = size(outimg, /dimens)
   outmask = bytarr(dims[0], dims[1])

   print, 'Computing image sigma...'
   djs_iterstat, outimg, sigma=sigval, sigrej=5

   if (skyfiltsize eq 0) then $
    djs_iterstat,outimg,median=skyimage $
    else skyimage = median(outimg, skyfiltsize)

   skysubimage = outimg - skyimage
   nzap = 0
   iter = 0
   nbad = 1

   while (iter LT maxiter and nbad GT 0) do begin
      iter = iter + 1
      print, 'Iteration', iter

      fmedimage = median(skysubimage, boxsize)

      crimage = skysubimage - fmedimage

      peaksimage = crimage GT nsigma*sigval

      if (NOT keyword_set(nofluxcompare)) then begin
         i = where(peaksimage NE 0)
         peaksimage[i] = (fmedimage[i] / crimage[i]) LT fluxratio
      endif

      peaksimage = smooth(float(peaksimage), 1+2*nrings, /edge) GT 0

      ibad = where(peaksimage NE 0, nbad)
      if (nbad GT 0) then begin
         outimg[ibad] = skyimage[ibad] + fmedimage[ibad]
         outmask[ibad] = 1
      endif

      print, 'Number zapped = ', nbad
      nzap = nzap + nbad

   endwhile

   if (size(outname, /tname) EQ 'STRING') then begin
      writefits, path+outname, outimg
   endif else begin
      outname = outimg
   endelse

   if (size(outmaskname, /tname) EQ 'STRING') then begin
      writefits, path+outmaskname, outmask
   endif else begin
      outmaskname = outmask
   endelse

end
;------------------------------------------------------------------------------
