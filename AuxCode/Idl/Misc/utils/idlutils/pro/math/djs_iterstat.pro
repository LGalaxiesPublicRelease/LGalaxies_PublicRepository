;------------------------------------------------------------------------------
;+
; NAME:
;   djs_iterstat
;
; PURPOSE:
;   Compute the mean, median and/or sigma of data with iterative sigma clipping.
;
; CALLING SEQUENCE:
;   djs_iterstat, image, [invvar=, sigrej=, maxiter=, $
;    mean=, median=, sigma=, mask=, newivar= ]
;
; INPUTS:
;   image:      Input data [N]
;
; OPTIONAL INPUTS:
;   invvar:     Inverse variance for each data point [N]; if set, then
;               the mean and rejection is computed using these errors
;   sigrej:     Sigma for rejection; default to 3.0
;   maxiter:    Maximum number of sigma rejection iterations; default to 10
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   mean:       Computed mean
;   median:     Computed median
;   sigma:      Computed sigma
;   mask:       Mask set to 1 for good points, and 0 for rejected points
;   newivar:    If INVVAR are set, then this is the formal invverse variance
;               of the returned MEAN
;
; PROCEDURES CALLED:
;
; COMMENTS:
;   Iteratively rejects outliers as determined by SIGREJ.  Stop
;   when one of the following conditions is met:
;   (1) The maximum number of iterations, as set by MAXITER, is reached.
;   (2) No new pixels are rejected, as compared to the previous iteration.
;   (3) At least 2 pixels remain from which to compute statistics.  If not,
;       then the returned values are based upon the previous iteration.
;
; BUGS:
;
; REVISION HISTORY:
;   16-Jun-1999  Written by David Schlegel, Princeton
;   11-Sep-2000  Speeded up by Hogg and Eisenstein
;   18-Sep-2000  Note change in MASK values to =1 for good (unrejected) points.
;-
;------------------------------------------------------------------------------
pro djs_iterstat, image, invvar=invvar, sigrej=sigrej, maxiter=maxiter, $
 mean=fmean, median=fmedian, sigma=fsig, mask=mask, newivar=newivar
 
   ; Need 1 parameter
   if N_params() LT 1 then begin
      print, 'Syntax - djs_iterstat, image, [invvar=, sigrej=, maxiter=, $'
      print, ' mean=, median=, sigma=, mask= ]'
      return
   endif

   if (NOT keyword_set(sigrej)) then sigrej = 3.0
   if (NOT keyword_set(maxiter)) then maxiter = 10

   ;----------
   ; Special cases of 0 or 1 data points

   ngood = n_elements(image)
   if (keyword_set(invvar)) then begin
      if (n_elements(invvar) NE ngood) then $
       message, 'Number of elements in IMAGE and INVVAR must agree!'
      if (min(invvar) LT 0) then $
       message, 'Invalid (negative) values for INVVAR!'
      if (total(invvar) GT 0) then begin
         qivar = 1B
         invsig = sqrt(invvar)
      endif
   endif
   if (ngood EQ 0) then begin
;      print, 'No data points'
      fmean = 0.0
      fmedian = 0.0
      fsig = 0.0
      mask = 0B
      newivar = 0.0
      return
   endif
   if (ngood EQ 1) then begin
;      print, 'Only 1 data point'
      fmean = image[0]
      fmedian = fmean
      fsig = 0.0
      mask = 1B
      if (keyword_set(qivar)) then newivar = invvar[0] $
       else newivar = 0.0
      return
   endif

   ;----------
   ; Compute the mean + stdev of the entire image.
   ; These values will be returned if there are fewer than 2 good points.

   if (keyword_set(qivar)) then begin
      mask = invvar GT 0
      fmean = total(image*invvar) / total(invvar)
   endif else begin
      mask = bytarr(ngood) + 1B
      fmean = total(image*mask) / ngood
   endelse
   fsig = sqrt(total((image-fmean)^2*mask) / (ngood-1))
   iiter = 1

   ;----------
   ; Iteratively compute the mean + stdev, updating the sigma-rejection
   ; thresholds each iteration.

   nlast = -1
   while (iiter LT maxiter AND nlast NE ngood AND ngood GE 2) do begin
      nlast = ngood
      if (keyword_set(qivar)) then begin
         mask = abs(image - fmean) * invsig LT sigrej AND invvar GT 0
      endif else begin
         mask = image GT (fmean - sigrej*fsig) $
          AND image LT (fmean + sigrej*fsig)
      endelse
      ngood = total(mask)

      if (ngood GE 2) then begin
         if (keyword_set(qivar)) then $
           fmean = total(image*invvar*mask) / total(invvar*mask) $
          else $
           fmean = total(image*mask) / ngood
         fsig = sqrt( total((image-fmean)^2*mask) / (ngood-1) )
         savemask = mask ; Save for computing the median using the same points
      endif

      iiter = iiter + 1
   endwhile

   if (arg_present(fmedian)) then begin
      if (keyword_set(savemask)) then $
       fmedian = median(image[where(savemask EQ 1)], /even) $
      else $
       fmedian = fmean
   endif

   if (arg_present(newivar)) then begin
      if (keyword_set(qivar)) then begin
         if (keyword_set(savemask)) then $
          newivar = total(invvar*savemask) $
         else $
          newivar = total(invvar)
      endif else begin
         newivar = 0.
      endelse
   endif

   return
end
;------------------------------------------------------------------------------
