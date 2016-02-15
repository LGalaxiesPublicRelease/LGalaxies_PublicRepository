;+
; BUGS:
;   - No correct code header!
;
; D. Finkbeiner - modified to compare probabilities of first and
;                 second peaks.  -DPF 4 Nov 2000
; set errflag if best value on edge - DPF 26 Nov 2000
;-
function offset_from_pairs, x1, y1, x2, y2, dmax=dmax, binsz=binsz, $
 minpeak=minpeak, errflag=errflag, bestsig=bestsig, verbose=verbose, $
 minprobratio=minprobratio

   if (not keyword_set(minprobratio)) then minprobratio= 1e-5
   sigcut = 5.0 ; A detection is sigcut sigma above the mean
   errflag = 0B
   if (NOT keyword_set(binsz)) then binsz = 1
   if (NOT keyword_set(minpeak)) then minpeak = 3.0
   xyshift = [0,0]

   num1 = n_elements(x1)
   num2 = n_elements(x2)

   if (num1 > num2) GT 32000 then $
    splog, 'WARNING:  Awful lot of stars...', num1, num2

   ;----------
   ; Construct an image and populate with the vector offsets for all
   ; distances between object positions X1,Y1 and X2,Y2.

;   nhalf = fix(2*dmax/binsz) + 1
   nhalf = fix(dmax/binsz)  ; I think this is what is intended - DPF
nhalf = nhalf > 2 ; ??? We need this to not crash further down - DJS
   nn = 2*nhalf + 1
   img = fltarr(nn,nn)

; Speed up change - 29 March 2001  - DPF
   xoff = fltarr(num1, num2)
   yoff = fltarr(num1, num2)
   for i=0L, num1-1 do begin
      xoff[i, *] = (x2 - x1[i])/float(binsz)
      yoff[i, *] = (y2 - y1[i])/float(binsz)
   endfor 
   xoff = reform(xoff, num1*num2) + (nhalf + 1)
   yoff = reform(yoff, num1*num2) + (nhalf + 1)

   ; Important: Only call populate_image once - calling overhead is large.
   populate_image, img, xoff, yoff, assign='cic'

   ;----------
   ; Find the peak value in this image smoothed 3x3

   mx = where(img EQ max(img))
   msk = img-img
   msk[mx] = 1.
   msk = smooth(msk, 5, /edge) EQ 0
   img2 = msk*img ; remove brightest pixel and those nearby

   std = stddev(img)

   nsig1 = (max(img)-mean(img))/std
   nsig2 = (max(img2)-mean(img))/std  ; use same mean and std
   prob1 = gauss_pdf(-double(nsig1))
   prob2 = gauss_pdf(-double(nsig2))
   if (keyword_set(verbose)) then begin
      splog, 'nsig1, nsig2, prob1, prob2'
      splog, nsig1, nsig2, prob1, prob2
   endif

   smimg = smooth(img, 3, /edge)
   maxval = max(smimg, imax)
   bestsig = (maxval-mean(smimg))/stddev(smimg)

   ; If prob1/prob2 is not really small then return
   IF prob1 GT minprobratio * prob2 THEN BEGIN 
      errflag = 1B
      splog, 'Maxval is only', bestsig, ' sigma above mean'
      splog, 'Odds of ', prob1/prob2, ' are not good enough'
      splog, xyshift
   ENDIF 

   if (maxval LE 0) then BEGIN 
       splog, 'WARNING: maxval LE 0!'
       return, xyshift
   endif
   xmax = imax MOD nn
   ymax = imax / nn

   ;----------
   ; Find the sub-pixel peak by looking at only the nearest 3x3 sub-image,
   ; of the peak, then computing an intensity-weighted centroid

   if (xmax GT 0 AND xmax LT nn-1 AND ymax GT 0 AND ymax LT nn-1) then begin
      subimg = img[(xmax-1)>0 : (xmax+1)<(nn-1), (ymax-1)>0 : (ymax+1)<(nn-1)]
      maxval = total(subimg)
      xmax = xmax + transpose(total(subimg,2)) # [-1,0,1] / maxval
      ymax = ymax + transpose(total(subimg,1)) # [-1,0,1] / maxval
   endif else begin 
      errflag = 1B
   endelse 

   if (maxval GT minpeak) then $
    xyshift = ([xmax,ymax] - nhalf - 1) * binsz

   return, xyshift
end
