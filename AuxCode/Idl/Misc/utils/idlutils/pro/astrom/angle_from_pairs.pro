function angle_from_pairs, x1, y1, x2, y2, dmax=dmax, binsz=binsz, $
 bestsig=bestsig, angrange=angrange, verbose=verbose

   sigcut = 5.0 ; A detection is sigcut sigma above the mean
   if (NOT keyword_set(binsz)) then binsz = 1
   xyshift = [0,0]

   num1 = n_elements(x1)
   num2 = n_elements(x2)

   if (num1 > num2) GT 32000 then $
    splog, 'WARNING:  Awful lot of stars...', num1, num2

   ;----------
   ; Construct an image and populate with the vector offsets for all
   ; distances between object positions X1,Y1 and X2,Y2.

   nhalf = fix(dmax/binsz)
   nn = 2*nhalf + 1

   xy1 = [[x1], [y1]]
   
   if NOT keyword_set(angrange) then  angrange = [-5, 5]
   dangle = 1.
   nang = (max(angrange)-min(angrange))/dangle

   theta = (0.5+findgen(nang))/nang * (angrange[1]-angrange[0])+angrange[0]
   pk = fltarr(nang)
   for k=0, nang-1 do begin 
      img = fltarr(nn,nn)
      angrad = theta[k] *!dtor
      mm = [[cos(angrad), sin(angrad)], [-sin(angrad), cos(angrad)]]
      xyr = xy1 # mm

; Speed up change - 29 March 2001  - DPF
      xoff = fltarr(num1, num2, /nozero)
      yoff = fltarr(num1, num2, /nozero)
      for i=0L, num1-1 do begin
         xoff[i, *] = (x2 - xyr[i, 0])/float(binsz)
         yoff[i, *] = (y2 - xyr[i, 1])/float(binsz)
      endfor 
      xoff = reform(xoff, num1*num2) + (nhalf + 1)
      yoff = reform(yoff, num1*num2) + (nhalf + 1)

      ; Important: Only call populate_image once - calling overhead is large.
      populate_image, img, xoff, yoff, assign='cic'

      pk[k] = max(img)
      if (keyword_set(verbose)) then $
       splog, 'Angle=', theta[k], '  Peak=', pk[k]
   endfor

   bestsig = (max(pk, ind)-mean(img))/stddev(img,/double)
   ang = theta[ind]

   return, ang
end
