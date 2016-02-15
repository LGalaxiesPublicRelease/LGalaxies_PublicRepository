pro find_simple, image, x, y, width=width, thresh=thresh, mask=mask

   if (NOT keyword_set(width)) then width = 5
   if (NOT keyword_set(thresh)) then thresh = 1e5

   dims = size(image, /dimens)
   x = 0
   y = 0

   ;----------
   ; Convolve the image with a boxcar.

   kern = fltarr(width,width) + 1.0
   cimg = convol(image, kern)

   ;----------
   ; For every pixel, ask if it is a local maximum by comparing it
   ; to its 8 neighbors.  Note this isn't done exactly correctly, since
   ; we allow the image to wrap when we do this comparisons.

   qmax = cimg GT shift(cimg,-1,0) $
      AND cimg GT shift(cimg,1,0) $
      AND cimg GT shift(cimg,0,-1) $
      AND cimg GT shift(cimg,0,1) $
      AND cimg GT shift(cimg,-1,-1) $
      AND cimg GT shift(cimg,-1,1) $
      AND cimg GT shift(cimg,1,-1) $
      AND cimg GT shift(cimg,1,1)

   ;----------
   ; Select all pixels that are local maxima and exceed the specified
   ; threshold.

   if (keyword_set(mask)) then $
    indx = where(cimg GT thresh AND qmax AND mask EQ 0) $
   else $
    indx = where(cimg GT thresh AND qmax)
   if (indx[0] EQ -1) then return
   x = float(indx MOD dims[0])
   y = float(indx / dims[0])

   ;----------
   ; Recenter the positions of these peaks using an intensity-weighted
   ; centroid.

   djs_photcen, x, y, image, calg='iweight', $
    cbox=width, cmaxiter=2, cmaxshift=0.5*width

   return
end
