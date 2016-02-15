;+
; NAME:
;	sshift
; PURPOSE: (one line)
;	Shift data using a damped sinc function for fractional part.
; DESCRIPTION:
;
;	This function will shift an array of data pointed to by x and
;	extending for n points.  The amount of the shift is given by shift.
;	The result of the operation is placed at xp.  A shift that is within
;	0.0001 of a whole number is treated to be that of the whole number.  If
;	the shift is by an integral number of pixels then the shift involves
;	reindexing the data, no interpolation is done.  If the shift is some
;	non-integral amount then the data is resampled using a damped sinc
;	function.
;
;	The sense of the shift is as follows: think of the array plotted on a
;	fixed scale.  A shift of 1 corresponds to shifting the data by one data
;	point to the right relative to the fixed scale, ie. x[3]=xp[4].
;
;	The data will fall off one end or another of the output vector as a
;	result of the shift.  However, this is not as significant as the edge
;	effect, the convolution is not complete for any data point within 10
;	points of the edge, so those points cannot be trusted.  The missing
;	points in the convolution are assumed to be equal to the end points.
;
; CATEGORY:
;       Numerical
; CALLING SEQUENCE:
;	xp = sshift(x,shift)
; INPUTS:
;	x     - Input data array to be shifted.
;	shift - Amount to shift (negative to left).
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;	Return value is the shifted array.
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
;	The input and output arrays cannot be the same.
; PROCEDURE:
; MODIFICATION HISTORY:
;	Adapted from Zodiac routine: shiftc/sshift
;	  Marc W. Buie, Lowell Observatory, 1992 October 2
;-

function sshift,x,shift

   EPS     = 1.0e-5 ; Smallest fractional shift allowed.
   DAMPFAC = 3.25   ; Damping factor for gaussian.
   NS      = 21     ; Number of points in the sinc convolution kernal.

   npts = n_elements(x)
   xp   = fltarr(npts)

; First, split the desired shift into a fractional and integer part.
   ishift = long(shift)
   fshift = shift - ishift

; Do the fractional shift first (if necessary).
   if ( ( abs(fshift) gt EPS ) and ( abs(fshift) lt 1.0-EPS) ) then begin

   ; Initialize the sinc array.
      y = fshift - (findgen(NS)-10)
      py = !pi * y
      sinc = exp( -y^2/DAMPFAC^2 ) * sin(py)/py
      sinc = rotate(sinc,2)

   ; Convolve the sinc array with the input data.  This is the shift.
      for point=0L,n_elements(x)-1 do begin
         lobe = (indgen(NS) - 10) + point
         vals = fltarr(NS)

         z = where(lobe lt 0,count)
         if (count ne 0) then vals[z] = x[0]

         z = where(lobe ge 0 and lobe lt npts,count)
         if (count ne 0) then vals[z] = x[lobe[z]]

         z = where(lobe ge npts,count)
         if (count ne 0) then vals[z] = x[npts-1]

         xp[point] = total(sinc*vals)

      endfor

; No fractional shift, just copy the data.
   endif else begin

      xp = x

   endelse

; Now perform the integer shift.
   left = -ishift
   right = npts-1 - ishift

   if (left lt 0) then begin
      xp[ishift:npts-1] = xp[0:npts-1-ishift]
      xp[0:ishift-1] = xp[ishift]
   endif else if (left gt 0) then begin
      xp[0:npts-1-left] = xp[left:npts-1]
      xp[npts-left:npts-1] = xp[npts-1-left]
   endif

   return,xp

end
