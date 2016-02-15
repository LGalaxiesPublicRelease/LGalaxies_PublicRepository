
;+
; NAME:
;   djs_maskinterp
;
; PURPOSE:
;   Interpolate over masked pixels in a vector, image or 3-D array.
;
; CALLING SEQUENCE:
;   ynew = djs_maskinterp( yval, mask, [ xval, iaxis=, /const ] )
;
; INPUTS:
;   yval       - Y values; 1-, 2-, or 3-dimensional
;   mask       - Mask values correspoding to YVAL; interpolate over all pixels
;                where MASK is not 0
;
; OPTIONAL INPUTS:
;   xval       - X (abscissa) values corresponding to YVAL; otherwise a
;                regular grid is assumed
;   iaxis      - Axis along which to interpolate if YVAL has more than one
;                dimension; required keyword in that case; dimensions are
;                0-indexed, so the X axis is IAXIS=0
;   const      - The default is to linearly interpolate beyond the endpoints
;                of good data.  Setting this keyword instead copied the
;                first (last) good points for the data beyond the first (last)
;                good points.
;
; OUTPUTS:
;   ynew       - Y values after linearly interpolating over masked pixels
;
; COMMENTS:
;   The IDL function INTERPOL is used for linear interpolation.
;
;   If no good points exist in a vector, then that vector is returned
;   unchanged.
;
; EXAMPLES:
;   Create a sin-wave function, and interpolate across points at the beginning
;   and in the middle of this function:
;     xval=findgen(100)/10
;     yval=sin(xval)
;     splot,xval,yval
;     mask=bytarr(100)
;     mask[0:10]=1
;     mask[40:60]=1
;     ynew = djs_maskinterp(yval, mask, xval)
;     plot,xval,yval
;     oplot,xval,ynew,ps=2
;
; BUGS:
;   This routine only supports 1-D, 2-D, and 3-D arrays.
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINS:
;   djs_maskinterp1()
;
; REVISION HISTORY:
;   27-Jan-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function djs_maskinterp1, yval, mask, xval, const=const

   ibad = where(mask NE 0, nbad)
   if (nbad EQ 0) then $
      return, yval

   igood = where(mask EQ 0, ngood)
   if (ngood EQ 0) then $
      return, yval
   if (ngood EQ 1) then $
      return, yval*0 + yval[igood[0]]

   ynew = yval
   ny = N_elements(yval)

   if (keyword_set(xval)) then begin
      ii = sort(xval) ; We must sort, since INTERPOL expects us to
      ibad = where(mask[ii] NE 0)
      igood = where(mask[ii] EQ 0)
      ynew[ii[ibad]] = interpol(yval[ii[igood]], xval[ii[igood]], xval[ii[ibad]])

      if (keyword_set(const)) then begin
         if (igood[0] NE 0) then $
          ynew[ii[ 0 : igood[0]-1 ]] = ynew[ ii[igood[0]] ]
         if (igood[ngood-1] NE ny-1) then $
          ynew[ii[ igood[ngood-1]+1 : ny-1 ]] = ynew[ ii[igood[ngood-1]] ]
      endif
   endif else begin
      ynew[ibad] = interpol(yval[igood], igood, ibad)

      if (keyword_set(const)) then begin
         if (igood[0] NE 0) then $
          ynew[0:igood[0]-1] = ynew[igood[0]]
         if (igood[ngood-1] NE ny-1) then $
          ynew[igood[ngood-1]+1:ny-1] = ynew[igood[ngood-1]]
      endif
   endelse

   return, ynew
end
;------------------------------------------------------------------------------
function djs_maskinterp, yval, mask, xval, iaxis=iaxis, _EXTRA=extra

   dims = size(yval, /dimens)
   ndim = size(yval, /n_dimens)

   if (size(mask,/n_dimens) NE ndim $
    OR total(size(mask,/dimens) NE dims) NE 0) then $
    message, 'MASK and YVAL are not the same dimensions'
   if (keyword_set(xval)) then $
    if (size(xval,/n_dimens) NE ndim $
     OR total(size(xval,/dimens) NE dims) NE 0) then $
     message, 'XVAL and YVAL are not the same dimensions'

   if (ndim EQ 1) then begin

      ynew = djs_maskinterp1(yval, mask, xval, _EXTRA=extra)

   endif else begin
      if (n_elements(iaxis) EQ 0) then $
       message, 'Must declare IAXIS if YVAL has more than 1 dimension'
      if (iaxis LT 0 OR iaxis GT ndim-1 OR iaxis-fix(iaxis) NE 0) then $
       message, 'IAXIS invalid'

      ynew = 0 * yval

      case ndim of 
      2 : begin
         if (keyword_set(xval)) then begin
            if (iaxis EQ 0) then begin
               for i=0, dims[1] do $
                ynew[*,i] = djs_maskinterp1(yval[*,i], mask[*,i], $
                 xval[*,i], _EXTRA=extra)
            endif else if (iaxis EQ 1) then begin
               for i=0, dims[0] do $
                ynew[i,*] = djs_maskinterp1(yval[i,*], mask[i,*], $
                 xval[i,*], _EXTRA=extra)
            endif
         endif else begin
            if (iaxis EQ 0) then begin
               for i=0, dims[1] - 1 do $
                ynew[*,i] = djs_maskinterp1(yval[*,i], mask[*,i], _EXTRA=extra)
            endif else if (iaxis EQ 1) then begin
               for i=0, dims[0] - 1 do $
                ynew[i,*] = djs_maskinterp1(yval[i,*], mask[i,*], _EXTRA=extra)
            endif
         endelse
      end
      3: begin
         if (keyword_set(xval)) then begin
            if (iaxis EQ 0) then begin
               for i=0, dims[1] - 1 do $
                for j=0, dims[2] - 1 do $
                 ynew[*,i,j] = djs_maskinterp1(yval[*,i,j], mask[*,i,j], $
                  xval[*,i,j], _EXTRA=extra)
            endif else if (iaxis EQ 1) then begin
               for i=0, dims[0] - 1 do $
                for j=0, dims[2] - 1 do $
                 ynew[i,*,j] = djs_maskinterp1(yval[i,*,j], mask[i,*,j], $
                  xval[i,*,j], _EXTRA=extra)
            endif else if (iaxis EQ 2) then begin
               for i=0, dims[0] - 1 do $
                for j=0, dims[1] - 1 do $
                 ynew[i,j,*] = djs_maskinterp1(yval[i,j,*], mask[i,j,*], $
                  xval[i,j,*], _EXTRA=extra)
            endif
         endif else begin
            if (iaxis EQ 0) then begin
               for i=0, dims[1] - 1 do $
                for j=0, dims[2] - 1 do $
                 ynew[*,i,j] = djs_maskinterp1(yval[*,i,j], mask[*,i,j], $
                  _EXTRA=extra)
            endif else if (iaxis EQ 1) then begin
               for i=0, dims[0] - 1 do $
                for j=0, dims[2] - 1 do $
                 ynew[i,*,j] = djs_maskinterp1(yval[i,*,j], mask[i,*,j], $
                  _EXTRA=extra)
            endif else if (iaxis EQ 2) then begin
               for i=0, dims[0] - 1 do $
                for j=0, dims[1] - 1 do $
                 ynew[i,j,*] = djs_maskinterp1(yval[i,j,*], mask[i,j,*], $
                  _EXTRA=extra)
            endif
         endelse
      end
      else: message, 'Unsupported number of dimensions'
      endcase
   endelse

   return, ynew
end
;------------------------------------------------------------------------------
