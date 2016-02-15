;+
; NAME:
;   traceset2xy
;
; PURPOSE:
;   Convert from a trace set to an array of x,y positions
;
; CALLING SEQUENCE:
;   traceset2xy, tset, xpos, ypos
;
; INPUTS:
;   tset       - Structure containing trace set
;
; OPTIONAL KEYWORDS:
;   xpos       - Input positions to evaluate YPOS; if not specified (or 0),
;                then generate an [NX,NTRACE] array of each pixel position
;
; OUTPUTS:
;   xpos       - X positions corresponding to YPOS
;   ypos       - Y centers as an [nx,nTrace] array
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   djs_laxisgen()
;   flegendre()
;   fpoly()
;
; REVISION HISTORY:
;   19-May-1999  Written by David Schlegel, Princeton.
;   01-Dec-2000  Handle scalar xpos correctly - D. Finkbeiner
;   10-Jul-2001  Added fpoly- S.Burles
;-
;------------------------------------------------------------------------------
pro traceset2xy, tset, xpos, ypos

   ; Need 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - traceset2xy, tset, xpos, ypos'
      return
   endif

   if (tset.func EQ 'legendre' OR $
       tset.func EQ 'chebyshev' OR $
       tset.func EQ 'poly') then begin

      ndim = size(tset.coeff, /n_dim)
      dims = size(tset.coeff, /dim)

      if (ndim EQ 1) then begin
         ncoeff = dims[0]
         nTrace = 1
      endif else if (ndim EQ 2) then begin
         ncoeff = dims[0]
         nTrace = dims[1]
      endif else begin
         message, 'TSET.COEFF contains invalid number of dimensions'
      endelse

      nx = long(tset.xmax - tset.xmin + 1)

      xmid = 0.5 * (tset.xmin + tset.xmax)
      xrange = tset.xmax - tset.xmin

      if (NOT keyword_set(xpos)) then $
        xpos = djs_laxisgen([nx, nTrace], iaxis=0) + tset.xmin

      ypos = xpos*0.0
      for iTrace=0, nTrace-1 do begin
         xvec = 2.0 * (xpos[*,iTrace]-xmid)/xrange
         if (tset.func EQ 'poly') then legarr = fpoly(xvec, ncoeff)
         if (tset.func EQ 'legendre') then legarr = flegendre(xvec, ncoeff)
         if (tset.func EQ 'chebyshev') then legarr = fchebyshev(xvec, ncoeff)
         ypos[*,iTrace] = legarr # [tset.coeff[*,iTrace]]
      endfor

   endif else begin
      error, 'Unknown function' + func
   endelse

   return
end
;------------------------------------------------------------------------------
