;-----------------------------------------------------------------------
;+
; NAME:
;   djs_laxisgen
;
; PURPOSE:
;   Return a longword integer array with the specified dimensions.
;   Each element of the array is set equal to its index number along
;   the dimension IAXIS.
;
; CALLING SEQUENCE:
;   result = djs_laxisgen( dimens, [ iaxis=iaxis ] )
;
; INPUT:
;   dimens:     Vector of the dimensions for the result.
;               Only up to 3 dimensions can be specified.
;   iaxis:      Axis number to use for indexing RESULT.  The first dimension
;               is axis number 0, the second 1, etc.  Default to 0.
;
; OUTPUTS:
;   result:     Output array
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Written by D. Schlegel, 7 Oct 1997, Durham
;   Modified 12 May 1998 to pass one vector with all dimensions.
;-
;-----------------------------------------------------------------------
function djs_laxisgen, dimens, iaxis=iaxis

   ; Need one parameter
   if N_params() LT 1 then begin
      print, 'Syntax - result = djs_laxisgen( dimens, [iaxis=iaxis] )'
      return, -1
   endif

   if (NOT keyword_set(iaxis)) then iaxis = 0

   ndimen = N_elements(dimens)
   naxis = long(dimens) ; convert to type LONG

   case ndimen of
   1 : $
      begin
         result = lindgen(naxis[0])
      end
   2 : $
      begin
         D1 = naxis[0]
         D2 = naxis[1]
         if (iaxis EQ 0) then result = lindgen(D1,D2) MOD D1 $
         else result = transpose( lindgen(D2,D1) MOD D2 )
      end
   3 : $
      begin
         nnax = N_elements(naxis)
         tax = [naxis[iaxis], naxis(where(indgen(nnax) NE iaxis))]
         result = lindgen( tax[0],tax[1],tax[2] ) MOD naxis[iaxis]
         if (iaxis EQ 0) then pvec = [0] $
         else pvec = [1+indgen(iaxis), 0]
         if (iaxis NE nnax-1) then pvec = [pvec, iaxis+1+indgen(nnax-iaxis-1)]
         result = transpose(result,pvec)
      end
   else : $
      begin
         print, ndimen, ' dimensions not supported'
         result = -1
      end
   endcase

   return, result
end 
;-----------------------------------------------------------------------
