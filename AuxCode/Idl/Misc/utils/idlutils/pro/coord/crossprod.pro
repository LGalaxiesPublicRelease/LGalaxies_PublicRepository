;+
; NAME:
;   crossprod
;
; PURPOSE:
;   compute cross product of two vectors or arrays of vectors
;
; CALLING SEQUENCE:
;   C=crossprod(A,B)
;
; INPUTS:
;   A, B  - either [3] or [Nvec,3]
;
; OUTPUT:
;   C     - either [3] or [Nvec,3]
;
; COMMENTS:
;   If A and B are BOTH row vectors, return a row vector.  
;    Otherwise, return [Nvec,3]
;
; REVISION HISTORY:
;   2003-May-13  Written by Finkbeiner & Schlegel, Princeton
;----------------------------------------------------------------------
function crossprod, A, B

  NvecA = n_elements(A)/3
  NvecB = n_elements(B)/3
  
  if (reverse(size(A, /dimen)))[0] NE 3 then $
   message, 'Dimensions of A are not [Nvec,3]'
  if (reverse(size(B, /dimen)))[0] NE 3 then $
   message, 'Dimensions of B are not [Nvec,3]'

  if NvecA NE NvecB and (NvecA NE 1 AND NvecB NE 1) then begin 
     message, 'inconsistent dimensions for A and B'
  endif 

  Nvec = NvecA > NvecB

  C = fltarr(Nvec, 3)*(A[0]*B[0]*0)    ; zero vector of type float or double

  if NvecA eq 1 then begin 
     A0 = A[0]
     A1 = A[1]
     A2 = A[2]
  endif else begin 
     A0 = A[*, 0]
     A1 = A[*, 1]
     A2 = A[*, 2]
  endelse

  if NvecB eq 1 then begin 
     B0 = B[0]
     B1 = B[1]
     B2 = B[2]
  endif else begin 
     B0 = B[*, 0]
     B1 = B[*, 1]
     B2 = B[*, 2]
  endelse

  C[*,0] =  A1*B2-A2*B1
  C[*,1] = -(A0*B2-A2*B0)
  C[*,2] =  A0*B1-A1*B0

  if size(A, /n_dimen) eq 1 and size(B, /n_dimen) eq 1 then C = transpose(C)

  return, C
end
