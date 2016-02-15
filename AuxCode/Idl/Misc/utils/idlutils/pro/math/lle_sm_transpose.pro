;+
; NAME:
;	  lle_sm_transpose
; PURPOSE: (one line)
;   transpose a sparse matrix for the LLE routines
; DESCRIPTION:
;   transposes the rather crappy sparse matrix format for NxM matrix:
;      .N - number of rows
;      .M - number of columns
;      .VALS - each nonzero value
;      .NINDX - row of each nonzero value 
;      .MINDX - column of each nonzero value 
;   but it handles nonsquare matrices
; CATEGORY:
;       Numerical
; CALLING SEQUENCE:
;   mattrans= lle_sm_transpose(mat)
; INPUTS:
;   mat - matrix in above format
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;   mattrans - transposed matrix in above format
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; MODIFICATION HISTORY:
;	  Blanton 2003-05-26 
;-
function lle_sm_transpose, a

at=a
at.n=a.m
at.m=a.n
msort=sort(a.nindx*a.m+a.mindx)
at.mindx=a.nindx[msort]
at.nindx=a.mindx[msort]
at.vals=a.vals[msort]

return, at

end
