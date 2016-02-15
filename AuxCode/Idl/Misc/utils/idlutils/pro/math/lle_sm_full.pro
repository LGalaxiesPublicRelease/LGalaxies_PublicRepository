;+
; NAME:
;	  lle_sm_full
; PURPOSE: (one line)
;   create full matrix from full matrix for LLE routines
; DESCRIPTION:
;   rather crappy sparse matrix format for NxM matrix:
;      .N - number of rows
;      .M - number of columns
;      .VALS - each nonzero value
;      .NINDX - row of each nonzero value 
;      .MINDX - column of each nonzero value 
;   but it handles nonsquare matrices
; CATEGORY:
;       Numerical
; CALLING SEQUENCE:
;   full_matrix= lle_sm_full(sparse_matrix)
; INPUTS:
;   sparse_matrix - matrix in above format
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;   full_matrix - complete NxM matrix 
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
;   always at double precision
; PROCEDURE:
; MODIFICATION HISTORY:
;	  Blanton 2003-05-26 
;-
function lle_sm_full, smarr

arr=dblarr(smarr.n,smarr.m)
arr[smarr.nindx,smarr.mindx]=smarr.vals
return,arr

end
