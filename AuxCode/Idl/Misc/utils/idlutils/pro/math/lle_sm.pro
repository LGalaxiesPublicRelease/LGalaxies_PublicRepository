;+
; NAME:
;	  lle_sm
; PURPOSE: (one line)
;   create sparse matrix from full matrix for LLE routines
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
;   sparse_matrix= lle_sm(full_matrix)
; INPUTS:
;   full_matrix - complete NxM matrix 
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;   sparse_matrix - full_matrix transformed into above format
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
;   always at double precision
; PROCEDURE:
; MODIFICATION HISTORY:
;	  Blanton 2003-05-26 
;-
function lle_sm, inarr

indx=where(inarr ne 0.,count)

outsm={n:0L, $
       m:0L, $
       vals:dblarr(count), $
       nindx:lonarr(count), $
       mindx:lonarr(count)}

ndim=(size(inarr))[0]
if(ndim eq 1) then begin
    outsm.n=1
    outsm.m=(size(inarr,/dimens))[0]
endif else begin
    outsm.n=(size(inarr,/dimens))[0]
    outsm.m=(size(inarr,/dimens))[1]
endelse
outsm.vals=inarr[indx]
outsm.mindx=indx/outsm.n
outsm.nindx=indx mod outsm.n

return, outsm

end
