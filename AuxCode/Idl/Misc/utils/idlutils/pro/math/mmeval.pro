;+
; NAME:
;   mmeval
; PURPOSE:
;   evaluate a matrix multiply sparsely
; CALLING SEQUENCE:
;   mmeval, cc, bb, aa
; INPUTS:
;   cc - sparse matrix struct (see below) [nz,ny]
;   bb - [nx, nz] matrix 
;   aa - [nx, ny] matrix 
; OUTPUTS:
;   cc - sparse matrix struct (see below) [nz,ny]
; COMMENTS:
;   The matrix multiply of bb.aa is evaluated specified by the sparse
;   matrix structure of cc.
;   The sparse matrix structure referred to above is:
;       .VAL[NVAL]      - actual values in matrix
;       .X[NVAL]        - columns for each value in matrix
;       .NX             - number of columns
;       .NY             - number of rows
;       .ROWSTART[NY]   - starting position of each row in VAL, X
;       .NXROW[NY]      - number of columns in each now 
;   This code is called by nmf_sparse.
; REVISION HISTORY:
;   2005-Feb-5  Written by Mike Blanton, NYU
;               Adapted from Matlab code of Sam Roweis
;
;----------------------------------------------------------------------
pro mmeval, cc, bb, aa

nx=(size(bb,/dim))[0]
ny=n_elements(aa)/nx 
nz=(size(bb,/dim))[1]

if(NOT keyword_set(soname)) then $
  soname=filepath('libmath.'+idlutils_so_ext(), $
                  root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')
val=cc.val
retval=call_external(soname, 'idl_mmeval', float(val), float(bb), $
                     float(aa), long(nx), long(ny), long(nz), $
                     long(cc.x), long(cc.rowstart), $
                     long(cc.nxrow))
cc.val=val

end
