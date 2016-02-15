;+
; NAME:
;   mmsparse
; PURPOSE:
;   multiply a regular matrix by a sparse matrix
; CALLING SEQUENCE:
;   mmsparse, cc, bb, aasparse 
; INPUTS:
;   bb - sparse matrix struct (see below) [nx, nz] of inverse var
;   aasparse - sparse matrix struct (see below) [nx, ny] 
; OUTPUTS:
;   cc - [nz, ny] result
; COMMENTS:
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
pro mmsparse, cc, bb, aasparse

nx=aasparse.nx
ny=aasparse.ny
nz=(size(bb,/dim))[1]

cc=fltarr(nz,ny)

if(NOT keyword_set(soname)) then $
  soname=filepath('libmath.'+idlutils_so_ext(), $
                  root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')
retval=call_external(soname, 'idl_mmsparse', float(cc), float(bb), $
                     long(nx), long(ny), long(nz), float(aasparse.val), $
                     long(aasparse.x), long(aasparse.rowstart), $
                     long(aasparse.nxrow))

end
