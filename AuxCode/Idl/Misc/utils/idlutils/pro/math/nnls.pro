;+
; NAME:
;   nnls
; PURPOSE:
;   non-negative least-square fitting routine
; COMMENTS:
;   See documentation in $IDLUTILS_DIR/src/math/nnls.f.
;-
pro nnls, a,mda,m,n,b,x,resnorm,w,zz,indx,mde
  a= reform(float(a),mda,n)
  mda= fix(mda)
  m= fix(m)
  n= fix(n)
  b= float(b)
  x= fltarr(n)
  resnorm= float(1.0)
  w= fltarr(n)
  zz= fltarr(m)
  indx= intarr(n)
  mde= fix(1)
  soname = filepath('libmath.'+idlutils_so_ext(), $
    root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')
  retval= call_external(soname,'idl_nnls', $
    a,mda,m,n,b,x,resnorm,w,zz,indx,mde)
end
