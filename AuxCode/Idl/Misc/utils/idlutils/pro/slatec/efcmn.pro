;+
; NAME:
;   efcmn
;
; PURPOSE:
;   Calculate a B-spline in the least-squares sense
;
; CALLING SEQUENCE:
;   
;   coeff = efcmn(x, y, invsig, nord, fullbkpt)
;
; INPUTS:
;   x          - data x values
;   y          - data y values
;   invsig     - inverse error array of y
;   nord       - Order of b-splines (default 4: cubic)
;   fullbkpt       - Breakpoint vector returned by efc
;
; RETURNS:
;   coeff      - B-spline coefficients calculated by efc
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;	This IDL proc mimics efc.f
;
; EXAMPLES:
;
;   x = findgen(100)
;   y = randomu(100,100)
;   invsig = fltarr(100) + 1.0
;   fullbkpt = [-3.0,-2.0,-1.0,findgen(12)*10.0,101.0,102.0,103.0]
;   nord = 4L
;
;   coeff = efcmn(x, y, invsig, nord, fullbkpt)
;
;   xfit = findgen(10)*10.0
;   yfit = slatec_bvalu(xfit, fullbkpt, coeff)
;
;
; PROCEDURES CALLED:
;
;   efc_idl in src/slatec/idlwrapper.c
;         which wraps to efc.o in libslatecidl.so
;
; REVISION HISTORY:
;   10-Mar-2000 Written by Scott Burles, FNAL
;-
;------------------------------------------------------------------------------
function efcmn, x, y, invsig, nordin, fullbkptin

      nord = long(nordin)
      nx = n_elements(x)
      fullbkpt = fullbkptin
      nbkpt= n_elements(fullbkpt)
      n = nbkpt - nord
      mdg = nbkpt + 1
      mdw = nbkpt - nord + 3
      mdata = max([nx,nbkpt])

      coeff = fltarr(n)
      bf = fltarr(nord*nord)
      g = fltarr(mdg, (nord+1))

      xtemp = fltarr(mdata)

      xmin = fullbkpt[nord-1]
      xmax = fullbkpt[n] 
      nordm1 = nord - 1
      nordp1 = nord + 1
      mdein = 1L

      ptemp = sort(x)

      xmin = min([xmin,x])
      xmax = max([xmax,x])

;
;	Make sure nord bkpts on each side lie outside [xmin,xmax]
;
	for i=0,nord-1 do fullbkpt[i] = min([fullbkpt[i],xmin])
	for i=n,nbkpt-1 do fullbkpt[i] = max([fullbkpt[i],xmax])

;
;     Initialize parameters of banded matrix processor, BNDACC( ).
;
      mt = 0L
      ip = 1L
      ir = 1L
      ileft = nord
      intseq = 1L
      xval = float(x[ptemp])
      yval = float(y[ptemp])
      invsigval = float(invsig[ptemp])

 
      for i=0L, nx-1 do begin

        if (xval[i] GE fullbkpt[ileft]) then begin
          test = call_external(getenv('IDLUTILS_DIR')+'/lib/libslatec.'+ $
                               idlutils_so_ext(), $
           'bndacc_idl', g, mdg, nord, ip, ir, mt, ileft-nordm1)

	  mt = 0L

	  while (ileft LE n AND xval[i] GE fullbkpt[ileft]) do begin
            ileft = ileft + 1L
          endwhile
        endif

;
;        Obtain B-spline function value.
;
        bf = bsplvn(fullbkpt, nord, xval[i], ileft - 1)

        irow = ir + mt
        mt = mt + 1
        g[irow-1,*] = [bf[*],yval[i]] * invsigval[i]

        if (irow EQ mdg-1) then begin 
           test = call_external(getenv('IDLUTILS_DIR')+'/lib/libslatec.'+ $
                                idlutils_so_ext(), $
                'bndacc_idl', g, mdg, nord, ip, ir, mt, ileft-nordm1)
           mt = 0L
        endif
    endfor
 
    test = call_external(getenv('IDLUTILS_DIR')+'/lib/libslatec.'+ $
                         idlutils_so_ext(), $
                'bndacc_idl', g, mdg, nord, ip, ir, 1L, ileft-nordm1)

    g[ir-1,*] = 0.0
    test = call_external(getenv('IDLUTILS_DIR')+'/lib/libslatec.'+ $
                         idlutils_so_ext(), $
                'bndacc_idl', g, mdg, nord, ip, ir, 1L, n+1)


    rnorm = 1.0
    test = call_external(getenv('IDLUTILS_DIR')+'/lib/libslatec.'+ $
                         idlutils_so_ext(), $
                'bndsol_idl', 1L, g, mdg, nord, ip, ir, coeff, n, rnorm)
   
    return, coeff
end 
