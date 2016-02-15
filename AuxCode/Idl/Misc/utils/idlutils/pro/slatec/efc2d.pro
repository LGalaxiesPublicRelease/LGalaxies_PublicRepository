;+
; NAME:
;   efc2d
;
; PURPOSE:
;   Calculate a B-spline in the least-squares sense 
;     based on two variables: x which is sorted and spans a large range
;				  where bkpts are required
;  		and 	      y which can be described with a low order
;				  polynomial	
;
; CALLING SEQUENCE:
;   
;   coeff = efc2d(x, y, z, invsig, npoly, nbkptord, fullbkpt)
;
; INPUTS:
;   x          - data x values
;   y          - data y values
;   z          - data z values
;   invsig     - inverse error array of y
;   npoly      - Order of polynomial (as a function of y)
;   nbkptord   - Order of b-splines (4 is cubic)
;   fullbkpt   - Breakpoint vector returned by efc
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
;	This IDL proc is an extension of efcmn
;
; EXAMPLES:
;
;   x = findgen(100)
;   y = randomu(200, 100, /normal)
;   zmodel = 10.0*sin(x/10.0) + y
;   z = zmodel + randomu(100,100,/normal)
;   invsig = fltarr(100) + 1.0
;   fullbkpt = [-3.0,-2.0,-1.0,findgen(11)*10.0,101.0,102.0,103.0]
;   npoly = 2L
;   nbkptord = 4L
;   coeff = efc2d(x, y, z, invsig, npoly, nbkptord, fullbkpt)
;
;   zfit = bvalu2d(x, y, fullbkpt, coeff)
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
function efc2d, x, y, z, invsig, npolyin, nbkptordin, fullbkptin

      nbkptord = long(nbkptordin)
      npoly= long(npolyin)
      nx = n_elements(x)
      fullbkpt = fullbkptin
      nbkpt= n_elements(fullbkpt)
      n = nbkpt - nbkptord
      mdg = (nbkpt+1)*npoly
      mdata = max([nx,nbkpt*npoly])

      nord = npoly*nbkptord
      coeff = fltarr(npoly*n)
      bf = fltarr(nord*nord)
      g = fltarr(mdg, (nord+1))

      xtemp = fltarr(mdata)

      xmin = fullbkpt[nbkptord-1]
      xmax = fullbkpt[n] 
      nordm1 = nbkptord - 1L
      mdein = 1L

      ptemp = sort(x)

      xmin = min([xmin,x])
      xmax = max([xmax,x])

;
;	Make sure nord bkpts on each side lie outside [xmin,xmax]
;
	for i=0,nbkptord-1 do fullbkpt[i] = min([fullbkpt[i],xmin])
	for i=n,nbkpt-1 do fullbkpt[i] = max([fullbkpt[i],xmax])

;
;     Initialize parameters of banded matrix processor, BNDACC( ).
;
      mt = 0L
      ip = 1L
      ir = 1L
      ileft = nbkptord
      intseq = 1L
      xval = float(x[ptemp])
      yval = float(y[ptemp])
      zval = float(z[ptemp])
      invsigval = float(invsig[ptemp])


      for i=0L, nx-1 do begin

        ypoly = yval[i] ^ findgen(npoly)

        if (xval[i] GE fullbkpt[ileft]) then begin
          test = call_external(getenv('IDLUTILS_DIR')+'/lib/libslatec.'+ $
                               idlutils_so_ext(), $
           'bndacc_idl', g, mdg, nord, ip, ir, mt, npoly*(ileft-nbkptord)+1)

	  mt = 0L

	  while (ileft LE n AND xval[i] GE fullbkpt[ileft]) do begin
            ileft = ileft + 1L
          endwhile
        endif

;
;        Obtain B-spline function value.
;
        bf = bsplvn(fullbkpt, nbkptord, xval[i], ileft - 1)
 
        bfall = ypoly # bf
   
        irow = ir + mt
        mt = mt + 1
        g[irow-1,*] = [bfall[*],zval[i]] * invsigval[i]

        if (irow EQ mdg-1) then begin 
           test = call_external(getenv('IDLUTILS_DIR')+'/lib/libslatec.'+ $
                                idlutils_so_ext(), $
             'bndacc_idl', g, mdg, nord, ip, ir, mt, npoly*(ileft-nbkptord)+1)
           mt = 0L
        endif
    endfor
 
    test = call_external(getenv('IDLUTILS_DIR')+'/lib/libslatec.'+ $
                         idlutils_so_ext(), $
             'bndacc_idl', g, mdg, nord, ip, ir, 1L, npoly*(ileft-nbkptord)+1)

    g[ir-1,*] = 0.0
    test = call_external(getenv('IDLUTILS_DIR')+'/lib/libslatec.'+ $
                         idlutils_so_ext(), $
                'bndacc_idl', g, mdg, nord, ip, ir, 1L, npoly*(n+1))

    rnorm = 1.0
    test = call_external(getenv('IDLUTILS_DIR')+'/lib/libslatec.'+ $
                         idlutils_so_ext(), $
                'bndsol_idl', 1L, g, mdg, nord, ip, ir, coeff, npoly*n, rnorm)
   
    return, reform(coeff,npoly,n)
end 
