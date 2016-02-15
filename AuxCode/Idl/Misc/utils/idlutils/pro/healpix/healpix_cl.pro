;+
; 
; A rapid power spectrum computation
; code. Doesn't try to do *anything*
; fancy or clever --- just straight
; spherical transforms and averaging.
;
; cl =  healpix_cl(data, [lmax=lmax], cross=cross_data)
;
; data === healpix array in ring format (can be multiple maps)
; lmax=lmax === lmax sent to healpix2alm
; cross === If you want to compute a cross power spectrum, 
;           then send the second healpix map here. Must be the
;           same resolution
;
; Nikhil Padmanabhan, Princeton 
; August 25, 2003
;-

function  healpix_cl, data, lmax=lmax, cross=cross

    if (NOT keyword_set(data)) then $ 
	message, 'ERROR: Must send in data'

    ; Determine the number of maps
    nmap = size(data, /n_dim) eq 1 ? 1 : (size(data, /dimens))[1]

    alm = healpix2alm(data, lmax=lmax)
    if (n_elements(cross) NE n_elements(data)) then begin
        alm = double(alm*conj(alm))  
    endif else begin
        nmap2 = size(cross, /n_dim) eq 1 ? 1 : (size(cross, /dimens))[1]
        if (nmap2 NE nmap) then $
          message,'ERROR : dimensions do not agree'
        blm = healpix2alm(cross, lmax=lmax)
        alm = double(alm*conj(blm))
    endelse

    ; Now loop over the different maps
    for i = 0L, nmap-1L do begin
       cl = total(alm[*,*,i], 2, /double)+total(alm[*,1:*,i],2,/double)
       cl = cl/(2.d0*dindgen(n_elements(cl)) + 1.d0)
       if (i EQ 0) then $
	  clarr = dblarr(n_elements(cl), nmap)
       clarr[*,i] = cl
    endfor

    return, clarr
 
; All done
end

