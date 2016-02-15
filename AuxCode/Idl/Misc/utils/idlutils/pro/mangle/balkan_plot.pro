;+
; NAME:
;   balkan_plot
; PURPOSE:
;   plot balkans of a set of polygons
; CALLING SEQUENCE:
;   balkan_plot, poly [, /over , /noplot, balkans=, weight=, $
;      poly=, shading=, verts= ]
; INPUTS:
;   poly - [np] mangle polygon structures
; OPTIONAL KEYWORDS:
;   /over - overplot previously existing plot
;   /noplot - don't actually plot
; OPTIONAL OUTPUTS:
;   poly - [np] polygons guessed from the input
;   balkans - [nb] balkanized polygons
;   weight - [nb] weight to add in each balkanized polygon
;   shading - [nb] shading values based on exposure time
;   verts - [nb] structure containing vertices 
;                     .NVERTS - number of vertices
;                     .RA - pointer to array of RA values
;                     .DEC - pointer to array of DEC values
; COMMENTS:
;   If you want better control over how things get plotted, call with
;   the /noplot option and have it return the "verts" and "shading"
;   outputs (or you can set the shading yourself based on the
;   "weight" output). Then you can do something like the following
;   inside your plotting routine:
;
;    loadct, 0
;    for i=0L, n_elements(verts)-1L do $
;      if(allverts[i].nverts gt 0) then $
;      polyfill, *verts[i].ra, *verts[i].dec, color=shading[i], $
;      noclip=0
;
; REVISION HISTORY:
;   2005-Jul-07   as exposure_plot for FITS files MRB (NYU)
;   2005-Dec-07   rewritten for more generality MRB (NYU)
;-
;------------------------------------------------------------------------------
pro balkan_plot, poly, over=over, balkans=balkans, weight=weight, $
                 noplot=noplot, shading=shading, verts=verts, $
                 minside=minside, dangle=dangle

if(NOT keyword_set(weight)) then $
  weight=fltarr(n_elements(poly))+1.

cmd = [ filepath('mrb_balkanize', root_dir=getenv('IDLUTILS_DIR'), $
                 subdir='bin'),'-X', '-ob' ]
spawn, cmd, /noshell, unit=unit
printf,unit,n_elements(poly)
for j=0L, n_elements(poly)-1L do begin
    printf,unit,n_elements(poly)
    printf,unit,lindgen(n_elements(poly))
endfor
write_mangle_polygons,'',poly,unit=unit
read_binary_polygons,'',balkans,unit=unit,/allow_doubles
free_lun,unit

;   c. for each balkan find total weight
xb=vmid(balkans)
bweight=fltarr(n_elements(balkans))
for i=0L, n_elements(balkans)-1L do begin
    for j=0L, n_elements(poly)-1L do begin
        inpoly=is_in_polygon(poly[j], xyz=xb[*,i])
        if(inpoly) then $
          bweight[i]=bweight[i]+weight[j]
    endfor
endfor

if(NOT keyword_set(noplot)) then begin
    verts=replicate({nverts:0L, $
                     ra:ptr_new(), $
                     dec:ptr_new()}, n_elements(balkans))
    for i=0L, n_elements(balkans)-1L do begin
        if(garea(balkans[i]) gt 0.) then begin
            gverts,balkans[i],angle=angle, verts=verts1, edges=edges, ends=ends, $
              dangle=dangle,minside=minside,/plot
            x_to_angles,verts1,ra,theta
            dec=90.-theta
            verts[i].nverts=n_elements(dec)
            verts[i].ra=ptr_new(reform(ra, n_elements(dec)))
            verts[i].dec=ptr_new(reform(dec, n_elements(dec)))
        endif
    endfor
    
    shading=230-byte(bweight/max(bweight)*220.)
    minra=361.
    maxra=-1.
    mindec=91.
    maxdec=-91.
    if(NOT keyword_set(over)) then begin
        for i=0L, n_elements(verts)-1L do begin
            if(verts[i].nverts gt 0) then begin
                minra=min([minra, *(verts[i].ra)])
                maxra=max([maxra, *(verts[i].ra)])
                mindec=min([mindec, *(verts[i].dec)])
                maxdec=max([maxdec, *(verts[i].dec)])
            endif
        endfor
        djs_plot, [0], [0], /nodata, xra=[minra, maxra], yra=[mindec, maxdec]
    endif
    loadct, 0
    for i=0L, n_elements(verts)-1L do $
      if(verts[i].nverts gt 0) then $
      polyfill, *verts[i].ra, *verts[i].dec, color=shading[i], $
      noclip=0
    if(NOT keyword_set(over)) then $
      djs_plot, [0], [0], /nodata, xra=[minra, maxra], yra=[mindec, maxdec], $
      /noerase
endif

end
;------------------------------------------------------------------------------
