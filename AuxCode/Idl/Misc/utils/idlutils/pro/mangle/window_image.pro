;+
; NAME:
;   window_image
; PURPOSE:
;   Given a set of mangle polygons, create a pixelized image of the 
;   window defined by the polygons. 
; CALLING SEQUENCE:
;   image=window_image(polygons [, nx=, ny=, /random])
; INPUTS:
;   polygons - list of polygons
; OPTIONAL INPUTS:
;   nx - number of ra pixels
;   ny - number of dec pixels
;   /random - fill image by assigning random points
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   30-Nov-2002  Written by MRB (NYU)
;-
;------------------------------------------------------------------------------
function window_image,polygons,nx=nx,ny=ny,random=random,nrandom=nrandom

if(n_elements(nx) eq 0) then nx=100L
if(n_elements(ny) eq 0) then ny=100L
if(n_elements(nrandom) eq 0) then nrandom=10000000L

xpos1=(-180.D)+(dindgen(nx)+0.5D)*360.D/double(nx)
ypos1=((-90.D)+(dindgen(ny)+0.5D)*180.D/double(ny))

if(NOT keyword_set(random)) then begin
    xpos=xpos1#replicate(1.,ny)
    ypos=transpose(ypos1#replicate(1.,nx))
    usexyz=angles_to_x(xpos,(90.D)-ypos)
    image=reform(is_in_window(usexyz,polygons),nx,ny)
endif else begin
    image=dblarr(nx,ny)
    write_mangle_polygons,'tmp_poly.ply',polygons
    for i=0L, nrandom/1000000L do begin 
        cmdstr='ransack -c'+strtrim(string(i),2)+' -r'+ $
          strtrim(string(1000000L),2)+' tmp_poly.ply tmp_wi_radec.dat'
        splog,cmdstr
;        spawn,cmdstr
        spawn, ['ransack', '-c'+strtrim(string(i),2), ' -r1000000', $
         'tmp_poly.ply', 'tmp_wi_radec.dat'], /noshell
        openr,unit,'tmp_wi_radec.dat',/get_lun
        dummy=''
        readf,unit,dummy
        inarr=dblarr(3,1000000L)
        readf,unit,inarr
        free_lun,unit
        inx=long(inarr[0,*]*double(nx)/360.)
        iny=long((inarr[1,*]+90.)*double(ny)/180.)
        image[inx,iny]=1L
    endfor
endelse 

return,image

end
