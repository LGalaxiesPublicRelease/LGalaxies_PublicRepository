;+
; NAME:
;   hogg_read_postscript
; PURPOSE:
;   read in a postscript file as an IDL data image
; INPUTS:
;   filename
; OPTIONAL INPUTS:
;   resolution   - in dots per inch; default 75
; KEYWORDS:
;   noantialias  - don't do Hogg's antialiasing trick
; OUTPUTS:
;   image        - [3,N,N] array of byte images for R,G,B planes
; BUGS:
;   Makes DUMB temporary files; I am sure there is a way to make smart ones.
; REVISION HISTORY:
;   2003-01-31  written - Hogg
;-
function hogg_read_postscript, filename,resolution=resolution, $
                               noantialias=noantialias

; set defaults
if NOT keyword_set(resolution) then resolution= 75
if keyword_set(noantialias) then aafactor= 1 else aafactor=2
resstr= strtrim(string(round(resolution*aafactor)),2)
prefix= '/tmp/tmp_hogg_read_postscript'

; bitmap postscript
splog, 'making bitmap'
cmd= 'echo showpage | ghostscript -dNOPAUSE -sDEVICE=ppm -r'+resstr+' -sOutputFile="'+prefix+'.ppm" '+filename
splog, cmd
spawn, cmd

; read bitmap
splog, 'reading bitmap'
image= read_image(prefix+'.ppm')
cmd= '/bin/rm -f '+prefix+'.ppm'
splog, cmd
spawn, cmd

; check rank
dims= size(image,/dimensions)
if n_elements(dims) EQ 2 then begin
    nx= dims[0]
    ny= dims[1]
    timage= image
    image= bytarr(3,(nx+dx),(ny+dy))+255
    for ii=0,2 do image[ii,*,*]= timage
    dims= size(image,/dimensions)
endif

; check dimensions
nx= dims[1]
ny= dims[2]
dx= nx MOD aafactor
dy= ny MOD aafactor
timage= image
image= bytarr(3,(nx+dx),(ny+dy))+255
image[*,0:nx-1,0:ny-1]= timage

; shrink!
if aafactor GT 1 then begin
splog, 'antialiasing image'
    nx= (nx+dx)/aafactor
    ny= (ny+dy)/aafactor
    timage= image
    image= bytarr(3,nx,ny)
    for ii=0,2 do image[ii,*,*]= ((rebin(timage[ii,*,*],1,nx,ny) > 0) < 255)
endif

return, image
end
