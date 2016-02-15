pro hogg_image_overlay_grid, hdr,overlay,factor=factor, $
                             ra_arrow=ra_arrow,dec_arrow=dec_arrow, $
                             name_arrow=name_arrow,length_arrow=length_arrow, $
                             _EXTRA=KeywordsForGrid
prefix= 'tmp_hogg_image_overlay_grid'
naxis1= round(sxpar(hdr,'NAXIS1'))
naxis2= round(sxpar(hdr,'NAXIS2'))
bangp= !P
bangx= !X
bangy= !Y
set_plot, 'PS'
xsize= 7.5
ysize= xsize*float(naxis2)/float(naxis1)
device, file=prefix+'.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color, bits=8
!P.MULTI= [0,1,1]
!P.POSITION= [0.,0.,1.,1.]
!X.MARGIN= [0,0]
!X.OMARGIN= [0,0]
!Y.MARGIN= !X.MARGIN
!Y.OMARGIN= !Y.OMARGIN
nw_overlay_range, naxis1,naxis2,xrange,yrange
xstyle= 5
ystyle= 5
plot, [0],[0],/nodata, $
  xstyle=xstyle,xrange=xrange, $
  ystyle=ystyle,yrange=yrange
nw_ad_grid, hdr,_EXTRA=KeywordsForGrid
if (keyword_set(ra_arrow)) then begin
    hogg_directions, sxpar(hdr,'CRVAL1'),sxpar(hdr,'CRVAL2'), $
      ra_arrow,dec_arrow,name_arrow,hdr,length=length_arrow
endif
device, /close
!P= bangp
!X= bangx
!Y= bangy
overlay1= 1.-hogg_image_overlay(prefix+'.ps',naxis1,naxis2,factor=factor)
if keyword_set(overlay) then overlay= overlay+overlay1 $
  else overlay= overlay1
return
end
