pro nw_overlay_range, NX,NY,xrange,yrange
xrange = [-.5,float(NX)-.5]-[0.000,0.001]*float(NX) ; empirical HACK
yrange = [-.5,float(NY)-.5]+[0.001,0.000]*float(NY) ; empirical HACK
return
end
