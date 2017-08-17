@read_script
plot,(g.bulgemass+g.diskmass+g.icm),g.burstmass,psym=1,/xlog,/ylog,xrange=[0.003,30],yrange=[0.001,1]
i0=where(g.type eq 0)
oplot,(g[i0].bulgemass+g[i0].diskmass+g[i0].icm),g[i0].burstmass,psym=1,color=2
i1=where(g.type eq 1)
oplot,(g[i1].bulgemass+g[i1].diskmass+g[i1].icm),g[i1].burstmass,psym=1,color=4
i2=where(g.type eq 2)
oplot,(g[i2].bulgemass+g[i2].diskmass+g[i2].icm),g[i2].burstmass,psym=1,color=3
oplot,[1e-4,100],[1e-4,100]
xyouts, 0.01, 0.8, 'Type 0', color=2
xyouts, 0.01, 0.4, 'Type 1', color=4
xyouts, 0.01, 0.2, 'Type 2', color=3

print,max(g.bulgemass,idmax)
print,g[idmax].sfh_burstmass
 
