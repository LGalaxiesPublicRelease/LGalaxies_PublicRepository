window,0
glog=alog10(Mdot_Msun*t_edd/(MBH_Msun-Mdot_Msun*263431803*3.16e7))
glog=glog[where(finite(glog))]
hist_plot,glog,color=2,binsize=0.1
flog=alog10(Mdot_Msun*t_edd/MBH_Msun)
flog=flog[where(finite(flog))]
hist_plot,flog,binsize=0.1,/overplot

window,1
plot,Mdot_Msun*3.16e7,Mdot_Msun*t_edd/MBH_Msun,psym=1,/xlog,xrange=[1e-3,10],/ylog,xtitle='Mdot/M!d!Mn!nyr!u-1!n',ytitle='Mdot/Mdot!dEdd!n'
