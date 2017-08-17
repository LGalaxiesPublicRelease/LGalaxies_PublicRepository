;=========================================================================
;
;  Script to read in two sets of L-galaxies snapshot data and compare
;  black hole masses
;
;-------------------------------------------------------------------------

File=5

;DirName = '../../output/BHGR9/'
DirName = './'
FileName = 'SA_z0.00'
ModelName = DirName + FileName
FirstFile = File
LastFile = File
@'../../output/snap_template.pro'

read_lgal,ModelName,Template,Gold $
         ,FirstFile=FirstFile,LastFile=LastFile ;,/swap_endian

;-------------------------------------------------------------------------------

DirName = '../../output/BHGR3/'
FileName = 'SA_z0.00'
ModelName = DirName + FileName
FirstFile = File
LastFile = File
@'../../output/snap_template.pro'

read_lgal,ModelName,Template,Gnew $
         ,FirstFile=FirstFile,LastFile=LastFile ;,/swap_endian

;-------------------------------------------------------------------------------

wset,0
plot,Gold.BlackHoleMass,gnew.BlackHoleMass,psym=1,/xlog,/ylog,xrange=[1e-8,1],yrange=[1e-8,1], $
     xtitle='BHMass/10!u10!nh!u-1!nM!d!mn!n',ytitle='BHMass/10!u10!nh!u-1!nM!d!mn!n'
xyouts,2e-8,0.2,'x-axis: BlackHoleGrowth=0',charsize=1.
xyouts,2e-8,0.05,'y-axis: BlackHoleGrowth=1',charsize=1.
xyouts,2e-8,0.0125,'  BlackHoleAccretionRate=9',charsize=1.
xyouts,2e-8,0.003125,'  BlackHoleSeedMass=1e-7',charsize=1.
oplot,[1e-8,1],[1e-8,1]
hard_jpg,'plot0.jpg'

wset,1
plot,Gold.QuasarAccretionRate,gnew.QuasarAccretionRate,psym=1,/xlog,/ylog,xrange=[1e-4,2],yrange=[1e-5,2], $
     xstyle=1,ystyle=1,xtitle='QAR/M!d!mn!nyr!u-1!n',ytitle='QAR/M!d!mn!nyr!u-1!n'
xyouts,0.00012,0.4,'x-axis: BlackHoleGrowth=0',charsize=1.
xyouts,0.00012,0.2,'y-axis: BlackHoleGrowth=1',charsize=1.
xyouts,0.00012,0.1,'  BlackHoleAccretionRate=9',charsize=1.
xyouts,0.00012,0.05,'  BlackHoleSeedMass=1e-7',charsize=1.
oplot,[1e-4,1],[1e-4,1]
hard_jpg,'plot1.jpg'

wset,2
plot,Gnew.BlackHoleMass,gnew.BlackHoleGas,psym=1,/xlog,/ylog,xrange=[1e-8,1],yrange=[1e-15,1], $
     xtitle='BHMass/10!u10!nh!u-1!nM!d!mn!n',ytitle='BHGas/10!u10!nh!u-1!nM!d!mn!n'
xyouts,2e-8,0.05,'BlackHoleGrowth=1',charsize=1.
xyouts,2e-8,0.005,'BlackHoleAccretionRate=9',charsize=1.
xyouts,2e-8,0.0005,'BlackHoleSeedMass=1e-7',charsize=1.
oplot,[1e-8,1],[1e-8,1]
hard_jpg,'plot2.jpg'

wset,3
plot,Gnew.BlackHoleMass-gnew.quasaraccretionrate*1.27e7*0.73/1e10,gnew.QuasarAccretionRate,psym=1,/xlog,/ylog,xrange=[1e-8,1], $
;plot,Gnew.BlackHoleMass,gnew.QuasarAccretionRate,psym=1,/xlog,/ylog,xrange=[1e-8,1], $
     xtitle='BHMass/10!u10!nh!u-1!nM!d!mn!n',ytitle='QAR/M!d!mn!nyr!u-1!n'
xyouts,2e-8,5,'BlackHoleGrowth=1',charsize=1.
xyouts,2e-8,0.5,'BlackHoleAccretionRate=9',charsize=1.
xyouts,2e-8,0.05,'BlackHoleSeedMass=1e-7',charsize=1.
Hubble_h=0.73
yr=3.16e7
tedd=0.45e9 ;in years
BlackHoleGrowthRate=9.
oplot,[1e-8,1],[1e-8,1]*BlackHoleGrowthRate*1e10/Hubble_h/tedd
hard_jpg,'plot3.jpg'

wset,4
c=3.e8                          ; m/s
yr=3.16e7                       ; s
epsilon=0.1                     ; Accretion (in)efficiency
M_sun=1.98e30                   ; kg 
L_bol_sun=3.846e26              ; kg.m^2/s^3
BoxSize = 500.
Hubble_h = 0.73
MaxTreeFiles = 512.
volume = ((BoxSize/Hubble_h)^3.0) * (Lastfile - Firstfile + 1) / MaxTreeFiles 
binsize=0.2
;Lbol=(Gold.QuasarAccretionRate+Gold.RadioAccretionRate)*epsilon/(1-epsilon)* (c^2/yr)*(M_sun/L_bol_sun)
Lbol=(Gold.QuasarAccretionRate)*epsilon/(1-epsilon)* c^2/yr*(M_sun/L_bol_sun)
counts=histogram(alog10(Lbol),locations=loc,min=8,max=14,binsize=binsize)
plot,loc,alog10(counts/binsize/Volume),psym=1,/nodata,$
     yrange=[-8,-2],ystyle=1, $
     ytitle='Log[dN!dQSO!n/dLog(L!dbol!n)/Mpc!u-3!n]', $
     xtitle='Log(L!Dbol!N/L!Dbol,!Mn!N)'
oplot,loc+0.5*binsize,alog10(counts/binsize/Volume),psym=1
;Lbol=(Gnew.QuasarAccretionRate+Gnew.RadioAccretionRate)*epsilon/(1-epsilon)* (c^2/yr)*(M_sun/L_bol_sun)
Lbol=(Gnew.QuasarAccretionRate)*epsilon/(1-epsilon)* c^2/yr*(M_sun/L_bol_sun)
counts=histogram(alog10(Lbol),locations=loc,min=8,max=14,binsize=binsize)
oplot,loc+0.5*binsize,alog10(counts/binsize/Volume),psym=4,color=2

xyouts,12,-2.6,'BlackHoleGrowth=0'
xyouts,12,-2.9,'BlackHoleGrowth=1',color=2
xyouts,12,-3.2,'BlackHoleSeedMass=1e-7',color=2
xyouts,12,-3.5,'BlackHoleAccretionRate=9',color=2

;------------------------------------------------------------------------------

wset,0
