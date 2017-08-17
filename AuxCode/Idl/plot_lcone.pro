;=========================================================================
  
pro plot_lcone,ps=ps

;=========================================================================
;
;  Example procedure to read in L-galaxies light-cone data
;  and create a mass function plot

; Inputs:
;   ps - keyword to cause the plot to be written to a postscript file.
;        (by default the plot will be drawn on the screen)
;
;-------------------------------------------------------------------------

compile_opt idl2

; Define which files you want to read in
;DirName = '../output/DLB07/' & format='DLB07_lcone'
DirName = '../output/HT09/' & format='HT09_lcone'
;DirName = '/export/virgo/SAM/lightcones/4times3/DLB07/' & format='DLB07_lcone'
;DirName = '/export/virgo/SAM/lightcones/4times3/HT09/' & format='HT09_lcone'
FileName = 'SA_z'
ModelName = DirName + FileName
FirstFile = 0
LastFile = 9
;LastFile = 511
KitzFileName='../input/kitz_mock5.dat'
;KitzFileName='/export/virgo/SAM/lightcones/kitz_mock5.dat'
   
;Nmock and Mmock are the number of replications along the x and y axis
Nmock=4
Mmock=3
   
;Box is the size of the simulation box in units of Mpc / h
BoxSize      = 500.0            ;Mpc/h
PartMass = 0.0860657            ;units 1e10 Msun/h
Hubble_h = 0.73
BaryonFrac = 0.17               ;0.18 in power spectrum
; Total volume in light-cone
TotalVolume=(500.0*500.0*500.0*hubble_h*hubble_h*hubble_h)/3.0

;-------------------------------------------------------------------------------

; read in data

read_lgal,ModelName,GalStruct $
         ,FileNr=FileNr,TreeNr=TreeNr,NGalsPerTree=NGalsPerTree $
         ,FirstFile=FirstFile,LastFile=LastFile,format=format ;,/swap_endian

;-------------------------------------------------------------------------------

print, ' Now defining useful quantities...'

Type    = GalStruct.Type
Mvir    = GalStruct.Mvir * 1.0e10 
Rvir    = GalStruct.Rvir * 1000.0 
Vvir    = GalStruct.Vvir
Vmax    = GalStruct.Vmax
x       = GalStruct.Pos[0]        
y       = GalStruct.Pos[1]
z       = GalStruct.Pos[2]
mockx   = GalStruct.MockPos[0]    
mocky   = GalStruct.MockPos[1]
mockz   = GalStruct.MockPos[2]
redshift= GalStruct.redshift
snapnum = GalStruct.SnapNum 
ra      = GalStruct.ra
dec     = GalStruct.dec
CentralMvir   = GalStruct.CentralMvir * 1.0e10
HotGas        = GalStruct.HotGas * 1.0e10
ColdGas       = GalStruct.ColdGas * 1.0e10
StellarMass   = GalStruct.StellarMass * 1.0e10
Mass          = Alog10(StellarMass) 
BulgeMass     = GalStruct.BulgeMass * 1.0e10
BlackHoleMass = GalStruct.BlackHoleMass * 1.0e10  
Sfr                = GalStruct.Sfr
SfrBulge           = GalStruct.SfrBulge
DiskRadius         = GalStruct.DiskRadius * 1000.0 
obsMagu = GalStruct.obsMagudust
obsMagg = GalStruct.obsMaggdust
obsMagr = GalStruct.obsMagrdust
obsMagi = GalStruct.obsMagidust
obsMagz = GalStruct.obsMagzdust
Mag_B = GalStruct.Mag_Bdust
Mag_I = GalStruct.Mag_Idust
Mag_V = GalStruct.Mag_Vdust
Mag_K = GalStruct.Mag_Kdust     

;----------------------------------------------------------------------------

; Example of how to create data for a plot, in this case the stellar
; mass function at a comoving distance of 5000-5500 Mpc/h

dist=sqrt(mockx*mockx+mocky*mocky+mockz*mockz)
sel_my=where(dist gt 5000.0 and dist lt 5500.0 and Alog10(StellarMass) gt 9.0)
mass1=mass[sel_my]

; I have no idea what the following is supposed to do.  It seems to be
; the volume in different distance slices, but why do those distances
; not lie in the range 5000-5500?  Oh well, it is only an example plot
a=(7534.-6849.)/10.0
distance=[6849.,6849.+a,6849.+2*a,6849.+3*a,6849.+4*a,6849.+5*a,$
          6849.+6*a,6849.+7*a,6849.+8*a,6849.+9*a,6849.+10.*a]
volume1=0.
for i=0,9 do $
   volume1=volume1+(tan(1./(Nmock*Nmock*Mmock*2.))*2.*distance[i])*$
           (tan(1./(Mmock*Nmock*Mmock*2.))*2.*distance[i])*$
           (distance[i+1]-distance[i])   
volume2=0.
for i=0,9 do $
   volume2=volume2+(tan(1./(Nmock*Nmock*Mmock*2.))*2.*distance[i+1])*$
           (tan(1./(Mmock*Nmock*Mmock*2.))*2.*distance[i+1])*$
           (distance[i+1]-distance[i])   
volume=(volume1+volume2)/2.

; Kitzbichler data for comparison
openr,3,KitzFileName
readf, 3, ndata
kitz_ra=fltarr(ndata)
kitz_dec=fltarr(ndata)
kitz_dis=fltarr(ndata)
kitz_z=fltarr(ndata)
kitz_mass=fltarr(ndata)
kitz_magk=fltarr(ndata)

loop=floor(ndata/40) 
j=0L
for k=0,39 do begin
   for i=0,loop-1 do begin         
      readf, 3, a, b, c, d, e, f
      kitz_ra[i+j]=a
      kitz_dec[i+j]=b
      kitz_dis[i+j]=c
      kitz_z[i+j]=d
      kitz_mass[i+j]=e  
      kitz_magk[i+j]=f
   endfor
   j=j+i
endfor
close,3

volume1=0.
for i=0,9 do $
   volume1=volume1+(tan(0.0244/2.)*2.*distance[i])^2*$
           (distance[i+1]-distance[i])   
volume2=0.
for i=0,9 do $
   volume2=volume2+(tan(0.0244/2.)*2.*distance[i+1])^2*$
           (distance[i+1]-distance[i])   
volume_kitz=(volume1+volume2)/2.
sel_kitz=where(kitz_dis gt 5000.0 and kitz_dis lt 5500.0)

;-------------------------------------------------------------------------------

; Create plot

start_plot

; Open file for all plots
if (keyword_set(ps)) then begin
    set_plot, 'PS'
    device, filename = 'mock.ps', /color, xsize = 15, ysize = 10, $
      xoffset=1, $
      yoffset=5
endif
  
xmin=11.5
xmax=9.0
ymin=0.000001
ymax=0.1

plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
      /ylog,xstyle = 1, ystyle = 1, xtitle = 'log!D10!N(M/M!D!9ng!3!N)', $
      ytitle = '!4U!3 (Mpc!U-3!N(log!D10!N(M!U-1!N))'
bin=0.1
hist=(histogram(mass1,locations=c,min=9.0,max=12.0,binsize=bin))
oplot,c+alog10(Hubble_h)+(bin/2.0),hist/(volume*bin), color=2

bin=0.1
hist=histogram(Alog10(kitz_mass[sel_kitz]*1e10),locations=c,$
               min=9.0,max=12.0,binsize=bin)
oplot,c+alog10(Hubble_h)+(bin/2.0),hist/(volume_kitz*bin), color=2,linestyle=2
oplot, [10.445,10.5], [0.0000046,0.0000046], linestyle=8,color=2
oplot, [10.445,10.5],  [0.00001,0.00001], linestyle=2,color=2
xyouts, 10.4, 0.000004, 'Our Light Cone'
xyouts, 10.4, 0.0000089, 'Kitzbichler & White'

if (keyword_set(ps)) then begin
    device, /close_file
    set_plot,'x'
endif
   
end
