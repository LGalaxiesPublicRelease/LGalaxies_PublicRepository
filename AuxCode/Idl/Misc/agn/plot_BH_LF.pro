; Define which files you want to read in

;FileName = 'SA_z0.99'
FileName = 'SA_z0.00'

;DirName = '/export/virgo/SAM/snaps/DLB07plusM05/'
;DirName = '/export/virgo/SAM/snaps/HT09plusM05/'
DirName = '~/eclipse/workspace/branch_metals/output/'
; Galaxy structure template
;@'/export/virgo/SAM/snaps/DLB07_template.pro'
;@'/export/virgo/SAM/snaps/HT09_template.pro'
@'~/eclipse/workspace/branch_metals/output/snap_template.pro'

ModelName = DirName + FileName

;FirstFile = 0
;LastFile = 100
;LastFile = 511
FirstFile = 5
LastFile = 5
   
; Define some simulation/SAM parameters.
; This should really be read in with the data but there is no way
; currently of doing this.
MaxTreeFiles = 512              ;Output split into this many files
BoxSize      = 500.0            ;Mpc/h
PartMass = 0.0860657            ;units 1e10 Msun/h
Hubble_h = 0.73
BaryonFrac = 0.17               ;0.18 in power spectrum

; Some constants
kev_to_K = 1.1605e7
volume = ((BoxSize/Hubble_h)^3.0) * (Lastfile - Firstfile + 1) / MaxTreeFiles 
log_hubble = alog10(Hubble_h)

;-------------------------------------------------------------------------------

; read in data

if n_elements(GalStruct) eq 0 then $
  read_lgal,ModelName,Template,GalStruct $
  ,FirstFile=FirstFile,LastFile=LastFile ;,/swap_endian

;-------------------------------------------------------------------------------

; Define simple arrays containing useful data from structure array
print, ' Now defining useful quantities...'

Type    = GalStruct.Type
BlackHoleMass = GalStruct.BlackHoleMass
QuasarAccretion=GalStruct.QuasarAccretionRate
RadioAccretion= GalStruct.RadioAccretionRate

;-------------------------------------------------------------------------------
; Parameters

Ngals=n_elements(Type)

c=3.e8                          ; m/s
c_km=3.e5                       ; km/s
yr=3.16e7                       ; s
t_edd=0.45*yr*1.e9              ; s
epsilon=0.1                     ; Accretion (in)efficiency
M_sun=1.98e30                   ; kg 
L_bol_sun=3.846e26              ; kg.m^2/s^3

iseed=1234321

;-----------------------------------------------------------------------------
;Variables

MBH_Msun=BlackHoleMass * 1.0e10 / Hubble_h
Mdot_Msun=(QuasarAccretion+RadioAccretion) / yr

L_edd_Msun=MBH_Msun/t_edd*c^2

Lbol=fltarr(Ngals)
Lbol1=fltarr(Ngals)
Lbol2=fltarr(Ngals)
Lbol3=fltarr(Ngals)

; This is the accretion rate compared to the Eddington rate
f=epsilon/(1-epsilon)*Mdot_Msun*t_edd/MBH_Msun
print,'Max f =',max(f)
n=randomu(iseed,Ngals) ; Have checked and does look random

Lbol=Mdot_Msun  * epsilon/(1-epsilon) * c^2 *(M_sun/L_bol_sun)
index1=where(n lt f)
Lbol1[index1]=L_edd_Msun[index1]* (M_sun/L_bol_sun)
index2=where(n lt 10*f)
Lbol2[index2]=0.1*L_edd_Msun[index2]* (M_sun/L_bol_sun)
index3=where(n lt 100*f)
Lbol3[index3]=0.01*L_edd_Msun[index3]* (M_sun/L_bol_sun)

LMBH=alog10(MBH_Msun)

Log_Lbol=alog10(Lbol)
Log_Lbol1=alog10(Lbol1)
Log_Lbol2=alog10(Lbol2)
Log_Lbol3=alog10(Lbol3)

;-------------------------------------------------------------------------------
;Create plot - Bolometric Luminosity Function

start_plot

window,win,xsize=700,ysize=600

binsize=0.2
counts=histogram(Log_Lbol,locations=loc,min=8,max=14,binsize=binsize)
counts1=histogram(Log_Lbol1,locations=loc,min=8,max=14,binsize=binsize)
counts2=histogram(Log_Lbol2,locations=loc,min=8,max=14,binsize=binsize)
counts3=histogram(Log_Lbol3,locations=loc,min=8,max=14,binsize=binsize)

plot,loc,alog10(counts/binsize/Volume),psym=1,/nodata,$
     yrange=[-8,-2],ystyle=1, $
     ytitle='Log[dN!dQSO!n/dLog(L!dbol!n)/Mpc!u-3!n]', $
     xtitle='Log(L!Dbol!N/L!Dbol,!Mn!N)'
oplot,loc+0.5*binsize,alog10(counts/binsize/Volume),psym=1
oplot,loc+0.5*binsize,alog10(counts1/binsize/Volume),psym=4,color=2
oplot,loc+0.5*binsize,alog10(counts2/binsize/Volume),psym=4,color=3
oplot,loc+0.5*binsize,alog10(counts3/binsize/Volume),psym=4,color=4

;xyouts,13,-3.,'f!dEdd!n='+string(f_edd,format='(f4.2)')
xyouts,12,-2.6,'Raw'
xyouts,12,-2.9,'f!dEdd!n=1',color=2
xyouts,12,-3.2,'f!dEdd!n=0.1',color=3
xyouts,12,-3.5,'f!dEdd!n=0.01',color=4

;------------------------------------------------------------------------------

