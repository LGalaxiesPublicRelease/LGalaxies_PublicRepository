;=========================================================================
  
pro plot_snap,ps=ps

;=========================================================================
;
;  Example procedure to read in L-galaxies snapshot data
;  and create a luminosity function plot

; Inputs:
;   ps - keyword to cause the plot to be written to a postscript file.
;        (by default the plot will be drawn on the screen)
;
;-------------------------------------------------------------------------

; Define which files you want to read in
;DirName = '../output/DLB07/' & format='DLB07_snap'
DirName = '../output/HT09/' & format='HT09_snap'
FileName = 'SA_z0.00'
;DirName = '/export/virgo/SAM/snaps/DLB07/Z000_63'
;DirName = '/export/virgo/SAM/snaps/HT09/Z000_63'
;FileName = 'SA_z0.00'
ModelName = DirName + FileName
;FirstFile = 0
;LastFile = 511
FirstFile = 5
LastFile = 5
   
; Decide whether to include dust or not
duston = 1

; Define some simulation/SAM parameters.
; This should really be read in with the data but there is no way
; currently of doing this.
MaxTreeFiles = 512              ;Output split into this many files
BoxSize      = 500.0            ;Mpc/h
PartMass = 0.0860657            ;units 1e10 Msun/h
Hubble_h = 0.73
BaryonFrac = 0.17               ;0.18 in power spectrum
Yield=0.03
Recycled=0.43

; Define some plotting parameters
dilute = 4.0 / (LastFile-FirstFile+1) ; for faster plotting
LFbinwidth = 0.25     ; for luminosity function
BARbinwidth = 0.1     ; for baryonic mass function

; Some constants
kev_to_K = 1.1605e7
volume = (BoxSize^3.0) * (Lastfile - Firstfile + 1) / MaxTreeFiles 
log_hubble = alog10(Hubble_h)

;-------------------------------------------------------------------------------

; read in data

read_lgal,ModelName,GalStruct $
         ,FileNr=FileNr,TreeNr=TreeNr,NGalsPerTree=NGalsPerTree $
         ,FirstFile=FirstFile,LastFile=LastFile,format=format ;,/swap_endian

;-------------------------------------------------------------------------------

; Define simple arrays containing useful data from structure array
print, ' Now defining useful quantities...'

Type    = GalStruct.Type
Mvir    = GalStruct.Mvir * 1.0e10
Rvir    = GalStruct.Rvir * 1000.0       ; kpc
Vvir    = GalStruct.Vvir
Vmax    = GalStruct.Vmax
Pos     = GalStruct.Pos

HaloIndex   = GalStruct.HaloIndex
CentralMvir = GalStruct.CentralMvir * 1.0e10

HotGas        = GalStruct.HotGas * 1.0e10
ColdGas       = GalStruct.ColdGas * 1.0e10
StellarMass   = GalStruct.StellarMass * 1.0e10
BulgeMass     = GalStruct.BulgeMass * 1.0e10
BlackHoleMass = GalStruct.BlackHoleMass * 1.0e10
Massage       = GalStruct.MassWeightAge

MetalsColdGas       = GalStruct.MetalsColdGas * 1.0e10
MetalsStellarMass   = GalStruct.MetalsStellarMass * 1.0e10
MetalsBulgeMass     = GalStruct.MetalsBulgeMass * 1.0e10
MetalsHotGas        = GalStruct.MetalsHotGas * 1.0e10

Sfr                = GalStruct.Sfr
SfrBulge           = GalStruct.SfrBulge
DiskRadius         = GalStruct.DiskRadius * 1000.0 ; kpc

;VBIK Magnitudes
if (duston eq 0) then begin
    Mag_B = GalStruct.Mag_B
    Mag_I = GalStruct.Mag_I
    Mag_V = GalStruct.Mag_V
    Mag_K = GalStruct.Mag_K
endif else begin
    Mag_B = GalStruct.Mag_Bdust
    Mag_I = GalStruct.Mag_Idust
    Mag_V = GalStruct.Mag_Vdust
    Mag_K = GalStruct.Mag_Kdust
endelse
Mag_I_unc = GalStruct.Mag_I
Mag_B_unc = GalStruct.Mag_B
Mag_V_unc = GalStruct.Mag_V
Mag_K_unc = GalStruct.Mag_K
Mag_Bb_unc = GalStruct.Mag_Bb
Mag_Bb = GalStruct.Mag_Bb
Mag_Ib = GalStruct.Mag_Ib
Mag_Vb = GalStruct.Mag_Vb
Mag_Kb = GalStruct.Mag_Kb

; Nullify orginal data arrays to save memory
FileNr=0
TreeNr=0
NGalsPerTree=0
GalStruct=0

;-------------------------------------------------------------------------------

; An example of how to create data for a plot, in this case
; luminosities

Ngals=n_elements(Type)

; Calculate luminosities
L_B  = fltarr(Ngals)
L_K  = fltarr(Ngals)
L_bj = fltarr(Ngals)
Mag_bj = fltarr(Ngals)
B_cut = 50.0
V_cut = 50.0
K_cut = 50.0
bj_cut = -16.0                  ; bj magnitude limit in survey
mbj_sun = 5.33                  ; magnitude of sun in bj band
mk_sun = 3.3                    ; magnitude of sun in K band
mb_sun = 5.48                   ; magnitude of sun in B band
L_bjmin = 10.0^((mbj_sun - bj_cut) / 2.5)
L_Bmin = 1.e6
L_Kmin = 1.e6
type_pivot = 0.8                ; type_pivot<0.8 ---> blue
                                ; type_pivot>0.8 ---> red

; Calculate luminosities and bj mag
for i=0L, Ngals-1 do begin
    if (Mag_B(i) lt B_cut) then L_B(i) = 10.0^((mb_sun - Mag_B(i)) / 2.5) $
    else L_B(i) = 0.0
    if (Mag_K(i) lt K_cut) then L_K(i) = 10.0^((mk_sun - Mag_K(i)) / 2.5) $
    else L_K(i) = 0.0
    if (Mag_B(i) lt B_cut and Mag_V(i) lt V_cut) then begin
        Mag_bj(i) = Mag_B(i) - 0.28 * (Mag_B(i) - Mag_V(i))
        L_bj(i) = 10.0^((mbj_sun - Mag_bj(i)) / 2.5)
    endif else begin
        Mag_bj(i) = 99.0
        L_bj(i) = 0.0
    endelse
endfor

;-------------------------------------------------------------------------------

; Create plot

start_plot

; Open file for all plots
if (keyword_set(ps)) then begin
    set_plot, 'PS'
    device, filename = 'plot.ps', /color, xsize = 15, ysize = 10, $
      xoffset=1, $
      yoffset=5
endif
  
xmin=-22.0
xmax=-16.0
ymin = 20.0
ymax = 3000.0
plot, findgen(10), /nodata, /ylog, $
      xrange = [xmin, xmax], xstyle=1, $
      yrange = [ymin, ymax], ystyle=1, $
      xtitle = 'M!dbJ!n', $
      ytitle = 'Number of Galaxies'

mag = Mag_b
bin=0.25
counts = histogram(mag, locations=c,min = -24.0, max = -15.0, binsize = bin)   
oplot, c-5*alog10(Hubble_h), counts/(5.0*0.25*0.7*0.7*0.7),color=2

Mstar = -19.99
alpha = -1.28
phistar = 650.0
M= [-21.875, -21.625, -21.375, -21.125, -20.875, -20.625, -20.375, $
    -20.125, -19.875, -19.625, -19.375, -19.125, -18.875, -18.625, $
    -18.375, -18.125, -17.875, -17.625, -17.375, -17.125, -16.875, $
    -16.625, -16.375, -16.125, -15.875]
xx = 10^(0.4*(Mstar-M))
SF = (2./5.)*alog(10.) * phistar * (xx^(alpha+1.)) * (exp(-xx))
oploterr, M,SF,sqrt(SF),3
;   oploterror, M,SF,sqrt(SF), color = 200, errcolor = 200, psym = 3
M = -(findgen(3000)/100)
xx = 10^(0.4*(Mstar-M))
SF = (2./5.)*alog(10.) * phistar * (xx^(alpha+1.)) * (exp(-xx))
oplot, M, SF,color=4

oplot, [-20.1,-19.6], [65.0,65.0], linestyle=8,color=4
xyouts, -19.5, 60.0, 'De Propris et al. 2002',color=4
oplot, [-20.1,-19.6], [42.0,42.0], linestyle=8,color=2
xyouts, -19.5, 40.0, 'Semi-Analytic Model',color=2

if (keyword_set(ps)) then begin
    device, /close_file
    set_plot,'x'
endif

;------------------------------------------------------------------------------


end
