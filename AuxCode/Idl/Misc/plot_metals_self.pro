;=========================================================================
;
;  plot_metals_self
;
;  Example script to read in L-galaxies snapshot data
;  with a star-formation history for the galaxies
;  and produce a plot of the fraction of hot gas metals that
;  come from self.
;  Assumes that all the galaxies are from a snap (ie all at the same time).

;  Parameters:
;   ps - keyword to cause the plot to be written to a postscript file.
;        (by default the plot will be drawn on the screen)
    PS = 0b
;   SNAP - which snapshot is it (affects axes)
    SNAP = 63
;    SNAP=41
;   Metals flag to read in and plot appropriate data
    METALS=0b
;   Whether want plot vs 0) StellarMass or 1) Mvir
    CHOICE=1b
;
;-------------------------------------------------------------------------

; Define which files you want to read in
DirName = '../output/'
@'../output/snap_msg_template.pro'
if SNAP eq 63 then FileName = 'SA_z0.00'
if SNAP eq 41 then FileName = 'SA_z0.99'
ModelName = DirName + FileName
FirstFile = 287
LastFile = 287
;FirstFile = 5
;LastFile = 5

;-------------------------------------------------------------------------------

; read in data

read_lgal,ModelName,Template,GalStruct $
         ,FileNr=FileNr,TreeNr=TreeNr,NGalsPerTree=NGalsPerTree $
         ,FirstFile=FirstFile,LastFile=LastFile;,/swap_endian

;-------------------------------------------------------------------------------

; Define simple arrays containing useful data from structure array
print, ' Now defining useful quantities...'

StellarMass = GalStruct.StellarMass * 1.0e10
Mvir = GalStruct.Mvir * 1.0e10
if METALS then $
   MetalsHotGas = GalStruct.MetalsHotGas[1] * 1.0e10 $
else $
   MetalsHotGas = GalStruct.MetalsHotGas * 1.0e10
if METALS then $
   MetalsHotGasSelf = GalStruct.MetalsHotGasSelf[1] * 1.0e10 $
else $
   MetalsHotGasSelf = GalStruct.MetalsHotGasSelf * 1.0e10
if METALS then Metals1a = GalStruct.MetalsStellarMass[0] * 1.0e10
if METALS then Metals2 = GalStruct.MetalsStellarMass[1] * 1.0e10
if METALS then MetalsAGB = GalStruct.MetalsStellarMass[2] * 1.0e10
 
; Nullify orginal data arrays to save memory
FileNr=0
TreeNr=0
NGalsPerTree=0
;GalStruct=0

;-------------------------------------------------------------------------------
  
; Create plot

start_plot

;-------------------------------------------------------------------------------

; Create plots

start_plot
symbols,1,0.4

; Open file for PS plot
if PS then set_plot, 'PS' else window,0
if PS then device, filename = 'metals_self.ps', /color, $
                   xsize = 15, ysize = 10, $
                   xoffset=1,  yoffset=5

index=where(MetalsHotGas gt 1.)
if CHOICE eq 1 then index=where(MetalsHotGas gt 1. and Mvir gt 1.)
x=StellarMass[index]
if CHOICE eq 1 then x=Mvir[index]
y=MetalsHotGasSelf[index]/MetalsHotGas[index]
xtitle='Stellar mass/h!u-1!nM!d!Mn!n'
if CHOICE then xtitle='M!d200!n/h!u-1!nM!d!Mn!n'
if CHOICE then xrange=[1e11,1e15] else xrange=[1e4,1e12]
if CHOICE then yrange=[1e-3,10.] else yrange=[1e-4,10.]

plot, x, y, /nodata, /xlog, /ylog, psym=3, $
  xtitle = xtitle, $
  xrange = xrange,xstyle=1, $
  ytitle = 'MetalsHotGasSelf / MetasHotGas', $
  yrange = yrange,ystyle=1
oplot, x, y, color=2, psym=8
oplot, [1e2,1e15],[1.,1.]

if PS then device, /close_file
if PS then set_plot,'x'
