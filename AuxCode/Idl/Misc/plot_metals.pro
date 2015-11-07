;=========================================================================
;
;  plot_metals
;
;  Example script to read in L-galaxies snapshot data
;  with a star-formation history for the galaxies
;  and produce a plot of stellar metallicities
;  Assumes that all the galaxies are from a snap (ie all at the same time).

;  Parameters:
;   ps - keyword to cause the plot to be written to a postscript file.
;        (by default the plot will be drawn on the screen)
    PS = 0b
;   SNAP - which snapshot is it (affects axes)
    SNAP = 63
;    SNAP=41
;   Metals flag to read in and plot appropriate data
;   Unforutately, don't seem to be allowed to treat @ as a
;   program statement, so can't bury it in an if..then..else construct.
;    METALS=0b
;    @'../output/snap_template.pro'
    METALS=0b
;    @'../output/snap_sfh_template.pro' $
    @'../output/snap_sfh_template.pro' $
;
;-------------------------------------------------------------------------

; Define which files you want to read in
DirName = '../output/'
if SNAP eq 63 then FileName = 'SA_z0.00'
if SNAP eq 41 then FileName = 'SA_z0.99'
ModelName = DirName + FileName
FirstFile = 287
LastFile = 287

;-------------------------------------------------------------------------------

; read in data

read_lgal,ModelName,Template,GalStruct $
         ,FileNr=FileNr,TreeNr=TreeNr,NGalsPerTree=NGalsPerTree $
         ,FirstFile=FirstFile,LastFile=LastFile ;,/swap_endian

;-------------------------------------------------------------------------------

; Define simple arrays containing useful data from structure array
print, ' Now defining useful quantities...'

;StellarMass = GalStruct.StellarMass * 1.0e10
StellarMass = (GalStruct.DiskMass+GalStruct.BulgeMass) * 1.0e10
;MetalsStellarMass = GalStruct.MetalsStellarMass * 1.0e10
if METALS then $
   ;MetalsStellarMass = GalStruct.MetalsStellarMass[1] * 1.0e10 $
   MetalsStellarMass = (GalStruct.MetalsDiskMass[1]+Galstruct.MetalsBulgeMass[1]) * 1.0e10 $
else $
   ;MetalsStellarMass = GalStruct.MetalsStellarMass * 1.0e10
   MetalsStellarMass = (GalStruct.MetalsDiskMass+Galstruct.MetalsBulgeMass) * 1.0e10
if METALS then Metals1a = GalStruct.MetalsStellarMass[0] * 1.0e10
if METALS then Metals2 = GalStruct.MetalsStellarMass[1] * 1.0e10
if METALS then MetalsAGB = GalStruct.MetalsStellarMass[2] * 1.0e10
 
; Nullify orginal data arrays to save memory
FileNr=0
TreeNr=0
NGalsPerTree=0
GalStruct=0

;-------------------------------------------------------------------------------
  
; Create plot

start_plot

;-------------------------------------------------------------------------------

; Create plots

start_plot

; Open file for PS plot
if PS then set_plot, 'PS' else window,0
if PS then device, filename = 'metals.ps', /color, $
                   xsize = 15, ysize = 10, $
                   xoffset=1,  yoffset=5

x=StellarMass
y=MetalsStellarMass

plot, x, y/x, /nodata, /xlog, /ylog, psym=3, $
  xtitle = 'Stellar mass/h!u-1!nM!d!Mn!n', $
  xrange = [1e2,1e12],xstyle=1, $
  ytitle = 'Stellar metallicity', $
  yrange = [1e-8,1],ystyle=1
oplot, x, y/x, color=2, psym=3

; Uncomment this to create a second plot on same axes
;y=MetalsIaStellarMass
;oplot, x, y, color=4, psym=3

; This will give a metal-ratio plot
;y=MetalsIaStellarMass/MetalsStellarMass
;plot, x, y, /xlog, psym=3, $
;  xtitle = 'Stellar mass/h!u-1!nM!dO!n', $
;  ytitle = 'MetalsIaStellarMass/MetalsStellarMass'

if PS then device, /close_file
if PS then set_plot,'x'

