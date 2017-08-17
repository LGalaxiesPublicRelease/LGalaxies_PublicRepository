;=========================================================================
;
;  Script to read in L-galaxies snapshot data
;  with a star-formation history for the galaxies
;  and produce a scatter plot of luminosity versus mass

;  Parameters:
;   ps - keyword to cause the plot to be written to a postscript file.
;        (by default the plot will be drawn on the screen)
    PS = 0b
;   SNAP - which snapshot is it (affects axes)
    SNAP = 63
;    SNAP=41
;   Unfortunaetly, have to set template here
@'../output/snap_sfh_template.pro'
;@'../output/snap_metals_template.pro'
;
;-------------------------------------------------------------------------

; Define which files you want to read in
DirName = '../output/'
if SNAP eq 63 then FileName = 'SA_z0.00'
if SNAP eq 41 then FileName = 'SA_z0.99'
ModelName = DirName + FileName
;FirstFile = 0
;LastFile = 511
FirstFile = 287
LastFile = 287

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

read_lgal,ModelName,Template,GalStruct $
         ,FileNr=FileNr,TreeNr=TreeNr,NGalsPerTree=NGalsPerTree $
         ,FirstFile=FirstFile,LastFile=LastFile ;,/swap_endian

;-------------------------------------------------------------------------------

; Define simple arrays containing useful data from structure array
print, ' Now defining useful quantities...'

;StellarMass = GalStruct.StellarMass * 1.0e10
StellarMass = (GalStruct.DiskMass+GalStruct.BulgeMass) * 1.0e10
ColdGas = GalStruct.ColdGas * 1.0e10
HotGas = GalStruct.HotGas * 1.0e10

; Nullify orginal data arrays to save memory
FileNr=0
TreeNr=0
NGalsPerTree=0
GalStruct=0

;-------------------------------------------------------------------------------

; Create plot

start_plot

; Open file for all plots
if PS then set_plot, 'PS'
if PS then device, filename = 'gas.ps', /color, $
                  xsize = 15, ysize = 10, $
                  xoffset=1,  yoffset=5
  
x=ColdGas
y=Hotgas

plot, x, y, /nodata, /xlog, /ylog, psym=3, $
  xtitle = 'Cold Gas mass/h!u-1!nM!d!Mn!n', $
  xrange = [1e6,1e12],xstyle=1, $
  ytitle = 'Hot Gas mass/h!u-1!nM!d!Mn!n', $
  yrange = [1e2,1e14],ystyle=1
oplot, x, y, color=2, psym=3

if PS then device, /close_file
if PS then set_plot,'x'

;------------------------------------------------------------------------------
