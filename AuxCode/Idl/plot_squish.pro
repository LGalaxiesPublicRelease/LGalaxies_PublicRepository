;=========================================================================
  
   PRO plot_squish,ps=ps

;=========================================================================
;
;  Example procedure to read in squished L-galaxies snapshot data 
;  and create metallicity plot
;
; Inputs:
;   ps - keyword to cause the plot to be written to a postscript file.
;        (by default the plot will be drawn on the screen)
;
;-------------------------------------------------------------------------

; Define which files you want to read in
DirName = '../output/DLB07/' & format='DLB07'
;DirName = '../output/HT09/' & format='HT09'
FileName = 'SA_z0.00'
ModelName = DirName + FileName

;-------------------------------------------------------------------------------

; read in data and convert to Msun/h

type=size(0.,/type)
StellarMass=read_squish(ModelName,'StellarMass',type)
MetalsStellarMass=read_squish(ModelName,'MetalsStellarMass',type)
StellarMass=1e10*StellarMass
MetalsStellarMass=1e10*MetalsStellarMass

;-------------------------------------------------------------------------------

; Create plot

start_plot

; Open file for all plots
if (keyword_set(ps)) then begin
    set_plot, 'PS'
    device, filename = 'metals.ps', /color, xsize = 15, ysize = 10, $
      xoffset=1, $
      yoffset=5
endif
  
x=StellarMass[where(StellarMass gt 1.e5)]
y=MetalsStellarMass[where(StellarMass gt 1.e5)]

plot, x, y/x, /nodata, /xlog, /ylog, psym=3, $
  xtitle = 'Stellar mass/h!u-1!nM!d!Mn!n', $
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

if (keyword_set(ps)) then begin
    device, /close_file
    set_plot,'x'
endif

;------------------------------------------------------------------------------


end
