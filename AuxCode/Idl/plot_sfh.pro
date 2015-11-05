;=========================================================================
;
;  plot_sfh
;
;  Example script to read in L-galaxies snapshot data
;  with a star-formation history for the galaxies
;  and produce an example plot of star-formation rate versus time.
;  Assumes that all the galaxies are from a snap (ie all at the same time).

;  Parameters:
;   ps - keyword to cause the plot to be written to a postscript file.
;        (by default the plot will be drawn on the screen)
    PS = 0b
;   SNAP - which snapshot is it (affects axes)
    SNAP = 63
;    SNAP=41
;
;-------------------------------------------------------------------------

; Define which files you want to read in
DirName = '../output/'
format='SFH_snap'
if SNAP eq 63 then FileName = 'SA_z0.00'
if SNAP eq 41 then FileName = 'SA_z0.99'
ModelName = DirName + FileName
FirstFile = 287
LastFile = 287

Volume=500^3/512 ; Volume in (Mpc/h)^3
   
;-------------------------------------------------------------------------------

; read in data

.comp read_lgal
read_lgal,ModelName,GalStruct $
         ,FileNr=FileNr,TreeNr=TreeNr,NGalsPerTree=NGalsPerTree $
         ,FirstFile=FirstFile,LastFile=LastFile,format=format ;,/swap_endian

;-------------------------------------------------------------------------------

; create data for plot

Ngals=n_elements(GalStruct)
Nbin=Galstruct[0].sfh_ibin+1

x=fltarr(Nbin)
for i = 0, Nbin-1 do $
   x[i]=0.5*(Galstruct[0].sfh_time[i]+Galstruct[0].sfh_time[i+1])

dt=fltarr(Nbin)
for i = 0, Nbin-1 do dt[i]=Galstruct[0].sfh_time[i]-Galstruct[0].sfh_time[i+1]
y=fltarr(Nbin)
for i = 0, Nbin-1 do y[i]=total(Galstruct.sfh_StellarMass[i])/dt[i]
y=y*1e10/Volume
ym=fltarr(Nbin)
for i = 0, Nbin-1 do ym[i]=total(Galstruct.sfh_MetalsStellarMass[i])/dt[i]
ym=ym*1e10/Volume

;-------------------------------------------------------------------------------

; Create plots

start_plot

; Open file for PS plot
if PS then set_plot, 'PS' else window,0
if PS then device, filename = 'sfh.ps', /color, $
                   xsize = 15, ysize = 10, $
                   xoffset=1,  yoffset=5

if SNAP eq 63 then xrange=[3e6,3e10]
if SNAP eq 41 then xrange=[1e6,1e10]
plot, [x,x], [y,ym], /nodata, /xlog, /ylog, $
      xrange = xrange, xstyle=1, $
      ;yrange = [ymin, ymax], ystyle=1, $
      xtitle = 'time/yr', $
      ytitle = 'SFR / M!d!Mn!n yr!u-1!n (h!u-1!n Mpc)^u3!n'
oplot, x, y, color=2
oplot, x, y, psym=2, color=2
oplot, x, ym, color=4, linestyle=1
oplot, x, ym, psym=2, color=4
if SNAP eq 63 then xyouts, 1e7, 0.025, 'StellarMass', color=2
if SNAP eq 63 then xyouts, 1e7, 0.0005, 'MetalsStellarMass', color=4
if SNAP eq 41 then xyouts, 3e8, 0.01, 'StellarMass', color=2
if SNAP eq 41 then xyouts, 3e8, 0.0001, 'MetalsStellarMass', color=4

if PS then device, /close_file
if PS then set_plot,'x'

if PS then set_plot, 'PS' else window,1
if PS then device, filename = 'sfh_met.ps', /color, $
                   xsize = 15, ysize = 10, $
                   xoffset=1,  yoffset=5
  
plot, x, ym/y, /nodata, /xlog, /ylog, $
      xrange = [3e6, 3e10], xstyle=1, $
      ;yrange = [0.001, 0.1], ystyle=1, $
      xtitle = 'time/yr', $
      ytitle = 'Metallicity'
oplot, x, ym/y
oplot, x, ym/y, psym=2

if PS then device, /close_file
if PS then set_plot,'x'

;------------------------------------------------------------------------------
