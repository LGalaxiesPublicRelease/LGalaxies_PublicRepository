;pro read_snap
;pro read_tree

;pro plot_kband_bband_smf
;pro chi_square_and_plot
;plot_obs_sfr_vmax_contour
;pro age_hist
;pro plot_obs_age_vmax
;pro ssfr_hist
;pro plot_obs_ssfr_vmax
;pro color_hist
;pro qi_obs_color_hist
;pro plot_UVJ_color
;pro plot_red_fraction_colorcut
;pro plot_smf_by_color
;pro plot_individual_obs_smf_by_color
;pro plot_gas_mass_function
;plot_mstar_metals
;pro plot_model_metals

;pro write_to_file
;pro write_to_file_with_errors
;pro plot_label
;pro get_slope
;pro make_float_array
;pro median_and_percentiles
;function make_log_array
;PRO COLORBAR
;function comdist
;pro plot_hist_2d
;pro invert_color_table
;PRO ROTATION

pro read_snap, structure, structure_definition, DirName, FileName, FirstFile, LastFile

   print,''
   print,''
   print,' doing file '+FileName
   print,''
   ModelName = DirName + FileName
   TotNTrees = 0L & TotNGals = 0L
   close, 1
   for fnr = FirstFile, LastFile do begin
       fname = strcompress(ModelName + '_' + string(fnr), /remove_all)
       openr, 1, fname;, /swap_endian 
       ;print,            fname
       Ntrees = 0L     & readu, 1, Ntrees
       NtotGals = 0L   & readu, 1, NtotGals
       TotNTrees = TotNTrees + Ntrees
       TotNGals =  TotNGals  + NtotGals      
       close, 1
    endfor  
   print, " Total number of Trees = ", ToTNTrees
   print, " Total number of galaxies = ", TotNGals
   structure = replicate(structure_definition, TotNGals) 
   offset = 0L
   for fnr = FirstFile, LastFile do begin
      print, ' file:', fnr
      fname = strcompress(ModelName+'_'+string(fnr), /remove_all)
      openr, 1, fname           ;, /swap_endian    
      Ntrees = 0L     & readu, 1, Ntrees   & print, ' Ntrees:', Ntrees
      NtotGals = 0L   & readu, 1, NtotGals & print, ' NtotGals:', NtotGals
      GalsPerTree = lonarr(Ntrees)
      if(NtotGals gt 0) then begin
      readu, 1, GalsPerTree
      GG = replicate(structure_definition, NtotGals)
      readu, 1, GG
      structure[offset:offset+NtotGals-1] = GG[*]     
      offset = offset + NtotGals  
   endif     
      close, 1
   endfor

end




pro read_tree, structure, structure_definition, DirName, FileName, FirstFile, LastFile
 
   print, ' Now reading the input files...'   
   ModelName = DirName + FileName
   TotNTrees = 0L & TotNGals = 0L
   close, 1
   for fnr = FirstFile, LastFile do begin
       fname = strcompress(ModelName + '_' + string(fnr), /remove_all)
       openr, 1, fname;, /swap_endian
       one = 0L & readu,1,one  
       nbytes = 0L & readu, 1, nbytes 
       ngals=0L & readu, 1, ngals   
       TotNGals =  TotNGals  + ngals   
       close, 1
    endfor
  
   print, " Total number of galaxies = ", TotNGals
  
   structure = replicate(structure_definition, TotNGals) 
   offset = 0L
   off = 0L
   for fnr = FirstFile, LastFile do begin
       print, ' file:', fnr
       fname = strcompress(ModelName+'_'+string(fnr), /remove_all)
       openr, 1, fname;, /swap_endian     
       one = 0L & readu,1,one
       nbytes = 0L & readu, 1, nbytes
       ngals=0L & readu, 1, ngals & print, ' NGals:', ngals
       nskip=nbytes/4-3
       ib=fltarr(nskip)
       readu,1,ib    
       GG = replicate(structure_definition, ngals)
       readu, 1, GG
       structure[offset:offset+ngals-1] = GG[*]
       offset = offset + ngals     
       close, 1
    endfor   
end



pro plot_kband_bband_smf,  G0, G1, G2, G3, G0_MRII, G1_MRII, G2_MRII, G3_MRII, $
                                 StellarMass_MR_0, StellarMass_MRII_0, StellarMass_MR_04, StellarMass_MRII_04, $
                                 StellarMass_MR_1, StellarMass_MRII_1, StellarMass_MR_2, StellarMass_MRII_2, $
                                 StellarMass_MR_3, StellarMass_MRII_3, total_chi2, $
                                 write_files, plot_after_write, file_to_write, $
                                 sample, sample_file, plot_obs, plot_only_new_obs, do_lums, $
                                 linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                                 linestyle_previous_model2, file_previous_model2, do_previous_model2, $
                                 prefix_previous_model1, prefix_previous_model2, prefix_this_model, $
                                 MRII,  volume_MR, volume_MRII, hubble_h, $
                                 Datadir, do_mcmc_regions, mcmc_folder, min_obs_err, redshift

  Hubble_h_WMAP7=0.7
  Hubble_h_WMAP1=0.73

  loadct,6, /SILENT

  if(plot_obs eq 0) then $
     joint_obs_symbols_size=1.2 $
  else joint_obs_symbols_size=0.7  

  yobs1=-0.8
  yobs12=-0.85
  yobs2=-1.05
  yobs22=-1.1
  yobs3=-1.3
  yobs32=-1.35
  yobs4=-1.55
  yobs42=-1.6
  yobs5=-1.8
  yobs52=-1.85
  yobs6=-2.05
  yobs62=-2.1
  yobs7=-2.3
  yobs72=-2.35
  yobs8=-2.55
  yobs82=-2.6
  
  if(MRII) then begin     
     yobs1=0.2
     yobs12=0.15
     yobs2=-0.05
     yobs22=-0.1
     yobs3=-0.3
     yobs32=-0.35
     yobs4=-0.55
     yobs42=-0.6
     yobs5=-0.8
     yobs52=-0.85
     yobs6=-1.05
     yobs62=-1.1
     yobs7=-1.3
     yobs72=-1.35
     yobs8=-1.55
     yobs82=-1.6
  endif
  

;-------------------------------------------------------------------------------
 


;******************************************************************************
;******************************************************************************
;******                                                            ************
;******               KBAND LUMINOSITY FUNCTION                    ************
;******                                                            ************
;******************************************************************************
;******************************************************************************


;NUMBER 0

  if(do_lums eq 1) then begin


  print,''
  print,''
  print,''
  print,'**********************'
  print,'*      Kband LF      *'
  print,'**********************'
  print,''

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                      ;;;
;;;    Kband LF Z=0      ;;;
;;;                      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
     G=G0   
     G_MRII=G0_MRII

     char_redshift=STRTRIM(number_formatter(redshift[0], DECIMAL=2),2)

     bin=0.25
     xmin=-25.9
     xmax=-17.0

     if(MRII) then xmax=-17.5
     ymin=-6.9
     ymax=-0.5
     if(MRII) then ymax =0.5
    
     plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
           xstyle = 1, ystyle = 1, $
           ytitle = '!4U!3 (h!U3!NMpc!U-3!Nmag!U-1!N)'
 

;JOINT OBS 
     symbols, 2,  joint_obs_symbols_size
     if (plot_obs eq 0) then begin
        oploterror, [-25.5,-25.5], [yobs2+0.04,yobs2+0.04], [0.15,0.15], $
                    psym=8, color = 200, errcolor=200, HATLENGTH = 80.0
        xyouts, -25.25, yobs22, 'Observations used in MCMC',charsize=1.2
     endif


 

;NEW MODEL
     chi_square_and_plot, mcmc_folder, 'kband', char_redshift, total_chi2, plot_obs, plot_only_new_obs,G.CentralMvir, G_MRII.CentralMvir, $
                          G.MagDust[9]-5.*Alog10(hubble_h), G_MRII.MagDust[9]-5.*Alog10(hubble_h), $
                          volume_MR, volume_MRII, Alog10(StellarMass_MR_0), Alog10(StellarMass_MRII_0),Hubble_h, $
                          write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                          linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                          linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII


;PREVIOUS MODELS  
     if(do_previous_model1 eq 1) then begin
        readcol,file_previous_model1+'_kband_z'+char_redshift+'.txt',mag,phi, /SILENT 
        oplot,mag,phi, color=70,linestyle=linestyle_previous_model1, thick=6
     endif
     if(do_previous_model2 eq 1) then begin
        readcol,file_previous_model2+'_kband_z'+char_redshift+'.txt',mag,phi, /SILENT 
        oplot,mag,phi, color=70,linestyle=linestyle_previous_model2, thick=6
     endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                              ;;;;
;;;;    INDIVIDUAL OBSERVATIONS   ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

     if (plot_obs eq 1) then begin

        ;JONES
        ;plus +0.05 to dont overplot with other errors
        ;WMAP1
        readcol,Datadir+'jones2006_kbandLF.dat',mag_k1,phi1,errorup1,errordown1, /SILENT  
        symbols, 32, 0.3
        oploterror, mag_k1+1.825, phi1, errorup1, color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mag_k1+1.825, phi1, errordown1, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0

        oploterror, [-25.5,-25.5], [yobs1,yobs1], [0.1,0.1], $
                    psym=8, color=0, errcolor=0, HATLENGTH = 80.0
        xyouts, -25.25, yobs12, 'Jones (2006)',charsize=0.9,charthick = 3   
 
        ;BELL
        ;WMAP1
        readcol,Datadir+'bell2003_kbandLF.dat',mag_k2,phi2,errordown2,errorup2, /SILENT
        symbols, 30, 0.25      
        oploterror, mag_k2+1.825, Alog10(phi2), Alog10(errorup2/phi2), color = 120, $
                    errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mag_k2+1.825, Alog10(phi2), Alog10(phi2/errordown2), color = 120, $
                    errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0
        oploterror, [-25.5,-25.5], [yobs2,yobs2], [0.1,0.1], $
                    psym=8, color=120, errcolor=120, HATLENGTH = 80.0 
        xyouts, -25.25, yobs22, 'Bell (2003)',charsize=0.9,charthick = 3    
 
;COLE
   ;WMAP1
        loadct,2, /SILENT
        symbols, 2, 0.7
        readcol,Datadir+'cole2001_kbandLF_rebin.dat',mag_k3,phi3,err3, /SILENT
        oploterror, mag_k3+1.825, Alog10(phi3), alog10((phi3+err3)/phi3), $ 
                    color=200, errcolor = 200, psym = 8,HATLENGTH = 80.0
        oploterror, [-25.5,-25.5], [yobs3,yobs3], [0.1,0.1], $
                    psym=8, color=200, errcolor=200, HATLENGTH = 80.0  
        xyouts, -25.25, yobs32, 'Cole et al. (2001)',charsize=0.9,charthick = 3    
        loadct,6, /SILENT
     endif  ;plot_obs



;check sample 
     if(sample eq 1) then begin
        readcol,sample_file+'1_z0.10.txt',kband,phi,err,model, /SILENT
        oplot, kband, Alog10(model), color=0, thick=6
     endif

 
     plot_label, xlog=0, ylog=0, type='label', label='K-band, z=0.0' , xmin, xmax, ymin, ymax, $
                 x_percentage=5.5, x2_percentage=0., y_percentage=4., $
                 charsize=1.3, charthick=4

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                      ;;;
;;;    Kband LF Z=1      ;;;
;;;                      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

     G=G1
     G_MRII=G1_MRII

     char_redshift=STRTRIM(number_formatter(redshift[2], DECIMAL=2),2)

     mag=G.MagDust[9]
     multiplot
  
     plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
           xstyle = 1, ystyle = 1         
    

;JONES - z=0
     readcol,Datadir+'jones2006_kbandLF.dat',mag_k1,phi1,errorup1,errordown1, /SILENT  
     symbols, 32, 0.3
     oplot, mag_k1+1.825, phi1,linestyle=1
 
;NEW MODEL
     chi_square_and_plot, mcmc_folder, 'kband', char_redshift, total_chi2, plot_obs,plot_only_new_obs, G.CentralMvir, G_MRII.CentralMvir, $
                          G.MagDust[9]-5.*Alog10(hubble_h), G_MRII.MagDust[9]-5.*Alog10(hubble_h), $
                          volume_MR, volume_MRII, Alog10(StellarMass_MR_1), Alog10(StellarMass_MRII_1),Hubble_h, $
                          write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                          linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                          linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII

;PREVIOUS MODELS 
     if(do_previous_model1 eq 1) then begin
        readcol,file_previous_model1+'_kband_z'+char_redshift+'.txt',mag,phi, /SILENT
        oplot,mag,phi, color=70,linestyle=linestyle_previous_model1, thick=6
     endif
     if(do_previous_model2 eq 1) then begin
        readcol,file_previous_model2+'_kband_z'+char_redshift+'.txt',mag,phi, /SILENT 
        oplot,mag,phi, color=70,linestyle=linestyle_previous_model2, thick=6
     endif

 
     plot_label, xlog=0, ylog=0, type='label', label='K-band, z=1.0' , xmin, xmax, ymin, ymax, $
                 x_percentage=5.5, x2_percentage=0., y_percentage=4., $
                 charsize=1.3, charthick=4
 

;check sample 
     if(sample eq 1) then begin
        readcol,sample_file+'1_z1.00.txt',kband,phi,err,model, /SILENT
        oplot, kband, Alog10(model), color=0, thick=6
     endif




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                              ;;;;
;;;;    INDIVIDUAL OBSERVATIONS   ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


     if (plot_obs eq 1) then begin


;DRORY
   ;h=0.65
        readcol,Datadir+'drory_kband_z0810.txt',mag,phi,errordown,errorup, /SILENT
        symbols, 2, 0.7
        oploterror, mag+1.872, Alog10(phi), Alog10((phi+errorup)/phi), color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        
        error=phi+errordown
        log_err=Alog10(error)
        sel=where(error eq 0.,n)
        if(n gt 0) then  log_err[sel] =0.
        
        oploterror, mag+1.872, Alog10(phi), Alog10(phi)-log_err, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0 
        
        
        oploterror, [-25.5,-25.5], [yobs3,yobs3], [0.1,0.1], $
                    psym=8, color=0, errcolor=0, HATLENGTH = 80.0
        xyouts, -25.25, yobs32, 'Drory (2003) - 0.8<z<1.0',charsize=0.9,charthick = 3   
        

;POZZETTI
        loadct,2, /SILENT
        readcol,Datadir+'pozzetti_kband_z07513.txt',mag,phi,errordown,errorup, /SILENT
        symbols, 30, 0.25  
        oploterror, mag+1.825-5.*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3), errorup, color = 200, $
                    errcolor = 200, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mag+1.825-5.*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3), errordown, color = 200, $
                    errcolor = 200, psym = 8, /loba,HATLENGTH = 80.0
        
        oploterror, [-25.5,-25.5], [yobs1,yobs1], [0.1,0.1], $
                    psym=8, color=200, errcolor=200, HATLENGTH = 80.0
        xyouts, -25.25, yobs12, 'Pozzetti (2003) - 0.75<z<1.3',charsize=0.9,charthick = 3   
        loadct,6, /SILENT

;CIRASUOLO 2010
  ;symbols, 2, 1.0 
        readcol,Datadir+'cirasuolo2010_new_z10125.txt',mag,phi,errordown, errorup, /SILENT 
        symbols, 1, 0.8  
        oploterror, mag-5.*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3), errorup, $
                    color = 0, errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mag-5.*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3), errordown, $
                    color = 0, errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
        oploterror, [-25.5,-25.5], [yobs2,yobs2], [0.1,0.1], $
                    psym=8,color = 0, errcolor=0, HATLENGTH = 80.0
        xyouts, -25.25, yobs22, 'Cirasuolo (2010) - 0.8<z<1.25',charsize=0.9,charthick = 3  

;Caputi2006
        readcol,Datadir+'caputi2006_kband_z10.txt',mag,phi,errordown,errorup, /SILENT
        symbols, 32, 0.3
        data_hubble=0.7
        oploterror, mag+1.825-5.*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3), errorup, color = 120, $
                    errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mag+1.825-5.*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3), errordown, color = 120, $
                    errcolor = 120, psym = 8, /loba, HATLENGTH = 80.0
        
        oploterror, [-25.5,-25.5], [yobs4,yobs4], [0.1,0.1], $
                    psym=8, color = 120, errcolor=120, HATLENGTH = 80.0
        xyouts, -25.25, yobs42, 'Caputi (2006) - z=1.0',charsize=0.9,charthick = 3   
        
     endif  ;plot_obs


;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                      ;;;
;;;    Kband LF Z=2      ;;;
;;;                      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

     G=G2
     G_MRII=G2_MRII

     char_redshift=STRTRIM(number_formatter(redshift[3], DECIMAL=2),2)

     mag=G.MagDust[9]  
     multiplot
 
     plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
           xstyle = 1, ystyle = 1, $
           xtitle = 'M!DK!N(AB)-5log!D10!N(h)',ytitle = '!4U!3 (h!U3!NMpc!U-3!Nmag!U-1!N)'
     

;JONES - z=0
     readcol,Datadir+'jones2006_kbandLF.dat',mag_k1,phi1,errorup1,errordown1, /SILENT  
     symbols, 32, 0.3
     oplot, mag_k1+1.825, phi1,linestyle=1
 
;NEW MODEL
     chi_square_and_plot, mcmc_folder, 'kband', char_redshift, total_chi2, plot_obs, plot_only_new_obs, G.CentralMvir, G_MRII.CentralMvir, $
                          G.MagDust[9]-5.*Alog10(hubble_h), G_MRII.MagDust[9]-5.*Alog10(hubble_h), $
                          volume_MR, volume_MRII, Alog10(StellarMass_MR_2), Alog10(StellarMass_MRII_2),Hubble_h, $
                          write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                          linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                          linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII


;PREVIOUS MODELS 
     if(do_previous_model1 eq 1) then begin
        readcol,file_previous_model1+'_kband_z'+char_redshift+'.txt',mag,phi, /SILENT 
        oplot,mag,phi, color=70,linestyle=linestyle_previous_model1, thick=6
     endif
     if(do_previous_model2 eq 1) then begin
        readcol,file_previous_model2+'_kband_z'+char_redshift+'.txt',mag,phi, /SILENT 
        oplot,mag,phi, color=70,linestyle=linestyle_previous_model2, thick=6
     endif


;check sample 
     if(sample eq 1) then begin
        readcol,sample_file+'1_z2.00.txt',kband,phi,err,model, /SILENT
        oplot, kband, Alog10(model), color=0, thick=6
     endif
     
     
     if (plot_obs eq 1) then begin


;CIRASUOLO 2010_new
        symbols, 1, 0.8 
        readcol,Datadir+'cirasuolo2010_new_z20_rebin.txt',mag,phi,err, /SILENT
        oploterror, mag-0.05, Alog10(phi), $
                    Alog10((phi+err)/phi),/hiba, psym=8, color=0, errcolor=0, HATLENGTH = 80.0  
        oploterror, mag-0.05, Alog10(phi), $
                    Alog10(phi/(phi-err)),/loba, psym=8, color=0, errcolor=0, HATLENGTH = 80.0 
        oploterror, [-25.5,-25.5], [yobs1,yobs1], [0.1,0.1], $
                    psym=8,color = 0, errcolor=0, HATLENGTH = 80.0
        xyouts, -25.25, yobs12, 'Cirasuolo (2010) - 1.75<z<2.5',charsize=0.9,charthick = 3   
 


;Caputi2006
        readcol,Datadir+'caputi2006_kband_z20.txt',mag,phi,errordown,errorup, /SILENT
        symbols, 32, 0.3
        data_hubble=0.7
        oploterror, mag+1.825-5.*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3), $
                    errorup, color = 120, errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mag+1.825-5.*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3), $
                    errordown, color = 120,errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0
        
        oploterror, [-25.5,-25.5], [yobs2,yobs2], [0.1,0.1], $
                    psym=8, color = 120, errcolor=120, HATLENGTH = 80.0
        xyouts, -25.25, yobs22, 'Caputi (2006) - z=2.0',charsize=0.9,charthick = 3   
        
     endif  ;plot_obs



     plot_label, xlog=0, ylog=0, type='label', label='K-band, z=2.0' , xmin, xmax, ymin, ymax, $
                 x_percentage=5.5, x2_percentage=0., y_percentage=4., $
                 charsize=1.3, charthick=4
 









;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                      ;;;
;;;    Kband LF Z=3      ;;;
;;;                      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

     G=G3
     G_MRII=G3_MRII

     char_redshift=STRTRIM(number_formatter(redshift[4], DECIMAL=2),2)

     bin=0.25   
     mag=G.MagDust[9]
  
     multiplot
     plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
           xstyle = 1, ystyle = 1, $
           xtitle = 'M!DK!N(AB)-5log!D10!N(h)'
  

;JONES - z=0
     readcol,Datadir+'jones2006_kbandLF.dat',mag_k1,phi1,errorup1,errordown1, /SILENT  
     symbols, 32, 0.3
     oplot, mag_k1+1.825, phi1,linestyle=1
 
;MODEL
     chi_square_and_plot, mcmc_folder, 'kband', char_redshift, total_chi2, plot_obs,  plot_only_new_obs,G.CentralMvir, G_MRII.CentralMvir, $
                          G.MagDust[9]-5.*Alog10(hubble_h), G_MRII.MagDust[9]-5.*Alog10(hubble_h), $
                          volume_MR, volume_MRII, Alog10(StellarMass_MR_3), Alog10(StellarMass_MRII_3),Hubble_h, $
                          write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                          linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                          linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII



;PREVIOUS MODELS 
     if(do_previous_model1 eq 1) then begin
        readcol,file_previous_model1+'_kband_z'+char_redshift+'.txt',mag,phi, /SILENT 
        oplot,mag,phi, color=70,linestyle=linestyle_previous_model1, thick=6
     endif
     if(do_previous_model2 eq 1) then begin
        readcol,file_previous_model2+'_kband_z'+char_redshift+'.txt',mag,phi, /SILENT 
        oplot,mag,phi, color=70,linestyle=linestyle_previous_model2, thick=6
     endif



;LABELS
;redshift
     plot_label, xlog=0, ylog=0, type='label', label='K-band, z=3.0' , xmin, xmax, ymin, ymax, $
                 x_percentage=5.5, x2_percentage=0., y_percentage=4., $
                 charsize=1.3, charthick=4
 
;observations
     if(plot_obs eq 1) then begin
        plot_label, xlog=0, ylog=0, type='label', label='Observations for MCMC' , xmin, xmax, ymin, ymax, $
                    x_percentage=4.0, x2_percentage=0., y_percentage=2.8 , $
                    charsize=1.2, charthick=3 
        plot_label, xlog=0, ylog=0, type='symbol', xmin, xmax, ymin, ymax, $
                    x_percentage=3.8, x2_percentage=0., y_percentage=3.0, $
                    color=200, errcolor=200, sym_num=2, sym_size=joint_obs_symbols_size, HATLENGTH = 80, err_size=0.15
     endif
;models

     plot_label, xlog=0, ylog=0, type='label', label=prefix_this_model , xmin, xmax, ymin, ymax, $
                 x_percentage=4.0, x2_percentage=0., y_percentage=2., $
                 charsize=1.1, charthick=3
     plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                 x_percentage=3.2, x2_percentage=3.9, y_percentage=2.15, $
                 linestyle=0, color=70, linethick=8
     
     if(do_previous_model2 eq 1) then begin
        plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model2, xmin, xmax, ymin, ymax, $
                    x_percentage=4.0, x2_percentage=0., y_percentage=1.3, $
                    charsize=1.1, charthick=3
        plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                    x_percentage=3.2, x2_percentage=3.9, y_percentage=1.45, $
                    linestyle=linestyle_previous_model2, color=70, linethick=8
     endif
     if(do_previous_model1 eq 1) then begin
        plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model1, xmin, xmax, ymin, ymax, $
                    x_percentage=4.0, x2_percentage=0., y_percentage=0.6, $
                    charsize=1.1, charthick=3
        plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                    x_percentage=3.2, x2_percentage=3.9, y_percentage=0.75, $
                    linestyle=linestyle_previous_model1, color=70, linethick=8
     endif
     


     if (plot_obs eq 1) then begin


;CIRASUOLO 2010_new
        symbols, 1, 0.8  
        readcol,Datadir+'cirasuolo2010_new_z30_rebin.txt',mag,phi,err, /SILENT
        oploterror, mag-0.05, Alog10(phi), $
                    Alog10((phi+err)/phi),/hiba, psym=8, color=0, errcolor=0, HATLENGTH = 80.0  
        oploterror, mag-0.05, Alog10(phi), $
                    Alog10(phi/(phi-err)),/loba, psym=8, color=0, errcolor=0, HATLENGTH = 80.0
        oploterror, [-25.5,-25.5], [yobs1,yobs1], [0.1,0.1], $
                    psym=8,color = 0, errcolor=0, HATLENGTH = 80.0
        xyouts, -25.25, yobs12, 'Cirasuolo (2010) - 2.5<z<3.5',charsize=0.9,charthick = 3   
 

;Caputi2006
        readcol,Datadir+'caputi2006_kband_z25.txt',mag,phi,errordown,errorup, /SILENT
        symbols, 32, 0.3 
        oploterror, mag+1.825-5.*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3), errorup, color = 120, errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mag+1.825-5.*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3), errordown, color = 120, errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0
        
        oploterror, [-25.5,-25.5], [yobs2,yobs2], [0.1,0.1], $
                    psym=8, color = 120, errcolor=120, HATLENGTH = 80.0
        xyouts, -25.25, yobs22, 'Caputi (2006) - z=2.5',charsize=0.9,charthick = 3   


;Saracco2006
   readcol,Datadir+'saracco_kband_z1940.txt',mag,phi,errordown,errorup, /SILENT
   symbols, 20, 0.3 
   oploterror, mag+1.825-5.*Alog10(Hubble_h_WMAP7), Alog10(phi/Hubble_h_WMAP7^3), Alog10((phi+errorup)/phi), color = 0, $
               errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
   oploterror, mag+1.825-5.*Alog10(Hubble_h_WMAP7), Alog10(phi/Hubble_h_WMAP7^3), Alog10(phi/(phi-errordown)), color = 0, $
               errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
  
   oploterror, [-25.5,-25.5], [yobs3,yobs3], [0.1,0.1], $
          psym=8, color = 0, errcolor=0, HATLENGTH = 80.0
   xyouts, -25.25, yobs32, 'Saracco (2006) - 1.9<z<4.0',charsize=0.9,charthick = 3   

endif ;plot_obs



;check sample 
     if(sample eq 1) then begin
        readcol,sample_file+'1_z3.00.txt',kband,phi,err,model, /SILENT
        oplot, kband, Alog10(model), color=0, thick=6
     endif
     
;******************************************************************************
;******************************************************************************
;******                                                            ************
;******               BBAND LUMINOSITY FUNCTION                    ************
;******                                                            ************
;******************************************************************************
;******************************************************************************

  print,''
  print,''
  print,''
  print,'**********************'
  print,'*      Bband LF      *'
  print,'**********************'
  print,''

  multiplot, /reset
  erase & multiplot, [2,2]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                      ;;;
;;;    Bband LF Z=0      ;;;
;;;                      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;bz0
  
  G=G0
  G_MRII=G0_MRII

  char_redshift=STRTRIM(number_formatter(redshift[0], DECIMAL=2),2)

  bin=0.5
  xmin=-25.0
  xmax=-16.1
if(MRII) then xmax=-14.1
  ymin=-5.9
  ymax=-0.5
  if(MRII) then ymax =0.5
 
 
  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, $
        ytitle = '!4U!3 (h!U3!NMpc!U-3!Nmag!U-1!N)'
 
  symbols, 2,  joint_obs_symbols_size
  if (plot_obs eq 0) then begin
      oploterror, [-24.5,-24.5], [yobs2+0.04,yobs2+0.04], [0.15,0.15], $
                  psym=8, color = 200, errcolor=200, HATLENGTH = 80.0
      xyouts, -24.25, yobs22, 'Observations used in MCMC',charsize=1.2
   endif

;check sample 
if(sample eq 1) then begin
  readcol,sample_file+'2_z0.10.txt',kband,phi,err,model, /SILENT
  oplot, kband, Alog10(model), color=0, thick=6
endif

;NEW MODEL
 mag_MR=G.MagDust[1]- 0.28 * (G.MagDust[1]-G.MagDust[2])
 mag_MRII=G_MRII.MagDust[1]- 0.28 * (G_MRII.MagDust[1]-G_MRII.MagDust[2])
 chi_square_and_plot, mcmc_folder, 'bband', char_redshift, total_chi2, plot_obs,  plot_only_new_obs,G.CentralMvir, G_MRII.CentralMvir, $
                      mag_MR-5.*Alog10(hubble_h), mag_MRII-5.*Alog10(hubble_h), $
                      volume_MR, volume_MRII, Alog10(StellarMass_MR_0), Alog10(StellarMass_MRII_0),Hubble_h, $
                      write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                      linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                      linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII

;PREVIOUS MODELS 
  if(do_previous_model1 eq 1) then begin
     readcol,file_previous_model1+'_bband_z'+char_redshift+'.txt',mag,phi, /SILENT 
     oplot,mag,phi, color=70,linestyle=linestyle_previous_model1, thick=6
  endif
  if(do_previous_model2 eq 1) then begin
     readcol,file_previous_model2+'_bband_z'+char_redshift+'.txt',mag,phi, /SILENT 
     oplot,mag,phi, color=70,linestyle=linestyle_previous_model2, thick=6
  endif


 if (plot_obs eq 1) then begin


;NORBERG

;if (MRII) then begin
  symbols, 30, 0.25
  readcol, Datadir+'norberg2002_bband_z0.txt',fullMagsbJ, fullbJ, fullErrbJ,exp, /SILENT
  oploterror, fullMagsbJ-0.083, Alog10(fullbJ*exp), $
              Alog10(fullbJ*exp+fullErrbJ*exp)-Alog10(fullbJ*exp), $
              color = 120, errcolor = 120, psym = 8,HATLENGTH = 80.0  

  symbols, 30, 0.25
   oploterror, [-24.5,-24.5], [yobs2,yobs2], [0.1,0.1], $
          psym=8, color=120, errcolor=120, HATLENGTH = 80.0
   xyouts, -24.25, yobs22, 'Norberg (2002)',charsize=0.9,charthick = 3   

;endif

;LABELS

  readcol,Datadir+'jones2006_bbandLF.dat',mag,phi,err, /SILENT
  symbols, 2, 0.7
 
  error=phi-err
  log_err=Alog10(error)
  sel=where(error eq 0.,n)
  if(n gt 0) then  log_err[sel] =0.
  
  oploterror, mag, Alog10(phi), Alog10(phi)-log_err, $ 
              color = 0, errcolor = 0, psym = 8,/loba,HATLENGTH = 80.0

  symbols, 2, 0.7
  oploterror, [-24.5,-24.5], [yobs1,yobs1], [0.1,0.1], $
              psym=8, color=0, errcolor=0, HATLENGTH = 80.0
  xyouts, -24.25, yobs12, 'Jones (2006)',charsize=0.9,charthick = 3   
  
endif 
 

 plot_label, xlog=0, ylog=0, type='label', label='b!Dj!N-band, z=0.0' , xmin, xmax, ymin, ymax, $
             x_percentage=5.5, x2_percentage=0., y_percentage=4., $
             charsize=1.3, charthick=4




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                      ;;;
;;;    Bband LF Z=1      ;;;
;;;                      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   G=G1
   G_MRII=G1_MRII
 
   char_redshift=STRTRIM(number_formatter(redshift[2], DECIMAL=2),2)

   multiplot
  
   plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
         xstyle = 1, ystyle = 1
 


;Z=0 LF Jones 2006  
 
  readcol,Datadir+'jones2006_bbandLF.dat',mag,phi,err, /SILENT
  symbols, 2, 0.7 
  oplot, mag, Alog10(phi), linestyle=1


;NEW MODEL 
 mag_MR=G.MagDust[1]
 mag_MRII=G_MRII.MagDust[1]
 chi_square_and_plot, mcmc_folder, 'bband', char_redshift, total_chi2, plot_obs,  plot_only_new_obs,G.CentralMvir, G_MRII.CentralMvir, $
                      mag_MR-5.*Alog10(hubble_h), mag_MRII-5.*Alog10(hubble_h), $
                      volume_MR, volume_MRII, Alog10(StellarMass_MR_1), Alog10(StellarMass_MRII_1),Hubble_h, $
                      write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                      linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                      linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII

;PREVIOUS MODELS 
  if(do_previous_model1 eq 1) then begin
     readcol,file_previous_model1+'_bband_z'+char_redshift+'.txt',mag,phi, /SILENT 
     oplot,mag,phi, color=70,linestyle=linestyle_previous_model1, thick=6
  endif
  if(do_previous_model2 eq 1) then begin
     readcol,file_previous_model2+'_bband_z'+char_redshift+'.txt',mag,phi, /SILENT 
     oplot,mag,phi, color=70,linestyle=linestyle_previous_model2, thick=6
  endif



 
  plot_label, xlog=0, ylog=0, type='label', label='B-band, z=1.0' , xmin, xmax, ymin, ymax, $
              x_percentage=5.5, x2_percentage=0., y_percentage=4., $
              charsize=1.3, charthick=4

;check sample 
if(sample eq 1) then begin
  readcol,sample_file+'2_z1.00.txt',kband,phi,err,model, /SILENT
  oplot, kband, Alog10(model), color=0, thick=6
endif



;INDIVIDUAL OBSERVATIONS

  if (plot_obs eq 1) then begin

  readcol,Datadir+'zucca_bband_z0810.txt',mag,phi,errordown,errorup, /SILENT
  symbols, 32, 0.3 
  hubble_h_not=1.043
  data_hubble=0.7
  phi=phi
  oploterror, mag-5.*Alog10(data_hubble), phi-Alog10(data_hubble^3), errorup, $
              color = 0, errcolor =0, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-5.*Alog10(data_hubble), phi-Alog10(data_hubble^3), errordown, $
              color = 0, errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0

  loadct,2, /SILENT
  readcol,Datadir+'ilbert_bband_z0810.txt',mag,phi,errordown,errorup, /SILENT
  symbols, 1, 0.7 
  oploterror, mag, phi, errorup, color = 200, $
                errcolor = 200, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag, phi, errordown, color = 200, $
                errcolor = 200, psym = 8, /loba,HATLENGTH = 80.0
  loadct,6, /SILENT

  loadct,2, /SILENT
  readcol,Datadir+'wilmer_bband_z0810.txt',mag,phi,errordown,errorup, /SILENT 
  symbols, 20, 0.3 
  data_hubble=0.7
  phi=phi-Alog10(data_hubble^3)
  oploterror, mag-0.0826-5.*Alog10(data_hubble), phi, errorup, color = 0, $
              errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-0.0826-5.*Alog10(data_hubble), phi, errordown, color = 0, $
              errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
   loadct,6, /SILENT


  readcol,Datadir+'giallongo_bband_z0710.txt',mag,phi,xerrdown, xerrup, errordown,errorup, /SILENT 
  symbols, 30, 0.25
  data_hubble=0.7
  phi=phi-Alog10(data_hubble^3)
  oploterror, mag-5.*Alog10(data_hubble), phi, xerrup,errorup, color = 120, $
                errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-5.*Alog10(data_hubble), phi,xerrdown,errordown, color = 120, $
                errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0
 
 



  symbols, 32, 0.3 
   oploterror, [-24.5,-24.5], [yobs1,yobs1], [0.1,0.1], $
               psym=8, color=0, errcolor=0, HATLENGTH = 80.0
   xyouts, -24.25, yobs12, 'Zucca (2009) - 0.8<z<1.0',charsize=0.9,charthick = 3   

   loadct,2, /SILENT
   symbols, 1, 0.7
   oploterror, [-24.5,-24.5], [yobs2,yobs2], [0.1,0.1], $
          psym=8, color=200, errcolor=200, HATLENGTH = 80.0
   xyouts, -24.25, yobs22, 'Ilbert (2005) - 0.8<z<1.0',charsize=0.9,charthick = 3   

   symbols, 22, 0.3 
   loadct,6, /SILENT
   
   symbols, 30, 0.25
   oploterror, [-24.5,-24.5], [yobs3,yobs3], [0.1,0.1], $
          psym=8, color=120, errcolor=120, HATLENGTH = 80.0
   xyouts, -24.25, yobs32, 'Giallongo (2005) - 0.7<z<1.0',charsize=0.9,charthick = 3   
 
   loadct,2, /SILENT
   symbols, 20, 0.3 
   oploterror, [-24.5,-24.5], [yobs4,yobs4], [0.1,0.1], $
          psym=8, color=0, errcolor=0, HATLENGTH = 80.0
   xyouts, -24.25, yobs42, 'Wilmer (2006) - 0.8<z<1.0',charsize=0.9,charthick = 3   
   loadct,6, /SILENT

endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                      ;;;
;;;    Bband LF Z=2      ;;;
;;;                      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   G=G2
   G_MRII=G2_MRII

   char_redshift=STRTRIM(number_formatter(redshift[3], DECIMAL=2),2)

  multiplot
 
  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, $
        xtitle = 'M!DB!N(AB)-5log!D10!N(h)',ytitle = '!4U!3 (h!U3!NMpc!U-3!Nmag!U-1!N)'
 


;Z=0 LF Jones 2006  
 
  readcol,Datadir+'jones2006_bbandLF.dat',mag,phi,err, /SILENT
  symbols, 2, 0.7 
  oplot, mag, Alog10(phi), linestyle=1

;NEW MODEL
 mag_MR=G.MagDust[1]
 mag_MRII=G_MRII.MagDust[1]
 chi_square_and_plot, mcmc_folder, 'bband', char_redshift, total_chi2, plot_obs,  plot_only_new_obs,G.CentralMvir, G_MRII.CentralMvir, $
                      mag_MR-5.*Alog10(hubble_h), mag_MRII-5.*Alog10(hubble_h), $
                      volume_MR, volume_MRII, Alog10(StellarMass_MR_2), Alog10(StellarMass_MRII_2),Hubble_h, $
                      write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                      linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                      linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII



;PREVIOUS MODELS 
  if(do_previous_model1 eq 1) then begin
     readcol,file_previous_model1+'_bband_z'+char_redshift+'.txt',mag,phi, /SILENT
     oplot,mag,phi, color=70,linestyle=linestyle_previous_model1, thick=6
  endif
  if(do_previous_model2 eq 1) then begin
     readcol,file_previous_model2+'_bband_z'+char_redshift+'.txt',mag,phi, /SILENT 
     oplot,mag,phi, color=70,linestyle=linestyle_previous_model2, thick=6
  endif


  plot_label, xlog=0, ylog=0, type='label', label='B-band, z=2.0' , xmin, xmax, ymin, ymax, $
              x_percentage=5.5, x2_percentage=0., y_percentage=4., $
              charsize=1.3, charthick=4
 
;check sample 
if(sample eq 1) then begin
  readcol,sample_file+'2_z2.00.txt',kband,phi,err,model, /SILENT
  oplot, kband, Alog10(model), color=0, thick=6
endif




  if (plot_obs eq 1) then begin

;INDIVIDUAL OBSERVATIONS
     loadct,2, /SILENT
  readcol,Datadir+'ilbert_bband_z1320.txt',mag,phi,errordown,errorup, /SILENT
  symbols, 1, 0.7 
 
  oploterror, mag, phi, errorup, color = 200, $
                errcolor = 200, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag, phi, errordown, color = 200, $
                errcolor = 200, psym = 8, /loba,HATLENGTH = 80.0
  loadct,6, /SILENT

  readcol,Datadir+'giallongo_bband_z1325.txt',mag,phi,xerrdown, xerrup, errordown,errorup, /SILENT
  symbols, 30, 0.25 
  oploterror, mag-5*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3), xerrup,errorup, color = 120, $
                errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-5*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3),xerrdown, errordown, color = 120, $
                errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0


  readcol,Datadir+'poli_bband_z1325.txt',mag,phi,xerrdown, xerrup,errordown,errorup, /SILENT
  symbols, 22, 0.3
  oploterror, mag-5*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3),xerrup, errorup, color = 0, $
                errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-5*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3), xerrdown,errordown, color = 0, $
                errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0


 loadct,2, /SILENT
   symbols, 1, 0.7 
   oploterror, [-24.5,-24.5], [yobs1,yobs1], [0.1,0.1], $
               psym=8, color=200, errcolor=200, HATLENGTH = 80.0
   xyouts, -24.25, yobs12, 'Ilbert (2005) - 1.3<z<2.0',charsize=0.9,charthick = 3   
loadct,6, /SILENT

   symbols, 30, 0.25
   oploterror, [-24.5,-24.5], [yobs2,yobs2], [0.1,0.1], $
          psym=8, color=120, errcolor=120, HATLENGTH = 80.0
   xyouts, -24.25, yobs22, 'Giallongo (2005) - 1.3<z<2.5',charsize=0.9,charthick = 3   


   symbols, 22, 0.3  
   oploterror, [-24.5,-24.5], [yobs3,yobs3], [0.1,0.1], $
          psym=8, color=0, errcolor=0, HATLENGTH = 80.0
   xyouts, -24.25, yobs32, 'Poli (2003) - 1.3<z<2.5',charsize=0.9,charthick = 3   
  

endif














;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                      ;;;
;;;    Bband LF Z=3      ;;;
;;;                      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  G=G3
  G_MRII=G3_MRII
 
  char_redshift=STRTRIM(number_formatter(redshift[4], DECIMAL=2),2)

  multiplot
  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, $
        xtitle = 'M!DB!N(AB)-5log!D10!N(h)'


;Z=0 LF Jones 2006  
 
  readcol,Datadir+'jones2006_bbandLF.dat',mag,phi,err, /SILENT
  symbols, 2, 0.7 
  oplot, mag, Alog10(phi), linestyle=1



;NEW MODEL
 mag_MR=G.MagDust[1]
 mag_MRII=G_MRII.MagDust[1]
 chi_square_and_plot, mcmc_folder, 'bband', char_redshift, total_chi2, plot_obs,  plot_only_new_obs,G.CentralMvir, G_MRII.CentralMvir, $
                      mag_MR-5.*Alog10(hubble_h), mag_MRII-5.*Alog10(hubble_h), $
                      volume_MR, volume_MRII, Alog10(StellarMass_MR_3), Alog10(StellarMass_MRII_3),Hubble_h, $
                      write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                      linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                      linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII

;PREVIOUS MODELS 
  if(do_previous_model1 eq 1) then begin
     readcol,file_previous_model1+'_bband_z'+char_redshift+'.txt',mag,phi, /SILENT 
     oplot,mag,phi, color=70,linestyle=linestyle_previous_model1, thick=6
  endif
  if(do_previous_model2 eq 1) then begin
     readcol,file_previous_model2+'_bband_z'+char_redshift+'.txt',mag,phi, /SILENT 
     oplot,mag,phi, color=70,linestyle=linestyle_previous_model2, thick=6
  endif





;LABELS
;redshift
   plot_label, xlog=0, ylog=0, type='label', label='B-band, z=3.0' , xmin, xmax, ymin, ymax, $
               x_percentage=5.5, x2_percentage=0., y_percentage=4., $
               charsize=1.3, charthick=4
 
;observations
   if(plot_obs eq 1) then begin
      plot_label, xlog=0, ylog=0, type='label', label='Observations for MCMC' , xmin, xmax, ymin, ymax, $
                  x_percentage=4.0, x2_percentage=0., y_percentage=2.8 , $
                  charsize=1.2, charthick=3 
      plot_label, xlog=0, ylog=0, type='symbol', xmin, xmax, ymin, ymax, $
                  x_percentage=3.8, x2_percentage=0., y_percentage=3.0, $
                  color=200, errcolor=200, sym_num=2, sym_size=joint_obs_symbols_size, HATLENGTH = 80, err_size=0.15
   endif
;models

   plot_label, xlog=0, ylog=0, type='label', label=prefix_this_model, xmin, xmax, ymin, ymax, $
               x_percentage=4.0, x2_percentage=0., y_percentage=2., $
               charsize=1.1, charthick=3
   plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
               x_percentage=3.2, x2_percentage=3.9, y_percentage=2.15, $
               linestyle=0, color=70, linethick=8
 
   if(do_previous_model2 eq 1) then begin
      plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model2, xmin, xmax, ymin, ymax, $
                  x_percentage=4.0, x2_percentage=0., y_percentage=1.3, $
                  charsize=1.1, charthick=3
      plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                  x_percentage=3.2, x2_percentage=3.9, y_percentage=1.45, $
                  linestyle=linestyle_previous_model2, color=70, linethick=8
   endif
   if(do_previous_model1 eq 1) then begin
      plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model1, xmin, xmax, ymin, ymax, $
                  x_percentage=4.0, x2_percentage=0., y_percentage=0.6, $
                  charsize=1.1, charthick=3
      plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                  x_percentage=3.2, x2_percentage=3.9, y_percentage=0.75, $
                  linestyle=linestyle_previous_model1, color=70, linethick=8
   endif



;check sample 
if(sample eq 1) then begin
  readcol,sample_file+'2_z3.00.txt',kband,phi,err,model, /SILENT
  oplot, kband, Alog10(model), color=0, thick=6
endif


 if (plot_obs eq 1) then begin

  loadct,2, /SILENT
  readcol,Datadir+'marchesini_bband_z2535.txt',mag,phi,xerrdown, xerrup,errordown,errorup, /SILENT
  symbols, 32, 0.3
  oploterror, mag-5.*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3), xerrup, errorup, color = 200, $
                errcolor = 200, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-5.*Alog10(Hubble_h_WMAP7), phi-Alog10(Hubble_h_WMAP7^3), xerrdown, errordown, color = 200, $
                errcolor = 200, psym = 8, /loba,HATLENGTH = 80.0
  loadct,6, /SILENT
 


  readcol,Datadir+'poli_bband_z2535.txt',mag,phi,xerrdown, xerrup,errordown,errorup, /SILENT
  symbols, 22, 0.3 
  oploterror, mag-5.*Alog10(Hubble_h_WMAP7), phi-3.*Alog10(Hubble_h_WMAP7),xerrup, errorup, color = 0, $
                errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-5.*Alog10(Hubble_h_WMAP7), phi-3.*Alog10(Hubble_h_WMAP7), xerrdown,errordown, color = 0, $
                errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0


  readcol,Datadir+'giallongo_bband_z2535.txt',mag,phi,xerrdown, xerrup, errordown,errorup, /SILENT
  symbols, 30, 0.25
  data_hubble=0.7 
  oploterror, mag-5.*Alog10(Hubble_h_WMAP7), phi-3.*Alog10(Hubble_h_WMAP7), xerrup,errorup, color = 120, $
                errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
  oploterror, mag-5.*Alog10(Hubble_h_WMAP7), phi-3.*Alog10(Hubble_h_WMAP7),xerrdown, errordown, color = 120, $
                errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0

  loadct,2, /SILENT
  symbols, 32, 0.3
  oploterror, [-24.5,-24.5], [yobs1,yobs1], [0.1,0.1], $
              psym=8, color=200, errcolor=200, HATLENGTH = 80.0
  xyouts, -24.25, yobs12, 'Marchesini (2007) - 2.5<z<3.5',charsize=0.9,charthick = 3   
  loadct,6, /SILENT

  symbols, 22, 0.3
  oploterror, [-24.5,-24.5], [yobs2,yobs2], [0.1,0.1], $
              psym=8, color=0, errcolor=0, HATLENGTH = 80.0
  xyouts, -24.25, yobs22, 'Poli (2003) - 2.5<z<3.5',charsize=0.9,charthick = 3   

  symbols, 30, 0.25
  oploterror, [-24.5,-24.5], [yobs3,yobs3], [0.1,0.1], $
         psym=8, color=120, errcolor=120, HATLENGTH = 80.0
  xyouts, -24.25, yobs32, 'Giallongo (2005) - 2.5<z<3.5',charsize=0.9,charthick = 3   

 
endif  ;plot_obs

endif   ;do_lums



;******************************************************************************
;******************************************************************************
;******                                                            ************
;******               Stellar Mass FUNCTION                        ************
;******                                                            ************
;******************************************************************************
;******************************************************************************

  print,''
  print,''
  print,''
  print,'*****************'
  print,'*      SMF      *'
  print,'*****************'
  print,''

  hubble_h_not=1.043
  multiplot, /reset
  erase & multiplot, [2,2]



;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                      ;;;
;;;    SMF Z=0      ;;;
;;;                      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  G=G0
  G_MRII=G0_MRII

  char_redshift=STRTRIM(number_formatter(redshift[0], DECIMAL=2),2)

  bin=0.25
  xmin=9.0
  xmax=12.4
  if (MRII eq 1) then  xmin=7.0
  ymin=-5.9
  ymax=-0.5
   if(MRII) then ymax =0.5
  
 
  plot, findgen(10), /nodata, xrange = [xmax,xmin], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, $
        ytitle = 'log!D10!N(!4U!3 [h!U3!NMpc!U-3!Nlog!D10!NM!U-1!N])'
 

;JOINT OBS 

   symbols, 2,  joint_obs_symbols_size
   if (plot_obs eq 0) then begin
      oploterror, [12.15,12.15], [yobs2+0.04,yobs2+0.04], [0.15,0.15], $
                  psym=8, color = 200, errcolor=200, HATLENGTH = 200.0
      xyouts, 12., yobs22, 'Observations used in MCMC'
   endif

   if(do_lums eq 0) then  total_chi2=0


;NEW MODEL
 chi_square_and_plot, mcmc_folder, 'smf', char_redshift, total_chi2, plot_obs,  plot_only_new_obs,G.CentralMvir, G_MRII.CentralMvir, $
                      Alog10(StellarMass_MR_0), Alog10(StellarMass_MRII_0), $
                      volume_MR, volume_MRII, Alog10(StellarMass_MR_0), Alog10(StellarMass_MRII_0),Hubble_h, $
                      write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                      linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                      linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII

;PREVIOUS MODELS 
  if(do_previous_model1 eq 1) then begin
     readcol,file_previous_model1+'_smf_z'+char_redshift+'.txt',mag,phi, /SILENT 
     oplot,mag,phi, color=70,linestyle=linestyle_previous_model1, thick=6
  endif
  if(do_previous_model2 eq 1) then begin
     readcol,file_previous_model2+'_smf_z'+char_redshift+'.txt',mag,phi, /SILENT 
     oplot,mag,phi, color=70,linestyle=linestyle_previous_model2, thick=6
  endif

;check sample 
if(sample eq 1) then begin
  readcol,sample_file+'0_z0.10.txt',mass,phi,err,model, /SILENT
  oplot, mass, Alog10(model), color=0, thick=6
endif


  if (plot_obs eq 1) then begin 
   
     data_hubble=0.7

     symbols, 30, 0.2
     Readcol,Datadir+'baldry2008_smf.txt',obsmass,num,fi,err,err2,err3, /SILENT 
     sel=where(obsmass lt 11.6)
     obsmass=obsmass[sel]
     num=num[sel]
     fi=fi[sel]
     err=err[sel]
     err2=err2[sel]
     err3=err3[sel]

     oploterror, obsmass+2.*Alog10(data_hubble), Alog10(fi)-3.*Alog10(data_hubble), (obsmass/obsmass)-1.0+0.05, Alog10(fi+err)-Alog10(fi), $ 
                 color = 120, errcolor = 120, psym = 8,/hiba,HATLENGTH = 80.0
     
     error=fi-err
     log_err=Alog10(error)
     sel=where(error eq 0.,n)
     if(n gt 0) then  log_err[sel] =0.
     
     oploterror, obsmass+2.*Alog10(data_hubble), Alog10(fi)-3.*Alog10(data_hubble), (obsmass/obsmass)-1.0+0.05, Alog10(fi)-log_err, $ 
                 color = 120, errcolor = 120, psym = 8,/loba,HATLENGTH = 80.0
     
 
     data_hubble=0.732
     symbols, 2, 0.7
     readcol,Datadir+'LiWhite2009_SMF.txt',cmass,mass2,mass3,cphi,cerror, /SILENT
     
     oploterror,mass3,Alog10(cphi),Alog10((cphi+cerror)/cphi), $
                psym=8,color = 0, errcolor = 0,HATLENGTH = 80.0
     
     
     symbols, 2, 0.7 
     oploterror, [12.2,12.2], [yobs1,yobs1], [0.1,0.1], $
                 psym=8, color=0, errcolor=0, HATLENGTH = 80.0
     xyouts, 12.1, yobs12, 'Li & White (2009)',charsize=0.9,charthick = 3   
     
     symbols, 30, 0.25
     oploterror, [12.2,12.2], [yobs2,yobs2], [0.1,0.1], $
                 psym=8, color=120, errcolor=120, HATLENGTH = 80.0
     xyouts, 12.1, yobs22, 'Baldry (2008)',charsize=0.9,charthick = 3   
     

     loadct,2, /SILENT
     symbols, 30, 0.25
     data_hubble=0.7
     readcol,Datadir+'baldry2012_smf.txt',obs_mass,obs_bin,obs_phi,obs_err,obs_num, /SILENT
     phi=obs_phi*1.e-3
     phi_up=(obs_phi+obs_err)*1.e-3
     phi_down=(obs_phi-obs_err)*1.e-3   
     oploterror,obs_mass+2.*Alog10(data_hubble),Alog10(phi)-3.*Alog10(data_hubble),Alog10(phi_up/phi), $
                psym=8,color = 200, errcolor = 200,HATLENGTH = 80.0,/hiba
     oploterror,obs_mass+2.*Alog10(data_hubble),Alog10(phi)-3.*Alog10(data_hubble),Alog10(phi/phi_down), $
                psym=8,color = 200, errcolor = 200,HATLENGTH = 80.0,/loba
      oploterror, [12.2,12.2], [yobs3,yobs3], [0.1,0.1], $
                 psym=8, color=200, errcolor=200, HATLENGTH = 80.0
     xyouts, 12.1, yobs32, 'Baldry (2012)',charsize=0.9,charthick = 3   
     
loadct,6, /SILENT
     


;overplot combined obs constraints
  close,1
  openr,1,mcmc_folder+'ObsConstraints/StellarMassFunction_z0.00.txt'
  Nbins = 0L & readf,1,Nbins
  data=FLTARR(4,Nbins)
  readf,1,data
  close,1  
  data_bin_low=data[0,*]
  data_bin_high=data[1,*]
  data_obs=data[2,*]
  data_err=data[3,*] 

  data_bin=data_bin_low+(data_bin_high-data_bin_low)/2.

  sel=where(data_err/data_obs lt min_obs_err,n)
  if(n gt 0) then  data_err[sel]=data_obs[sel]*min_obs_err

  if(data_bin[1]-data_bin[0] gt 0.) then begin 
     Bin = data_bin[1]-data_bin[0]
  endif else begin 
     Bin = data_bin[0]-data_bin[1]
  endelse

;overplot OBS constraints  
  if(plot_obs eq 0) then $
     joint_obs_symbols_size=1.0 $
  else joint_obs_symbols_size=0.7

  if(plot_only_new_obs eq 1) then joint_obs_symbols_size=1.0

  symbols, 2,  joint_obs_symbols_size
  oploterror, data_bin, Alog10(data_obs), Alog10((data_obs+data_err)/data_obs), $
              /hiba, psym=8, color=200, errcolor=200, HATLENGTH = 200.0  , ERRTHICK=4
  oploterror, data_bin, Alog10(data_obs), Alog10(data_obs/(data_obs-data_err)), $
              /loba, psym=8, color=200, errcolor=200, HATLENGTH = 200.0 , ERRTHICK=4

  endif



  plot_label, xlog=0, ylog=0, type='label', label='SMF, z=0.0' , xmin, xmax, ymin, ymax, $
              x_percentage=4.5, x2_percentage=0., y_percentage=4.



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                      ;;;
;;;    SMF Z=1           ;;;
;;;                      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  G=G1
  G_MRII=G1_MRII

  char_redshift=STRTRIM(number_formatter(redshift[2], DECIMAL=2),2)

  multiplot
 
  plot, findgen(10), /nodata, xrange = [xmax,xmin], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1
 
  


;Li&White - z=0
  data_hubble=0.732

  symbols, 2, 0.7
  readcol,Datadir+'LiWhite2009_SMF.txt',cmass,mass2,mass3,cphi,cerror, /SILENT
  oplot,mass3,Alog10(cphi),linestyle=1


;NEW MODEL
   chi_square_and_plot, mcmc_folder, 'smf', char_redshift, total_chi2, plot_obs,  plot_only_new_obs, G.CentralMvir, G_MRII.CentralMvir, $
                        Alog10(StellarMass_MR_1), $
                        Alog10(StellarMass_MRII_1), $
                        volume_MR, volume_MRII, Alog10(StellarMass_MR_1), Alog10(StellarMass_MRII_1),Hubble_h, $
                        write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                        linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                        linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII

;PREVIOUS MODELS 
   if(do_previous_model1 eq 1) then begin
      readcol,file_previous_model1+'_smf_z'+char_redshift+'.txt',mag,phi, /SILENT 
      oplot,mag,phi, color=70,linestyle=linestyle_previous_model1, thick=6
   endif
   if(do_previous_model2 eq 1) then begin
      readcol,file_previous_model2+'_smf_z'+char_redshift+'.txt',mag,phi, /SILENT 
      oplot,mag,phi, color=70,linestyle=linestyle_previous_model2, thick=6
   endif

;check sample 
if(sample eq 1) then begin
  readcol,sample_file+'0_z1.00.txt',kband,phi,err,model, /SILENT
  oplot, kband, Alog10(model), color=0, thick=6
endif


 if (plot_obs eq 1) then begin    

;ILBERT2010
    loadct,2, /SILENT
    readcol,Datadir+'ilbert2010_smf_z0810.txt',mass,phi,errordown,errorup, /SILENT
    symbols, 20, 0.25
    data_hubble=0.7 
    oploterror, mass-0.14+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),errorup, color = 200, $
                errcolor = 200, psym = 8, /hiba,HATLENGTH = 80.0
    oploterror, mass-0.14+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),errordown, color = 200, $
                errcolor = 200, psym = 8, /loba,HATLENGTH = 80.0
    oploterror, [12.2,12.2], [yobs1,yobs1], [0.1,0.1], $
                psym=8, color=200, errcolor=200, HATLENGTH = 80.0
    xyouts, 12.1, yobs12, 'Ilbert (2010) - 0.8<z<1.0',charsize=0.9,charthick = 3   

    readcol,Datadir+'ilbert2010_smf_z1012.txt',mass,phi,errordown,errorup, /SILENT
    symbols, 30, 0.25   
    oploterror, mass-0.14+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),errorup, color = 200, $
                errcolor = 200, psym = 8, /hiba,HATLENGTH = 80.0
    oploterror, mass-0.14+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),errordown, color = 200, $
                errcolor = 200, psym = 8, /loba,HATLENGTH = 80.0
   
    oploterror, [12.2,12.2], [yobs2,yobs2], [0.1,0.1], $
                psym=8, color=200, errcolor=200, HATLENGTH = 80.0
    xyouts, 12.1, yobs22, 'Ilbert (2010) - 1.0<z<1.2',charsize=0.9,charthick = 3   
    loadct,6, /SILENT
    

;TOMCZAK 2013
    readcol,Datadir+'tomczak2013_smf_z0.75_1.00.txt',obs_mass, obs_phi, obs_error_up, obs_error_down, /SILENT
    symbols, 22, 0.3
    data_hubble=0.7 
    oploterror, obs_mass-0.14+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_up, color = 70, $
                errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
    oploterror, obs_mass-0.14+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_down, color = 70, $
                errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
    oploterror, [12.2,12.2], [yobs3,yobs3], [0.1,0.1], $
                psym=8, color=70, errcolor=70, HATLENGTH = 80.0
    xyouts, 12.1, yobs32, 'Tomczak (2014) - 0.75<z<1.0',charsize=0.9,charthick = 3 


    readcol,Datadir+'tomczak2013_smf_z1.00_1.25.txt',obs_mass, obs_phi, obs_error_up, obs_error_down, /SILENT
    symbols, 32, 0.3
    data_hubble=0.7  
    oploterror, obs_mass-0.14+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_up, color = 70, $
                errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
    oploterror, obs_mass-0.14+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_down, color = 70, $
                errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
    oploterror, [12.2,12.2], [yobs4,yobs4], [0.1,0.1], $
                psym=8, color=70, errcolor=70, HATLENGTH = 80.0
    xyouts, 12.1, yobs42, 'Tomczak (2014) - 1.0<z<1.25',charsize=0.9,charthick = 3 


   
;ILBERT2013
    readcol,Datadir+'/ilbert2013/MF_Vmax_All_z0.8_1.1.dat',mass,phi,phi_down,phi_up, /SILENT
    min_mass=min(mass)
    readcol,Datadir+'/ilbert2013/MF_Vmax_All_z0.8_1.1.dat',mass,phi,err_down,err_up, /SILENT
    symbols, 31, 0.25
    data_hubble=0.7  
    sel=where(mass gt min_mass)
    oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
    oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0 
    
    oploterror, [12.2,12.2], [yobs5,yobs5], [0.1,0.1], $
                psym=8, color=120, errcolor=120, HATLENGTH = 80.0
    xyouts, 12.1, yobs52, 'Ilbert (2013) - 0.8<z<1.1',charsize=0.9,charthick = 3   
    
;MUZZIN2013
    readcol,Datadir+'muzzin2013_maraston_z0.5_1.0.txt',mass,x_err, phi,error_up,error_down, $
            passive_phi, passive_error_up, passive_error_down, $
            active_phi, active_error_up, active_error_down, /SILENT
    symbols, 1, 0.7
    data_hubble=0.7  
    oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),error_up, color = 0, $
                errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
    oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),error_down, color = 0, $
                errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0

    oploterror, [12.2,12.2], [yobs6,yobs6], [0.1,0.1], $
                psym=8, color=0, errcolor=0, HATLENGTH = 80.0
    xyouts, 12.1, yobs62, 'Muzzin (2013) - 0.5<z<1.0',charsize=0.9,charthick = 3   
    


    readcol,Datadir+'muzzin2013_maraston_z1.0_1.5.txt',mass,x_err, phi,error_up,error_down, $
            passive_phi, passive_error_up, passive_error_down, $
            active_phi, active_error_up, active_error_down, /SILENT
    symbols, 2, 0.7
    data_hubble=0.7  
    oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),error_up, color = 0, $
                errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
    oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),error_down, color = 0, $
                errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
    
    oploterror, [12.2,12.2], [yobs7,yobs7], [0.1,0.1], $
                psym=8, color=0, errcolor=0, HATLENGTH = 80.0
    xyouts, 12.1, yobs72, 'Muzzin (2013) - 1.0<z<1.5',charsize=0.9,charthick = 3   
    

 


;overplot combined obs constraints
  close,1
  openr,1,mcmc_folder+'ObsConstraints/StellarMassFunction_z1.00.txt'
  Nbins = 0L & readf,1,Nbins
  data=FLTARR(4,Nbins)
  readf,1,data
  close,1  
  data_bin_low=data[0,*]
  data_bin_high=data[1,*]
  data_obs=data[2,*]
  data_err=data[3,*] 

  data_bin=data_bin_low+(data_bin_high-data_bin_low)/2.

  sel=where(data_err/data_obs lt min_obs_err,n)
  if(n gt 0) then  data_err[sel]=data_obs[sel]*min_obs_err

  if(data_bin[1]-data_bin[0] gt 0.) then begin 
     Bin = data_bin[1]-data_bin[0]
  endif else begin 
     Bin = data_bin[0]-data_bin[1]
  endelse

;overplot OBS constraints  
  if(plot_obs eq 0) then $
     joint_obs_symbols_size=1.0 $
  else joint_obs_symbols_size=0.7

  if(plot_only_new_obs eq 1) then joint_obs_symbols_size=1.0

  symbols, 2,  joint_obs_symbols_size
  oploterror, data_bin, Alog10(data_obs), Alog10((data_obs+data_err)/data_obs), $
              /hiba, psym=8, color=200, errcolor=200, HATLENGTH = 200.0  , ERRTHICK=4
  oploterror, data_bin, Alog10(data_obs), Alog10(data_obs/(data_obs-data_err)), $
              /loba, psym=8, color=200, errcolor=200, HATLENGTH = 200.0 , ERRTHICK=4


    
 endif
 

   plot_label, xlog=0, ylog=0, type='label', label='SMF, z=1.0' , xmin, xmax, ymin, ymax, $
               x_percentage=4.5, x2_percentage=0., y_percentage=4.




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                      ;;;
;;;     SMF Z=2          ;;;
;;;                      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  G=G2
  G_MRII=G2_MRII

  char_redshift=STRTRIM(number_formatter(redshift[3], DECIMAL=2),2)

  multiplot
 
  plot, findgen(10), /nodata, xrange = [xmax,xmin], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, xtitle = 'log!D10!N(M!D*!N[h!U-2!NM!D!9ng!3!N])', $
        ytitle = 'log!D10!N(!4U!3 [h!U3!NMpc!U-3!Nlog!D10!NM!U-1!N])'
 

;Li&White - z=0
  symbols, 2, 0.7
  data_hubble=0.732
  readcol,Datadir+'LiWhite2009_SMF.txt',cmass,mass2,mass3,cphi,cerror, /SILENT
  oplot,mass3,Alog10(cphi),linestyle=1

;NEW MODEL
  chi_square_and_plot, mcmc_folder, 'smf', char_redshift, total_chi2, plot_obs,  plot_only_new_obs,G.CentralMvir, G_MRII.CentralMvir, $
                       Alog10(StellarMass_MR_2), Alog10(StellarMass_MRII_2), $
                       volume_MR, volume_MRII, Alog10(StellarMass_MR_2), Alog10(StellarMass_MRII_2),Hubble_h, $
                       write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                       linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                       linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII


;PREVIOUS MODELS 
  if(do_previous_model1 eq 1) then begin
     readcol,file_previous_model1+'_smf_z'+char_redshift+'.txt',mag,phi, /SILENT 
     oplot,mag,phi, color=70,linestyle=linestyle_previous_model1, thick=6
  endif
  if(do_previous_model2 eq 1) then begin
     readcol,file_previous_model2+'_smf_z'+char_redshift+'.txt',mag,phi, /SILENT 
     oplot,mag,phi, color=70,linestyle=linestyle_previous_model2, thick=6
  endif


;check sample 
  if(sample eq 1) then begin
     readcol,sample_file+'0_z2.00.txt',kband,phi,err,model, /SILENT
     oplot, kband, Alog10(model), color=0, thick=6
  endif


;INDIVIDUAL OBSERVATIONS

  if (plot_obs eq 1) then begin    

     ;SANCHEZ 2011
     readcol,Datadir+'sanchez2011_smf_z1620.txt',mass,phi,errordown,errorup, /SILENT
     symbols, 30, 0.3
     loadct,2, /SILENT
     data_hubble=0.7   
     oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble), errorup, color = 200, $
                 errcolor = 200, psym = 8, /hiba,HATLENGTH = 80.0
     oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble), errordown, color = 200, $
                 errcolor = 200, psym = 8, /loba,HATLENGTH = 80.0   
     oploterror, [12.2,12.2], [yobs1,yobs1], [0.1,0.1], $
                 psym=8, color=200, errcolor=200, HATLENGTH = 80.0
     xyouts, 12.1, yobs12, 'Sanchez (2011) - 1.6<z<2.0',charsize=0.9,charthick = 3   
     loadct,6, /SILENT


     ;TOMCZAK2013
     readcol,Datadir+'tomczak2013_smf_z1.50_2.00.txt',obs_mass, obs_phi, obs_error_up, obs_error_down, /SILENT
     symbols, 22, 0.3
     data_hubble=0.7     
     oploterror, obs_mass-0.14+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_up, color = 70, $
                 errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
     oploterror, obs_mass-0.14+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_down, color = 70, $
                 errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
     oploterror, [12.2,12.2], [yobs2,yobs2], [0.1,0.1], $
                psym=8, color=70, errcolor=70, HATLENGTH = 80.0
     xyouts, 12.1, yobs22, 'Tomczak (2014) - 0.75<z<1.0',charsize=0.9,charthick = 3 

     readcol,Datadir+'tomczak2013_smf_z2.00_2.50.txt',obs_mass, obs_phi, obs_error_up, obs_error_down, /SILENT
     symbols, 32, 0.3
     data_hubble=0.7   
     oploterror, obs_mass-0.14+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_up, color = 70, $
                 errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
     oploterror, obs_mass-0.14+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_down, color = 70, $
                 errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
     oploterror, [12.2,12.2], [yobs3,yobs3], [0.1,0.1], $
                 psym=8, color=70, errcolor=70, HATLENGTH = 80.0
     xyouts, 12.1, yobs32, 'Tomczak (2014) - 0.75<z<1.0',charsize=0.9,charthick = 3 
     
     ;ILBERT2013
     readcol,Datadir+'/ilbert2013/MF_Vmax_All_z1.5_2.0.dat',mass,phi,phi_down,phi_up, /SILENT
     min_mass=min(mass)
     readcol,Datadir+'/ilbert2013/MF_Vmax_All_z1.5_2.0.dat',mass,phi,err_down,err_up, /SILENT
     symbols, 21, 0.3
     data_hubble=0.7       
     sel=where(mass gt min_mass)
     oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                 errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
     oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                 errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0 
     oploterror, [12.2,12.2], [yobs4,yobs4], [0.1,0.1], $
                 psym=8, color=120, errcolor=120, HATLENGTH = 80.0
     xyouts, 12.1, yobs42, 'Ilbert (2013) - 1.5<z<2.0',charsize=0.9,charthick = 3   

     readcol,Datadir+'/ilbert2013/MF_Vmax_All_z2.0_2.5.dat',mass,phi,phi_down,phi_up, /SILENT
     min_mass=min(mass)
     readcol,Datadir+'/ilbert2013/MF_Vmax_All_z2.0_2.5.dat',mass,phi,err_down,err_up, /SILENT
     symbols, 31, 0.3
     data_hubble=0.7       
     oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                 errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
     oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                 errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0  
     oploterror, [12.2,12.2], [yobs5,yobs5], [0.1,0.1], $
                 psym=8, color=120, errcolor=120, HATLENGTH = 80.0
     xyouts, 12.1, yobs52, 'Ilbert (2013) - 2.0<z<2.5',charsize=0.9,charthick = 3   


     ;MUZZIN2013
     readcol,Datadir+'muzzin2013_maraston_z1.5_2.0.txt',mass,x_err, phi,error_up,error_down, $
             passive_phi, passive_error_up, passive_error_down, $
             active_phi, active_error_up, active_error_down, /SILENT
     symbols, 1, 0.7
     data_hubble=0.7     
     oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),error_up, color = 0, $
                 errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
     oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),error_down, color = 0, $
                 errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
     oploterror, [12.2,12.2], [yobs6,yobs6], [0.1,0.1], $
                 psym=8, color=0, errcolor=0, HATLENGTH = 80.0
     xyouts, 12.1, yobs62, 'Muzzin (2013) - 1.5<z<2.0',charsize=0.9,charthick = 3   
     
     readcol,Datadir+'muzzin2013_maraston_z2.0_2.5.txt',mass,x_err, phi,error_up,error_down, $
             passive_phi, passive_error_up, passive_error_down, $
             active_phi, active_error_up, active_error_down, /SILENT
     symbols, 2, 0.7
     data_hubble=0.7     
     oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),error_up, color = 0, $
                 errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
     oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),error_down, color = 0, $
                 errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0   
     oploterror, [12.2,12.2], [yobs7,yobs7], [0.1,0.1], $
                 psym=8, color=0, errcolor=0, HATLENGTH = 80.0
     xyouts, 12.1, yobs72, 'Muzzin (2013) - 2.0<z<2.5',charsize=0.9,charthick = 3   




;overplot combined obs constraints
  close,1
  openr,1,mcmc_folder+'ObsConstraints/StellarMassFunction_z2.00.txt'
  Nbins = 0L & readf,1,Nbins
  data=FLTARR(4,Nbins)
  readf,1,data
  close,1  
  data_bin_low=data[0,*]
  data_bin_high=data[1,*]
  data_obs=data[2,*]
  data_err=data[3,*]
  data_bin=data_bin_low+(data_bin_high-data_bin_low)/2.
 
  sel=where(data_err/data_obs lt min_obs_err,n)
  if(n gt 0) then  data_err[sel]=data_obs[sel]*min_obs_err

  if(data_bin[1]-data_bin[0] gt 0.) then begin 
     Bin = data_bin[1]-data_bin[0]
  endif else begin 
     Bin = data_bin[0]-data_bin[1]
  endelse

;overplot OBS constraints  
  if(plot_obs eq 0) then $
     joint_obs_symbols_size=1.0 $
  else joint_obs_symbols_size=0.7

  if(plot_only_new_obs eq 1) then joint_obs_symbols_size=1.0

  symbols, 2,  joint_obs_symbols_size
  oploterror, data_bin, Alog10(data_obs), Alog10((data_obs+data_err)/data_obs), $
              /hiba, psym=8, color=200, errcolor=200, HATLENGTH = 200.0  , ERRTHICK=4
  oploterror, data_bin, Alog10(data_obs), Alog10(data_obs/(data_obs-data_err)), $
              /loba, psym=8, color=200, errcolor=200, HATLENGTH = 200.0 , ERRTHICK=4


  endif


;OVERPLOT CURRENT MODEL
  ;readcol,file_to_write+'_smf_z2.00.txt',mag,phi, /SILENT 
  ;oplot,mag,phi, color=70,linestyle=0, thick=6

  plot_label, xlog=0, ylog=0, type='label', label='SMF, z=2.0' , xmin, xmax, ymin, ymax, $
               x_percentage=4.5, x2_percentage=0., y_percentage=4.






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                      ;;;
;;;       SMF Z=3        ;;;
;;;                      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  G=G3
  G_MRII=G3_MRII

  char_redshift=STRTRIM(number_formatter(redshift[4], DECIMAL=2),2)

  multiplot
 
  
   plot, findgen(10), /nodata, xrange = [xmax,xmin], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, xtitle = 'log!D10!N(M!D*!N[h!U-2!NM!D!9ng!3!N])'
  

;Li&White - z=0
  symbols, 2, 0.7
  data_hubble=0.732
  readcol,Datadir+'LiWhite2009_SMF.txt',cmass,mass2,mass3,cphi,cerror, /SILENT
  oplot,mass3,Alog10(cphi),linestyle=1


;NEW MODEL
   chi_square_and_plot, mcmc_folder, 'smf', char_redshift, total_chi2, plot_obs,  plot_only_new_obs,G.CentralMvir, G_MRII.CentralMvir, $
                        Alog10(StellarMass_MR_3), Alog10(StellarMass_MRII_3), $
                        volume_MR, volume_MRII, Alog10(StellarMass_MR_3), Alog10(StellarMass_MRII_3),Hubble_h, $
                        write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                        linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                        linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII


;PREVIOUS MODELS 
  if(do_previous_model1 eq 1) then begin
     readcol,file_previous_model1+'_smf_z'+char_redshift+'.txt',mag,phi, /SILENT 
     oplot,mag,phi, color=70,linestyle=linestyle_previous_model1, thick=6
  endif
  if(do_previous_model2 eq 1) then begin
     readcol,file_previous_model2+'_smf_z'+char_redshift+'.txt',mag,phi, /SILENT 
     oplot,mag,phi, color=70,linestyle=linestyle_previous_model2, thick=6
  endif


   if (plot_obs eq 1) then begin    
      ;SANCHEZ2011
      readcol,Datadir+'sanchez2011_smf_z2530.txt',mass,phi,errordown,errorup, /SILENT
      loadct,2, /SILENT
      symbols, 30, 0.3
      data_hubble=0.7    
      oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble), errorup, color = 200, $
                  errcolor = 200, psym = 8, /hiba,HATLENGTH = 80.0
      oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble), errordown, color = 200, $
                  errcolor = 200, psym = 8, /loba,HATLENGTH = 80.0
      oploterror, [12.2,12.2], [yobs1,yobs1], [0.1,0.1], $
                  psym=8, color=200, errcolor=200, HATLENGTH = 80.0
      xyouts, 12.1, yobs12, 'Sanchez (2011) - 2.5<z<3.0',charsize=0.9,charthick = 3   
      loadct,6, /SILENT

      ;MARCHESINI2009    
      readcol,Datadir+'marchesini_smf_z2030.txt',mass,phi,errup,errdown, $
              errorup_poi,errordown_poi, errorup_zran, errordown_zran, $
              errorup_cv, errordown_cv, errorup_zsys, errordown_zsys, $
              errorup_sed, errordown_sed, /SILENT      
      errorup=sqrt(errorup_poi^2+errorup_zran^2+errorup_cv^2)           ;+errorup_sed^2)
      errordown=sqrt(errordown_poi^2+errordown_zran^2+errordown_cv^2)   ;+errordown_sed^2)
      data_hubble=0.7
      symbols, 22, 0.3    
      oploterror, mass-0.14+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble), errorup, color = 70, $
                  errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
      oploterror, mass-0.14+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble), errordown, color = 70, $
                  errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0      
      oploterror, [12.2,12.2], [yobs2,yobs2], [0.1,0.1], $
                  psym=8, color=70, errcolor=70, HATLENGTH = 80.0
      xyouts, 12.1, yobs22, 'Marchesini (2009) - 2.0<z<3.0',charsize=0.9,charthick = 3   
      
      readcol,Datadir+'marchesini_smf_z3040.txt',mass,phi,errup,errdown, $
              errorup_poi,errordown_poi, errorup_zran, errordown_zran, $
              errorup_cv, errordown_cv, errorup_zsys, errordown_zsys, $
              errorup_sed, errordown_sed, /SILENT      
      errorup=sqrt(errorup_poi^2+errorup_zran^2+errorup_cv^2)           ;+errorup_sed^2)
      errordown=sqrt(errordown_poi^2+errordown_zran^2+errordown_cv^2)   ;+errordown_sed^2)      
      symbols, 32, 0.3    
      oploterror, mass-0.14+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble), errorup, color = 70, $
                  errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
      oploterror, mass-0.14+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble), errordown, color = 70, $
                  errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0     
      oploterror, [12.2,12.2], [yobs3,yobs3], [0.1,0.1], $
                  psym=8, color=70, errcolor=70, HATLENGTH = 80.0
      xyouts, 12.1, yobs32, 'Marchesini (2009) - 3.0<z<4.0',charsize=0.9,charthick = 3   
     

      ;MARCHESINI2010
      readcol,Datadir+'marchesini2010_smf_z3040.txt',mass,phi,errup,errdown, $
              errorup_poi,errordown_poi, errorup_zran, errordown_zran, $
              errorup_cv, errordown_cv, errorup_sed, errordown_sed, /SILENT      
      errorup=sqrt(errorup_poi^2+errorup_zran^2+errorup_cv^2)           ;+errorup_sed^2)
      errordown=sqrt(errordown_poi^2+errordown_zran^2+errordown_cv^2)   ;+errordown_sed^2)         
      symbols, 20, 0.3     
      oploterror, mass-0.14+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble), errorup, color = FSC_COLOR("orange"), $
                  errcolor = FSC_COLOR("orange"), psym = 8, /hiba,HATLENGTH = 80.0
      oploterror, mass-0.14+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble), errordown, color = FSC_COLOR("orange"), $
                  errcolor = FSC_COLOR("orange"), psym = 8, /loba,HATLENGTH = 80.0      
      oploterror, [12.2,12.2], [yobs4,yobs4], [0.1,0.1], $
                  psym=8, color=FSC_COLOR("orange"), errcolor=FSC_COLOR("orange"), HATLENGTH = 80.0
      xyouts, 12.1, yobs42, 'Marchesini (2010) - 3.0<z<4.0',charsize=0.9,charthick = 3   
      

      ;ILBERT2013
      readcol,Datadir+'/ilbert2013/MF_Vmax_All_z2.5_3.0.dat',mass,phi,phi_down,phi_up, /SILENT
      min_mass=min(mass)
      readcol,Datadir+'/ilbert2013/MF_Vmax_All_z2.5_3.0.dat',mass,phi,err_down,err_up, /SILENT
      symbols, 21, 0.25
      data_hubble=0.7        
      sel=where(mass gt min_mass)
      oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                  errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
      oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                  errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0 
      oploterror, [12.2,12.2], [yobs5,yobs5], [0.1,0.1], $
                  psym=8, color=120, errcolor=120, HATLENGTH = 80.0
      xyouts, 12.1, yobs52, 'Ilbert (2013) - 2.5<z<3.0',charsize=0.9,charthick = 3   
        
      readcol,Datadir+'/ilbert2013/MF_Vmax_All_z3.0_4.0.dat',mass,phi,phi_down,phi_up, /SILENT
      min_mass=min(mass)
      readcol,Datadir+'/ilbert2013/MF_Vmax_All_z3.0_4.0.dat',mass,phi,err_down,err_up, /SILENT
      symbols, 31, 0.25
      data_hubble=0.7       
      oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                  errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
      oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                  errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0  
      oploterror, [12.2,12.2], [yobs6,yobs6], [0.1,0.1], $
                  psym=8, color=120, errcolor=120, HATLENGTH = 80.0
      xyouts, 12.1, yobs62, 'Ilbert (2013) - 3.0<z<4.0',charsize=0.9,charthick = 3   




;MUZZIN2013
     readcol,Datadir+'muzzin2013_maraston_z2.5_3.0.txt',mass,x_err, phi,error_up,error_down, $
             passive_phi, passive_error_up, passive_error_down, $
             active_phi, active_error_up, active_error_down, /SILENT
     symbols, 1, 0.7
     data_hubble=0.7       
     oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),error_up, color = 0, $
                 errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
     oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),error_down, color = 0, $
                 errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
     
     readcol,Datadir+'muzzin2013_maraston_z3.0_4.0.txt',mass,x_err, phi,error_up,error_down, $
             passive_phi, passive_error_up, passive_error_down, $
             active_phi, active_error_up, active_error_down, /SILENT
     symbols, 2, 0.7
     data_hubble=0.7   
     oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),error_up, color = 0, $
                 errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
     oploterror, mass+2.*Alog10(data_hubble), phi-3.*Alog10(data_hubble),error_down, color = 0, $
                 errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0

        
        symbols, 1, 0.7
        oploterror, [12.2,12.2], [yobs7,yobs7], [0.1,0.1], $
                    psym=8, color=0, errcolor=0, HATLENGTH = 80.0
        xyouts, 12.1, yobs72, 'Muzzin (2013) - 2.5<z<3.0',charsize=0.9,charthick = 3   
        
        symbols, 2, 0.7
        oploterror, [12.2,12.2], [yobs8,yobs8], [0.1,0.1], $
                    psym=8, color=0, errcolor=0, HATLENGTH = 80.0
        xyouts, 12.1, yobs82, 'Muzzin (2013) - 3.0<z<4.0',charsize=0.9,charthick = 3   




;overplot combined obs constraints
  close,1
  openr,1,mcmc_folder+'ObsConstraints/StellarMassFunction_z3.00.txt'
  Nbins = 0L & readf,1,Nbins
  data=FLTARR(4,Nbins)
  readf,1,data
  close,1  
  data_bin_low=data[0,*]
  data_bin_high=data[1,*]
  data_obs=data[2,*]
  data_err=data[3,*]
  data_bin=data_bin_low+(data_bin_high-data_bin_low)/2.
 
  sel=where(data_err/data_obs lt min_obs_err,n)
  if(n gt 0) then  data_err[sel]=data_obs[sel]*min_obs_err

  if(data_bin[1]-data_bin[0] gt 0.) then begin 
     Bin = data_bin[1]-data_bin[0]
  endif else begin 
     Bin = data_bin[0]-data_bin[1]
  endelse

;overplot OBS constraints  
  if(plot_obs eq 0) then $
     joint_obs_symbols_size=1.0 $
  else joint_obs_symbols_size=0.7

  if(plot_only_new_obs eq 1) then joint_obs_symbols_size=1.0

  symbols, 2,  joint_obs_symbols_size
  oploterror, data_bin, Alog10(data_obs), Alog10((data_obs+data_err)/data_obs), $
              /hiba, psym=8, color=200, errcolor=200, HATLENGTH = 200.0  , ERRTHICK=4
  oploterror, data_bin, Alog10(data_obs), Alog10(data_obs/(data_obs-data_err)), $
              /loba, psym=8, color=200, errcolor=200, HATLENGTH = 200.0 , ERRTHICK=4


  endif


;LABELS
;redshift
   plot_label, xlog=0, ylog=0, type='label', label='SMF, z=3.0' , xmin, xmax, ymin, ymax, $
               x_percentage=4.5, x2_percentage=0., y_percentage=4.
 
;observations
   if(plot_obs eq 1) then begin
      plot_label, xlog=0, ylog=0, type='label', label='Observations for MCMC' , xmin, xmax, ymin, ymax, $
                  x_percentage=6.0, x2_percentage=0., y_percentage=2.8 

      plot_label, xlog=0, ylog=0, type='symbol', xmin, xmax, ymin, ymax, $
                  x_percentage=6.2, x2_percentage=0., y_percentage=3.0, $
                  color=200, errcolor=200, sym_num=2, sym_size=joint_obs_symbols_size, HATLENGTH = 80, err_size=0.15
   endif
;models

   plot_label, xlog=0, ylog=0, type='label', label=prefix_this_model , xmin, xmax, ymin, ymax, $
               x_percentage=6.0-1., x2_percentage=0., y_percentage=2., $
                  charsize=1.25, charthick=4 
   plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
               x_percentage=6.9-1., x2_percentage=6.2-1., y_percentage=2.15, $
               linestyle=0, color=70, linethick=8
 
   if(do_previous_model2 eq 1) then begin
      plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model2, xmin, xmax, ymin, ymax, $
                  x_percentage=6.0-1., x2_percentage=0., y_percentage=1.3, $
                  charsize=1.25, charthick=4
      plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                  x_percentage=6.9-1., x2_percentage=6.2-1., y_percentage=1.45, $
                  linestyle=linestyle_previous_model2, color=70, linethick=8
   endif
   if(do_previous_model1 eq 1) then begin
      plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model1, xmin, xmax, ymin, ymax, $
                  x_percentage=6.0-1., x2_percentage=0., y_percentage=0.6, $
                  charsize=1.25, charthick=4
      plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                  x_percentage=6.9-1., x2_percentage=6.2-1., y_percentage=0.75, $
                  linestyle=linestyle_previous_model1, color=70, linethick=8
   endif


  

;check sample 
if(sample eq 1) then begin
  readcol,sample_file+'0_z3.00.txt',kband,phi,err,model, /SILENT
  oplot, kband, Alog10(model), color=0, thick=6
endif


end




pro chi_square_and_plot, mcmc_folder, property_name, redshift, total_chi2, plot_obs,  plot_only_new_obs, CentralMvir_MR, CentralMvir_MRII, $
                         property_MR, property_MRII, $
                         volume_MR, volume_MRII, StellarMass_MR, StellarMass_MRII,Hubble_h, $
                         write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                         linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                         linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII


  if(property_name eq 'smf') then observation_name='StellarMassFunction'
  if(property_name eq 'kband') then observation_name='KBandLF'
  if(property_name eq 'bband') then observation_name='BBandLF'
  if(property_name eq 'smf_red') then observation_name='StellarMassFunctionRed'
  if(property_name eq 'smf_blue') then observation_name='StellarMassFunctionBlue'

  property_MR_no_err=0
  property_MRII_no_err=0
 
  
  close,1
  openr,1,mcmc_folder+'ObsConstraints/'+observation_name+'_z'+redshift+'.txt'
  Nbins = 0L & readf,1,Nbins
  data=FLTARR(4,Nbins)
  readf,1,data
  close,1 
  data_bin_low=data[0,*]
  data_bin_high=data[1,*]
  data_obs=data[2,*]
  data_err=data[3,*]

  data_bin=data_bin_low+(data_bin_high-data_bin_low)/2.

  sel=where(data_err/data_obs lt min_obs_err,n)
  if(n gt 0) then  data_err[sel]=data_obs[sel]*min_obs_err

  if(data_bin[1]-data_bin[0] gt 0.) then begin 
     Bin = data_bin[1]-data_bin[0]
  endif else begin 
     Bin = data_bin[0]-data_bin[1]
  endelse

  if(do_mcmc_regions eq 1) then $
     plot_mcmc_regions, mcmc_folder, property_name, redshift, data_bin 


;overplot OBS constraints  
  if(plot_obs eq 0) then $
     joint_obs_symbols_size=1.0 $
  else joint_obs_symbols_size=0.7

  if(plot_only_new_obs eq 1) then joint_obs_symbols_size=1.0

  data_hubble=0.7
  symbols, 2,  joint_obs_symbols_size

  oploterror, data_bin, Alog10(data_obs), Alog10((data_obs+data_err)/data_obs), $
              /hiba, psym=8, color=200, errcolor=200, HATLENGTH = 200.0  , ERRTHICK=4
  oploterror, data_bin, Alog10(data_obs), Alog10(data_obs/(data_obs-data_err)), $
              /loba, psym=8, color=200, errcolor=200, HATLENGTH = 200.0 , ERRTHICK=4


  if(MRII eq 1) then begin
     sel_MR=where(StellarMass_MR gt 9.5, N_MR)
     sel_MRII=where(StellarMass_MRII lt 9.5, N_MRII)
  endif else begin
     sel_MR=where(StellarMass_MR gt -9.5, N_MR)
     sel_MRII=where(StellarMass_MRII lt -9.5, N_MRII)
  endelse
 
  model_MR=property_MR[sel_MR]
  weight_MR=model_MR/model_MR/(volume_MR)
  if( N_MRII gt 0) then begin
     model_MRII=property_MRII[sel_MRII]
     weight_MRII=model_MRII/model_MRII/(volume_MRII)
  endif else begin
     model_MRII=0.
     weight_MRII=0.
  endelse

;COMPUTE CHI_2 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  hist=fltarr(Nbins)

  for i=0,Nbins-1 do begin

     hist[i]=0

     sel=where(model_MR ge data_bin[i]-Bin/2. and model_MR le data_bin[i]+Bin/2.,n)
     if(n gt 0) then hist[i]+=total(weight_MR[sel])

     sel=where(model_MRII ge data_bin[i]-Bin/2. and model_MRII le data_bin[i]+Bin/2.,n)
     if(n gt 0) then hist[i]+=total(weight_MRII[sel])

     hist[i]/=Bin

  endfor

  chsq=0
  for i=0,Nbins-1 do begin        
     chsq += ((hist[i]-data_obs[i])*(hist[i]-data_obs[i]))/(data_err[i]*data_err[i])  
     ;print,data_bin[i], ((hist[i]-data_obs[i])*(hist[i]-data_obs[i]))/(data_err[i]*data_err[i])  
  endfor

  print,'chi^2 for '+property_name+' at z='+redshift+' ->',chsq
  print,'Like for '+property_name+' at z='+redshift+' ->',exp(-chsq/2.)
  total_chi2+=chsq
  print, 'total chi^2 =0',total_chi2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;DO PLOTS DOWN TO LOW MASSES (BELLOW OBSERVATIONS)

  new_Nbins=50

  if(data_bin[1]-data_bin[0] gt 0.) then begin 
     Bin = data_bin[1]-data_bin[0]
  endif else begin 
     Bin = data_bin[0]-data_bin[1]
  endelse

  bin_start=min(data_bin)-5.0
  new_data_bin=findgen(new_Nbins)*Bin+bin_start


  hist=fltarr(new_Nbins)

  for i=0,new_Nbins-1 do begin

     hist[i]=0

     sel=where(model_MR ge new_data_bin[i]-Bin/2. and model_MR le new_data_bin[i]+Bin/2.,n)
     if(n gt 0) then hist[i]+=total(weight_MR[sel])

     sel=where(model_MRII ge new_data_bin[i]-Bin/2. and model_MRII le new_data_bin[i]+Bin/2.,n)
     if(n gt 0) then hist[i]+=total(weight_MRII[sel])

     hist[i]/=Bin

  endfor

;***PLOT***
  color=70
  if(property_name eq 'smf_blue') then color =200
  oplot,new_data_bin,Alog10(hist), color=color, thick=6


;WRITE TO FILE
  file=file_to_write+'_'+property_name+'_'+'z'+ STRTRIM(number_formatter(redshift, DECIMAL=2),2)+'.txt'
  write_to_file, file, new_data_bin,Alog10(hist), write_files, plot_after_write

  if(property_name eq 'smf' and redshift ne '0.00' and do_previous_model2 eq 0) then begin
     plot_smf_no_err, property_MR_no_err,property_MRII_no_err , N_MR, N_MRII, $
                      sel_MR,sel_MRII,weight_MR,weight_MRII,data_bin,Bin,Nbins, $
                      write_files, plot_after_write, file_to_write, property_name, redshift
  endif



end



pro plot_obs_sfr_vmax_contour, dr7_gal_final, xmin, xmax, xbin, ymin, ymax, ybin
 
  hubble_h_WMAP7=0.7
  Nbinsx=(xmax-xmin)/xbin
  Nbinsy=(ymax-ymin)/ybin
  obs_mstar = dr7_gal_final.jarle_median_mass+Alog10(hubble_h_WMAP7^2)  
  obs_sfr=dr7_gal_final.median_sfr+Alog10(hubble_h_WMAP7^2)
 
  mag=dr7_gal_final.jarle_mag[2] ;absolute mag
  zz=dr7_gal_final.z 
  max_d=10^((17.6-mag)/5.+1.)/1.e6 ;apparent
  weight = 1./max_d^3

  nii=total(weight)
  indx=0
  indy=0
      
  frac=fltarr(Nbinsx,Nbinsy) 
  mstar_c=indgen(Nbinsx)*xbin+xbin/2.+xmin
  sfr_c=indgen(Nbinsy)*ybin+ybin/2.+ymin

  xind=0


  for xid=xmin,xmax-xbin,xbin do begin 
   yind=0  
   ;print,xid,xid+xbin
   sel=where(obs_mstar gt xid and obs_mstar lt xid+xbin, NN)
   nii=total(weight[NN])
    for yid=ymin,ymax-ybin,ybin do begin     
        ww=where(obs_mstar gt xid and obs_mstar lt xid+xbin and $
                 obs_sfr gt yid and obs_sfr lt yid+ybin,n)    
        if(n gt 1) then $
           frac[xind,yind]=total(weight[ww]) $
        else $          
         frac[xind,yind]=0.      
        yind++
     endfor
  
     xind++
   endfor
 
  sel=where(frac gt 0)
  min=min(frac[sel])
  max=max(frac[sel])

  Nbins=11
  
  levels=make_log_array(min, max, Nbins)
  sel=[5,7,9,10]
  Contour, (frac>0.),mstar_c,sfr_c, levels=levels[sel], /overplot, $
           xtitle ='Mass', ytitle = 'Density', xrange=[xmin,xmax],yrange=[ymin,ymax]
 
end




pro age_hist, dr4_gal_final, ibin, mass, type, age, min_mass, max_mass, redshift=redshift, $
              xmin, xmax, ymin, ymax, bin, Datadir,label_x,label_y, $
              write_files, plot_after_write, file_to_write, $
              ssfr_and_age_hist_multiple_plots, min_z, max_z, $
              linestyle_previous_model1, file_previous_model1, do_previous_model1, $
              linestyle_previous_model2, file_previous_model2, do_previous_model2, $
              prefix_previous_model1, prefix_previous_model2, prefix_this_model


  char_redshift=STRTRIM(number_formatter(redshift, DECIMAL=2),2)

  ;help,dr4_gal_final,/struct

  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, $
        xtitle = label_x, ytitle = label_y

;overplot axis
 ; if(min_mass ge 10.0) then begin
 ;    axis_Nbins=1000 
 ;    axis_bin=(xmax-xmin)/(axis_Nbins*1.)
 ;    x=findgen(axis_Nbins)*((xmin+(axis_Nbins*axis_bin))-(xmin+(axis_Nbins*0.0)))/(axis_Nbins*1.)+xmin
 ;    y=x*1./x*1.-1.0+ymax-0.001
  ;oplot,x,y,thick=9
;endif

  sel=where(type lt 3 and age gt 0. and mass gt 10^min_mass and mass lt 10^max_mass,n_gals)
  if(n_gals gt 1) then begin
     hist=histogram(age[sel],locations=c,min=xmin-2.0,max=xmax+2.0,binsize=bin)  
     oplot,c+bin/2.0,hist*1.0/(total(hist*1.0)), color=70, thick=6,psym=10
 
     file=file_to_write+'_age_hist_'+STRTRIM(number_formatter(min_mass, DECIMAL=1),2)+ $
          '_'+STRTRIM(number_formatter(max_mass, DECIMAL=1),2)+'_z'+char_redshift+'.txt'
     write_to_file, file, c+bin/2.0,hist*1.0/(total(hist*1.0)), write_files, plot_after_write
  endif


;PREVIOUS MODELS
  if(do_previous_model1) then begin
     file=file_previous_model1+'_age_hist_'+STRTRIM(number_formatter(min_mass, DECIMAL=1),2)+ $
          '_'+STRTRIM(number_formatter(max_mass, DECIMAL=1),2)+'_z'+char_redshift+'.txt'
     readcol,file,ssfr,frac, /SILENT
     oplot,ssfr,frac,psym=10,color=70,linestyle=linestyle_previous_model1,  thick=6
  endif
  if(do_previous_model2) then begin
     file=file_previous_model2+'_age_hist_'+STRTRIM(number_formatter(min_mass, DECIMAL=1),2)+ $
          '_'+STRTRIM(number_formatter(max_mass, DECIMAL=1),2)+'_z'+char_redshift+'.txt'
     readcol,file,ssfr,frac, /SILENT
     oplot,ssfr,frac,psym=10,color=70,linestyle=linestyle_previous_model2,  thick=6
  endif

  plot_obs_age_vmax, dr4_gal_final, min_mass, max_mass, xmin, xmax,bin, min_z, max_z
 
;LABELS
  label= STRTRIM(number_formatter(min_mass, DECIMAL=1),2) + $
         '<log!D10!N(M!D*!N[h!U-2!NM!D!9ng!3!N])<' + $
         STRTRIM(number_formatter(max_mass, DECIMAL=1),2)
  plot_label, xlog=0, ylog=0, type='label', label=label , xmin, xmax, ymin, ymax, $
              x_percentage=0.5, x2_percentage=0., y_percentage=9.0, $
              charthick=3., charsize=1.1
 
  if(ibin eq 3) then begin
     plot_label, xlog=0, ylog=0, type='label', label='SDSS/DR7' , xmin, xmax, ymin, ymax, $
                 x_percentage=1.5, x2_percentage=0., y_percentage=8.0, $
                 charthick=3., charsize=1.1 
     plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                 x_percentage=0.4, x2_percentage=1.3, y_percentage=8.1, $
                 linestyle=0, color=0, linethick=8 

     plot_label, xlog=0, ylog=0, type='label', label=prefix_this_model , xmin, xmax, ymin, ymax, $
              x_percentage=1.5, x2_percentage=0., y_percentage=7.2, $
              charthick=3., charsize=1.1
     plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.4, x2_percentage=1.3, y_percentage=7.3, $
              linestyle=0, color=70, linethick=8  
     if(do_previous_model2) then begin
        plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model2 , xmin, xmax, ymin, ymax, $
                    x_percentage=1.5, x2_percentage=0., y_percentage=6.4, $
                    charthick=3., charsize=1.1
        plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                    x_percentage=0.4, x2_percentage=1.3, y_percentage=6.5, $
                    linestyle=linestyle_previous_model2, color=70, linethick=8
     endif
     if(do_previous_model1) then begin
        plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model1, xmin, xmax, ymin, ymax, $
                    x_percentage=1.5, x2_percentage=0., y_percentage=5.6, $
                    charthick=3., charsize=1.1
        plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                    x_percentage=0.4, x2_percentage=1.3, y_percentage=5.7, $
                    linestyle=linestyle_previous_model1, color=70, linethick=8
     endif

     
  endif


;re-plot axis
  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, $
        xtitle = label_x, ytitle = label_y


end




pro plot_obs_age_vmax, dr4_gal_final, min_mass, max_mass, xmin, xmax,bin, min_z, max_z
 
  Hubble_h_WMAP7=0.7

  Nbins=(xmax-xmin)/bin+1
  mstar = dr4_gal_final.mass+alog10(Hubble_h_WMAP7^2)  
  obsout={r:fltarr(Nbins),frac:fltarr(Nbins)}
  
  obs_age=10^(dr4_gal_final.age)/1.e9
 
  ii=where(mstar gt min_mass and mstar le max_mass and $
           dr4_gal_final.z gt min_z and dr4_gal_final.z lt max_z,nii)
  
  frac=fltarr(Nbins)
  tmp_age=obs_age[ii] 
  tmpsample = dr4_gal_final[ii]
  
  lumdist=(1.+dr4_gal_final[ii].z)*comdist(dr4_gal_final[ii].z,0.7,0.27)/0.7*1.e6
  abs_mag=dr4_gal_final[ii].kcor_mag[1]-5.*(Alog10(lumdist)-1.)

  app_mag=dr4_gal_final[ii].kcor_mag[1] 
  max_d=10^((17.6-abs_mag)/5.+1.)/1.e6
 
  weight = double(1./max_d^3)


;WEIGHTS HAVE BEEN TESTED BY SELECTING SMALLER INTERVALS OF MASS LIMITED SAMPLES

   sel=where(abs_mag ne 0.)
   if(min_mass gt 11.4) then sel=where(abs_mag ne 0. and tmp_age gt 2.0)
   weight=weight[sel]
   tmp_age=tmp_age[sel]


   nii=total(weight[where(weight gt 0)])
   ind=0
    
   for cid =xmin,xmax-bin,bin do begin
      ww=where(tmp_age gt cid and tmp_age lt cid+bin and weight gt 0,nww)   
      if nww eq 0 then begin
         ind++
         continue
      endif
      frac[ind]=total(weight[ww])/nii    
      ind++
   endfor
   
   obsout.r=indgen(Nbins)*bin+bin/2.+xmin
   obsout.frac=frac
 
   oplot,obsout.r,obsout.frac,thick=6,color=0,psym=10
 
end


pro ssfr_hist, dr7_gal_final, ibin, mass, ssfr, type, min_mass, max_mass, redshift=redshift, $
               xmin, xmax, ymin, ymax, bin, Datadir,label_x,label_y, $
               write_files, plot_after_write, file_to_write, $
               ssfr_and_age_hist_multiple_plots, min_z, max_z, $
               linestyle_previous_model1, file_previous_model1, do_previous_model1, $
               linestyle_previous_model2, file_previous_model2, do_previous_model2, $
               prefix_previous_model1, prefix_previous_model2, prefix_this_model

  char_redshift=STRTRIM(number_formatter(redshift, DECIMAL=2),2)
 
  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, $
        xtitle = label_x, ytitle = label_y

;Replot axis
if(min_mass eq 9.0 or min_mass eq 11.0) then begin
  axis_Nbins=1000
  axis_bin=(ymax-ymin)/(axis_Nbins*1.)
  y=findgen(axis_Nbins)*((ymin+(axis_Nbins*axis_bin))-(ymin+(axis_Nbins*0.0)))/(axis_Nbins*1.)+ymin
  x=y*1./y*1.-1.0+xmin+0.0001
  ;oplot,x,y,thick=9
endif
if(min_mass ge 10.0) then begin
  axis_Nbins=1000 
  axis_bin=(xmax-xmin)/(axis_Nbins*1.)
  x=findgen(axis_Nbins)*((xmin+(axis_Nbins*axis_bin))-(xmin+(axis_Nbins*0.0)))/(axis_Nbins*1.)+xmin
  y=x*1./x*1.-1.0+ymax-0.0005
  ;oplot,x,y,thick=9
endif


;NEW MODEL

  sel=where(type lt 3 and mass gt 10^min_mass and mass lt 10^max_mass,n_gals)
  ssfr=ssfr[sel]



;fit to observations of low SSFR
    ;SSFR=-0.3*M -8.6
  mean_bin_mass=(min_mass+max_mass)/2.

  ;if(mean_bin_mass gt 9.6) then begin
     slope=-0.3
     b=-8.6
     width=0.5
 
     sel=where(ssfr lt 10^(-12.),n)
     ssfr[sel]=randomn(1l,n_elements(ssfr[sel]))*width*10^(slope*mean_bin_mass+b) + 10^(slope*mean_bin_mass+b)


     if(n_gals gt 1) then begin
        hist=histogram(Alog10(ssfr),locations=c,min=xmin-2.0,max=xmax+2.0,binsize=bin)  
        oplot,c+bin/2.0,hist*1.0/(total(hist*1.0)), color=70, thick=6,psym=10
    

        file=file_to_write+'_ssfr_hist_'+STRTRIM(number_formatter(min_mass, DECIMAL=1),2)+ $
             '_'+STRTRIM(number_formatter(max_mass, DECIMAL=1),2)+'_z'+char_redshift+'.txt'
        write_to_file, file, c+bin/2.0,hist*1.0/(total(hist*1.0)), write_files, plot_after_write
     endif

;PREVIOUS MODELS
  if(do_previous_model1) then begin
     file=file_previous_model1+'_ssfr_hist_'+STRTRIM(number_formatter(min_mass, DECIMAL=1),2)+ $
          '_'+STRTRIM(number_formatter(max_mass, DECIMAL=1),2)+'_z'+char_redshift+'.txt'
     readcol,file,ssfr,frac, /SILENT
     oplot,ssfr,frac,psym=10,color=70,linestyle=linestyle_previous_model1,  thick=6
  endif
  if(do_previous_model2) then begin
     file=file_previous_model2+'_ssfr_hist_'+STRTRIM(number_formatter(min_mass, DECIMAL=1),2)+ $
          '_'+STRTRIM(number_formatter(max_mass, DECIMAL=1),2)+'_z'+char_redshift+'.txt'
     readcol,file,ssfr,frac, /SILENT
     oplot,ssfr,frac,psym=10,color=70,linestyle=linestyle_previous_model2,  thick=6
  endif



;OBSERVATIONS - DR7

  plot_obs_ssfr_vmax, dr7_gal_final, rob, min_mass, max_mass, xmin, xmax, bin, min_z, max_z


;LABELS
  label= STRTRIM(number_formatter(min_mass, DECIMAL=1),2) + $
         '<log!D10!N(M!D*!N[h!U-2!NM!D!9ng!3!N])<' + $
         STRTRIM(number_formatter(max_mass, DECIMAL=1),2)
  plot_label, xlog=0, ylog=0, type='label', label=label , xmin, xmax, ymin, ymax, $
              x_percentage=0.5, x2_percentage=0., y_percentage=9.0, $
              charthick=3., charsize=1.1
 
  if(ibin eq 3) then begin
     plot_label, xlog=0, ylog=0, type='label', label='SDSS/DR7' , xmin, xmax, ymin, ymax, $
                 x_percentage=1.5, x2_percentage=0., y_percentage=8.0, $
                 charthick=3., charsize=1.1 
     plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                 x_percentage=0.4, x2_percentage=1.3, y_percentage=8.1, $
                 linestyle=0, color=0, linethick=8 

     plot_label, xlog=0, ylog=0, type='label', label=prefix_this_model, xmin, xmax, ymin, ymax, $
              x_percentage=1.5, x2_percentage=0., y_percentage=7.2, $
              charthick=3., charsize=1.1
     plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.4, x2_percentage=1.3, y_percentage=7.3, $
              linestyle=0, color=70, linethick=8  
     if(do_previous_model2) then begin
        plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model2 , xmin, xmax, ymin, ymax, $
                    x_percentage=1.5, x2_percentage=0., y_percentage=6.4, $
                    charthick=3., charsize=1.1
        plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                    x_percentage=0.4, x2_percentage=1.3, y_percentage=6.5, $
                    linestyle=linestyle_previous_model2, color=70, linethick=8
     endif
     if(do_previous_model1) then begin
        plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model1, xmin, xmax, ymin, ymax, $
                    x_percentage=1.5, x2_percentage=0., y_percentage=5.6, $
                    charthick=3., charsize=1.1
        plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                    x_percentage=0.4, x2_percentage=1.3, y_percentage=5.7, $
                    linestyle=linestyle_previous_model1, color=70, linethick=8
     endif

     
  endif


;re-plot axis
 plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, $
        xtitle = label_x, ytitle = label_y


end






pro plot_obs_ssfr_vmax, dr7_gal_final, rob, min_mass, max_mass, xmin, xmax,bin, min_z, max_z
 
;help,dr7_gal_final, /struct
   hubble_h_WMAP7=0.7
   Nbins=(xmax-xmin)/bin+1
   mstar = dr7_gal_final.jarle_median_mass+Alog10(hubble_h_WMAP7^2)  
   obsout={r:fltarr(Nbins),frac:fltarr(Nbins)} 
   obs_ssfr=dr7_gal_final.median_ssfr

   ii=where(mstar gt min_mass and mstar le max_mass and dr7_gal_final.z gt min_z and dr7_gal_final.z lt max_z,nii)
    
   frac=fltarr(Nbins)
   tmp_ssfr=obs_ssfr[ii]
   tmpsample = dr7_gal_final[ii]  
   mag=dr7_gal_final[ii].jarle_mag[2] ;absolute mag
   zz=dr7_gal_final[ii].z

   max_d=10^((17.6-mag)/5.+1.)/1.e6 ;apparent
   weight = 1./max_d^3
 
   nii=total(weight)
   ind=0
      
   for cid =xmin,xmax-bin,bin do begin
      ww=where(tmp_ssfr gt cid and tmp_ssfr lt cid+bin,nww)
      if nww eq 0 then begin
         ind++
         continue
      endif
      frac[ind]=total(weight[ww])/nii
      ind++
   endfor
   
   obsout.r=indgen(Nbins)*bin+bin/2.+xmin
   obsout.frac=frac
 
   oplot,obsout.r,obsout.frac,thick=6,color=0,psym=10
 
end


pro color_hist, dr7_gal_final, ibin, sam_mass, u_i, type, min_mass, max_mass, redshift=redshift, $
                xmin, xmax, ymin, ymax, bin, Datadir, $
                label_x, label_y, obs, min_z, max_z, $                
                write_files, plot_after_write, file_to_write, $
                linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                linestyle_previous_model2, file_previous_model2, do_previous_model2, $
                prefix_previous_model1, prefix_previous_model2, prefix_this_model  
  

  plot, findgen(10), /nodata, xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xthick=5, ythick=5, xstyle = 1, ystyle = 1, $ 
        xtitle = label_x, ytitle = label_y


  char_redshift=STRTRIM(number_formatter(redshift, DECIMAL=2),2)

;OBSERVATIONS
  qi_obs_color_hist, dr7_gal_final, Datadir, obs, min_mass, max_mass, min_z, max_z



;NEW MODEL

   sel=where(sam_mass gt min_mass and sam_mass lt max_mass,n)
   color=u_i[sel] 
   if(n gt 0) then begin
      counts=histogram(color,locations=c,min=xmin-bin,max=xmax+bin,binsize=bin)  
      oplot, c+bin/2., (counts*1.0)/(total(counts)*1.0), color=70, thick=6
   endif

   if(write_files eq 1 and n gt 0) then begin
      file=file_to_write+'_color_hist_'+STRTRIM(number_formatter(min_mass, DECIMAL=1),2)+ $
           '_'+STRTRIM(number_formatter(max_mass, DECIMAL=1),2)+'_z'+char_redshift+'.txt'
      write_to_file, file, c+bin/2., (counts*1.0)/(total(counts)*1.0), write_files, plot_after_write
   endif
 

;PREVIOUS_MODELS 
  if(do_previous_model1 eq 1) then begin
     file=file_previous_model1+'_color_hist_'+STRTRIM(number_formatter(min_mass, DECIMAL=1),2)+ $
          '_'+STRTRIM(number_formatter(max_mass, DECIMAL=1),2)+'_z'+char_redshift+'.txt'
     readcol,file,color,frac, /SILENT 
     oplot,color,frac, color=70,linestyle=linestyle_previous_model1, thick=6
  endif
  if(do_previous_model2 eq 1) then begin
     file=file_previous_model2+'_color_hist_'+STRTRIM(number_formatter(min_mass, DECIMAL=1),2)+ $
          '_'+STRTRIM(number_formatter(max_mass, DECIMAL=1),2)+'_z'+char_redshift+'.txt'
     readcol,file,color,frac, /SILENT
     oplot,color,frac, color=70,linestyle=linestyle_previous_model2, thick=6
  endif
 



;LABELS
  label= STRTRIM(number_formatter(min_mass, DECIMAL=1),2) + $
         '<log!D10!N(M!D*!N[h!U-2!NM!D!9ng!3!N])<' + $
         STRTRIM(number_formatter(max_mass, DECIMAL=1),2)
  plot_label, xlog=0, ylog=0, type='label', label=label , xmin, xmax, ymin, ymax, $
              x_percentage=0.5, x2_percentage=0., y_percentage=9.0, $
              charthick=3., charsize=1.1
 
  if(ibin eq 3) then begin
     plot_label, xlog=0, ylog=0, type='label', label='SDSS/DR7' , xmin, xmax, ymin, ymax, $
                 x_percentage=1.5, x2_percentage=0., y_percentage=8.0, $
                 charthick=3., charsize=1.1 
     plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                 x_percentage=0.4, x2_percentage=1.3, y_percentage=8.1, $
                 linestyle=0, color=0, linethick=8 

     plot_label, xlog=0, ylog=0, type='label', label=prefix_this_model , xmin, xmax, ymin, ymax, $
              x_percentage=1.5, x2_percentage=0., y_percentage=7.2, $
              charthick=3., charsize=1.1
     plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.4, x2_percentage=1.3, y_percentage=7.3, $
              linestyle=0, color=70, linethick=8  
     if(do_previous_model2) then begin
        plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model2 , xmin, xmax, ymin, ymax, $
                    x_percentage=1.5, x2_percentage=0., y_percentage=6.4, $
                    charthick=3., charsize=1.1
        plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                    x_percentage=0.4, x2_percentage=1.3, y_percentage=6.5, $
                    linestyle=linestyle_previous_model2, color=70, linethick=8
     endif
     if(do_previous_model1) then begin
        plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model1, xmin, xmax, ymin, ymax, $
                    x_percentage=1.5, x2_percentage=0., y_percentage=5.6, $
                    charthick=3., charsize=1.1
        plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                    x_percentage=0.4, x2_percentage=1.3, y_percentage=5.7, $
                    linestyle=linestyle_previous_model1, color=70, linethick=8
     endif

     
  endif

;re-plot axis 
  plot, findgen(10), /nodata, xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xthick=5, ythick=5, xstyle = 1, ystyle = 1, $ 
        xtitle = label_x, ytitle = label_y


 
end

pro qi_obs_color_hist, dr7_gal_final, Datadir, obs, min_mass, max_mass, min_z, max_z
 

  h=0.7

  mstar = obs.mstar             ;-alog10(h)*2.  
  obsout={r:fltarr(30),frac:fltarr(30)}
  
   ;obscolor=obs.jarle_mag[0]-obs.jarle_mag[3]
  obscolor=obs.mag[0]-obs.mag[3]

  ii=where(mstar gt min_mass and mstar le max_mass $
           and obs.z ge min_z and obs.z le max_z,nii)
    
  ;print,'N=',nii
  
  weight = fltarr(nii)+1
  frac=fltarr(30)
  tmpcolor=obscolor[ii]
  tmpsample = obs[ii]

   ;different weight then other histograms because there is different
   ;information. are identical
  ii= where(tmpsample.zs[1] le max_z,nii)
         
  if nii ne 0 then $
     weight[ii] = ( (comdis2(max_z,0.25,0.75)^3 - comdis2(min_z,0.25,0.75)^3) $
                    /(comdis2(tmpsample[ii].zs[1],0.25,0.75)^3-comdis2(min_z,0.25,0.75)^3)) 
  
   nii=total(weight)
   ind=0
      
   for cid =1.,3.9,0.1 do begin
      ww=where(tmpcolor gt cid and tmpcolor lt cid+0.1,nww)
      if nww eq 0 then begin
         ind++
         continue
      endif
      frac[ind]=total(weight[ww])/nii
      ind++
   endfor
   
   obsout.r=indgen(30)*0.1+0.05+1.
   obsout.frac=frac
   oplot,obsout.r,obsout.frac,thick=6,color=0

   h=0.7
   xmax=4.0
   xmin=0.0
   bin=0.1
   Nbins=(xmax-xmin)/bin+1
   mstar = dr7_gal_final.jarle_median_mass  
   obsout={r:fltarr(Nbins),frac:fltarr(Nbins)} 
   obs_color=dr7_gal_final.jarle_mag[0]-dr7_gal_final.jarle_mag[3]
  
   ii=where(mstar gt min_mass and mstar le max_mass and dr7_gal_final.z gt min_z and dr7_gal_final.z lt max_z,nii)
    

   frac=fltarr(Nbins)
   tmp_color=obs_color[ii]
   tmpsample = dr7_gal_final[ii] 
   mag=dr7_gal_final[ii].jarle_mag[2] ;absolute mag
   zz=dr7_gal_final[ii].z


   max_d=10^((17.6-mag)/5.+1.)/1.e6 ;apparent
   weight = 1./max_d^3

   nii=total(weight)
   ind=0
      
   for cid =xmin,xmax-bin,bin do begin
      ww=where(tmp_color gt cid and tmp_color lt cid+bin,nww)
      if nww eq 0 then begin
         ind++
         continue
      endif
      frac[ind]=total(weight[ww])/nii
      ind++
   endfor
   
   obsout.r=indgen(Nbins)*bin+bin/2.+xmin
   obsout.frac=frac
 
end




pro plot_UVJ_color,  G, i, hubble_h, redshift, $
                     file_to_write, write_files, plot_after_write, Datadir, $
                     linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                     linestyle_previous_model2, file_previous_model2, do_previous_model2, $
                     prefix_previous_model1, prefix_previous_model2, $
                     slope_red_fraction, offset_red_fraction, minimum_y_red_fraction

  
  color_UV=G.MagDust[0]-G.MagDust[2]
  color_VJ=G.MagDust[2]-G.MagDust[7]
  
  xmin=-0.4
  xmax=2.5
  ymin=-0.5
  ymax=2.4
  
  
  ;slope=0.88
  ;if(redshift lt 1.0) then offset=0.69 else offset=0.59   
  ;interception_x=(1.3-offset)/slope
  ;sel=where((color_VJ lt interception_x and color_UV lt 1.3) or $
  ;          (color_VJ gt interception_x and color_UV lt color_VJ*slope+offset))
  
  if(i eq 0) then ylabel = 'U-V' else ylabel = ''
  
       
  levels=findgen(30)*0.145+1.
  plot_hist_2d,color_VJ,color_UV,xmin = xmin-1.0,xmax=xmax+1.0, ymin = ymin-1.0, ymax=ymax+1.0, $
               xstyle = 1, ystyle = 1, xtitle ='V-J', $
               ytitle =  ylabel,xbin=0.1,ybin=0.05,/cont,/log,$                 
               levels=levels,xrange=[xmin,xmax],yrange=[ymin,ymax],xs=1,ys=1,/fill
    
  if(i eq 0) then begin   
     COLORBAR, BOTTOM=0, VERTICAL=0, POSITION=[0.12, 0.31, 0.24, 0.36], $
               MAX=0.,Min=alog10(10^min(levels)/10^max(levels)),DIVISIONS=4, FORMAT='(F6.1)'
  
     ;needed to re-set keywords
     plot_hist_2d,color_VJ,color_UV,xmin = xmin-1.0,xmax=xmax+1.0, ymin = ymin-1.0, ymax=ymax+1.0, $
                  xstyle = 1, ystyle = 1, xtitle ='V-J', $
                  ytitle =  ylabel,xbin=0.1,ybin=0.05,/cont,/log,$                
                  levels=levels,xrange=[xmin,xmax],yrange=[ymin,ymax],xs=1,ys=1,/fill  
  endif

;DATA

   N=1000.
   slope=0.88
   if(redshift lt 1.0) then offset=0.69 else offset=0.59

   a=findgen(N)*(xmax-xmin)/N+xmin
  
   cut2=findgen(N)/findgen(N)+0.3 
   sel2=where(a lt (1.3-offset)/slope)
   oplot,a[sel2],cut2[sel2],color=0   
 
   cut3=a*slope+offset
   sel3=where(a gt (1.3-offset)/slope) 
   oplot,a[sel3],cut3[sel3],color=0




;MODEL 
   if(redshift lt 0.8) then begin
      get_slope, 0.5, 1.35, 2.5, 1.9, slope, offset  
      print,'new UVJ slope z=0.4:',slope, offset        
   endif
 
   if(redshift gt 0.8 and redshift lt 1.1) then begin
      get_slope, 0.5, 1.33, 2.5, 1.93, slope, offset  
      print,'new UVJ slope z=1.0:',slope, offset    
   endif


   if(redshift gt 1.9 and redshift lt 2.1) then begin     
      get_slope, 0.8, 1.2, 2.5, 1.8, slope, offset 
      print,'new UVJ slope z=2.0:',slope, offset 
   endif

   if(redshift gt 2.1) then begin     
      get_slope, 0.8, 1.1, 2.5, 1.7, slope, offset 
      print,'new UVJ slope z=3.0:',slope, offset 
   endif


   if (redshift lt 2.1) then minimum_y=1.3 else minimum_y=1.3

   N=1000.
   a=findgen(N)*(xmax-xmin)/N+xmin  
   cut2=findgen(N)/findgen(N)+minimum_y-1.  
   sel2=where(a lt (minimum_y-offset)/slope)
   oplot,a[sel2],cut2[sel2],color=0 ,linestyle=2  
 
   cut3=a*slope+offset
   sel3=where(a gt (minimum_y-offset)/slope,n)
   oplot,a[sel3],cut3[sel3],color=0,linestyle=2
  

   char_redshift=STRTRIM(number_formatter(redshift, DECIMAL=1),2)

   label='z='+char_redshift
   plot_label, xlog=0, ylog=0, type='label', label=label, xmin, xmax, ymin, ymax, $
               x_percentage=7.5, x2_percentage=0., y_percentage=1.0, $
               charthick=3., charsize=1.1
   
   if(i eq 3) then begin
      plot_label, xlog=0, ylog=0, type='label', label='Observational cut' , xmin, xmax, ymin, ymax, $
                  x_percentage=1.5, x2_percentage=0., y_percentage=9.0, $
                  charthick=3., charsize=1.0
      plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                  x_percentage=0.4, x2_percentage=1.3, y_percentage=9.2, $
                  linestyle=0, color=linecolor, linethick=8
      
      plot_label, xlog=0, ylog=0, type='label', label='BestFit model cut' , xmin, xmax, ymin, ymax, $
                  x_percentage=1.5, x2_percentage=0., y_percentage=8.2, $
                  charthick=3., charsize=1.0
      plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                  x_percentage=0.4, x2_percentage=1.3, y_percentage=8.4, $
                  linestyle=2, color=linecolor, linethick=8
   endif
end



pro plot_red_fraction_colorcut, i_z, StellarMass_MR, StellarMass_MRII, G_MR, G_MRII, sel_MR, sel_MRII, model_err_red_frac, $
                       redshift=redshift, hubble_h, write_files, plot_after_write, file_to_write, $
                       linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                       linestyle_previous_model2, file_previous_model2, do_previous_model2, $
                       prefix_previous_model1, prefix_previous_model2, sample, sample_file, mcmc_folder,  prefix_this_model

  char_redshift=STRTRIM(number_formatter(redshift, DECIMAL=2),2)
  
  bin=0.25
  xmin=12.0
  xmax=7.5
  ymin=0.0
  ymax=1.2
 
  if(redshift lt 0.2) then  ytitle =  '!4U!3!Dred!N/!4U!3!Dtotal!N' else ytitle =  ''

  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, xtitle ='log!D10!N(M!D*!N[h!U-2!NM!D!9ng!3!N])', $
        ytitle =  ytitle
  
;MODEL
  combined_frac=fltarr(((xmin+bin)-(xmax-bin))/bin+1)

;MR
   hist_red=histogram(Alog10(StellarMass_MR[sel_MR]),locations=c,min= xmax-bin,max=xmin+bin,binsize=bin)   
   hist_total=histogram(Alog10(StellarMass_MR),locations=c,min=xmax-bin,max=xmin+bin,binsize=bin) 

   frac=hist_red*1./hist_total*1.     
   sel=where(c gt 9.5)
   combined_frac[sel]=frac[sel]

;MRII  
   hist_red=histogram(Alog10(StellarMass_MRII[sel_MRII]),locations=c,min= xmax-bin,max=xmin+bin,binsize=bin)   
   hist_total=histogram(Alog10(StellarMass_MRII),locations=c,min=xmax-bin,max=xmin+bin,binsize=bin) 

   frac=hist_red*1./hist_total*1.    
   sel=where(c le 9.5)
   combined_frac[sel]=frac[sel]

;plot
   symbols,2,0.7    
   oplot,c+bin/2.0+0.035,combined_frac, color=70
 
   write_to_file, file_to_write+'_redfrac_colorcut_z'+char_redshift+'.txt', c+bin/2.0, combined_frac, write_files, plot_after_write


;PREVIOUS MODELS
   previous_models_redshift=[0.1,0.4,1.,2.,3.0]
   previous_models_char_redshift=STRTRIM(number_formatter(previous_models_redshift[i_z], DECIMAL=2),2)
   if(do_previous_model1 eq 1) then begin
      readcol,file_previous_model1+'_redfrac_colorcut_z'+previous_models_char_redshift+'.txt',mass,frac, /SILENT
      symbols,1,0.7
      sel=where(mass lt 11.4)
      oplot,mass[sel],frac[sel],linestyle=linestyle_previous_model1,color=70  
   endif
 
   if(do_previous_model2 eq 1) then begin
      readcol,file_previous_model2+'_redfrac_colorcut_z'+char_redshift+'.txt',mass,frac, /SILENT
      symbols,1,0.7
      sel=where(mass lt 11.6)
      oplot,mass[sel],frac[sel],linestyle=linestyle_previous_model2,color=70     
   endif


;OBSERVATIONS
   observation_name='RedFraction'
   close,1
   openr,1,mcmc_folder+'ObsConstraints/'+observation_name+'_z'+char_redshift+'.txt'
   Nbins = 0L & readf,1,Nbins
   data=FLTARR(4,Nbins)
   readf,1,data
   close,1 
   data_bin_low=data[0,*]
   data_bin_high=data[1,*]
   data_obs=data[2,*]
   data_err=data[3,*] 

   oploterror,data_bin_low+(data_bin_high-data_bin_low)/2., data_obs, data_err, $
              psym=8,color=200, errcolor=200, HATLENGTH = 80.0


;sample
 if(sample eq 1) then begin
    readcol,sample_file+'10_z'+char_redshift+'.txt',mass,phi,err,model, /SILENT
    oplot, mass,model, color=0, thick=3
 endif


;LABELS
  plot_label, xlog=0, ylog=0, type='label', label='z='+char_redshift , xmin, xmax, ymin, ymax, $
                 x_percentage=7.0, x2_percentage=0., y_percentage=5.5, $
                 charthick=3., charsize=1.1

  if(i_z eq 4) then begin    
     plot_label, xlog=0, ylog=0, type='label', label='Observations used in MCMC' , xmin, xmax, ymin, ymax, $
                 x_percentage=1.1, x2_percentage=0., y_percentage=9.0, $
                 charthick=2.5, charsize= 0.95
     plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
              x_percentage=0.7, x2_percentage=0., y_percentage=9.2, $
              color=200, errcolor=200, sym_num=1, sym_size=0.6, HATLENGTH = 100, err_size=0.025

     plot_label, xlog=0, ylog=0, type='label', label=prefix_this_model , xmin, xmax, ymin, ymax, $
                 x_percentage=1.5, x2_percentage=0., y_percentage=8.3, $
                 charthick=2.5, charsize=0.95
     plot_label, xlog=0, ylog=0, type='line', label=label, xmin, xmax, ymin, ymax, $
                 x_percentage=0.4, x2_percentage=1.3, y_percentage=8.4, $
                 linestyle=0, color=70, linethick=8
  
     if(do_previous_model2) then begin
        plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model2 , xmin, xmax, ymin, ymax, $
                    x_percentage=1.5, x2_percentage=0., y_percentage=7.5, $
                    charthick=2.5, charsize=0.95
        plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                    x_percentage=0.4, x2_percentage=1.3, y_percentage=7.6, $
                    linestyle=linestyle_previous_model2, color=70, linethick=8
     endif

     if(do_previous_model1) then begin
        plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model1, xmin, xmax, ymin, ymax, $
                    x_percentage=1.5, x2_percentage=0., y_percentage=6.7, $
                    charthick=2.5, charsize=0.95
        plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                    x_percentage=0.4, x2_percentage=1.3, y_percentage=6.8, $
                    linestyle=linestyle_previous_model1, color=70, linethick=8
     endif
    
  endif

end



pro plot_smf_by_color, G0_MR, G04_MR, G1_MR, G2_MR, G3_MR,G0_MRII, G04_MRII, G1_MRII, G2_MRII, G3_MRII, $
                       StellarMass_MR_0, StellarMass_MRII_0, StellarMass_MR_04, StellarMass_MRII_04, $
                       StellarMass_MR_1, StellarMass_MRII_1, StellarMass_MR_2, StellarMass_MRII_2, $
                       StellarMass_MR_3, StellarMass_MRII_3, Volume_MR, Volume_MRII, $
                       hubble_h, redshift, slope_red_fraction,offset_red_fraction, minimum_y_red_fraction, $
                       file_to_write, write_files, plot_after_write, Datadir, sample, sample_file, $
                       linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                       linestyle_previous_model2, file_previous_model2, do_previous_model2, $
                       prefix_previous_model1, prefix_previous_model2, do_mcmc_regions, prefix_this_model, $
                       plot_mcmc_regions, mcmc_folder, plot_obs, total_chi2, plot_only_new_obs, min_obs_err, MRII

  multiplot, /reset
  erase & multiplot, [5,2]

  bin=0.25
  xmin=12.5
  xmax=8.0
  ymin=-6
  ymax=0.5

  if(plot_obs eq 0) then  ymax=-0.3


  type=['passive','active']
  type_label=['red','blue']
 
  color=[70,200]
  for type_index=0,1 do begin
     for i_z=0,4 do begin  
     
        if(type_index gt 0 or (type_index eq 0 and i_z ne 0)) then multiplot
      
        if(type_index eq 1) then x_title = 'log!D10!N(M!D*!N[h!U-2!NM!D!9ng!3!N)' else x_title = ''
        if(i_z eq 0) then y_title = 'log!D10!N(!4U!3[h!U3!NMpc!U-3!Nlog!D10!NM!U-1!N])' else y_title=''

        plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
              xstyle = 1, ystyle = 1, xtitle = x_title, ytitle = y_title
        
        char_redshift=STRTRIM(number_formatter(redshift[i_z], DECIMAL=2),2)
        
        if(redshift[i_z] eq 0.0) then begin
           G_MR=G0_MR
           G_MRII=G0_MRII
           stellarmass_MR=Alog10(stellarmass_MR_0)
           stellarmass_MRII=Alog10(stellarmass_MRII_0)
        endif
        if(redshift[i_z] eq 0.4) then begin
           G_MR=G04_MR
           G_MRII=G04_MRII
           stellarmass_MR=Alog10(stellarmass_MR_04)
           stellarmass_MRII=Alog10(stellarmass_MRII_04)
        endif
        if(redshift[i_z] eq 1.0) then begin
           G_MR=G1_MR
           G_MRII=G1_MRII
           stellarmass_MR=Alog10(stellarmass_MR_1)
           stellarmass_MRII=Alog10(stellarmass_MRII_1)
        endif
        if(redshift[i_z] eq 2.0) then begin
           G_MR=G2_MR
           G_MRII=G2_MRII
           stellarmass_MR=Alog10(stellarmass_MR_2)
           stellarmass_MRII=Alog10(stellarmass_MRII_2)
        endif        
        if(redshift[i_z] eq 3.0) then begin
           G_MR=G3_MR
           G_MRII=G3_MRII
           stellarmass_MR=Alog10(stellarmass_MR_3)
           stellarmass_MRII=Alog10(stellarmass_MRII_3)
        endif
      
        slope=slope_red_fraction
        offset=offset_red_fraction
        minimum_y=minimum_y_red_fraction

        color_UV_MR=G_MR.MagDust[0]-G_MR.MagDust[2]
        color_VJ_MR=G_MR.MagDust[2]-G_MR.MagDust[7]
        color_UV_MRII=G_MRII.MagDust[0]-G_MRII.MagDust[2]
        color_VJ_MRII=G_MRII.MagDust[2]-G_MRII.MagDust[7]

        color_MR=G_MR.MagDust[15]-G_MR.MagDust[17]
        color_MRII=G_MRII.MagDust[15]-G_MRII.MagDust[17]
        rband_MR=G_MR.MagDust[17]-5.*alog10(hubble_h)
        rband_MRII=G_MRII.MagDust[17]-5.*alog10(hubble_h)

        if(i_z eq 0) then begin ;Z=0
           if(type_index eq 0) then begin ;red 
              sel_MR = where(color_MR gt (offset[i_z]-slope[i_z]*tanh((rband_MR+18.07)/1.09)))
              sel_MRII = where(color_MRII gt (offset[i_z]-slope[i_z]*tanh((rband_MRII+18.07)/1.09)))
          
              chi_square_and_plot, mcmc_folder, 'smf_red', char_redshift, total_chi2, plot_obs, plot_only_new_obs, $
                                   G_MR[sel_MR].CentralMvir, G_MRII[sel_MRII].CentralMvir, stellarmass_MR[sel_MR], stellarmass_MRII[sel_MRII], $
                                   volume_MR, volume_MRII, stellarmass_MR[sel_MR], stellarmass_MRII[sel_MRII],Hubble_h, $
                                   write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                                   linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                                   linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII
              
           endif 
           if(type_index eq 1) then begin ;blue
              sel_MR = where(color_MR lt (offset[i_z]-slope[i_z]*tanh((rband_MR+18.07)/1.09)))
              sel_MRII = where(color_MRII lt (offset[i_z]-slope[i_z]*tanh((rband_MRII+18.07)/1.09)))
          
              chi_square_and_plot, mcmc_folder, 'smf_blue', char_redshift, total_chi2, plot_obs, plot_only_new_obs, $
                                   G_MR[sel_MR].CentralMvir, G_MRII[sel_MRII].CentralMvir, stellarmass_MR[sel_MR], stellarmass_MRII[sel_MRII], $
                                   volume_MR, volume_MRII, stellarmass_MR[sel_MR], stellarmass_MRII[sel_MRII],Hubble_h, $
                                   write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                                   linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                                   linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII
  
           endif

        endif else begin        ; HIGH-Z
           
           if(type_index eq 0) then begin ;red        
              sel_MR = where( (color_VJ_MR lt (minimum_y[i_z]-offset[i_z])/slope[i_z] and color_UV_MR gt minimum_y[i_z]) or $                    
                              (color_VJ_MR gt (minimum_y[i_z]-offset[i_z])/slope[i_z] and color_UV_MR gt (color_VJ_MR*slope[i_z] + offset[i_z])) )   
              sel_MRII = where( (color_VJ_MRII lt (minimum_y[i_z]-offset[i_z])/slope[i_z] and color_UV_MRII gt minimum_y[i_z]) or $              
                          (color_VJ_MRII gt (minimum_y[i_z]-offset[i_z])/slope[i_z] and color_UV_MRII gt (color_VJ_MRII*slope[i_z] + offset[i_z])) )
           
              chi_square_and_plot, mcmc_folder, 'smf_red', char_redshift, total_chi2, plot_obs, plot_only_new_obs, $
                                   G_MR[sel_MR].CentralMvir, G_MRII[sel_MRII].CentralMvir, stellarmass_MR[sel_MR], stellarmass_MRII[sel_MRII], $
                                   volume_MR, volume_MRII, stellarmass_MR[sel_MR], stellarmass_MRII[sel_MRII],Hubble_h, $
                                   write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                                   linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                                   linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII
              
           endif 
           if(type_index eq 1) then begin ;blue
              sel_MR = where( (color_VJ_MR lt (minimum_y[i_z]-offset[i_z])/slope[i_z] and color_UV_MR lt minimum_y[i_z]) or $         
                              (color_VJ_MR gt (minimum_y[i_z]-offset[i_z])/slope[i_z] and color_UV_MR lt (color_VJ_MR*slope[i_z] + offset[i_z])) ) 
              sel_MRII = where( (color_VJ_MRII lt (minimum_y[i_z]-offset[i_z])/slope[i_z] and color_UV_MRII lt minimum_y[i_z]) or $       
                       (color_VJ_MRII gt (minimum_y[i_z]-offset[i_z])/slope[i_z] and color_UV_MRII lt (color_VJ_MRII*slope[i_z] + offset[i_z])) )  
       
              chi_square_and_plot, mcmc_folder, 'smf_blue', char_redshift, total_chi2, plot_obs, plot_only_new_obs, $
                                   G_MR[sel_MR].CentralMvir, G_MRII[sel_MRII].CentralMvir, stellarmass_MR[sel_MR], stellarmass_MRII[sel_MRII], $
                                   volume_MR, volume_MRII, stellarmass_MR[sel_MR], stellarmass_MRII[sel_MRII],Hubble_h, $
                                   write_files, plot_after_write, file_to_write, do_mcmc_regions, $
                                   linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                                   linestyle_previous_model2, file_previous_model2, do_previous_model2, min_obs_err, MRII
              
              
           endif
        endelse   
        
        mass_MR=stellarmass_MR[sel_MR]
        mass_MRII=stellarmass_MRII[sel_MRII]

        old_redshift=[0.1,0.4,1.0,2.0,3.0,4.0]
        old_char_redshift=STRTRIM(number_formatter(old_redshift[i_z], DECIMAL=2),2)

        if(do_previous_model1) then begin
           file=file_previous_model1+'_smf_'+type_label[type_index]+'_z'+old_char_redshift+'.txt'
           readcol,file,mass,phi, /SILENT
           oplot,mass,phi,color=color[type_index],linestyle=linestyle_previous_model1,  thick=6
        endif
        if(do_previous_model2) then begin
           file=file_previous_model2+'_smf_'+type_label[type_index]+'_z'+char_redshift+'.txt'
           readcol,file,mass,phi, /SILENT
           oplot,mass,phi,color=color[type_index],linestyle=linestyle_previous_model2,  thick=6
        endif


;sample
        if(sample eq 1 and type_index eq 0 and i_z ne 1) then begin
           readcol,sample_file+'8_z'+char_redshift+'.txt',mass,phi,err,model, /SILENT
           oplot, mass,Alog10(model), color=0, thick=6
        endif
        if(sample eq 1 and type_index eq 1 and i_z ne 1) then begin
           readcol,sample_file+'9_z'+char_redshift+'.txt',mass,phi,err,model, /SILENT
           oplot, mass,Alog10(model), color=0, thick=6
        endif

        char_redshift_new=STRTRIM(number_formatter(redshift[i_z], DECIMAL=1),2)
        label=type_label[type_index]+', z='+char_redshift_new
        plot_label, xlog=0, ylog=0, type='label', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=4.0, x2_percentage=0., y_percentage=1.0, $
                    charthick=3., charsize=1.1

;OBSERVATIONS  
        if(plot_obs eq 1) then begin
       
           yobs =[9.4,8.9,8.4,7.9,7.4,6.9]
           yobs2=[9.5,9.0,8.5,8.0,7.5,7.0]
           
           plot_individual_obs_smf_by_color, Datadir, type_index,i_z,yobs,yobs2, xmin, xmax, ymin, ymax
           
           
        endif else begin        ;LABELS FOR COMBINED OBSERVATIONS
           
           if(i_z eq 0 and type_index eq 0) then begin  
              plot_label, xlog=0, ylog=0, type='label', label='Observations used in MCMC' , xmin, xmax, ymin, ymax, $
                          x_percentage=1.1, x2_percentage=0., y_percentage=9.0, $
                          charthick=3., charsize=1.0
              if(plot_obs eq 0) then $
                 joint_obs_symbols_size=1.2 $
              else joint_obs_symbols_size=0.7
              plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                          x_percentage=0.7, x2_percentage=0., y_percentage=9.2, $
                          color=200, errcolor=200, sym_num=2, sym_size=joint_obs_symbols_size, HATLENGTH = 100, err_size=0.15
           endif
           


           if(i_z eq 4) then begin  
              
              if(type_index eq 0) then linecolor=70 else linecolor=200
              
              plot_label, xlog=0, ylog=0, type='label', label=prefix_this_model , xmin, xmax, ymin, ymax, $
                          x_percentage=1.5, x2_percentage=0., y_percentage=9.2, $
                          charthick=3., charsize=1.0
              plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                          x_percentage=0.4, x2_percentage=1.3, y_percentage=9.4, $
                          linestyle=0, color=linecolor, linethick=8
              
              if(do_previous_model2) then begin
                 plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model2 , xmin, xmax, ymin, ymax, $
                             x_percentage=1.5, x2_percentage=0., y_percentage=8.6, $
                             charthick=3., charsize=1.0
                 plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                             x_percentage=0.4, x2_percentage=1.3, y_percentage=8.8, $
                             linestyle=linestyle_previous_model2, color=linecolor, linethick=8
              endif
              
              if(do_previous_model1) then begin
                 plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model1, xmin, xmax, ymin, ymax, $
                             x_percentage=1.5, x2_percentage=0., y_percentage=8.0, $
                             charthick=3., charsize=1.0
                 plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                             x_percentage=0.4, x2_percentage=1.3, y_percentage=8.2, $
                             linestyle=linestyle_previous_model1, color=linecolor, linethick=8
              endif 
           endif           
           
        endelse
  

     endfor ;i_z=0,4
  endfor ;type_index=0,1
 
  multiplot, /reset

end ;plot_smf_by_color



pro plot_individual_obs_smf_by_color, Datadir, type_index,i, yobs, yobs2, xmin, xmax, ymin, ymax

  if(type_index eq 0) then begin ; passive

     if(i eq 0) then begin      ;z=0.0

        ;BELL2003
        readcol,Datadir+'bell2003_smf_bycolor_red_z0.txt',mass, p_phi, p_errorup, p_errordown, $
                a_phi, a_errorup, a_errordown, /SILENT                                                
        symbols, 30, 0.225
        h=0.7
        oploterror, mass-Alog10(1.6), p_phi,p_errorup, color = 120, $
                    errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass-Alog10(1.6), p_phi,p_errordown, color = 120, $
                    errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0
 
        plot_label, xlog=0, ylog=0, type='label', label='Bell (2003)', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[0], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[0], $
                    color=120, errcolor=120, sym_num=30, sym_size=0.225, HATLENGTH = 100, err_size=0.1
        
        ;Baldry2004
        readcol,Datadir+'baldry2004_smf_passive_z0.txt',mass,phi,errordown,errorup, /SILENT
        symbols, 2, 0.7
        h=0.7
        oploterror, mass-Alog10(1.6)+2.*Alog10(h), Alog10(phi)-3.*Alog10(h),Alog10((phi+errorup)/phi), color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass-Alog10(1.6)+2.*Alog10(h), Alog10(phi)-3.*Alog10(h),Alog10(phi/(phi-errordown)), color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0

        plot_label, xlog=0, ylog=0, type='label', label='Baldry (2004)', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[1], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[1], $
                    color=0, errcolor=0, sym_num=2, sym_size=0.7, HATLENGTH = 100, err_size=0.1
        
     endif

     if(i eq 1) then begin ;z=0.5 
        readcol,Datadir+'muzzin2013_maraston_z0.2_0.5.txt',mass,x_err, phi,error_up,error_down, $
                passive_phi, passive_error_up, passive_error_down, $
                active_phi, active_error_up, active_error_down, /SILENT
        symbols, 2, 0.7
        data_hubble=0.7
        oploterror, mass+2.*Alog10(data_hubble), passive_phi-3.*Alog10(data_hubble),passive_error_up, color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass+2.*Alog10(data_hubble), passive_phi-3.*Alog10(data_hubble),passive_error_down, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0  
        plot_label, xlog=0, ylog=0, type='label', label='Muzzin (2013) - 0.2<z<0.5', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[0], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[0], $
                    color=0, errcolor=0, sym_num=2, sym_size=0.7, HATLENGTH = 100, err_size=0.1
  
        readcol,Datadir+'/ilbert2013/MF_Vmax_red_z0.2_0.5.dat',mass,phi,phi_down,phi_up, /SILENT
        min_mass=min(mass)
        readcol,Datadir+'/ilbert2013/MF_Vmax_red_z0.2_0.5.dat',mass,phi,err_down,err_up, /SILENT
        symbols, 31, 0.3
        data_hubble=0.7
        sel=where(mass gt min_mass)
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                    errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                    errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0 
        plot_label, xlog=0, ylog=0, type='label', label='Ilbert (2013) - 0.2<z<0.5', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[1], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[1], $
                    color=120, errcolor=120, sym_num=31, sym_size=0.25, HATLENGTH = 100, err_size=0.1             
     endif

     if(i eq 2) then begin ;z=1.0
        readcol,Datadir+'muzzin2013_maraston_z0.5_1.0.txt',mass,x_err, phi,error_up,error_down, $
                passive_phi, passive_error_up, passive_error_down, $
                active_phi, active_error_up, active_error_down, /SILENT
        symbols, 1, 0.7
        data_hubble=0.7
        oploterror, mass+2.*Alog10(data_hubble), passive_phi-3.*Alog10(data_hubble),passive_error_up, color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass+2.*Alog10(data_hubble), passive_phi-3.*Alog10(data_hubble),passive_error_down, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Muzzin (2013) - 0.5<z<1.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[2], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[2], $
                    color=0, errcolor=0, sym_num=1, sym_size=0.7, HATLENGTH = 100, err_size=0.1      


        readcol,Datadir+'muzzin2013_maraston_z1.0_1.5.txt',mass,x_err, phi,error_up,error_down, $
                passive_phi, passive_error_up, passive_error_down, $
                active_phi, active_error_up, active_error_down, /SILENT
        symbols, 2, 0.7
        data_hubble=0.7
        oploterror, mass+2.*Alog10(data_hubble), passive_phi-3.*Alog10(data_hubble),passive_error_up, color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass+2.*Alog10(data_hubble), passive_phi-3.*Alog10(data_hubble),passive_error_down, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Muzzin (2013) - 1.0<z<1.5', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[3], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[3], $
                    color=0, errcolor=0, sym_num=2, sym_size=0.7, HATLENGTH = 100, err_size=0.1      

     
        readcol,Datadir+'/ilbert2013/MF_Vmax_red_z0.8_1.1.dat',mass,phi,phi_down,phi_up, /SILENT
        min_mass=min(mass)
        readcol,Datadir+'/ilbert2013/MF_Vmax_red_z0.8_1.1.dat',mass,phi,err_down,err_up, /SILENT
        symbols, 31, 0.3
        data_hubble=0.7
        sel=where(mass gt min_mass)
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                    errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                    errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0 
        plot_label, xlog=0, ylog=0, type='label', label='Ilbert (2013) - 0.8<z<1.1', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[4], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[4], $
                    color=120, errcolor=120, sym_num=31, sym_size=0.25, HATLENGTH = 100, err_size=0.1      



        readcol,Datadir+'tomczak2013_red_smf_z0.75_1.00.txt',obs_mass, obs_phi, obs_error_up, obs_error_down, /SILENT
        symbols, 22, 0.35
        data_hubble=0.7
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_up, color = 70, $
                    errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_down, color = 70, $
                    errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Tomczak (2014) - 0.75<z<1.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[0], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[0], $
                    color=70, errcolor=70, sym_num=22, sym_size=0.35, HATLENGTH = 100, err_size=0.1      

        readcol,Datadir+'tomczak2013_red_smf_z1.00_1.25.txt',obs_mass, obs_phi, obs_error_up, obs_error_down, /SILENT
        symbols, 32, 0.35
        data_hubble=0.7
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_up, color = 70, $
                    errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_down, color = 70, $
                    errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Tomczak (2014) - 1.0<z<1.25', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[1], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[1], $
                    color=70, errcolor=70, sym_num=32, sym_size=0.35, HATLENGTH = 100, err_size=0.1      


     endif

     if(i eq 3) then begin ;z=2.0

        readcol,Datadir+'muzzin2013_maraston_z1.5_2.0.txt',mass,x_err, phi,error_up,error_down, $
                passive_phi, passive_error_up, passive_error_down, $
                active_phi, active_error_up, active_error_down, /SILENT
        symbols, 1, 0.7
        data_hubble=0.7
        oploterror, mass+2.*Alog10(data_hubble), passive_phi-3.*Alog10(data_hubble),passive_error_up, color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass+2.*Alog10(data_hubble), passive_phi-3.*Alog10(data_hubble),passive_error_down, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Muzzin (2013) - 1.5<z<2.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[2], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[2], $
                    color=0, errcolor=0, sym_num=1, sym_size=0.7, HATLENGTH = 100, err_size=0.1  

        readcol,Datadir+'muzzin2013_maraston_z2.0_2.5.txt',mass,x_err, phi,error_up,error_down, $
                passive_phi, passive_error_up, passive_error_down, $
                active_phi, active_error_up, active_error_down, /SILENT
        symbols, 2, 0.7
        data_hubble=0.7
        oploterror, mass+2.*Alog10(data_hubble), passive_phi-3.*Alog10(data_hubble),passive_error_up, color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass+2.*Alog10(data_hubble), passive_phi-3.*Alog10(data_hubble),passive_error_down, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Muzzin (2013) - 2.0<z<2.5', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[3], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[3], $
                    color=0, errcolor=0, sym_num=2, sym_size=0.7, HATLENGTH = 100, err_size=0.1  

        readcol,Datadir+'/ilbert2013/MF_Vmax_red_z1.5_2.0.dat',mass,phi,phi_down,phi_up, /SILENT
        min_mass=min(mass)
        readcol,Datadir+'/ilbert2013/MF_Vmax_red_z1.5_2.0.dat',mass,phi,err_down,err_up, /SILENT
        symbols, 21, 0.3
        data_hubble=0.7
        sel=where(mass gt min_mass)
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                    errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                    errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0 
        plot_label, xlog=0, ylog=0, type='label', label='Ilbert (2013) - 1.5<z<2.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[4], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[4], $
                    color=120, errcolor=120, sym_num=21, sym_size=0.25, HATLENGTH = 100, err_size=0.1  

        readcol,Datadir+'/ilbert2013/MF_Vmax_red_z2.0_2.5.dat',mass,phi,phi_down,phi_up, /SILENT
        min_mass=min(mass)
        readcol,Datadir+'/ilbert2013/MF_Vmax_red_z2.0_2.5.dat',mass,phi,err_down,err_up, /SILENT
        symbols, 31, 0.3
        data_hubble=0.7
        sel=where(mass gt min_mass)
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                    errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                    errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0 
        plot_label, xlog=0, ylog=0, type='label', label='Ilbert (2013) - 1.5<z<2.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[5], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[5], $
                    color=120, errcolor=120, sym_num=31, sym_size=0.25, HATLENGTH = 100, err_size=0.1  


        readcol,Datadir+'tomczak2013_red_smf_z1.50_2.00.txt',obs_mass, obs_phi, obs_error_up, obs_error_down, /SILENT
        symbols, 22, 0.35
        data_hubble=0.7
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_up, color = 70, $
                    errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_down, color = 70, $
                    errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Tomczak (2014) - 1.5<z<2.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[0], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[0], $
                    color=70, errcolor=70, sym_num=22, sym_size=0.35, HATLENGTH = 100, err_size=0.1

        readcol,Datadir+'tomczak2013_red_smf_z2.00_2.50.txt',obs_mass, obs_phi, obs_error_up, obs_error_down, /SILENT
        symbols, 32, 0.35
        data_hubble=0.7
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_up, color = 70, $
                    errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_down, color = 70, $
                    errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
         plot_label, xlog=0, ylog=0, type='label', label='Tomczak (2014) - 2.0<z<2.5', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[1], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[1], $
                    color=70, errcolor=70, sym_num=32, sym_size=0.35, HATLENGTH = 100, err_size=0.1


     endif

  
     if(i eq 4) then begin      ;z=3.0
        readcol,Datadir+'muzzin2013_maraston_z2.5_3.0.txt',mass,x_err, phi,error_up,error_down, $
                passive_phi, passive_error_up, passive_error_down, $
                active_phi, active_error_up, active_error_down, /SILENT
        symbols, 1, 0.7
        data_hubble=0.7
        oploterror, mass+2.*Alog10(data_hubble), passive_phi-3.*Alog10(data_hubble),passive_error_up, color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass+2.*Alog10(data_hubble), passive_phi-3.*Alog10(data_hubble),passive_error_down, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Muzzin (2013) - 2.5<z<3.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[0], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[0], $
                    color=0, errcolor=0, sym_num=1, sym_size=0.7, HATLENGTH = 100, err_size=0.1  

        readcol,Datadir+'muzzin2013_maraston_z3.0_4.0.txt',mass,x_err, phi,error_up,error_down, $
                passive_phi, passive_error_up, passive_error_down, $
                active_phi, active_error_up, active_error_down, /SILENT
        symbols, 2, 0.7
        data_hubble=0.7
        oploterror, mass+2.*Alog10(data_hubble), passive_phi-3.*Alog10(data_hubble),passive_error_up, color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass+2.*Alog10(data_hubble), passive_phi-3.*Alog10(data_hubble),passive_error_down, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Muzzin (2013) - 3.0<z<4.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[1], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[1], $
                    color=0, errcolor=0, sym_num=2, sym_size=0.7, HATLENGTH = 100, err_size=0.1  
      

        readcol,Datadir+'/ilbert2013/MF_Vmax_red_z2.5_3.0.dat',mass,phi,phi_down,phi_up, /SILENT
        min_mass=min(mass)
        readcol,Datadir+'/ilbert2013/MF_Vmax_red_z2.5_3.0.dat',mass,phi,err_down,err_up, /SILENT
        symbols, 31, 0.3
        data_hubble=0.7
        sel=where(mass gt min_mass)
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                    errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                    errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0 
        plot_label, xlog=0, ylog=0, type='label', label='Ilbert (2013) - 3.0<z<4.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[2], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[2], $
                    color=120, errcolor=120, sym_num=31, sym_size=0.3, HATLENGTH = 100, err_size=0.1  
     endif

  endif




 if(type_index eq 1) then begin ;active

    if(i eq 0) then begin     
       readcol,Datadir+'bell2003_smf_bycolor_blue_z0.txt',mass, p_phi, p_errorup, p_errordown, $
               a_phi, a_errorup, a_errordown, /SILENT                                                
       symbols, 30, 0.225
       h=0.7
       oploterror, mass-Alog10(1.6), a_phi,a_errorup, color = 120, $
                   errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
       oploterror, mass-Alog10(1.6), a_phi,a_errordown, color = 120, $
                   errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0
       plot_label, xlog=0, ylog=0, type='label', label='Bell (2003)', xmin, xmax, ymin, ymax, $
                   x_percentage=1.0, x2_percentage=0., y_percentage=yobs[0], $
                   charthick=3., charsize=0.9 
       plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                   x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[0], $
                   color=120, errcolor=120, sym_num=30, sym_size=0.225, HATLENGTH = 100, err_size=0.1
        

       readcol,Datadir+'baldry2004_smf_active_z0.txt',mass,phi,errordown,errorup, /SILENT
       symbols, 2, 0.7
       data_hubble=0.7
       oploterror, mass-Alog10(1.6)+2.*Alog10(data_hubble), Alog10(phi)-3.*Alog10(data_hubble),Alog10((phi+errorup)/phi), color = 0, $
                   errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
       oploterror, mass-Alog10(1.6)+2.*Alog10(data_hubble), Alog10(phi)-3.*Alog10(data_hubble),Alog10(phi/(phi-errordown)), color = 0, $
                   errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0       
        plot_label, xlog=0, ylog=0, type='label', label='Baldry (2004)', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[1], $
                    charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[1], $
                    color=0, errcolor=0, sym_num=2, sym_size=0.7, HATLENGTH = 100, err_size=0.1
        
    endif


  
      
     if(i eq 1) then begin         
        readcol,Datadir+'/ilbert2013/MF_Vmax_blue_z0.2_0.5.dat',mass,phi,phi_down,phi_up, /SILENT
        min_mass=min(mass)
        readcol,Datadir+'/ilbert2013/MF_Vmax_blue_z0.2_0.5.dat',mass,phi,err_down,err_up, /SILENT
        symbols, 31, 0.3
        data_hubble=0.7
        sel=where(mass gt min_mass)
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                    errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                    errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0 
        plot_label, xlog=0, ylog=0, type='label', label='Ilbert (2013) - 0.2<z<0.5', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[1], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[1], $
                    color=120, errcolor=120, sym_num=2, sym_size=0.7, HATLENGTH = 100, err_size=0.1
  
        readcol,Datadir+'muzzin2013_maraston_z0.2_0.5.txt',mass,x_err, phi,error_up,error_down, $
                passive_phi, passive_error_up, passive_error_down, $
                active_phi, active_error_up, active_error_down, /SILENT
        symbols, 2, 0.7
        data_hubble=0.7
        oploterror, mass+2.*Alog10(data_hubble), active_phi-3.*Alog10(data_hubble),active_error_up, color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass+2.*Alog10(data_hubble), active_phi-3.*Alog10(data_hubble),active_error_down, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Muzzin (2013) - 0.2<z<0.5', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[0], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[0], $
                    color=0, errcolor=0, sym_num=31, sym_size=0.25, HATLENGTH = 100, err_size=0.1  

     endif


     if(i eq 2) then begin
        readcol,Datadir+'muzzin2013_maraston_z0.5_1.0.txt',mass,x_err, phi,error_up,error_down, $
                active_phi, active_error_up, active_error_down, $
                active_phi, active_error_up, active_error_down, /SILENT
        symbols, 1, 0.7
        data_hubble=0.7
        oploterror, mass+2.*Alog10(data_hubble), active_phi-3.*Alog10(data_hubble),active_error_up, color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass+2.*Alog10(data_hubble), active_phi-3.*Alog10(data_hubble),active_error_down, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Muzzin (2013) - 0.5<z<1.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[2], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[2], $
                    color=0, errcolor=0, sym_num=1, sym_size=0.7, HATLENGTH = 100, err_size=0.1      

        readcol,Datadir+'muzzin2013_maraston_z1.0_1.5.txt',mass,x_err, phi,error_up,error_down, $
                active_phi, active_error_up, active_error_down, $
                active_phi, active_error_up, active_error_down, /SILENT
        symbols, 2, 0.7
        data_hubble=0.7
        oploterror, mass+2.*Alog10(data_hubble), active_phi-3.*Alog10(data_hubble),active_error_up, color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass+2.*Alog10(data_hubble), active_phi-3.*Alog10(data_hubble),active_error_down, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Muzzin (2013) - 1.0<z<1.5', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[3], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[3], $
                    color=0, errcolor=0, sym_num=2, sym_size=0.7, HATLENGTH = 100, err_size=0.1      

        readcol,Datadir+'/ilbert2013/MF_Vmax_blue_z0.8_1.1.dat',mass,phi,phi_down,phi_up, /SILENT
        min_mass=min(mass)
        readcol,Datadir+'/ilbert2013/MF_Vmax_blue_z0.8_1.1.dat',mass,phi,err_down,err_up, /SILENT
        symbols, 31, 0.3
        data_hubble=0.7
        sel=where(mass gt min_mass)
        oploterror, mass[sel]-0.14+2*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                    errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                    errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0      
        plot_label, xlog=0, ylog=0, type='label', label='Ilbert (2013) - 0.8<z<1.1', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[4], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[4], $
                    color=120, errcolor=120, sym_num=31, sym_size=0.25, HATLENGTH = 100, err_size=0.1      

        readcol,Datadir+'tomczak2013_blue_smf_z0.75_1.00.txt',obs_mass, obs_phi, obs_error_up, obs_error_down, /SILENT
        symbols, 22, 0.35
        data_hubble=0.7
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_up, color = 70, $
                    errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_down, color = 70, $
                    errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Tomczak (2014) - 0.75<z<1.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[0], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[0], $
                    color=70, errcolor=70, sym_num=22, sym_size=0.35, HATLENGTH = 100, err_size=0.1      

        readcol,Datadir+'tomczak2013_blue_smf_z1.00_1.25.txt',obs_mass, obs_phi, obs_error_up, obs_error_down, /SILENT
        symbols, 32, 0.35
        data_hubble=0.7
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_up, color = 70, $
                    errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_down, color = 70, $
                    errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Tomczak (2014) - 1.0<z<1.25', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[1], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[1], $
                    color=70, errcolor=70, sym_num=32, sym_size=0.35, HATLENGTH = 100, err_size=0.1      

     endif


     if(i eq 3) then begin
        readcol,Datadir+'muzzin2013_maraston_z1.5_2.0.txt',mass,x_err, phi,error_up,error_down, $
                active_phi, active_error_up, active_error_down, $
                active_phi, active_error_up, active_error_down, /SILENT
        symbols, 1, 0.7
        data_hubble=0.7
        oploterror, mass+2.*Alog10(data_hubble), active_phi-3.*Alog10(data_hubble),active_error_up, color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass+2.*Alog10(data_hubble), active_phi-3.*Alog10(data_hubble),active_error_down, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Muzzin (2013) - 1.5<z<2.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[2], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[2], $
                    color=0, errcolor=0, sym_num=1, sym_size=0.7, HATLENGTH = 100, err_size=0.1

        readcol,Datadir+'muzzin2013_maraston_z2.0_2.5.txt',mass,x_err, phi,error_up,error_down, $
                active_phi, active_error_up, active_error_down, $
                active_phi, active_error_up, active_error_down, /SILENT
        symbols, 2, 0.7
        data_hubble=0.7
        oploterror, mass+2.*Alog10(data_hubble), active_phi-3.*Alog10(data_hubble),active_error_up, color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass+2.*Alog10(data_hubble), active_phi-3.*Alog10(data_hubble),active_error_down, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Muzzin (2013) - 2.0<z<2.5', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[3], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[3], $
                    color=0, errcolor=0, sym_num=2, sym_size=0.7, HATLENGTH = 100, err_size=0.1  
      
        readcol,Datadir+'/ilbert2013/MF_Vmax_blue_z1.5_2.0.dat',mass,phi,phi_down,phi_up, /SILENT
        min_mass=min(mass)
        readcol,Datadir+'/ilbert2013/MF_Vmax_blue_z1.5_2.0.dat',mass,phi,err_down,err_up, /SILENT
        symbols, 21, 0.3
        data_hubble=0.7
        sel=where(mass gt min_mass)
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                    errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                    errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0 
        plot_label, xlog=0, ylog=0, type='label', label='Ilbert (2013) - 1.5<z<2.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[4], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[4], $
                    color=120, errcolor=120, sym_num=21, sym_size=0.25, HATLENGTH = 100, err_size=0.1  

        readcol,Datadir+'/ilbert2013/MF_Vmax_blue_z2.0_2.5.dat',mass,phi,phi_down,phi_up, /SILENT
        min_mass=min(mass)
        readcol,Datadir+'/ilbert2013/MF_Vmax_blue_z2.0_2.5.dat',mass,phi,err_down,err_up, /SILENT
        symbols, 31, 0.3
        data_hubble=0.7
        sel=where(mass gt min_mass)
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                    errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                    errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Ilbert (2013) - 1.5<z<2.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[5], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[5], $
                    color=120, errcolor=120, sym_num=31, sym_size=0.25, HATLENGTH = 100, err_size=0.1 

        readcol,Datadir+'tomczak2013_blue_smf_z1.50_2.00.txt',obs_mass, obs_phi, obs_error_up, obs_error_down, /SILENT
        symbols, 22, 0.35
        data_hubble=0.7
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_up, color = 70, $
                    errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_down, color = 70, $
                    errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Tomczak (2014) - 1.5<z<2.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[0], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[0], $
                    color=70, errcolor=70, sym_num=22, sym_size=0.35, HATLENGTH = 100, err_size=0.1

        readcol,Datadir+'tomczak2013_blue_smf_z2.00_2.50.txt',obs_mass, obs_phi, obs_error_up, obs_error_down, /SILENT
        symbols, 32, 0.35
        data_hubble=0.7
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_up, color = 70, $
                    errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, obs_mass+2.*Alog10(data_hubble), obs_phi-3.*Alog10(data_hubble),obs_error_down, color = 70, $
                    errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
         plot_label, xlog=0, ylog=0, type='label', label='Tomczak (2014) - 2.0<z<2.5', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[1], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[1], $
                    color=70, errcolor=70, sym_num=32, sym_size=0.35, HATLENGTH = 100, err_size=0.1
     endif


    if(i eq 4) then begin
        readcol,Datadir+'muzzin2013_maraston_z2.5_3.0.txt',mass,x_err, phi,error_up,error_down, $
                active_phi, active_error_up, active_error_down, $
                active_phi, active_error_up, active_error_down, /SILENT
        symbols, 1, 0.7
        data_hubble=0.7
        oploterror, mass+2.*Alog10(data_hubble), active_phi-3.*Alog10(data_hubble),active_error_up, color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass+2.*Alog10(data_hubble), active_phi-3.*Alog10(data_hubble),active_error_down, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Muzzin (2013) - 2.5<z<3.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[0], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[0], $
                    color=0, errcolor=0, sym_num=1, sym_size=0.7, HATLENGTH = 100, err_size=0.1  

        readcol,Datadir+'muzzin2013_maraston_z3.0_4.0.txt',mass,x_err, phi,error_up,error_down, $
                active_phi, active_error_up, active_error_down, $
                active_phi, active_error_up, active_error_down, /SILENT
        symbols, 2, 0.7
        data_hubble=0.7
        oploterror, mass+2.*Alog10(data_hubble), active_phi-3.*Alog10(data_hubble),active_error_up, color = 0, $
                    errcolor = 0, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass+2.*Alog10(data_hubble), active_phi-3.*Alog10(data_hubble),active_error_down, color = 0, $
                    errcolor = 0, psym = 8, /loba,HATLENGTH = 80.0
        plot_label, xlog=0, ylog=0, type='label', label='Muzzin (2013) - 3.0<z<4.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[1], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[1], $
                    color=0, errcolor=0, sym_num=2, sym_size=0.7, HATLENGTH = 100, err_size=0.1  

        readcol,Datadir+'/ilbert2013/MF_Vmax_blue_z2.5_3.0.dat',mass,phi,phi_down,phi_up, /SILENT
        min_mass=min(mass)
        readcol,Datadir+'/ilbert2013/MF_Vmax_blue_z2.5_3.0.dat',mass,phi,err_down,err_up, /SILENT
        symbols, 21, 0.3
        data_hubble=0.7
        sel=where(mass gt min_mass)
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                    errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                    errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0 
        plot_label, xlog=0, ylog=0, type='label', label='Ilbert (2013) - 2.5<z<3.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[2], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[2], $
                    color=120, errcolor=120, sym_num=21, sym_size=0.3, HATLENGTH = 100, err_size=0.1  

        readcol,Datadir+'/ilbert2013/MF_Vmax_blue_z3.0_4.0.dat',mass,phi,phi_down,phi_up, /SILENT
        min_mass=min(mass)
        readcol,Datadir+'/ilbert2013/MF_Vmax_blue_z3.0_4.0.dat',mass,phi,err_down,err_up, /SILENT
        symbols, 31, 0.3
        data_hubble=0.7
        sel=where(mass gt min_mass)
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble),err_up[sel], color = 120, $
                    errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
        oploterror, mass[sel]-0.14+2.*Alog10(data_hubble), phi[sel]-3.*Alog10(data_hubble), err_down[sel], color = 120, $
                    errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0 
        plot_label, xlog=0, ylog=0, type='label', label='Ilbert (2013) - 3.0<z<4.0', xmin, xmax, ymin, ymax, $
                    x_percentage=1.0, x2_percentage=0., y_percentage=yobs[3], charthick=3., charsize=0.9 
        plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                    x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[3], $
                    color=120, errcolor=120, sym_num=31, sym_size=0.3, HATLENGTH = 100, err_size=0.1  
     endif

  endif


end ;plot_individual_obs_smf_by_color




pro plot_gas_mass_function, G_MR, volume_MR, G_MRII, volume_MRII, Hubble_h, hubble_h_WMAP7, Datadir, sample, sample_file, $                  
                            write_files, plot_after_write, file_to_write, $
                            linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                            linestyle_previous_model2, file_previous_model2, do_previous_model2, $
                            prefix_previous_model1, prefix_previous_model2, prefix_this_model, MRII

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                      ;;;
;;;    COLD GAS Z=0      ;;;
;;;                      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  multiplot, /reset
  !p.multi=0
  loadct,6, /SILENT

 
  bin=0.25
  xmin=11.50
  xmax=8.0
  ymin=-6.
  ymax=0.
  
  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, $
        xtitle = 'log!D10!N(M!DHI!N[h!U-2!NM!D!9ng!3!N])', ytitle = 'log!D10!N (!4U!3 [h!U3!NMpc!U-3!Nlog!D10!NM!U-1!N])'
  

;NEW MODEL
   hist_final=fltarr(n_elements(c)) 

   if(MRII eq 1) then mass_cut=9.2 else mass_cut=-9.2

   sel=where(c gt mass_cut)  
   hist_final[sel]=Alog10(hist_MR[sel]/(volume_MR*bin))

   if (MRII eq 1) then begin
      sel=where(c lt mass_cut)
      hist_final[sel]=Alog10(hist_MRII[sel]/(volume_MRII*bin))
   endif

   oplot,c+bin/2.0,hist_final, color=70, thick=6

   file=file_to_write+'_coldgas_MF.txt'
   write_to_file, file, c+bin/2.0,hist_final, write_files, plot_after_write


;previous models
   if(do_previous_model1) then begin
      file=file_previous_model1+'_coldgas_MF.txt'
      readcol,file,mass,phi, /SILENT
      oplot,mass,phi,color=70,linestyle=linestyle_previous_model1,  thick=6
   endif

   if(do_previous_model2) then begin
      file=file_previous_model2+'_coldgas_MF.txt'
      readcol,file,mass,phi, /SILENT
      oplot,mass,phi,color=70,linestyle=linestyle_previous_model2,  thick=6
   endif



;labels
   plot_label, xlog=0, ylog=0, type='label', label=prefix_this_model , xmin, xmax, ymin, ymax, $
               x_percentage=3.5, x2_percentage=0., y_percentage=2., $
                  charsize=1.25, charthick=4
   plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
               x_percentage=2.7, x2_percentage=3.3, y_percentage=2.15, $
               linestyle=0, color=70, linethick=8
 
   if(do_previous_model2 eq 1) then begin
      plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model2, xmin, xmax, ymin, ymax, $
                  x_percentage=3.5, x2_percentage=0., y_percentage=1.3, $
                  charsize=1.25, charthick=4
      plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                  x_percentage=2.7, x2_percentage=3.3, y_percentage=1.45, $
                  linestyle=linestyle_previous_model2, color=70, linethick=8
   endif
   if(do_previous_model1 eq 1) then begin
      plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model1, xmin, xmax, ymin, ymax, $
                  x_percentage=3.5, x2_percentage=0., y_percentage=0.6, $
                  charsize=1.25, charthick=4
      plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                  x_percentage=2.7, x2_percentage=3.3, y_percentage=0.75, $
                  linestyle=linestyle_previous_model1, color=70, linethick=8
   endif





;OBSERVATIONS
   h=0.75
   readcol,Datadir+'zwaan2005_GMF.txt',mass,phi,errordown,errorup, /SILENT  
   mass=mass  
   symbols, 1, 0.7   
   oploterror, mass, phi-3.*Alog10(h), errorup, color = 200, $
               errcolor = 200, psym = 8, /hiba,HATLENGTH = 80.0
   oploterror, mass, phi-3.*Alog10(h), errordown, color = 200, $
               errcolor = 200, psym = 8, /loba,HATLENGTH = 80.0

   h=0.7
   readcol,Datadir+'haynes2011_gmf.txt',mass,phi,errordown,errorup, /SILENT  
   symbols, 2, 0.75  
   oploterror, mass+2.*Alog10(hubble_h_WMAP7), phi-3.*Alog10(hubble_h_WMAP7), errorup, color = 120, $
               errcolor = 120, psym = 8, /hiba,HATLENGTH = 80.0
   oploterror, mass+2.*Alog10(hubble_h_WMAP7), phi-3.*Alog10(hubble_h_WMAP7), errordown, color = 120, $
               errcolor = 120, psym = 8, /loba,HATLENGTH = 80.0



;LABELS
 
   plot_label, xlog=0, ylog=0, type='label', label='M!DHI!N=0.54M!Dcold!N', xmin, xmax, ymin, ymax, $
               x_percentage=6.5, x2_percentage=0., y_percentage=4.6, $
               charsize=1.25, charthick=4


   xyouts, 9.8, -4.0, 'z=0.0',charsize=1.3,charthick=4

   symbols, 2, 0.75
   oploterror, [11.3,11.3], [-0.38,-0.38], [0.1,0.1], $
               psym=8, color=120, errcolor=120, HATLENGTH = 80.0
   xyouts, 11.2, -0.5, 'Haynes (2011)'


   symbols, 1, 0.7
   oploterror, [11.3,11.3], [-0.78,-0.78], [0.1,0.1], $
               psym=8, color=200, errcolor=200, HATLENGTH = 80.0
   xyouts, 11.2, -0.9, 'Zwaan (2005)'

 

if(sample eq 1) then begin
  readcol,sample_file+'15_z0.10.txt',coldgas,phi,err,model, /SILENT  
  oplot, coldgas, model, color=0, thick=6
endif
 
end





pro plot_mstar_metals, G, redshift=redshift , hubble_h, $
                       write_files, plot_after_write, file_to_write, $
                       linestyle_previous_model1, file_previous_model1, do_previous_model1, prefix_previous_model1, $
                       linestyle_previous_model2, file_previous_model2, do_previous_model2, prefix_previous_model2, prefix_this_model


  char_redshift=STRTRIM(number_formatter(redshift, DECIMAL=2),2)

  xmin=9.0
  xmax=12.0
  ymin=-1.5
  ymax=0.5

  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
	 xstyle = 1, ystyle = 1, xtitle = 'log!D10!N(M!D*!N[h!U-2!NM!D!9ng!3!N])', $
	 ytitle =  'log!D10!N(Z/Z!D!9ng!3!N)'
 
  sam_mass=Alog10(G.StellarMass*hubble_h*1.e10)   
  metal=Alog10((G.MetalsStellarMass)/((G.StellarMass)*0.02)) 

  N=30
  median=fltarr(N)
  pc16=fltarr(N)
  pc84=fltarr(N)
  bin=0.2
  xmass=findgen(N)*bin+8.8

  for i=0,N-1 do begin
     sel=where(sam_mass gt xmass[i] and sam_mass lt xmass[i]+bin,n)
     if(n gt 0) then begin    
        metal1=metal[sel]
        median[i]=median(metal1)
        pc16[i] = metal1((SORT(metal1))(16*N_ELEMENTS(metal1)/100))
        pc84[i] = metal1((SORT(metal1))(84*N_ELEMENTS(metal1)/100))  
     endif   
  endfor

  xmass=xmass+bin/2.0

  sel=where(xmass lt 11.8)
  oplot,xmass[sel],median[sel],color=70,thick=6
  oplot,xmass[sel],pc16[sel],color=70,thick=6
  oplot,xmass[sel],pc84[sel],color=70,thick=6


  file=file_to_write+'_metals_z'+char_redshift+'.txt'
  write_to_file_with_errors, file, xmass[sel], median[sel], pc84[sel]-median[sel], median[sel]-pc16[sel], write_files, plot_after_write


;oplot,sam_mass+alog10(hubble_h),metal,psym=3

;observations from GALLAZI
  obsp50=[-0.60,-0.61,-0.65,-0.61,-0.52,-0.41,-0.23,-0.11,-0.01,0.04,0.07,0.10,0.12,0.13,0.14,0.15]
  obsp16=[-1.11,-1.07,-1.10,-1.03,-0.97,-0.90,-0.80,-0.65,-0.41,-0.24,-0.14,-0.09,-0.06,-0.04,-0.03,-0.03]
  obsp84=[-0.00,-0.00,-0.05,-0.01,0.05,0.09,0.14,0.17,0.20,0.22,0.24,0.25,0.26,0.28,0.29,0.30]
  
  symbols, 30, 0.25
  oplot,xmass,obsp50,psym=8,color=200
  oplot,xmass,obsp16,color=200
  oplot,xmass,obsp84,color=200


;previous models
  if(do_previous_model1 eq 1) then begin
     readcol,file_previous_model1+'_metals_median_z'+char_redshift+'.txt',mass,phi, /SILENT   
     oplot,mass,phi, color=70,linestyle=linestyle_previous_model1, thick=6
     ;readcol,file_previous_model1+'_metals_pc16_z'+char_redshift+'.txt',mass,phi, /SILENT   
     ;oplot,mass,phi, color=70,linestyle=linestyle_previous_model1, thick=6
     ;readcol,file_previous_model1+'_metals_pc84_z'+char_redshift+'.txt',mass,phi, /SILENT   
     ;oplot,mass,phi, color=70,linestyle=linestyle_previous_model1, thick=6
  endif
  if(do_previous_model2 eq 1) then begin
     readcol,file_previous_model2+'_metals_z'+char_redshift+'.txt',mass,phi, /SILENT   
     oplot,mass,phi, color=70,linestyle=linestyle_previous_model2, thick=6   
  endif
 


  ;loadct,0, /SILENT
  plot_label, xlog=0, ylog=0, type='label', label=prefix_this_model+' - z=0.1' , xmin, xmax, ymin, ymax, $
              x_percentage=3.7, x2_percentage=0., y_percentage=2.5, $
              charsize=1.1, charthick=3
  plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=2.9, x2_percentage=3.5, y_percentage=2.55,  $
              linestyle=0, color=70, linethick=8
  ;loadct,6, /SILENT

  loadct,2, /SILENT
  plot_label, xlog=0, ylog=0, type='label', label=prefix_this_model+' - z=3' , xmin, xmax, ymin, ymax, $
              x_percentage=3.7, x2_percentage=0., y_percentage=2.0, $
              charsize=1.1, charthick=3
  plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=2.9, x2_percentage=3.5, y_percentage=2.1, $
              linestyle=0, color=200, linethick=8
   loadct,6, /SILENT

  if(do_previous_model2 eq 1) then begin
     plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model2, xmin, xmax, ymin, ymax, $
                 x_percentage=3.7, x2_percentage=0., y_percentage=1.5, $
                 charsize=1.1, charthick=3
     plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                 x_percentage=2.9, x2_percentage=3.5, y_percentage=1.6, $
                 linestyle=linestyle_previous_model2, color=70, linethick=8
  endif
  if(do_previous_model1 eq 1) then begin
     plot_label, xlog=0, ylog=0, type='label', label=prefix_previous_model1, xmin, xmax, ymin, ymax, $
                 x_percentage=3.7, x2_percentage=0., y_percentage=1.0, $
                 charsize=1.1, charthick=3
     plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                 x_percentage=2.9, x2_percentage=3.5, y_percentage=1.1, $
                 linestyle=linestyle_previous_model1, color=70, linethick=8
  endif



  label= 'Gallazzi (2005)'
  x=xmin+(xmax-xmin)*3.7/10.
  y=ymin+(ymax-ymin)*(0.5)/10.
  xyouts, x, y,label,charsize=1.1, charthick=3
   
  x=xmin+(xmax-xmin)*3.3/10.
  y=ymin+(ymax-ymin)*(0.6)/10.
  symbols, 30, 0.25
  oplot, [x,x], [y,y], psym=8, color=200

 
 

end




pro plot_model_metals,  G, redshift=redshift , hubble_h, $
                        write_files, plot_after_write, file_to_write, plot_color=plot_color
                                      
  sam_mass=Alog10(G.StellarMass*hubble_h*1.e10)  
  metal=Alog10((G.MetalsStellarMass)/((G.StellarMass)*0.02)) 
   
  N=30
  median=fltarr(N)
  pc16=fltarr(N)
  pc84=fltarr(N)
  bin=0.2
  xmass=findgen(N)*bin+8.8

  for i=0,N-1 do begin
     sel=where(sam_mass gt xmass[i] and sam_mass lt xmass[i]+bin,n)
     if(n gt 0) then begin      
        metal1=metal[sel]
        median[i]=median(metal1)
        pc16[i] = metal1((SORT(metal1))(16*N_ELEMENTS(metal1)/100))
        pc84[i] = metal1((SORT(metal1))(84*N_ELEMENTS(metal1)/100))  
     endif   
  endfor

  xmass=xmass+bin/2.0

  sel=where(xmass lt 11.8) 
  loadct,2, /SILENT
  oplot,xmass[sel],median[sel],color=plot_color,thick=6
  loadct,6, /SILENT

  char_redshift=STRTRIM(number_formatter(redshift, DECIMAL=2),2)
  file=file_to_write+'_metals_z'+char_redshift+'.txt'
 
  write_to_file_with_errors, file, xmass[sel], median[sel], pc84[sel]-median[sel],median[sel]-pc16[sel], write_files, plot_after_write


end

















;;;;;;;;;;;;;;;;;;;
;;
;; MISC ROUTINES
;;
;;;;;;;;;;;;;;;;;;;

pro write_to_file, file, x, y, write_files, plot_after_write

  if(write_files eq 1) then begin
     close,1    
     openw,1,file
     n=n_elements(x)
     for i=0,n-1 do printf, 1,x[i],y[i]
     close,1
  endif

  if(plot_after_write eq 1) then begin
   readcol,file,x,y, /SILENT  
   oplot,x,y,color=0,thick=6
  endif

end

pro write_to_file_with_errors, file, x, y, y_err_up, y_err_down, write_files, plot_after_write

  if(write_files eq 1) then begin
     close,1    
     openw,1,file
     n=n_elements(x)
     for i=0,n-1 do printf, 1, x[i],y[i],y_err_up[i],y_err_down[i]
     close,1
  endif

  if(plot_after_write eq 1) then begin
   readcol,file,x,y,y_err_up,y_err_down, /SILENT  
   oploterror,x,y,y_err_up,/hiba,color=0,thick=6
   oploterror,x,y,y_err_down,/loba,color=0,thick=6
  endif

end



pro plot_label, xlog=xlog, ylog=ylog, type=type, label=label, xmin, xmax, ymin, ymax, $
                x_percentage=x_percentage, x2_percentage=x2_percentage, $
                y_percentage=y_percentage, linestyle=linestyle, linethick=linethick, $
                color=color, charthick=charthick, charsize=charsize, errcolor=errcolor, $
                sym_num=sym_num, sym_size=sym_size,HATLENGTH = HATLENGTH, err_size=err_size

  if(STRCMP(type,'symbol',6)) then symbols, sym_num, sym_size

  if(xlog eq 1 and ylog eq 1) then begin  
     if(STRCMP(type,'label',5)) then begin
        x=10^(Alog10(xmin)+(Alog10(xmax)-Alog10(xmin))*x_percentage/10.)
        y=10^(Alog10(ymin)+(Alog10(ymax)-Alog10(ymin))*y_percentage/10.)
        xyouts, x, y,label, charthick=charthick, charsize=charsize
     endif 
     
     if(STRCMP(type,'line',4)) then begin
        x1=10^(Alog10(xmin)+(Alog10(xmax)-Alog10(xmin))*x_percentage/10.)
        x2=10^(Alog10(xmin)+(Alog10(xmax)-Alog10(xmin))*x2_percentage/10.)
        y=10^(Alog10(ymin)+(Alog10(ymax)-Alog10(ymin))*y_percentage/10.)
        oplot, [x1,x2], [y,y],Color=color,thick=linethick,linestyle=linestyle
     endif 

     if(STRCMP(type,'symbol',6)) then begin
        x=10^(Alog10(xmin)+(Alog10(xmax)-Alog10(xmin))*x_percentage/10.)
        y=10^(Alog10(ymin)+(Alog10(ymax)-Alog10(ymin))*y_percentage/10.)        
        oploterror, [x,x], [y,y], [err_size,err_size], $
                    psym=8, color=color, errcolor=errcolor, HATLENGTH = HATLENGTH
     endif
  endif

 if(xlog eq 0 and ylog eq 1) then begin  
     if(STRCMP(type,'label',5)) then begin
        x=xmin+(xmax-xmin)*x_percentage/10.
        y=10^(Alog10(ymin)+(Alog10(ymax)-Alog10(ymin))*y_percentage/10.)
        xyouts, x, y,label, charthick=charthick, charsize=charsize
     endif 
     
     if(STRCMP(type,'line',4)) then begin
        x1=xmin+(xmax-xmin)*x_percentage/10.
        x2=xmin+(xmax-xmin)*x2_percentage/10.
        y=10^(Alog10(ymin)+(Alog10(ymax)-Alog10(ymin))*y_percentage/10.)
        oplot, [x1,x2], [y,y],Color=color,thick=linethick,linestyle=linestyle
     endif 

     if(STRCMP(type,'symbol',6)) then begin
        x=xmin+(xmax-xmin)*x_percentage/10.
        y=10^(Alog10(ymin)+(Alog10(ymax)-Alog10(ymin))*y_percentage/10.)        
        oploterror, [x,x], [y,y], [err_size,err_size], $
                    psym=8, color=color, errcolor=errcolor, HATLENGTH = HATLENGTH
     endif
  endif

 if(xlog eq 1 and ylog eq 0) then begin  
     if(STRCMP(type,'label',5)) then begin
        x=10^(Alog10(xmin)+(Alog10(xmax)-Alog10(xmin))*x_percentage/10.)
        y=ymin+(ymax-ymin)*y_percentage/10.
        xyouts, x, y,label, charthick=charthick, charsize=charsize
     endif 
     
     if(STRCMP(type,'line',4)) then begin
        x1=10^(Alog10(xmin)+(Alog10(xmax)-Alog10(xmin))*x_percentage/10.)
        x2=10^(Alog10(xmin)+(Alog10(xmax)-Alog10(xmin))*x2_percentage/10.)
        y=ymin+(ymax-ymin)*y_percentage/10.
        oplot, [x1,x2], [y,y],Color=color,thick=linethick,linestyle=linestyle
     endif 

     if(STRCMP(type,'symbol',6)) then begin
        x=10^(Alog10(xmin)+(Alog10(xmax)-Alog10(xmin))*x_percentage/10.)
        y=ymin+(ymax-ymin)*y_percentage/10.        
        oploterror, [x,x], [y,y], [err_size,err_size], $
                    psym=8, color=color, errcolor=errcolor, HATLENGTH = HATLENGTH
     endif
  endif

  if(xlog eq 0 and ylog eq 0) then begin
     if(STRCMP(type,'label',5)) then begin
        x=xmin+(xmax-xmin)*x_percentage/10.
        y=ymin+(ymax-ymin)*y_percentage/10.
        xyouts, x, y,label, charthick=charthick, charsize=charsize
     endif 
     
     if(STRCMP(type,'line',4)) then begin
        x1=xmin+(xmax-xmin)*x_percentage/10.
        x2=xmin+(xmax-xmin)*x2_percentage/10.
        y=ymin+(ymax-ymin)*y_percentage/10.
        oplot, [x1,x2], [y,y],Color=color,thick=linethick,linestyle=linestyle
     endif 

     if(STRCMP(type,'symbol',6)) then begin
        x=xmin+(xmax-xmin)*x_percentage/10.
        y=ymin+(ymax-ymin)*y_percentage/10.        
        oploterror, [x,x], [y,y], [err_size,err_size], $
                    psym=8, color=color, errcolor=errcolor, HATLENGTH = HATLENGTH
     endif
  endif

end


pro get_slope, x1, y1, x2, y2, slope, b
  slope=(y2-y1)/(x2-x1)
  b=y1-(y2-y1)/(x2-x1)*x1 
end

pro make_float_array, arr, min, max, bin
  arr=findgen((max-min)/bin)*bin+min+bin/2.0
end

pro median_and_percentiles, bin, xmin, xmax, x_variable, y_variable, x_binned, median, pc16, pc84 
 
  min_x=min(x_variable[where(x_variable gt -1.e20)])
  if(min_x gt xmin) then xmin=min_x
 
  Nbins=(xmax-xmin)/bin+1

  median=fltarr(Nbins)
  pc16=fltarr(Nbins) 
  pc84=fltarr(Nbins)  
  x_min=xmin-bin/2.

  for i=0,Nbins-1 do begin  

     sel=where(x_variable gt x_min+(i*bin) and x_variable lt x_min+((i+1.0)*bin),n) 
     if(n gt 0) then begin
        x_variable_sel=x_variable[sel]
        y_variable_sel=y_variable[sel]
        median[i]=median(y_variable_sel)
        pc16[i] = y_variable_sel((SORT(y_variable_sel))(16*N_ELEMENTS(y_variable_sel)/100))      
        pc84[i] = y_variable_sel((SORT(y_variable_sel))(84*N_ELEMENTS(y_variable_sel)/100))  
     endif  

  endfor

  x_binned=findgen(Nbins)*((x_min+(Nbins*bin))-(x_min+(Nbins*0.0)))/(Nbins*1.)+x_min+bin/2.
 
end


function make_log_array, min, max, Nbins

  log_min=Alog10(min)
  log_max=Alog10(max)
  log_bin=(log_max-log_min)*1./Nbins*1.

  new_arr=findgen(Nbins)/(Nbins)*(Nbins)*log_bin+log_min

  return, 10^new_arr
end



;+
; NAME:
;   COLORBAR
;
; PURPOSE:
;       The purpose of this routine is to add a color bar to the current
;       graphics window.
;
; CATEGORY:
;       Graphics, Widgets.
;
; CALLING SEQUENCE:
;       COLORBAR
;
; INPUTS:
;       None.
;
; KEYWORD PARAMETERS:
;
;       BOTTOM: The lowest color index of the colors to be loaded in
;                 the bar.
;
;       CHARSIZE: The character size of the color bar annotations. Default is 1.0.
;
;       COLOR:    The color index of the bar outline and characters. Default
;                 is ncolors - 1 + bottom.
;
;       DIVISIONS: The number of divisions to divide the bar into. There will
;                 be (divisions + 1) annotations. The default is 2.
;
;       FORMAT:   The format of the bar annotations. Default is '(F6.2)'.
;
;       MAX:      The maximum data value for the bar annotation. Default is
;                 NCOLORS-1.
;
;       MIN:      The minimum data value for the bar annotation. Default is 0.
;
;       NCOLORS:  This is the number of colors in the color bar.
;
;       POSITION: A four-element array of normalized coordinates in the same
;                 form as the POSITION keyword on a plot. Default is
;                 [0.88, 0.15, 0.95, 0.95] for a vertical bar and
;                 [0.15, 0.88, 0.95, 0.95] for a horizontal bar.
;
;       PSCOLOR:  This keyword is only applied if the output is being sent to
;                 a PostScript file. It indicates that the PostScript device
;                 is configured for color output. If this keyword is set, then
;                 the annotation is drawn in the color specified by the COLOR
;                 keyword. If the keyword is not set, the annotation is drawn
;                 in the color specified by the !P.COLOR system variable
;                 (usually this will be the color black). In general, this
;                 gives better looking output on non-color or gray-scale
;                 printers. If you are not specifically setting the annotation
;                 color (with the COLOR keyword), it will probably
;                 be better NOT to set this keyword either, even if you
;                 are outputting to a color PostScript printer.
;
;       RIGHT:    This puts the labels on the right-hand side of a vertical
;                 color bar. It applies only to vertical color bars.
;
;       TITLE:    This is title for the color bar. The default is to have
;                 no title.
;
;       TOP:      This puts the labels on top of the bar rather than under it.
;                 The keyword only applies if a horizontal color bar is rendered.
;
;       VERTICAL: Setting this keyword give a vertical color bar. The default
;                 is a horizontal color bar.
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;       Color bar is drawn in the current graphics window.
;
; RESTRICTIONS:
;       The number of colors available on the display device (not the
;       PostScript device) is used unless the NCOLORS keyword is used.
;
; EXAMPLE:
;       To display a horizontal color bar above a contour plot, type:
;
;       LOADCT, 5, NCOLORS=100
;       CONTOUR, DIST(31,41), POSITION=[0.15, 0.15, 0.95, 0.75], $
;          C_COLORS=INDGEN(25)*4, NLEVELS=25
;       COLORBAR, NCOLORS=100
;
; MODIFICATION HISTORY:
;       Written by: David Fanning, 10 JUNE 96.
;       10/27/96: Added the ability to send output to PostScript. DWF
;       11/4/96: Substantially rewritten to go to screen or PostScript
;           file without having to know much about the PostScript device
;           or even what the current graphics device is. DWF
;       1/27/97: Added the RIGHT and TOP keywords. Also modified the
;            way the TITLE keyword works. DWF
;       7/15/97: Fixed a problem some machines have with plots that have
;            no valid data range in them. DWF
;       12/5/98: Fixed a problem in how the colorbar image is created that
;            seemed to tickle a bug in some versions of IDL. DWF.
;       1/12/99: Fixed a problem caused by RSI fixing a bug in IDL 5.2. Sigh... DWF.
;-

PRO COLORBAR, BOTTOM=bottom, CHARSIZE=charsize, COLOR=color, DIVISIONS=divisions, $
   FORMAT=format, POSITION=position, MAX=max, MIN=min, NCOLORS=ncolors, $
   PSCOLOR=pscolor, TITLE=title, VERTICAL=vertical, TOP=top, RIGHT=right

   ; Is the PostScript device selected?

postScriptDevice = (!D.NAME EQ 'PS')

  ; Check and define keywords.

IF N_ELEMENTS(ncolors) EQ 0 THEN BEGIN

   ; Most display devices to not use the 256 colors available to
   ; the PostScript device. This presents a problem when writing
   ; general-purpose programs that can be output to the display or
   ; to the PostScript device. This problem is especially bothersome
   ; if you don't specify the number of colors you are using in the
   ; program. One way to work around this problem is to make the
   ; default number of colors the same for the display device and for
   ; the PostScript device. Then, the colors you see in PostScript are
   ; identical to the colors you see on your display. Here is one way to
   ; do it.

   IF postScriptDevice THEN BEGIN
      oldDevice = !D.NAME

         ; What kind of computer are we using? SET_PLOT to appropriate
         ; display device.

      thisOS = !VERSION.OS_FAMILY
      thisOS = STRMID(thisOS, 0, 3)
      thisOS = STRUPCASE(thisOS)
      CASE thisOS of
         'MAC': SET_PLOT, thisOS
         'WIN': SET_PLOT, thisOS
         ELSE: SET_PLOT, 'X'
      ENDCASE

         ; Open a window (to make sure !D.N_COLORS is accurate).

      WINDOW, /FREE, /PIXMAP, XSIZE=10, YSIZE=10
      WDELETE, !D.WINDOW

         ; Here is how many colors we should use.

      ncolors = !D.N_COLORS < 256
      SET_PLOT, oldDevice
    ENDIF ELSE ncolors = !D.N_COLORS < 256
ENDIF
IF N_ELEMENTS(bottom) EQ 0 THEN bottom = 0B
IF N_ELEMENTS(charsize) EQ 0 THEN charsize = 1.0
IF N_ELEMENTS(format) EQ 0 THEN format = '(F8.2)'
IF N_ELEMENTS(color) EQ 0 THEN color = ncolors - 1 + bottom
IF N_ELEMENTS(min) EQ 0 THEN min = 0.0
IF N_ELEMENTS(max) EQ 0 THEN max = FLOAT(ncolors) - 1
IF N_ELEMENTS(divisions) EQ 0 THEN divisions = 2
IF N_ELEMENTS(title) EQ 0 THEN title = ''
pscolor = KEYWORD_SET(pscolor)

IF KEYWORD_SET(vertical) THEN BEGIN
   bar = REPLICATE(1B,10) # BINDGEN(ncolors)
   IF N_ELEMENTS(position) EQ 0 THEN position = [0.88, 0.15, 0.95, 0.95]
ENDIF ELSE BEGIN
   bar = BINDGEN(ncolors) # REPLICATE(1B, 10)
   IF N_ELEMENTS(position) EQ 0 THEN position = [0.15, 0.88, 0.95, 0.95]
ENDELSE

   ; Scale the color bar.

 bar = BYTSCL(bar, TOP=ncolors-1) + bottom

   ; Get starting locations in DEVICE coordinates.

xstart = position(0) * !D.X_VSIZE
ystart = position(1) * !D.Y_VSIZE

   ; Get the size of the bar in DEVICE coordinates.

xsize = (position(2) - position(0)) * !D.X_VSIZE
ysize = (position(3) - position(1)) * !D.Y_VSIZE

   ; For PostScript output only, draw the annotation in !P.COLOR
   ; unless "pscolor" is set. This makes better output on grayscale
   ; printers.

IF postScriptDevice AND (pscolor NE 1) THEN BEGIN
   oldcolor = color
   color = !P.COLOR
ENDIF

   ; Display the color bar in the window. Sizing is
   ; different for PostScript and regular display.

IF postScriptDevice THEN BEGIN

   TV, bar, xstart, ystart, XSIZE=xsize, YSIZE=ysize

ENDIF ELSE BEGIN

   bar = CONGRID(bar, CEIL(xsize), CEIL(ysize), /INTERP)
   TV, bar, xstart, ystart

ENDELSE

   ; Annotate the color bar.

IF KEYWORD_SET(vertical) THEN BEGIN

   IF KEYWORD_SET(right) THEN BEGIN

      PLOT, [min,max], [min,max], /NODATA, XTICKS=1, $
         YTICKS=divisions, XSTYLE=1, YSTYLE=9, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT='(A1)', XTICKFORMAT='(A1)', YTICKLEN=0.1 , $
         YRANGE=[min, max], YTITLE=title

      AXIS, YAXIS=1, YRANGE=[min, max], YTICKFORMAT=format, YTICKS=divisions, $
         YTICKLEN=0.1, YSTYLE=1, COLOR=color, CHARSIZE=charsize

   ENDIF ELSE BEGIN

      PLOT, [min,max], [min,max], /NODATA, XTICKS=1, $
         YTICKS=divisions, XSTYLE=1, YSTYLE=9, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT=format, XTICKFORMAT='(A1)', YTICKLEN=0.1 , $
         YRANGE=[min, max]

      AXIS, YAXIS=1, YRANGE=[min, max], YTICKFORMAT='(A1)', YTICKS=divisions, $
         YTICKLEN=0.1, YTITLE=title, YSTYLE=1, COLOR=color, CHARSIZE=charsize

   ENDELSE

ENDIF ELSE BEGIN

   IF KEYWORD_SET(top) THEN BEGIN

      PLOT, [min,max], [min,max], /NODATA, XTICKS=divisions, $
         YTICKS=1, XSTYLE=9, YSTYLE=1, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT='(A1)', XTICKFORMAT='(A1)', XTICKLEN=0.1, $
         XRANGE=[min, max], XTITLE=title

      AXIS, XTICKS=divisions, XSTYLE=1, COLOR=color, CHARSIZE=charsize, $
         XTICKFORMAT=format, XTICKLEN=0.1, XRANGE=[min, max], XAXIS=1

   ENDIF ELSE BEGIN

      PLOT, [min,max], [min,max], /NODATA, XTICKS=divisions, $
         YTICKS=1, XSTYLE=1, YSTYLE=1, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT='(A1)', XTICKFORMAT=format, XTICKLEN=0.1, $
         XRANGE=[min, max], TITLE=title

    ENDELSE

ENDELSE

   ; Restore color variable if changed for PostScript.

IF postScriptDevice AND (pscolor NE 1) THEN color = oldcolor

END



;returnin comoving distance in Mpc/h
function comdist, redshift, Hubble_h, Omega_m
 
  H0 = Hubble_h*100. 
  WM = Omega_m                       ;Omega(matter)                
  WV = 1.0 - WM - 0.4165/(H0*H0)     ;Omega(lambda)
  WR = 4.165E-5/(Hubble_h*Hubble_h)  ;Omega(radiation)
  WK = 1-WM-WR-WV                    ;Omega curvaturve = 1-Omega(total)
 
  c = 299792.458 ;velocity of light in km/sec 
  DCMR = 0.0     ;comoving radial distance in units of c/H0  
  a = 1.0        ;//# 1/(1+z), the scale factor of the Universe
  az = 0.5       ;//# 1/(1+z(object)) 
  az = 1.0/(1+1.0*redshift)
 
  n=1000   ;//# number of points in integrals
 

;do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
  for i=0, n  do begin    
      a = az+(1-az)*(i+0.5)/n
      adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))    
      DCMR = DCMR + 1./(a*adot);   
  endfor

 
  DCMR = (1.-az)*DCMR/n 
  DCMR_Mpc = (c/H0)*DCMR*Hubble_h

 return, DCMR_Mpc

end



pro plot_hist_2d,x,y,$
                 xbin=xbin,ybin=ybin,$
                 xmin=xmin,xmax=xmax,$
                 ymin=ymin,ymax=ymax,$
                 log=log,COLORS=COLORS,_extra=extra,hist=hist,$
                 xarray=xarray,yarray=yarray,$
                 point=point,psym=psym,pcol=pcol, $
                 levels=levels, cont=cont, frac=frac, weight=weight


; NAME: plot_hist_2d
;
; PURPOSE: plots a 2D histogram of data pairs with grey-scale, contour + points plot
;
;
; CATEGORY: graphics
;
;
; CALLING SEQUENCE:plot_hist_2d,x,y                 
;                 [xbin=xbin,ybin=ybin,$
;                  xmin=xmin,xmax=xmax,$
;                  ymin=ymin,ymax=ymax,$
;                  log=log,_extra=extra,hist=hist,$
;                  xarray=xarray,yarray=yarray,$
;                  point=point,psym=psym,pcol=pcol,cont=cont, frac=frac]
;
;  
;
; INPUTS: x,y : coordinates of points to be binned
;  
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;          xbin: size of bins in x axis (default 1)
;          ybin: size of bins in y axis (default 1)
;           log: if set then the histogram is logged before plotting
;       _ extra: any keywords accepted by contour
;          hist: returns the histogram
;        xarray: returns the x-axis for the histogram
;        yarray: returns the y-axis for the histogram
;         point: if set then points are plotted when histogram is below
;                point threshold
;          psym: symbol to be used for points
;          pcol: colour to be used for points
;          cont: if set then contour plot rather than grey-scale is used
;          frac: if cont and frac are set then contour levels are plotted
;                where density is at frac fractions of the total number of points
;




  xgood=where(finite(x) eq 1)
  ygood=where(finite(y) eq 1)

  if not keyword_set(xbin) then xbin=1
  if not keyword_set(ybin) then ybin=1
  if not keyword_set(xmin) then xmin=min(x[xgood])
  if not keyword_set(xmax) then xmax=max(x[xgood])
  if not keyword_set(ymin) then ymin=min(y[ygood])
  if not keyword_set(ymax) then yax=max(y[ygood])  
  if not keyword_set(psym) then psym=3
  if not keyword_set(pcol) then pcol=0  
  if not keyword_set(weight) then weight=1 

  hist=hist_2d(x,y,min1=xmin,min2=ymin,max1=xmax,max2=ymax,bin1=xbin,bin2=ybin)*weight


  sz=size(hist)
  xarray=findgen(sz[1])*xbin+xmin
  yarray=findgen(sz[2])*ybin+ymin


;----------------------------------------------------------------------
; plot contour
;----------------------------------------------------------------------


  if keyword_set(cont) then begin 

     if keyword_set(log) then begin
        ;print,'Plotting Log contour'
        if keyword_set(COLORS) then $
           contour,alog10(hist>1),xarray+xbin/2.0,yarray+ybin/2.0,_extra=extra,C_COLORS=COLORS,xst=1,yst=1,levels=levels $
        else $
           contour,alog10(hist>1),xarray+xbin/2.0,yarray+ybin/2.0,_extra=extra,xst=1,yst=1,levels=levels      
     endif else BEGIN
        ;print,'Plotting non-Log contour'
        if keyword_set(COLORS) then $
           contour,hist,xarray+xbin/2.0,yarray+ybin/2.0,C_COLORS=COLORS,_extra=extra,levels=levels $
        else $
           contour,hist,xarray+xbin/2.0,yarray+ybin/2.0,_extra=extra,levels=levels 
     endelse

  endif else begin
     
     ; a lot of this code taken from icplot
     ; Check if device has scalable pixels
     if ( !d.flags and 1 ) then scale = 1 else scale = 0

     ; Set the number of reserved pens
     if not keyword_set( respen ) then respen = 0

     message, 'this bit doesnt work properly - icplot needs fixing'
     if keyword_set(log) then begin
        icplot,alog10(hist>1),xax=xarray,yax=yarray,/pix,/fill,_extra=extra
     endif else begin
        icplot,hist,xax=xarray,yax=yarray,/pix ;,level=10.^(findgen(10)/3.+1),/fill,_extra=extra
     endelse

  endelse



;----------------------------------------------------------------------

  ix=fix((x-xmin)/xbin)
  iy=fix((y-ymin)/ybin)    
  den=hist(ix,iy)


;----------------------------------------------------------------------
; pick contour with levels set at given fractions of the total number
; of points
;----------------------------------------------------------------------
  
  if keyword_set(frac) and keyword_set(cont) then begin

     histsm=smooth(float(hist),3)
     densm=histsm(ix,iy)
     order=sort(densm)
     n=round(n_elements(x)*frac)
     den_frac=interpol(densm[order],findgen(n_elements(x)),n)    
 
     contour,histsm,xarray,yarray,/over,levels=den_frac,c_colo=3,xs=1,ys=1
     
  endif

;----------------------------------------------------------------------
; plot up points in sparse areas
;----------------------------------------------------------------------

if keyword_set(point) then begin     
    use=where(den le point,nuse)
    if nuse gt 0 then oplot,x[use]+xbin/2.0*0.0,y[use],psym=psym,col=pcol
 endif

 ;contour,alog10(hist>1),xarray+xbin/2.0,yarray+ybin/2.0,_extra=extra,xst=1,yst=1

end


pro invert_color_table

  tvlct, r_orig, g_orig, b_orig, /get
  r_new=fltarr(256)
  g_new=fltarr(256)
  b_new=fltarr(256)
  
  for i=0, 255, 1 do begin
     r_new[i]=r_orig[255-i]
     g_new[i]=g_orig[255-i]
     b_new[i]=b_orig[255-i]       
  endfor
      
  tvlct, r_new, g_new, b_new
end


PRO ROTATION,X,Y,DEG,NX,NY
ang=deg*!dtor

;convert to polar coordinates for rotation
r = sqrt(x*x + y*y)
theta = r*0.
;get angle in for loop so that zero radii will be left as zero angle
for i = 0,n_elements(r)-1 do $
if r[i] ne 0 then theta[i] = atan(y[i],x[i])  ;range from -pi to +pi
;
;add rotation angle
theta = theta + ang
;
;convert back to rectangular coordinates, now rotated
nx = r * cos(theta)
ny = r * sin(theta)
;
return
end
