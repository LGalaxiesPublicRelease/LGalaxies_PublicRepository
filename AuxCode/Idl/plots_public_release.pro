;=========================================================================

;Henriques2014
@LGalaxy_Henriques2014a
@LightGalaxy_Henriques2014a
@LGalaxy_tree_Henriques2014a

@procedures
;=========================================================================

  PRO sam, NN

  Gstruct={LGalaxy_Henriques2014a}; - paper - final code used for DB (no time on sfh)
  ;Gstruct={LightGalaxy_Henriques2014a};
  ;Gstruct={LGalaxy_tree_Henriques2014a}
 

  LIGHT=0

  ;FileName = 'SA_galtree'
  FileName = ''

  !Path = Expand_Path( '+./') 
 
   Datadir = './data/'
   DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/MR/'
   DirName_MRII = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/MR/'  

   mcmc_folder= '/net/bootes/export/data1/Workspace/PR_GitHub_Hen15/MCMC/'


   IF (not keyword_set(NN)) THEN NN = 0  
   FirstFile = 40
   LastFile = FirstFile + NN - 1
 
   WMAP1=0
   WMAP7=0 
   PLANCK=1
 
   MRII=0
 
   read=[1,1,1,1,1,0]

   default_plot_square=0
   reduced_plot_square=0
   extended_plot_square=1


   plot_kband_bband_vband_smf=1
   stellar_mass_vs_sfr=1     ;evolution contour plots
  
   color_vs_magnitude=1
   color_ssfr_age_hist=1

   UVJ_color=1
   red_fraction_colorcut=1
   stellar_mass_function_by_color=1

   morphology=1
   bhbm=1
   gas_mass_function=1
   mstar_metals=1
 


   prefix_this_model='This Work'
 
  
;WRITE_FILES
   write_files=1
   plot_after_write=0  
   file_to_write=Datadir+'Henriques2014c'
 
  
;PREVIOUS MODELS
  do_previous_model1=1
  file_previous_model1=DataDir+'Guo2013a_m05'
  prefix_previous_model1='Guo2013a - WMAP7'
  linestyle_previous_model1=1

  ;do_previous_model2=1
  ;file_previous_model2=DataDir+'Henriques2013a'
  ;prefix_previous_model2='Henriques2013a - WMAP7'
  ;linestyle_previous_model2=2

  do_previous_model2=1
  file_previous_model2=DataDir+'Henriques2015a'
  prefix_previous_model2='Henriques2015a - PLANCK'
  linestyle_previous_model2=2
   



;***********************************************
;*                                            ;*
;*   plot_kband_bband_vband_smf               ;*
;*                                            ;*
;***********************************************
;READ SAMPLE                                  ;*
   sample=0                                   ;*
   ;dir='/galformod/scratch/bmh20/Workspace/Development_Branch/'
   ;dir='/galformod/scratch/bmh20/Workspace/Henriques2014a/'
   ;dir='/galformod/scratch/bmh20/Workspace/nifty/mcmc_working/'
   dir='/net/bootes/export/data1/Workspace/PublicRelease_Henriques2015/'
   sample_file=dir+'output/mcmc_plus_obs'     ;* 
                                              ;*
;PLOT OBSERVATIONS                            ;*
   plot_obs=0                                  ;*
   plot_only_new_obs=0                        ;*
                                              ;*
;MISCELANEOUS                                 ;*
   do_lums=1  
   min_obs_err=0.1
   total_chi2=0.
                                              ;*
;2SGIMA REGIONS FROM MCMC                     ;*
   do_mcmc_regions=0                          ;*   
                                              ;*
;***********************************************

 slope_red_fraction=[0.1,0.275, 0.3, 0.35,0.35]
 offset_red_fraction=[1.93,1.213, 1.18,0.92,0.82]
 minimum_y_red_fraction=[0.0,1.3,1.3,1.3,1.3]


    

;***********************************************
;*                                            ;*
;*           plot_red_fraction                ;*
;*                                            ;*
;***********************************************
                                              ;*
model_err_red_frac=0.025                      ;*
                                              ;*
;***********************************************



;***********************************************
;*                                            ;*
;*           ssfr_and_age_hist                ;*
;*                                            ;*
;***********************************************
                                              ;*
ssfr_and_age_hist_multiple_plots=1            ;*
;min_z=0.005                                  ;*
;max_z=0.05                                   ;*
;min_z=0.01                                    ;*
;max_z=0.15                                    ;*
min_z=0.005                                    ;*
max_z=0.2    

;min_z=0.0                                    ;*
;max_z=0.06    
                                              ;*
;***********************************************




 
   if(WMAP1 eq 1) then begin
      BoxSize_MR    = 500.      ;full MR 
      BoxSize_MRII  = 100.      ;full MRII
      Hubble_h      = 0.732
   endif
   if(WMAP7 eq 1) then begin
      BoxSize_MR    = 521.555   ;full MR 
      ;BoxSize_MR    = 452.986
      BoxSize_MRII  = 104.3     ;full MRII 
      Hubble_h      = 0.704
   endif  
   if(PLANCK eq 1) then begin   
      BoxSize_MR    = 500.* 0.960558    ;full MR 
      BoxSize_MRII  = 100.* 0.960558    ;full MRII        
      Hubble_h      = 0.673
      Omega_M=0.315 
      Omega_Lambda=0.683
   endif
  

   MaxTreeFiles = 512        ;full MR  
   volume_MR = (BoxSize_MR^3.0) * (Lastfile - Firstfile + 1) / MaxTreeFiles 
   volume_MRII = (BoxSize_MRII^3.0) * (Lastfile - Firstfile + 1) / MaxTreeFiles 
 
   Hubble_h_wmap1 = 0.732
   hubble_h_wmap7 = 0.704 
 
   if(WMAP1 eq 1) then begin
      redshift=[0.0,0.4,1.0,2.0,3.0]
      FileName0 = 'SA_z0.00'
      ;FileName0 = 'SA_z0.09'
      FileName04 = 'SA_z0.41' 
      FileName1 = 'SA_z0.99'
      FileName2 = 'SA_z2.07'
      FileName3 = 'SA_z3.06'
      if(LIGHT eq 1) then begin 
         FileName0 = 'Light_SA_z0.00'
         ;FileName0 = 'Light_SA_z0.09'
         FileName04 = 'Light_SA_z0.51'
         FileName1 = 'Light_SA_z0.99'
         FileName2 = 'Light_SA_z2.07'
         FileName3 = 'Light_SA_z3.06'
      endif 
   endif

   if(WMAP7 eq 1) then begin
      ;redshift=[0.1,0.5,1.0,2.0,3.0]
      redshift=[0.0,0.4,1.0,2.0,3.0]
      FileName0 = 'SA_z0.00'    
      ;FileName0 = 'SA_z0.08'    
      FileName04 = 'SA_z0.39' 
      FileName1 = 'SA_z1.02'
      FileName2 = 'SA_z1.92'
      FileName3 = 'SA_z2.92'     
      if(LIGHT eq 1) then begin 
         FileName0 = 'Light_SA_z0.00'
         ;FileName0 = 'Light_SA_z0.08'
         FileName04 = 'Light_SA_z0.39'      
         FileName1 = 'Light_SA_z1.02'
         FileName2 = 'Light_SA_z1.92'
         FileName3 = 'Light_SA_z2.92'       
      endif 
   endif

   if(PLANCK eq 1) then begin    
         redshift=[0.0,0.4,1.0,2.0,3.0]    
         FileName0 = 'SA_z0.00'        
         ;FileName0 = 'SA_z0.11'
         FileName04 = 'SA_z0.40'   
         FileName1 = 'SA_z1.04'
         FileName2 = 'SA_z2.07'
         FileName3 = 'SA_z3.11' 
         FileName4 = 'SA_z3.95' 
         if(LIGHT eq 1) then begin  
            FileName0 = 'Light_SA_z0.00'
            ;FileName0 = 'Light_SA_z0.11'          
            FileName04 = 'Light_SA_z0.40'     
            FileName1 = 'Light_SA_z1.04'
            FileName2 = 'Light_SA_z2.07'
            FileName3 = 'Light_SA_z3.11'
            FileName4 = 'Light_SA_z3.95'
         endif  
      endif
  

      ;treeoutput
      if  (FileName eq 'SA_galtree') then begin  
 
         read_tree, tree_struct, Gstruct, DirName_MR, FileName, FirstFile, LastFile
         print,'tree read'
         if(convert_to_nifty eq 1) then convert_to_nifty_struct, tree_struct, DirName_MR, hubble_h
         print,'conversion done'

         G0_MR =tree_struct[where(tree_struct.snapnum eq 58)]
         ;G0_MR =tree_struct[where(tree_struct.snapnum eq 54)]  ;z=0.1
         G04_MR=tree_struct[where(tree_struct.snapnum eq 47)]
         G1_MR =tree_struct[where(tree_struct.snapnum eq 38)]
         G2_MR =tree_struct[where(tree_struct.snapnum eq 30)]
         G3_MR =tree_struct[where(tree_struct.snapnum eq 25)]
      
         if(MRII eq 1) then begin
            read_tree, tree_struct, Gstruct, DirName_MRII, FileName, FirstFile, LastFile
            
            ;read_tree, tree_struct, Gstruct, DirName_MR, FileName, FirstFile, LastFile
            G0_MRII =tree_struct[where(tree_struct.snapnum eq 62)]
            ;G0_MRII =tree_struct[where(tree_struct.snapnum eq 58)] ;z=0.1
            G04_MRII=tree_struct[where(tree_struct.snapnum eq 51)]
            G1_MRII =tree_struct[where(tree_struct.snapnum eq 42)]
            G2_MRII =tree_struct[where(tree_struct.snapnum eq 34)]
            G3_MRII =tree_struct[where(tree_struct.snapnum eq 29)]
         endif else begin
            G0_MRII = G0_MR
            G04_MRII = G0_MR
            G1_MRII = G0_MR
            G2_MRII = G0_MR
            G3_MRII = G0_MR
         endelse

;SNAP OUTPUT
      endif else begin

  ;read_snap_merged_files,G0_MR, Gstruct, DirName_MR, FileName0, FirstFile, LastFile

  ;Get LGalaxy structure from LGalaxy.pro 
  
         if(read[0] eq 1) then $       
            read_snap, G0_MR, Gstruct, DirName_MR, FileName0, FirstFile, LastFile         
         if(read[1] eq 1) then $       
         read_snap, G04_MR, Gstruct, DirName_MR, FileName04, FirstFile, LastFile
         if(read[2] eq 1) then $
         read_snap, G1_MR, Gstruct, DirName_MR, FileName1, FirstFile, LastFile
         if(read[3] eq 1) then $
         read_snap, G2_MR, Gstruct, DirName_MR, FileName2, FirstFile, LastFile
         if(read[4] eq 1) then $
         read_snap, G3_MR, Gstruct, DirName_MR, FileName3, FirstFile, LastFile
         if(read[5] eq 1) then $
            read_snap, G4_MR, Gstruct, DirName_MR, FileName4, FirstFile, LastFile
 
         if(MRII) then begin
            if(read[0] eq 1) then $        
               read_snap, G0_MRII, Gstruct, DirName_MRII, FileName0, FirstFile, LastFile        
            if(read[1] eq 1) then $         
            read_snap, G04_MRII, Gstruct, DirName_MRII, FileName04, FirstFile, LastFile
            if(read[2] eq 1) then $
            read_snap, G1_MRII, Gstruct, DirName_MRII, FileName1, FirstFile, LastFile
            if(read[3] eq 1) then $
            read_snap, G2_MRII, Gstruct, DirName_MRII, FileName2, FirstFile, LastFile
            if(read[4] eq 1) then $
            read_snap, G3_MRII, Gstruct, DirName_MRII, FileName3, FirstFile, LastFile 
            if(read[5] eq 1) then $
            read_snap, G4_MRII, Gstruct, DirName_MRII, FileName4, FirstFile, LastFile
         endif else begin
            if(read[0] eq 1) then G_MRII=G0_MR
            if(read[0] eq 1) then G0_MRII=G0_MR
            if(read[1] eq 1) then G04_MRII=G04_MR
            if(read[2] eq 1) then G1_MRII=G1_MR
            if(read[3] eq 1) then G2_MRII=G2_MR
            if(read[4] eq 1) then G3_MRII=G3_MR
            if(read[5] eq 1) then G4_MRII=G4_MR          
         endelse


      endelse ;snap output


   if(read[0] eq 1) then begin     
      stellarmass_MR=(G0_MR.BulgeMass+G0_MR.DiskMass)*1.e10*hubble_h
      stellarmass_MRII=(G0_MRII.BulgeMass+G0_MRII.DiskMass)*1.e10*hubble_h      
      StellarMass_MR_0=10^(Alog10(stellarmass_MR)+randomn(1l,n_elements(stellarmass_MR))*0.08*(1+redshift[0]))
      StellarMass_MRII_0=10^(Alog10(stellarmass_MRII)+randomn(1l,n_elements(stellarmass_MRII))*0.08*(1+redshift[0]))
   endif
   if(read[1] eq 1) then begin
      stellarmass_MR=(G04_MR.BulgeMass+G04_MR.DiskMass)*1.e10*hubble_h
      stellarmass_MRII=(G04_MRII.BulgeMass+G04_MRII.DiskMass)*1.e10*hubble_h  
      StellarMass_MR_04=10^(Alog10(stellarmass_MR)+randomn(1l,n_elements(stellarmass_MR))*0.08*(1+redshift[1]))
      StellarMass_MRII_04=10^(Alog10(stellarmass_MRII)+randomn(1l,n_elements(stellarmass_MRII))*0.08*(1+redshift[1]))
   endif
   if(read[2] eq 1) then begin
      stellarmass_MR=(G1_MR.BulgeMass+G1_MR.DiskMass)*1.e10*hubble_h
      stellarmass_MRII=(G1_MRII.BulgeMass+G1_MRII.DiskMass)*1.e10*hubble_h  
      StellarMass_MR_1=10^(Alog10(stellarmass_MR)+randomn(1l,n_elements(stellarmass_MR))*0.08*(1+redshift[2]))
      StellarMass_MRII_1=10^(Alog10(stellarmass_MRII)+randomn(1l,n_elements(stellarmass_MRII))*0.08*(1+redshift[2]))    
   endif
   if(read[3] eq 1) then begin
      stellarmass_MR=(G2_MR.BulgeMass+G2_MR.DiskMass)*1.e10*hubble_h
      stellarmass_MRII=(G2_MRII.BulgeMass+G2_MRII.DiskMass)*1.e10*hubble_h  
      StellarMass_MR_2=10^(Alog10(stellarmass_MR)+randomn(1l,n_elements(stellarmass_MR))*0.08*(1+redshift[3]))
      StellarMass_MRII_2=10^(Alog10(stellarmass_MRII)+randomn(1l,n_elements(stellarmass_MRII))*0.08*(1+redshift[3]))  
   endif
   if(read[4] eq 1) then begin
      stellarmass_MR=(G3_MR.BulgeMass+G3_MR.DiskMass)*1.e10*hubble_h
      stellarmass_MRII=(G3_MRII.BulgeMass+G3_MRII.DiskMass)*1.e10*hubble_h  
      StellarMass_MR_3=10^(Alog10(stellarmass_MR)+randomn(1l,n_elements(stellarmass_MR))*0.08*(1+redshift[4]))
      StellarMass_MRII_3=10^(Alog10(stellarmass_MRII)+randomn(1l,n_elements(stellarmass_MRII))*0.08*(1+redshift[4]))  
   endif
   if(read[5] eq 1) then begin
      stellarmass_MR=(G4_MR.BulgeMass+G4_MR.DiskMass)*1.e10*hubble_h
      stellarmass_MRII=(G4_MRII.BulgeMass+G4_MRII.DiskMass)*1.e10*hubble_h  
      StellarMass_MR_4=10^(Alog10(stellarmass_MR)+randomn(1l,n_elements(stellarmass_MR))*0.08*(1+redshift[5]))
      StellarMass_MRII_4=10^(Alog10(stellarmass_MRII)+randomn(1l,n_elements(stellarmass_MRII))*0.08*(1+redshift[5]))  
   endif


   print,''
   print,''
   print,''
   print,' ******************************************'
   print,' *                                        *'
   print,' *        All catalogs loaded             *' 
   print,' *         Starting Analysis              *' 
   print,' *                                        *'
   print,' ******************************************'
   print,''
   print,''

;-------------------------------------------------------------------------------
 

;********************************************************************************
;********************************************************************************
;********************************************************************************
;********************************************************************************
;********************************************************************************
;********************************************************************************

  ; Open ps file for all plots
   loadct,6, /SILENT
   p_charsize = 1.4
   p_charthick = 4
   p_thick = 4
   x_thick = 4
   y_thick = 4
   !p.charsize = p_charsize
   !p.charthick = p_charthick
   !p.thick = p_thick
   !x.thick = x_thick
   !y.thick = y_thick
   set_plot, 'PS'
    
   device, filename = './read_gal.ps', xsize = 20, ysize = 20, /color, xoffset=1, yoffset=5
 
   if(default_plot_square eq 1) then $
      device, filename = './read_gal.ps', xsize = 25, ysize = 20, /color, xoffset=1, yoffset=5   
   
   if(reduced_plot_square eq 1) then $
      device, filename = './read_gal.ps', xsize = 18, ysize = 14, /color, xoffset=1, yoffset=5
   
   if(extended_plot_square eq 1) then $
      device, filename = './read_gal.ps', xsize = 30, ysize = 20, /color, xoffset=1, yoffset=5
   



;KBAND+BBAND+VBAND+SMF
if(plot_kband_bband_vband_smf eq 1) then begin
  multiplot, /reset 
  erase & multiplot, [2,2] 
 
  plot_kband_bband_smf, G0_MR, G1_MR, G2_MR, G3_MR, G0_MRII, G1_MRII, G2_MRII, G3_MRII, $ 
                        StellarMass_MR_0, StellarMass_MRII_0, StellarMass_MR_04, StellarMass_MRII_04, $
                        StellarMass_MR_1, StellarMass_MRII_1, StellarMass_MR_2, StellarMass_MRII_2, $
                        StellarMass_MR_3, StellarMass_MRII_3, total_chi2, $
                        write_files, plot_after_write, file_to_write, $
                        sample, sample_file, plot_obs, plot_only_new_obs, do_lums, $
                        linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                        linestyle_previous_model2, file_previous_model2, do_previous_model2, $
                        prefix_previous_model1, prefix_previous_model2, prefix_this_model, $
                        MRII, volume_MR, volume_MRII, hubble_h, $
                        Datadir, do_mcmc_regions, mcmc_folder, min_obs_err, redshift
 
  multiplot, /reset 
endif







if(stellar_mass_vs_sfr eq 1) then begin

   print,''
   print,''
   print,''
   print,'*************************'
   print,'*      SFR vs Mass      *'
   print,'*************************'
   print,''

   multiplot,/reset
   erase & multiplot, [5,1]

   yobs =[9.3,8.7,8.1,7.5,6.9,6.3]
   yobs2=[9.45,8.85,8.25,7.65,7.15,6.45]

   xmin=8.6
   xmax=11.9
   ymin=-3.0
   ymax=4.5  
   bin=0.25


   N_redshifts=5
   for i=0, N_redshifts-1 do char_redshift=STRTRIM(number_formatter(redshift, DECIMAL=2),2)

   get_slope, 8.0, -1.8, 12.6, 2.7, slope0,b0   
   get_slope, 8.0, -0.7, 12.6, 3.8, slope1,b1     
   get_slope, 8.0, -0.4, 12.0, 3.5, slope2,b2

   obs_slope_elbaz2007 =[0.77, -99.0, 0.9, -99.0, -99.0]
   obs_offset_elbaz2007=[Alog10(8.7)-(0.77*11.), -99.0, Alog10(7.2)-9, -99.0, -99.0]
   obs_offset_low_elbaz2007=[Alog10(5.0)-(0.77*11.), -99.0, Alog10(3.6)-9, -99.0, -99.0]
   obs_offset_high_elbaz2007=[Alog10(16.1)-(0.77*11.), -99.0, Alog10(14.4)-9, -99.0, -99.0]
 
   obs_slope_daddi2007=[-99.0, -99.0, -99.0, 0.9, -99.0]
   obs_offset_daddi2007=[-99.0, -99.0, -99.0, -7.6, -99.0]

   readcol,Datadir+'karim2011_sfr_mass_sf.txt', karim_low_mass_limit, karim_high_mass_limit, karim_medium_mass, $
           karim_low_z_limit, karim_high_z_limit, karim_medium_z, peak_flux_14, upper_peak_flux_14, lower_peak_flux_14, $
           total_flux_14, upper_total_flux_14, lower_total_flux_14, back_noise, upper_back_noise, lower_back_noise, lum_14, $
           upper_lum_14, lower_lum_14, karim_sfr, karim_sfr_error_up, karim_sfr_error_down, /SILENT


    make_float_array, obs_x, xmin, xmax, 0.01

    for i_redshift=0, N_redshifts-1 do begin

       if(i_redshift gt 0) then multiplot
       if(i_redshift eq 0) then y_label='log!D10!N(SFR[h!U-2!NM!D!9ng!3!Nyr!U-1!N])' else y_label=''

       if(i_redshift eq 0) then begin
          Sfr=G0_MR.Sfr*hubble_h^2
          Mass=StellarMass_MR_0
          type=G0_MR.type
       endif
       if(i_redshift eq 1) then begin
          Sfr=G04_MR.Sfr*hubble_h^2
          Mass=StellarMass_MR_04
      endif
       if(i_redshift eq 2) then begin
          Sfr=G1_MR.Sfr*hubble_h^2
          Mass=StellarMass_MR_1
       endif
       if(i_redshift eq 3) then begin
          Sfr=G2_MR.Sfr*hubble_h^2
          Mass=StellarMass_MR_2
       endif
       if(i_redshift eq 4) then begin
          Sfr=G3_MR.Sfr*hubble_h^2
          Mass=StellarMass_MR_3
       endif

       Sfr_new=Sfr

       ;assign threshold SFR to galaxies with SFR=0
       slope1=-0.3
       b1=-8.6
       width1=0.5
       
       slope2=-0.5
       b2=-6.5
       width2=1.0     

       sel=where(sfr/mass lt 10^(-15.) and Alog10(mass) ge 10.0,n)
       if(n gt 0) then $
          sfr_new[sel]=randomn(1l,n_elements(sfr[sel]))*width1*10^(slope1*Alog10(mass)+b1) + 10^(slope1*Alog10(mass)+b1)
       sel=where(sfr/mass lt 10^(-15.) and Alog10(mass) lt 10.0,n)
       if(n gt 0) then $
          sfr_new[sel]=randomn(1l,n_elements(sfr[sel]))*width2*10^(slope2*Alog10(mass)+b2) + 10^(slope2*Alog10(mass)+b2)


       ;PLOT MODEL
       XTICKS=2
       XTICKV = [xmin+0.5,xmin+1.5,xmin+2.5]
       XTICKNAME = ['9.0','10.0','11.0']
       XMINOR=10
       loadct,0, /SILENT    
       levels=findgen(15)*0.2+1.5
     
       plot_hist_2d,Alog10(Mass),Alog10(SFR_new),xmin = xmin-1.0,xmax=xmax+1.0, ymin = ymin-1.0, ymax=ymax+1.0, $
                    xstyle = 1, ystyle = 1, xtitle ='log!D10!N(M!D*!N[h!U-2!NM!D!9ng!3!N])', ytitle = y_label, $
                    xbin=0.2,ybin=0.2,/cont,/log, XTICKS = XTICKS, XTICKV=XTICKV, $
                    XTICKNAME = XTICKNAME,  XMINOR=XMINOR, $
                    ;COLORS=[0,20,40,60,80,100,120,140,160,180,200,220,240,260,280], $
                    levels=levels, xrange=[xmin,xmax],yrange=[ymin,ymax],xs=1,ys=1,/fill
 
       if(i_redshift eq 4) then begin   
          COLORBAR, BOTTOM=0, VERTICAL=0, POSITION=[0.87, 0.38, 0.96, 0.43], $
                    MAX=0.,Min=alog10(10^min(levels)/10^max(levels)),DIVISIONS=3, FORMAT='(F6.1)'
  
          ;needed to re-set keywords
          plot_hist_2d,Alog10(Mass),Alog10(SFR_new),xmin = xmin-1.0,xmax=xmax+1.0, ymin = ymin-1.0, ymax=ymax+1.0, $
                       xstyle = 1, ystyle = 1, xtitle ='log!D10!N(M!D*!N[h!U-2!NM!D!9ng!3!N])', ytitle = y_label, $
                       xbin=0.2,ybin=0.2,/cont,/log, XTICKS = XTICKS, XTICKV=XTICKV, $
                       XTICKNAME = XTICKNAME,  XMINOR=XMINOR, $        
                       levels=levels, xrange=[xmin,xmax],yrange=[ymin,ymax],xs=1,ys=1,/fill
 
       endif   
     
       ;PLOT OBSERVATIONS
       loadct,6, /SILENT
       if(i_redshift eq 0) then begin 

          ;salim SDSS
          file=Datadir+'dr7_gal_final.fit' 
          close,1
          openr,1,file
          dr7_gal_final = MRDFITS(1,1)
          close,1

          bin=0.25
          sel=where(dr7_gal_final.z gt min_z and dr7_gal_final.z lt max_z,nii)
          plot_obs_sfr_vmax_contour, dr7_gal_final[sel], 7.5, 12.0, bin, ymin, ymax, bin

        
          oplot,obs_x+2.*Alog10(hubble_h_WMAP7), $
                obs_x*obs_slope_elbaz2007[i_redshift]+obs_offset_elbaz2007[i_redshift]+2.*Alog10(hubble_h_WMAP7),color=30
          oplot,obs_x+2.*Alog10(hubble_h_WMAP7), $
                obs_x*obs_slope_elbaz2007[i_redshift]+obs_offset_low_elbaz2007[i_redshift]+2.*Alog10(hubble_h_WMAP7),color=30, linestyle=2
          oplot,obs_x+2.*Alog10(hubble_h_WMAP7), $
                obs_x*obs_slope_elbaz2007[i_redshift]+obs_offset_high_elbaz2007[i_redshift]+2.*Alog10(hubble_h_WMAP7),color=30, linestyle=2

          plot_label, xlog=0, ylog=0, type='label', label= 'Elbaz et al. (2007)', xmin, xmax, ymin, ymax, $
                      x_percentage=1.5, x2_percentage=0., y_percentage=yobs[0], $
                      charthick=4., charsize=0.95
          plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                      x_percentage=0.4, x2_percentage=1.3, y_percentage=yobs2[0], $
                      linestyle=0, color=30, linethick=8

          plot_label, xlog=0, ylog=0, type='label', label= 'SDSS/DR7', xmin, xmax, ymin, ymax, $
                      x_percentage=1.5, x2_percentage=0., y_percentage=yobs[1], $
                      charthick=4., charsize=0.95
          plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                      x_percentage=0.4, x2_percentage=1.3, y_percentage=yobs2[1], $
                      linestyle=0, color=0, linethick=8

       endif
               
       log_karim_sfr_error_up=Alog10((karim_sfr+karim_sfr_error_up)/karim_sfr)
       log_karim_sfr_error_down=Alog10(karim_sfr/(karim_sfr-karim_sfr_error_down))

       if(i_redshift eq 1) then begin        
          sel=where(karim_low_z_limit eq 0.2 and karim_medium_mass gt 8.8)        
          symbols,1,1.0
          oploterror, karim_medium_mass[sel]+2.*Alog10(hubble_h_WMAP7), $
                      Alog10(karim_sfr[sel])+2.*Alog10(hubble_h_WMAP7), log_karim_sfr_error_up, $
                      color = 70, errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
          oploterror, karim_medium_mass[sel]+2.*Alog10(hubble_h_WMAP7), $
                      Alog10(karim_sfr[sel])+2.*Alog10(hubble_h_WMAP7), log_karim_sfr_error_down, $
                      color = 70, errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
          plot_label, xlog=0, ylog=0, type='label', label='Karim (2011) - 0.2<z<0.4', xmin, xmax, ymin, ymax, $
                      x_percentage=1.0, x2_percentage=0., y_percentage=yobs[0], charthick=4., charsize=0.95
          plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                      x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[0], $
                      color=70, errcolor=70, sym_num=1, sym_size=1.0, HATLENGTH = 100, err_size=0.1

          sel=where(karim_low_z_limit eq 0.4 and karim_medium_mass gt 8.9)
          symbols,2,1.0
          oploterror, karim_medium_mass[sel]+2.*Alog10(hubble_h_WMAP7), $
                      Alog10(karim_sfr[sel])+2.*Alog10(hubble_h_WMAP7), log_karim_sfr_error_up, $
                      color = 70, errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
          oploterror, karim_medium_mass[sel]+2.*Alog10(hubble_h_WMAP7), $
                      Alog10(karim_sfr[sel])+2.*Alog10(hubble_h_WMAP7), log_karim_sfr_error_down, $
                      color = 70, errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
          plot_label, xlog=0, ylog=0, type='label', label='Karim (2011) - 0.4<z<0.6', xmin, xmax, ymin, ymax, $
                      x_percentage=1.0, x2_percentage=0., y_percentage=yobs[1], charthick=4., charsize= 0.95
          plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                      x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[1], $
                      color=70, errcolor=70, sym_num=2, sym_size=1.0, HATLENGTH = 100, err_size=0.1    
      endif

       if(i_redshift eq 2) then begin
          oplot,obs_x+2.*Alog10(hubble_h_WMAP7), $
                obs_x*obs_slope_elbaz2007[i_redshift]+obs_offset_elbaz2007[i_redshift]+2.*Alog10(hubble_h_WMAP7),color=30
       
          sel=where(karim_low_z_limit eq 0.8 and karim_medium_mass gt 9.1)
          symbols,1,1.0
          oploterror, karim_medium_mass[sel]+2.*Alog10(hubble_h_WMAP7), $
                      Alog10(karim_sfr[sel])+2.*Alog10(hubble_h_WMAP7), log_karim_sfr_error_up, $
                      color = 70, errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
          oploterror, karim_medium_mass[sel]+2.*Alog10(hubble_h_WMAP7), $
                      Alog10(karim_sfr[sel])+2.*Alog10(hubble_h_WMAP7), log_karim_sfr_error_down, $
                      color = 70, errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0      
          plot_label, xlog=0, ylog=0, type='label', label='Karim (2011) - 0.8<z<1.0', xmin, xmax, ymin, ymax, $
                      x_percentage=1.0, x2_percentage=0., y_percentage=yobs[0], charthick=4., charsize=0.95
          plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                      x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[0], $
                      color=70, errcolor=70, sym_num=1, sym_size=1.0, HATLENGTH = 100, err_size=0.1

          sel=where(karim_low_z_limit eq 1.0 and karim_medium_mass gt 9.3)
          symbols,2,1.0
          oploterror, karim_medium_mass[sel]+2.*Alog10(hubble_h_WMAP7), $
                      Alog10(karim_sfr[sel])+2.*Alog10(hubble_h_WMAP7), log_karim_sfr_error_up, $
                      color = 70, errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
          oploterror, karim_medium_mass[sel]+2.*Alog10(hubble_h_WMAP7), $
                      Alog10(karim_sfr[sel])+2.*Alog10(hubble_h_WMAP7), log_karim_sfr_error_down, $
                      color = 70, errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0  
          plot_label, xlog=0, ylog=0, type='label', label='Karim (2011) - 1.0<z<1.2', xmin, xmax, ymin, ymax, $
                      x_percentage=1.0, x2_percentage=0., y_percentage=yobs[1], charthick=4., charsize= 0.95
          plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                      x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[1], $
                      color=70, errcolor=70, sym_num=2, sym_size=1.0, HATLENGTH = 100, err_size=0.1

          plot_label, xlog=0, ylog=0, type='label', label= 'Elbaz et al. (2007)', xmin, xmax, ymin, ymax, $
                      x_percentage=1.5, x2_percentage=0., y_percentage=yobs[4], $
                      charthick=4., charsize=0.95
          plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                      x_percentage=0.4, x2_percentage=1.3, y_percentage=yobs2[4]-0.1, $
                      linestyle=0, color=30, linethick=8
 

   ;whitaker2013
          symbols,20,0.3
          readcol,Datadir+'whitaker2013_mass_vs_sfr_z0.5_1.0.txt', log_whitaker_mass, log_whitaker_sfr, log_whitaker_sfr_err, /SILENT
          log_obs_mass=log_whitaker_mass+2.*Alog10(hubble_h_WMAP7)
          log_obs_sfr=log_whitaker_sfr+2.*Alog10(hubble_h_WMAP7)
          oploterror, log_obs_mass, log_obs_sfr, log_whitaker_sfr_err, color = 200, errcolor = 200, psym = 8, HATLENGTH = 80.0

          symbols,30,0.3
          readcol,Datadir+'whitaker2013_mass_vs_sfr_z1.0_1.5.txt', log_whitaker_mass, log_whitaker_sfr, log_whitaker_sfr_err, /SILENT
          log_obs_mass=log_whitaker_mass+2.*Alog10(hubble_h_WMAP7)
          log_obs_sfr=log_whitaker_sfr+2.*Alog10(hubble_h_WMAP7)
          oploterror, log_obs_mass, log_obs_sfr, log_whitaker_sfr_err, color = 200, errcolor = 200, psym = 8, HATLENGTH = 80.0
       
          plot_label, xlog=0, ylog=0, type='label', label='Whitaker (2014) - 0.5<z<1.0', xmin, xmax, ymin, ymax, $
                      x_percentage=1.0, x2_percentage=0., y_percentage=yobs[2], charthick=4., charsize= 0.95
          plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                      x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[2], $
                      color=200, errcolor=200, sym_num=20, sym_size=0.3, HATLENGTH = 100, err_size=0.1
          plot_label, xlog=0, ylog=0, type='label', label='Whitaker (2014) - 1.0<z<1.5', xmin, xmax, ymin, ymax, $
                      x_percentage=1.0, x2_percentage=0., y_percentage=yobs[3], charthick=4., charsize= 0.95
          plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                      x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[3], $
                      color=200, errcolor=200, sym_num=30, sym_size=0.3, HATLENGTH = 100, err_size=0.1

      endif

       if(i_redshift eq 3) then begin         
          oplot,obs_x+2.*Alog10(hubble_h_WMAP7), $
                obs_x*obs_slope_daddi2007[i_redshift]+obs_offset_daddi2007[i_redshift]+2.*Alog10(hubble_h_WMAP7),color=120

          sel=where(karim_low_z_limit eq 1.6 and karim_medium_mass gt 9.6)
          symbols,1,1.0
          oploterror, karim_medium_mass[sel]+2.*Alog10(hubble_h_WMAP7), $
                      Alog10(karim_sfr[sel])+2.*Alog10(hubble_h_WMAP7), log_karim_sfr_error_up, $
                      color = 70, errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
          oploterror, karim_medium_mass[sel]+2.*Alog10(hubble_h_WMAP7), $
                      Alog10(karim_sfr[sel])+2.*Alog10(hubble_h_WMAP7), log_karim_sfr_error_down, $
                      color = 70, errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
          plot_label, xlog=0, ylog=0, type='label', label='Karim (2011) - 1.6<z<2.0', xmin, xmax, ymin, ymax, $
                      x_percentage=1.0, x2_percentage=0., y_percentage=yobs[0], charthick=4., charsize= 0.95
          plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                      x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[0], $
                      color=70, errcolor=70, sym_num=1, sym_size=1.0, HATLENGTH = 100, err_size=0.1

          sel=where(karim_low_z_limit eq 2.0 and karim_medium_mass gt 9.8)
          symbols,2,1.0
          oploterror, karim_medium_mass[sel]+2.*Alog10(hubble_h_WMAP7), $
                      Alog10(karim_sfr[sel])+2.*Alog10(hubble_h_WMAP7), log_karim_sfr_error_up, $
                      color = 70, errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
          oploterror, karim_medium_mass[sel]+2.*Alog10(hubble_h_WMAP7), $
                      Alog10(karim_sfr[sel])+2.*Alog10(hubble_h_WMAP7), log_karim_sfr_error_down, $
                      color = 70, errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0
          plot_label, xlog=0, ylog=0, type='label', label='Karim (2011) - 2.0<z<2.5', xmin, xmax, ymin, ymax, $
                      x_percentage=1.0, x2_percentage=0., y_percentage=yobs[1], charthick=4., charsize= 0.95
          plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                      x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[1], $
                      color=70, errcolor=70, sym_num=2, sym_size=1.0, HATLENGTH = 100, err_size=0.1

          plot_label, xlog=0, ylog=0, type='label', label='Daddi et al. (2007)', xmin, xmax, ymin, ymax, $
                      x_percentage=1.5, x2_percentage=0., y_percentage=yobs[4], charthick=4., charsize= 0.95                
          plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                      x_percentage=0.4, x2_percentage=1.3, y_percentage=yobs2[4]-0.1, $
                      linestyle=0, color=120, linethick=8                 
  

          symbols,20,0.3
          readcol,Datadir+'whitaker2013_mass_vs_sfr_z1.5_2.0.txt', log_whitaker_mass, log_whitaker_sfr, log_whitaker_sfr_err, /SILENT
          log_obs_mass=log_whitaker_mass+2.*Alog10(hubble_h_WMAP7)
          log_obs_sfr=log_whitaker_sfr+2.*Alog10(hubble_h_WMAP7)
          oploterror, log_obs_mass, log_obs_sfr, log_whitaker_sfr_err, color = 200, errcolor = 200, psym = 8, HATLENGTH = 80.0

          symbols,30,0.3
          readcol,Datadir+'whitaker2013_mass_vs_sfr_z2.0_2.5.txt', log_whitaker_mass, log_whitaker_sfr, log_whitaker_sfr_err, /SILENT
          log_obs_mass=log_whitaker_mass+2.*Alog10(hubble_h_WMAP7)
          log_obs_sfr=log_whitaker_sfr+2.*Alog10(hubble_h_WMAP7)
          oploterror, log_obs_mass, log_obs_sfr, log_whitaker_sfr_err, color = 200, errcolor = 200, psym = 8, HATLENGTH = 80.0
       

          plot_label, xlog=0, ylog=0, type='label', label='Whitaker (2014) - 1.5<z<2.0', xmin, xmax, ymin, ymax, $
                      x_percentage=1.0, x2_percentage=0., y_percentage=yobs[2], charthick=4., charsize= 0.95
          plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                      x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[2], $
                      color=200, errcolor=200, sym_num=20, sym_size=0.3, HATLENGTH = 100, err_size=0.1
          plot_label, xlog=0, ylog=0, type='label', label='Whitaker (2014) - 2.0<z<2.5', xmin, xmax, ymin, ymax, $
                      x_percentage=1.0, x2_percentage=0., y_percentage=yobs[3], charthick=4., charsize= 0.95
          plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                      x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[3], $
                      color=200, errcolor=200, sym_num=30, sym_size=0.3, HATLENGTH = 100, err_size=0.1

     endif

       if(i_redshift eq 4) then begin    
          sel=where(karim_low_z_limit eq 2.5 and karim_medium_mass gt 10.0)
          symbols,1,1.0
          oploterror, karim_medium_mass[sel]+2.*Alog10(hubble_h_WMAP7), $
                      Alog10(karim_sfr[sel])+2.*Alog10(hubble_h_WMAP7), log_karim_sfr_error_up, $
                      color = 70, errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
          oploterror, karim_medium_mass[sel]+2.*Alog10(hubble_h_WMAP7), $
                      Alog10(karim_sfr[sel])+2.*Alog10(hubble_h_WMAP7), log_karim_sfr_error_down, $
                      color = 70, errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0  
          plot_label, xlog=0, ylog=0, type='label', label='Karim (2011) - 2.5<z<3.0', xmin, xmax, ymin, ymax, $
                      x_percentage=1.0, x2_percentage=0., y_percentage=yobs[0], charthick=4., charsize= 0.95
          plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
                      x_percentage=0.6, x2_percentage=0., y_percentage=yobs2[0], $
                      color=70, errcolor=70, sym_num=1, sym_size=1.0, HATLENGTH = 100, err_size=0.1
       endif

       plot_label, xlog=0, ylog=0, type='label', label= 'z='+char_redshift[i_redshift], xmin, xmax, ymin, ymax, $
                      x_percentage=6.8, x2_percentage=0., y_percentage=0.5, $
                      charthick=5., charsize=1.4


    endfor

    multiplot, /reset 

 endif



if(color_vs_magnitude eq 1) then begin

   print,''
   print,''
   print,''
   print,'*************************'
   print,'*  Colour vs Magnitude  *'
   print,'*************************'
   print,''

   multiplot, /reset
 
   !p.multi=0

   G=G0_MRII
      
   color=G.MagDust[15]-G.MagDust[17]
   mass=Alog10(StellarMass_MRII_0)
   xmin=-23.0
   xmax=-14.0
   ymin=0.7
   ymax=2.7

   loadct,0, /SILENT
 
   plot_hist_2d,G.MagDust[17],color,xmin = xmin-1.0,xmax=xmax+1.0, ymin = ymin-1.0, ymax=ymax+1.0, $
                xstyle = 1, ystyle = 1, xtitle ='r', ytitle ='u-r',$
                xbin=0.25,ybin=0.05,/cont,/log, levels=findgen(30)*0.1+0.5, $
                xrange=[xmin,xmax],yrange=[ymin,ymax],xs=1,ys=1,/fill
  
   loadct,6, /SILENT  

   ;baldry cut at z=0 u-r vs r
   x=findgen(1000.)*(xmax-xmin)/1000.+xmin
   y=(offset_red_fraction[0]-slope_red_fraction[0]*tanh((x+18.07)/1.09))
   oplot,x,y,color=200,thick=8

endif


  ;COLOR+SSFR+AGE HIST
if( color_ssfr_age_hist eq 1) then begin

   Nbins=8
   lower_mass_limits=8.0+indgen(Nbins)*0.5
   upper_mass_limits=8.5+indgen(Nbins)*0.5
 


;AGE 

   print,''
   print,''
   print,''
   print,'************************'
   print,'*       Age Hist       *'
   print,'************************'
   print,''
  

   multiplot, /reset 
   erase & multiplot, [4,2]

   ymin=0.
   ymax=0.5
   xmin=0.0
   xmax=11.5    
   bin=0.5  

;DR4 - GALLAZZI - AGES
   file=Datadir+'dr4_gal_final.fit' 
   close,1
   openr,1,file
   dr4_gal_final = MRDFITS(1,1)
   close,1
   dr4_gal_final=dr4_gal_final[where(dr4_gal_final.mass gt lower_mass_limits[0])]


   for i=0,Nbins-1 do begin
        
      if(i gt 0) then multiplot    

      ;AXIS LABELS
      if(i eq 0 or i eq 4) then label_y='fraction' else label_y=''
      if(i gt 3) then label_x='r-Weighted Age (Gyr)' else label_x=''   
    
      if(lower_mass_limits[i] lt 9.5) then begin ;MRII
         mass=StellarMass_MRII_0  
         type=G0_MRII.type
         age= G0_MRII.rbandWeightAge        
      endif else begin ;MR
         mass=StellarMass_MR_0  
         type=G0_MR.type
         age= G0_MR.rbandWeightAge      
      endelse 

      age_hist, dr4_gal_final, i, mass, type, age, $                 
                lower_mass_limits[i], upper_mass_limits[i], redshift=0.1, $
                xmin, xmax, ymin, ymax, bin, Datadir,label_x,  label_y , $
                write_files, plot_after_write, file_to_write, $
                ssfr_and_age_hist_multiple_plots, min_z, max_z, $
                linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                linestyle_previous_model2, file_previous_model2, do_previous_model2, $
                prefix_previous_model1, prefix_previous_model2, prefix_this_model
     
   endfor




;SSFR

   multiplot, /reset
   erase & multiplot, [4,2]

   ymin=0.
   ymax=0.42
   xmin=-13.5
   xmax=-8.5
   bin=0.1

   print,''
   print,''
   print,''
   print,'*************************'
   print,'*       SSFR Hist       *'
   print,'*************************'
   print,''
  

;DR7
   file=Datadir+'dr7_gal_final.fit' 
   close,1
   openr,1,file
   dr7_gal_final = MRDFITS(1,1)
   close,1  
   dr7_gal_final=dr7_gal_final[where(dr7_gal_final.jarle_median_mass   gt lower_mass_limits[0])]
 
   ;lower limit doesn't matter
   z_min=[0.005,0.005,0.005, 0.005,0.005,0.005, 0.005,0.005]
   z_max=[0.15,0.15,0.15, 0.15,0.15,0.15, 0.25,0.25]
 

   for i=0,Nbins-1 do begin
     
      if(i gt 0) then multiplot

      ;AXIS LABELS
      if(i eq 0 or i eq 4) then label_y='fraction' else label_y=''
      if(i gt 3) then label_x='log!D10!N(SSFR)[yr!U-1!N]' else label_x=''
              
      if(lower_mass_limits[i] lt 9.5) then begin 
         mass=StellarMass_MRII_0
         ssfr=G0_MRII.sfr*hubble_h^2/(StellarMass_MRII_0)
         type=G0_MRII.type  
      endif else begin 
         mass=StellarMass_MR_0
         ssfr=G0_MR.sfr*hubble_h^2/(StellarMass_MR_0)
         type=G0_MR.type  
      endelse
     
      ssfr_hist, dr7_gal_final, i, mass, ssfr, type,$
                 lower_mass_limits[i], upper_mass_limits[i], redshift=0.1, $
                 xmin, xmax, ymin, ymax, bin, Datadir,label_x,label_y, $
                 write_files, plot_after_write, file_to_write , $
                 ssfr_and_age_hist_multiple_plots, z_min[i], z_max[i], $
                 linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                 linestyle_previous_model2, file_previous_model2, do_previous_model2, $
                 prefix_previous_model1, prefix_previous_model2, prefix_this_model
     
   endfor


;COLOR

   print,''
   print,''
   print,''
   print,'*************************'
   print,'*      Colour Hist      *'
   print,'*************************'
   print,''
 

   multiplot, /reset
   erase & multiplot, [4,2]
   
   xmin=1.0
   xmax=4.0  
   ymin=0.0  
   ymax=0.75
   bin=0.1

   ;load SDSS DATA
   str={ids:lonarr(3), ra:0.0, dec:0.0, z:0.0, p:0.0, rmags:fltarr(2), $
        zs:fltarr(2), mag:fltarr(5), mstar:0.0, magdump:fltarr(2), $
        Radius:fltarr(2), grcolor:fltarr(2), c:0.0, D4000:0.0, mu:0.0}
   
   N=533731L
   print,'reading SDSS color_hist data'
   openr,1,Datadir+'reference_lss_dr72bbright0.final.dat'
   obs=replicate(str,n)
   readf,1,obs
   close,1
   
   z_min=[0.01,0.01,0.01,0.03,0.03,0.03,0.16,0.24]
   z_max=[0.025,0.025,0.025,0.05,0.05,0.05,0.19,0.28]
 
   for i=0, Nbins-1 do begin 
            
      if(i gt 0) then multiplot
     
      ;AXIS LABELS
      if(i eq 0 or i eq 4) then  label_y='fraction' else label_y=''
      if(i gt 3) then label_x='u-i' else label_x=''
        
      if(lower_mass_limits[i] lt 9.5) then begin 
         G=G0_MRII 
         sam_mass=Alog10(StellarMass_MRII_0)
      endif else begin
         G=G0_MR
         sam_mass=Alog10(StellarMass_MR_0)
      endelse

      u_i=G.Magdust[15]- G.Magdust[18] 
      type=G.Type
    
      color_hist, dr7_gal_final,i, sam_mass, u_i, type, lower_mass_limits[i], upper_mass_limits[i], redshift=0.1, $
                  xmin, xmax, ymin, ymax, bin, Datadir, $
                  label_x,label_y, obs, z_min[i], z_max[i], $                  
                  write_files, plot_after_write, file_to_write, $
                  linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                  linestyle_previous_model2, file_previous_model2, do_previous_model2, $
                  prefix_previous_model1, prefix_previous_model2, prefix_this_model      
   endfor
  
 multiplot, /reset 

endif




if(UVJ_color eq 1) then begin

   print,''
   print,''
   print,''
   print,'************************'
   print,'*      UVJ colour      *'
   print,'************************'
   print,''
 

;change colors for the contour plot
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

   multiplot, /reset
   erase & multiplot, [4,1]

  
   redshift_new=[0.4,1.,2.,3.] 
   for i=0, 3 do begin

      if(redshift_new[i] eq 0.4) then G=G04_MR[where(G04_MR.MagDust[0] lt -5.0)]
      if(redshift_new[i] eq 1.0) then G=G1_MR[where(G1_MR.MagDust[0] lt -5.0)]
      if(redshift_new[i] eq 2.0) then G=G2_MR[where(G2_MR.MagDust[0] lt -5.0)]
      if(redshift_new[i] eq 3.0) then G=G3_MR[where(G3_MR.MagDust[0] lt -5.0)]
      
      if(i gt 0) then multiplot
      
      plot_UVJ_color,  G, i, hubble_h, redshift_new[i], $
                       file_to_write, write_files, plot_after_write, Datadir, $
                       linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                       linestyle_previous_model2, file_previous_model2, do_previous_model2, $
                       prefix_previous_model1, prefix_previous_model2, $
                       slope_red_fraction, offset_red_fraction, minimum_y_red_fraction
      
   endfor

   loadct,6, /SILENT
   multiplot, /reset
endif
 

;RED_FRACTION_COLORCUT
if(red_fraction_colorcut eq 1) then begin

   print,''
   print,''
   print,''
   print,'**********************'
   print,'*    Red Fraction    *'
   print,'**********************'
   print,''

   multiplot, /reset
   erase & multiplot, [5,1]


   for i_z=0,4 do begin
      if(i_z eq 0) then begin 
         mass_MR=StellarMass_MR_0
         mass_MRII=StellarMass_MRII_0
         G_MR=G0_MR
         G_MRII=G0_MRII
      endif
      if(i_z eq 1) then begin 
         mass_MR=StellarMass_MR_04
         mass_MRII=StellarMass_MRII_04
         G_MR=G04_MR
         G_MRII=G04_MRII
      endif
      if(i_z eq 2) then begin 
         mass_MR=StellarMass_MR_1
         mass_MRII=StellarMass_MRII_1
         G_MR=G1_MR
         G_MRII=G1_MRII
      endif
      if(i_z eq 3) then begin 
         mass_MR=StellarMass_MR_2
         mass_MRII=StellarMass_MRII_2
         G_MR=G2_MR
         G_MRII=G2_MRII
      endif
      if(i_z eq 4) then begin 
         mass_MR=StellarMass_MR_3
         mass_MRII=StellarMass_MRII_3
         G_MR=G3_MR
         G_MRII=G3_MRII
      endif

     
      slope=slope_red_fraction
      offset=offset_red_fraction
      minimum_y=minimum_y_red_fraction

      ;cut for z=0 on u-r vs r
      if(i_z eq 0) then  begin  
         color_MR=G_MR.MagDust[15]-G_MR.MagDust[17]
         color_MRII=G_MRII.MagDust[15]-G_MRII.MagDust[17]
         rband_MR=G_MR.MagDust[17]-5.*alog10(hubble_h)
         rband_MRII=G_MRII.MagDust[17]-5.*alog10(hubble_h)

         sel_MR = where(color_MR gt (offset[i_z]-slope[i_z]*tanh((rband_MR+18.07)/1.09)))
         sel_MRII = where(color_MRII gt (offset[i_z]-slope[i_z]*tanh((rband_MRII+18.07)/1.09)))
        
      ;cut for z>0 on UVJ
      endif else begin         
         color_UV_MR=G_MR.MagDust[0]-G_MR.MagDust[2]
         color_VJ_MR=G_MR.MagDust[2]-G_MR.MagDust[7]
         color_UV_MRII=G_MRII.MagDust[0]-G_MRII.MagDust[2]
         color_VJ_MRII=G_MRII.MagDust[2]-G_MRII.MagDust[7]
         
         sel_MR = where( (color_VJ_MR lt (minimum_y[i_z]-offset[i_z])/slope[i_z] and color_UV_MR gt minimum_y[i_z]) or $      
             (color_VJ_MR gt (minimum_y[i_z]-offset[i_z])/slope[i_z] and color_UV_MR gt (color_VJ_MR*slope[i_z] + offset[i_z])) ) 
         sel_MRII = where( (color_VJ_MRII lt (minimum_y[i_z]-offset[i_z])/slope[i_z] and color_UV_MRII gt minimum_y[i_z]) or $  
            (color_VJ_MRII gt (minimum_y[i_z]-offset[i_z])/slope[i_z] and color_UV_MRII gt (color_VJ_MRII*slope[i_z] + offset[i_z])) )
      endelse

      if(i_z gt 0) then  multiplot
      plot_red_fraction_colorcut, i_z, mass_MR, mass_MRII, G_MR, G_MRII, sel_MR, sel_MRII, model_err_red_frac, $
                         redshift=redshift[i_z], hubble_h, write_files, plot_after_write, file_to_write, $
                         linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                         linestyle_previous_model2, file_previous_model2, do_previous_model2, $
                         prefix_previous_model1, prefix_previous_model2, sample, sample_file, mcmc_folder, prefix_this_model


      ;mcmc_regions
      if(do_mcmc_regions eq 1) then begin
         char_redshift=STRTRIM(number_formatter(redshift[i_z], DECIMAL=2),2)
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
         data_bin=data_bin_low+(data_bin_high-data_bin_low)/2. 
                
         plot_mcmc_regions, mcmc_folder, observation_name, char_redshift, data_bin 
      endif

   endfor ;i_z=0,4

endif



if(stellar_mass_function_by_color eq 1) then begin

   print,''
   print,''
   print,''
   print,'***********************'
   print,'*    SMF by Colour    *'
   print,'***********************'
   print,''


plot_smf_by_color, G0_MR, G04_MR, G1_MR, G2_MR, G3_MR,G0_MRII, G04_MRII,G1_MRII, G2_MRII, G3_MRII, $
                   StellarMass_MR_0, StellarMass_MRII_0, StellarMass_MR_04, StellarMass_MRII_04, $
                   StellarMass_MR_1, StellarMass_MRII_1, StellarMass_MR_2, StellarMass_MRII_2, $
                    StellarMass_MR_3, StellarMass_MRII_3, Volume_MR, Volume_MRII, $
                    hubble_h, redshift, slope_red_fraction,offset_red_fraction, minimum_y_red_fraction, $
                    file_to_write, write_files, plot_after_write, Datadir, sample, sample_file, $
                    linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                    linestyle_previous_model2, file_previous_model2, do_previous_model2, $
                    prefix_previous_model1, prefix_previous_model2, do_mcmc_regions, prefix_this_model, $
                    plot_mcmc_regions, mcmc_folder, plot_obs ,  total_chi2, plot_only_new_obs, min_obs_err, MRII

endif


if(morphology eq 1) then begin

   print,''
   print,''
   print,''
   print,'**********************'
   print,'*     Morphology     *'
   print,'**********************'
   print,''

   multiplot, /reset

   bin=0.5
   xmin=8.0
   xmax=12.5
   ymin=0.0
   ymax=1.0
      
   plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
         xstyle = 1, ystyle = 1, xtitle = 'log!D10!N(M!D*!N[h!U-2!NM!D!9ng!3!N])', ytitle = 'fraction'
   
   symbols,2,0.75
   readcol,Datadir+'conselice2006_bulge_fract.txt',mass,frac,err, /SILENT
   oploterror,mass+alog10(0.7^2),frac,err, HATLENGTH = 100,color=70,errcolor=70,psym=8
   
   readcol,Datadir+'conselice2006_disk_fract.txt',mass,frac,err, /SILENT
   oploterror,mass+alog10(0.7^2),frac,err, HATLENGTH = 100,color=200,errcolor=200,psym=8
   
   readcol,Datadir+'conselice2006_irr_fract.txt',mass,frac,err, /SILENT
   oploterror,mass+alog10(0.7^2),frac,err, HATLENGTH = 100,color=120,errcolor=120,psym=8
   
   
   N_bins=(xmax-xmin)/bin
   mass_range = 8.+bin/2.+indgen(N_bins)*bin
   bulge=fltarr(N_bins)
   intermediate=fltarr(N_bins)
   disk=fltarr(N_bins)


   G=G0_MR
   mass=G.BulgeMass+G.DiskMass

   for ii=0, N_bins-1 do begin
      sel=where(alog10(StellarMass_MR_0) gt mass_range[ii]-bin/2. and $
                alog10(StellarMass_MR_0) lt mass_range[ii]+bin/2. and $
                G0_MR.BulgeMass/mass gt 0.7, n_bulge)

      sel=where(alog10(StellarMass_MR_0) gt mass_range[ii]-bin/2. and $
                alog10(StellarMass_MR_0) lt mass_range[ii]+bin/2. and $
                G0_MR.BulgeMass/mass lt 0.7 and G0_MR.BulgeMass/mass gt 0.01, n_intermediate)

      sel=where(alog10(StellarMass_MR_0) gt mass_range[ii]-bin/2. and $
                alog10(StellarMass_MR_0) lt mass_range[ii]+bin/2. and $
                G0_MR.BulgeMass/mass lt 0.01, n_disk)

      bulge[ii]=n_bulge*1./(n_bulge+n_intermediate+n_disk)
      intermediate[ii]=n_intermediate*1./(n_bulge+n_intermediate+n_disk)
      disk[ii]=n_disk*1./(n_bulge+n_intermediate+n_disk)   
   endfor
   oplot,mass_range,bulge,color=70,linestyle=2
   oplot,mass_range,intermediate,color=200,linestyle=2
   oplot,mass_range,disk,color=120,linestyle=2

   file=file_to_write+'_MR_morphology_bulge_z0.00.txt'
   write_to_file, file, mass_range,bulge, write_files, plot_after_write
   file=file_to_write+'_MR_morphology_intermediate_z0.00.txt'
   write_to_file, file, mass_range,intermediate, write_files, plot_after_write
   file=file_to_write+'_MR_morphology_disk_z0.00.txt'
   write_to_file, file, mass_range,disk, write_files, plot_after_write



   if(MRII eq 1) then begin

      mass=G0_MRII.BulgeMass+G0_MRII.DiskMass

      for ii=0, N_bins-1 do begin
         sel=where(alog10(StellarMass_MRII_0) gt mass_range[ii]-bin/2. and $
                   alog10(StellarMass_MRII_0) lt mass_range[ii]+bin/2. and $
                   G0_MRII.BulgeMass/mass gt 0.7, n_bulge)
         
         sel=where(alog10(StellarMass_MRII_0) gt mass_range[ii]-bin/2. and $
                   alog10(StellarMass_MRII_0) lt mass_range[ii]+bin/2. and $
                   G0_MRII.BulgeMass/mass lt 0.7 and G0_MRII.BulgeMass/Mass gt 0.01, n_intermediate)
         
         sel=where(alog10(StellarMass_MRII_0) gt mass_range[ii]-bin/2. and $
                   alog10(StellarMass_MRII_0) lt mass_range[ii]+bin/2. and $
                   G0_MRII.BulgeMass/Mass lt 0.01, n_disk)
         
         bulge[ii]=n_bulge*1./(n_bulge+n_intermediate+n_disk)
         intermediate[ii]=n_intermediate*1./(n_bulge+n_intermediate+n_disk)
         disk[ii]=n_disk*1./(n_bulge+n_intermediate+n_disk)         
      endfor

      oplot,mass_range,bulge,color=70
      oplot,mass_range,intermediate,color=200
      oplot,mass_range,disk,color=120

      file=file_to_write+'_MRII_morphology_bulge_z0.00.txt'
      write_to_file, file, mass_range,bulge, write_files, plot_after_write
      file=file_to_write+'_MRII_morphology_intermediate_z0.00.txt'
      write_to_file, file, mass_range,intermediate, write_files, plot_after_write
      file=file_to_write+'_MRII_morphology_disk_z0.00.txt'
      write_to_file, file, mass_range,disk, write_files, plot_after_write

   endif ;MRII

endif





if(bhbm eq 1) then begin

   print,''
   print,''
   print,''
   print,'****************'
   print,'*     BHBM     *'
   print,'****************'
   print,''

   loadct,0, /SILENT

   xmin=8.5
   xmax=12.5
   ymin=5.0
   ymax=10.5
   bin=0.1
 
;z=0
 
;MODEL
   bulgemass=StellarMass_MR_0
   BulgeMass=G0_MR.BulgeMass*hubble_h*1.e10
   ssfr=G0_MR.sfr/(StellarMass_MR_0)
   bh_mass=G0_MR.BlackHoleMass*1.e10
   
   levels=findgen(30)*0.1+1.5
   plot_hist_2d,Alog10(BulgeMass),Alog10(bh_mass),xmin = xmin-1.0,xmax=xmax+1.0, ymin = ymin-1.0, ymax=ymax+1.0, $
                xstyle = 1, ystyle = 1, xtitle ='log!D10!N(M!DBulge!N[h!U-2!NM!D!9ng!3!N])', $
                ytitle ='log!D10!N(M!DBH!N[h!U-1!NM!D!9ng!3!N])',$
                xbin=0.1,ybin=0.05,/cont,point=10.^1.5,/log,$
                levels=levels,xrange=[xmin,xmax],yrange=[ymin,ymax],xs=1,ys=1,/fill
  
  
   loadct,6, /SILENT


;observations
   readcol,Datadir+'mcconnel2012_BHBM.dat', con_gal_name, con_distance, con_MBH, con_MBH_lower, con_MBH_upper, $
           con_method, con_sigma, con_sigma_lower, con_sigma_upper, con_log_LV, con_error_log_LV, $
           con_log_L_36, error_log_L_36, con_Mbulge, con_radius_of_influence, con_Morphology, con_Profile, $
           con_Reff_Vband, con_Reff_iband, con_Reff_36um, FORMAT = 'A,f,f,f,f,A,f,f,f,f,f,f,f,f,f,A,A,f,f,f', /SILENT 

   symbols,2,1.0
   oploterror, Alog10(con_Mbulge*hubble_h_WMAP7^2), Alog10(con_MBH*hubble_h_WMAP7), $
               fltarr(n_elements(con_Mbulge))+0.24, Alog10(con_MBH_upper/con_MBH), $
               /hiba, HATLENGTH=80., psym=8, color=120, errcolor=120
   oploterror, Alog10(con_Mbulge*hubble_h_WMAP7^2), Alog10(con_MBH*hubble_h_WMAP7), $
               fltarr(n_elements(con_Mbulge))+0.24, Alog10(con_MBH/con_MBH_lower), $
               /loba, HATLENGTH=80., psym=8, color=120, errcolor=120

   plot_label, xlog=0, ylog=0, type='label', label='McConnell & Ma (2013)' , xmin, xmax, ymin, ymax, $
               x_percentage=0.7, x2_percentage=0., y_percentage=9.075, $
               charthick=3., charsize=1.0   
   plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
               x_percentage=0.5, x2_percentage=0., y_percentage=9.2, $
               color=120, errcolor=120, sym_num=2, sym_size=1.0, HATLENGTH = 100, err_size=0.1
  
   ;COLORBAR
   loadct,0, /SILENT
   COLORBAR, BOTTOM=0, VERTICAL=0, POSITION=[0.67, 0.24, 0.9, 0.28], $
             MAX=0.,Min=alog10(10^min(levels)/10^max(levels)),DIVISIONS=3, FORMAT='(F6.1)'
   loadct,6, /SILENT

endif




 if(gas_mass_function eq 1) then begin 

    print,''
    print,''
    print,''
    print,'*********************'
    print,'*        GMF        *'
    print,'*********************'
    print,''

    plot_gas_mass_function, G0_MR, Volume_MR, G0_MRII, Volume_MRII, Hubble_h, hubble_h_WMAP7, Datadir, sample, sample_file, $   
                            write_files, plot_after_write, file_to_write, $
                            linestyle_previous_model1, file_previous_model1, do_previous_model1, $
                            linestyle_previous_model2, file_previous_model2, do_previous_model2, $
                            prefix_previous_model1, prefix_previous_model2, prefix_this_model, MRII
    
 endif


if(mstar_metals eq 1) then begin 

   print,''
   print,''
   print,''
   print,'**********************'
   print,'*   Metals vs Mass   *'
   print,'**********************'
   print,''

   G=G0_MR
   plot_mstar_metals, G, redshift=0.1, hubble_h, $
                      write_files, plot_after_write, file_to_write, $
                      linestyle_previous_model1, file_previous_model1, do_previous_model1, prefix_previous_model1, $
                      linestyle_previous_model2, file_previous_model2, do_previous_model2, prefix_previous_model2, prefix_this_model

   plot_model_metals,  G3_MR, redshift=redshift[4] , hubble_h, $
                        write_files, plot_after_write, file_to_write, plot_color=200
endif







 print,''
  print,''
  print,''
  print,'**************************'
  print,'*   ISM Metals vs Mass   *'
  print,'**************************'
  print,''

  DTDchoice = 'BiModal' ;'PowerLaw' ;'Gaussian'

  theta = findgen(30)/29.*360.
  usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill
  loadct, 39, /silent

  ;Atomic Weights:
  H_aw = 1.008
  He_aw = 4.003
  C_aw = 12.01
  N_aw = 14.01
  O_aw = 16.0
  Ne_aw = 20.18
  Mg_aw = 24.31
  Si_aw = 28.09
  S_aw = 32.07
  Ca_aw = 40.08
  Fe_aw = 55.84

  GMZR_M_ALL_noCuts = alog10(StellarMass_MR_0) ;alog10(StellarMass_MR)
  GMZR_SFR_ALL_noCuts = alog10(G0_MR.Sfr*hubble_h^2)
  GMZR_Z_ALL_noCuts =
alog10((G0_MR.ColdGas_elements[4]/G0_MR.ColdGas_elements[0])*(H_aw/O_aw))+12.0

  ww = where(GMZR_SFR_ALL_noCuts ge -2.0 AND GMZR_SFR_ALL_noCuts le 1.6
AND GMZR_M_ALL_noCuts ge 8.5 AND GMZR_M_ALL_noCuts le 11.5)
  GMZR_M_ALL = GMZR_M_ALL_noCuts(ww)
  GMZR_SFR_ALL = GMZR_SFR_ALL_noCuts(ww)
  GMZR_Z_ALL = GMZR_Z_ALL_noCuts(ww)

  GMZR_M_ALL = GMZR_M_ALL(WHERE((FINITE(GMZR_Z_ALL, /NAN, SIGN=0) ne 1)))
  GMZR_SFR_ALL = GMZR_SFR_ALL(WHERE((FINITE(GMZR_Z_ALL, /NAN, SIGN=0) ne 1)))
  GMZR_Z_ALL = GMZR_Z_ALL(WHERE((FINITE(GMZR_Z_ALL, /NAN, SIGN=0) ne 1)))

  print, "Total number of z=0 star-forming galaxies (for MZRs): ",
n_elements(GMZR_Z_ALL)

  xlimits = [9.0, 11.5]
  ylimits = [7.9, 9.7]
  plot, findgen(10), /nodata, xtitle = 'log(M!D*!N)
[M!D!9'+string(110B)+'!3!N]', ytitle = '12 + log(O/H)!Dcold!N', $
    xrange = xlimits, yrange = ylimits, xstyle = 1, ystyle = 1,
CHARSIZE=1.5, CHARTHICK=4., XTHICK=5., YTHICK=5.

  countNum = 8 ;20

  ;Bin galaxies by log(M*) and Z:
  Mmin = xlimits[0]
  Mmax = xlimits[1]
  Mbinno = 64.
  Mbinwidth = (Mmax-Mmin) / Mbinno
  Zmin = ylimits[0]
  Zmax = ylimits[1]
  Zbinno = 64.
  Zbinwidth = (Zmax-Zmin) / Zbinno

  AveM = fltarr(Mbinno*Zbinno)
  AveZ = fltarr(Mbinno*Zbinno)
  AveSFR = fltarr(Mbinno*Zbinno)
  count = intarr(Mbinno*Zbinno)
  deleg = intarr(n_elements(GMZR_M_ALL))
  k=0
  FOR i=Mmin, Mmax-Mbinwidth, Mbinwidth DO BEGIN
    FOR j=Zmin, Zmax-Zbinwidth, Zbinwidth DO BEGIN
      w99 = WHERE(GMZR_M_ALL ge i and GMZR_M_ALL lt i+Mbinwidth and
GMZR_Z_ALL ge j and GMZR_Z_ALL lt j+Zbinwidth)
      IF(w99[0] ne -1) THEN BEGIN
        IF(n_elements(w99) ge countNum) THEN deleg(w99) = 1
        AveM[k] = i + (0.5*Mbinwidth)
        AveZ[k] = j + (0.5*Zbinwidth)
        AveSFR[k] = AVG(GMZR_SFR_ALL(w99(WHERE(GMZR_SFR_ALL(w99) gt
-100.0))))
        count[k] = n_elements(w99)
      ENDIF
      k+=1
    ENDFOR
  ENDFOR
  ww2 = WHERE(count ge countNum)

  ;Plot low-res points:
  usersym, [-1.1,1.1,1.1,-1.1], [1.1,1.1,-1.1,-1.1], /fill
  ;ctload, 0, /REVERSE, CLIP=[50,200]
  ctload, 0, CLIP=[50,200]
  plots,AveM(ww2), AveZ(ww2), psym=8, color=bytscl(count(ww2))

  ;Plot other points:
  loadct, 13, /silent
  usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill
  ww99 = WHERE(deleg eq 0)
  oplot, GMZR_M_ALL(ww99), GMZR_Z_ALL(ww99), psym = 8, symsize = ssiz,
color=fsc_color("grey")

  ;Plot mean, and 1-sigma spread lines:
  projmin=MIN(GMZR_M_ALL)
  projmax=MAX(GMZR_M_ALL)
  dim=number_formatter((projmax-projmin)/0.1, DECIMALS=0)
  xx = findgen(dim)/10. + projmin
  lower1sigmag2 = fltarr(dim)
  upper1sigmag2 = fltarr(dim)
  meanZg2 = fltarr(dim)
  j=0
  FOR i=projmin,projmax,0.1 DO BEGIN
    w1s = WHERE(GMZR_M_ALL ge i-0.05 AND GMZR_M_ALL lt i+0.05)
    IF(w1s[0] ne -1) THEN BEGIN
      upper1sigmag2[j] =
GMZR_Z_ALL(w1s((SORT(GMZR_Z_ALL(w1s)))(84*N_ELEMENTS(GMZR_Z_ALL(w1s))/100)))
      meanZg2[j] = AVG(GMZR_Z_ALL(w1s))
    ENDIF
    j=j+1
  ENDFOR
  www = WHERE(xx ge 8.5 and xx le 11.3)
  oplot, xx(www), upper1sigmag2(www), thick=6., linestyle = 2,
color=fsc_color("black")
  oplot, xx(www), meanZg2(www), thick=6., linestyle = 0,
color=fsc_color("black")

  ;Plot Qi Guo's old gas MZR relation:
  xxOLD =
[8.50016,8.60016,8.70016,8.80016,8.90016,9.00016,9.10016,9.20016,9.30016,9.40016,9.50016,9.60016,9.70016,9.80016,9.90016,10.0002,10.1002,10.2002,10.3002,10.4002,10.5002,10.6002,10.7002,10.8002,10.9002,11.0002,11.1002]
  oplot, xxOLD,
[8.29328,8.31744,8.36549,8.38975,8.45409,8.48258,8.52472,8.57333,8.61011,8.64945,8.69886,8.73729,8.78996,8.82114,8.87590,8.91173,8.94738,8.98213,9.01795,9.05720,9.07654,9.10225,9.14133,9.16229,9.16034,9.16464,9.15175]
$
    , thick=6., linestyle = 0, color=fsc_color("red")
  oplot, xxOLD,
[8.15828,8.17994,8.23356,8.25732,8.31912,8.34660,8.38584,8.45410,8.47954,8.52805,8.58798,8.62263,8.67632,8.70672,8.76378,8.79798,8.82560,8.87950,8.90626,8.94199,8.97523,8.99093,9.02006,9.03025,9.04985,9.02321,9.01312]
$
    , thick=6., linestyle = 2, color=fsc_color("red")
  oplot, xxOLD,
[8.42827,8.45493,8.49743,8.52218,8.58906,8.61855,8.66362,8.69257,8.74068,8.77085,8.80974,8.85195,8.90359,8.93555,8.98802,9.02548,9.06917,9.08476,9.12964,9.17242,9.17786,9.21357,9.26259,9.29433,9.27082,9.30607,9.29037]
$
    , thick=6., linestyle = 2, color=fsc_color("red")

  ;Yates et al (2012a):
  xxx = findgen(100)/(100.0/(11.2374-8.73742)) + 8.73742
  yyy = 26.6864 - 6.6399*xxx + 0.768653*xxx*xxx - 0.0282147*xxx*xxx*xxx
  oplot, xxx, yyy, thick = 6, linestyle = 0, color = fsc_color("orange")

  xxxSigma = [8.73742, 8.83742, 8.93742, 9.03742, 9.13742, 9.23742,
9.33742, 9.43742, 9.53742, 9.63742, 9.73742, 9.83742, 9.93742, 10.0374,
10.1374, 10.2374, 10.3374, 10.4374, 10.5374, 10.6374, 10.7374, 10.8374,
10.9374, 11.0374, 11.1374, 11.2374]
  UpperSigma = [8.76605, 8.75496, 8.76407, 8.78068, 8.79891, 8.81812,
8.84087, 8.87536, 8.90295, 8.93559, 8.96448, 8.99686, 9.02141, 9.04769,
9.06818, 9.08794, 9.10189, 9.11680, 9.13128, 9.14483, 9.15107, 9.16058,
9.17166, 9.17798, 9.18027, 9.19381]
  LowerSigma = [8.41567, 8.46743, 8.52952, 8.55702, 8.57257, 8.59747,
8.62097, 8.64273, 8.67724, 8.70625, 8.74066, 8.77924, 8.81401, 8.85304,
8.88790, 8.91583, 8.93949, 8.96325, 8.97943, 8.98876, 9.00417, 9.01550,
9.02068, 9.02914, 9.02327, 9.02569]
  oplot, xxxSigma, UpperSigma, thick = 6, linestyle = 2, color =
fsc_color("orange")
  oplot, xxxSigma, LowerSigma, thick = 6, linestyle = 2, color =
fsc_color("orange")

  legend, ['L-Galaxies (w/ new GCE)', 'L-Galaxies (w/o new GCE)',
'SDSS-DR7'], box = 0, linestyle = [-99,-99,-99], textcolor =
[fsc_color("black"),fsc_color("red"),fsc_color("orange")], thick =
[6,6,6], charsize=1.5, charthick=4., position = [11.5,8.1+0.1], /right
;/bottom ;, linsize = 0.75
  IF (DTDchoice eq 'BiModal') THEN  legend, ['BM'], box = 0, linestyle =
[-99], textcolor = [fsc_color("black")], charsize=1.5, charthick=4.,
/top, /left
  IF (DTDchoice eq 'PowerLaw') THEN  legend, ['PL'], box = 0, linestyle =
[-99], textcolor = [fsc_color("black")], charsize=1.5, charthick=4.,
/top, /left
  IF (DTDchoice eq 'Gaussian') THEN  legend, ['NG'], box = 0, linestyle =
[-99], textcolor = [fsc_color("black")], charsize=1.5, charthick=4.,
/top, /left


  axis, xaxis=0, xrange = xlimits, xstyle = 1, XTHICK=5., CHARSIZE=1.5
  axis, yaxis=0, yrange = ylimits, ystyle = 1, YTHICK=5., CHARSIZE=1.5
  axis, 0, ylimits[1], xaxis=1, xrange = xlimits, xstyle = 1, XTHICK=5.,
CHARSIZE=1.5, XTICKNAME=[' ',' ',' ',' ',' ',' ',' ',' ',' ']
  axis, xlimits[1], yaxis=1, yrange = ylimits, ystyle = 1, YTHICK=5.,
CHARSIZE=1.5, YTICKNAME=[' ',' ',' ',' ',' ',' ',' ',' ']

;------------------

  print,''
  print,''
  print,''
  print,'******************************'
  print,'*   Stellar [a/Fe] vs Mass   *'
  print,'******************************'
  print,''

  Wind = 0 ;Turn on if you ran with METALRICHWIND on.

  ;Solar mass ratios (from Grevesse, Noels & Sauval 1996):
  ;H_mf_GNS = ??
  HeH_mf_GNS = (4.003/1.008)*10.0^(10.99-12.)
  CH_mf_GNS = (12.01/1.008)*10.0^(8.55-12.)
  NH_mf_GNS = (14.01/1.008)*10.0^(7.97-12.)
  OH_mf_GNS = (16.00/1.008)*10.0^(8.87-12.)
  NeH_mf_GNS = (20.18/1.008)*10.0^(8.08-12.)
  MgH_mf_GNS = (24.31/1.008)*10.0^(7.58-12.)
  SiH_mf_GNS = (28.09/1.008)*10.0^(7.55-12.)
  SH_mf_GNS = (32.07/1.008)*10.0^(7.33-12.)
  CaH_mf_GNS = (40.08/1.008)*10.0^(6.36-12.)
  FeH_mf_GNS = (55.84/1.008)*10.0^(7.5-12.)
  aFe_mf_GNS =
AVG([alog10(MgH_mf_GNS/FeH_mf_GNS),alog10(SiH_mf_GNS/FeH_mf_GNS),alog10(SH_mf_GNS/FeH_mf_GNS),alog10(CaH_mf_GNS/FeH_mf_GNS)])
  ZH_GNS = 0.0169 ;The calculated value of Z/(H+He) from Grevesse & Sauval
(1998) (see Asplund et al. 2009). NB: Grevesse & Noels (1993) found
Z/(H+He)=0.0179


  M_ALL = alog10(StellarMass_MR_0) ;alog10(StellarMass_MR)
  OFe_ALL =
alog10((G0_MR.DiskMass_elements[4]+G0_MR.BulgeMass_elements[4])/(G0_MR.DiskMass_elements[10]+G0_MR.BulgeMass_elements[10]))
- alog10(OH_mf_GNS/FeH_mf_GNS)

  ;Basic elliptical sample:
  w1 = WHERE((G0_MR.BulgeMass/(G0_MR.BulgeMass + G0_MR.DiskMass) ge 0.7) $
    AND (G0_MR.magdust[1]-G0_MR.magdust[2] ge
(0.051*alog10(G0_MR.BulgeMass + G0_MR.DiskMass)) + 0.14) $
;Genevieve's g-r cut
    AND (FINITE(G0_MR.DiskMass_elements[*], /NAN, SIGN=0) ne 1) $
    AND (FINITE(G0_MR.BulgeMass_elements[*], /NAN, SIGN=0) ne 1))

  M_ALL = M_ALL(w1)
  OFe_ALL = OFe_ALL(w1)

xlimits = [9.5, 11.8] ;[9.5, 12.0]
ylimits = [0.05, 0.6] ;[0.05, 0.45]
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N)
[M!D!9'+string(110B)+'!3!N]', ytitle = '[O/Fe]!Dbulge+disc!N', $
  xrange = xlimits, yrange = ylimits, xstyle = 1, ystyle = 1,
CHARSIZE=1.5, CHARTHICK=4., XTHICK=5., YTHICK=5.,
POSITION=[0.15,0.1,0.95,0.6]

;;Do shaded contours for mass-age-selected points (95th percentile):
binno = 40
xbinwidth = (xlimits[1]-xlimits[0])/(binno-1)
ybinwidth = (ylimits[1]-ylimits[0])/(binno-1)

;Do contours for all points:
;hist2d, M_ALL, OFe_ALL, hist, xlimits, ylimits, binno, binno
hist=hist_2d(M_ALL,OFe_ALL,min1=xlimits[0],min2=ylimits[0],max1=xlimits[1],max2=ylimits[1],bin1=xbinwidth,bin2=ybinwidth)
;*weight
contour, hist, findgen(binno)/(binno/(xlimits[1]-xlimits[0])) +
xlimits[0], findgen(binno)/(binno/(ylimits[1]-ylimits[0])) + ylimits[0],
levels=[((100.-95.)/100.)*MAX(hist),((100.-68.)/100.)*MAX(hist)], $
;nlevels = 3, $
  xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle =
'[O/Fe]!Ddisc+bulge!N', CHARSIZE=1.5, CHARTHICK=4., XTHICK=5.,
YTHICK=5., $
  thick = 8., c_color = fsc_color("dodger blue"), c_linestyle = 0,
/CLOSED, /OVERPLOT

usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill

;Plot Jonas' fit using MPA M* directly:
slopeJ3 = 0.051308367
interceptJ3 = -0.32725757
OnesigmaJ3 = 0.092184387
logMstar3 = findgen(100+1)/(float(100)/(11.5-9.5)) + 9.5
;robsarray(100,9.5,11.5)
yJ3 = (slopeJ3*logMstar3) + interceptJ3
oplot, logMstar3, yJ3, linestyle = 0, thick = 15., color = fsc_color("black")
oplot, logMstar3, yJ3, linestyle = 0, thick = 10., color =
fsc_color("orange")
oplot, logMstar3, yJ3-OnesigmaJ3, linestyle = 1, thick = 15., color =
fsc_color("black")
oplot, logMstar3, yJ3-OnesigmaJ3, linestyle = 1, thick = 10., color =
fsc_color("orange")
oplot, logMstar3, yJ3+OnesigmaJ3, linestyle = 1, thick = 15., color =
fsc_color("black")
oplot, logMstar3, yJ3+OnesigmaJ3, linestyle = 1, thick = 10., color =
fsc_color("orange")

IF (DTDchoice eq 'Gaussian') THEN BEGIN
  polyfill, [9.6,9.82,9.82,9.6,9.6], [0.35,0.35,0.42,0.42,0.35], /fill,
color=fsc_color("white")
  IF (Wind eq 0) THEN polyfill, [11.67,11.85,11.85,11.67,11.67],
[0.405,0.405,0.45,0.45,0.405], /fill, color=fsc_color("white")
  IF (Wind eq 1) THEN polyfill, [11.075,11.375,11.375,11.075,11.075],
[0.39,0.39,0.46,0.46,0.39], /fill, color=fsc_color("white")
ENDIF

IF (DTDchoice eq 'BiModal' AND Wind ne 1) THEN  legend, ['BM'], box = 0,
linestyle = [-99], textcolor = [fsc_color("black")], charsize=1.5,
charthick=4., /top, /right
IF (DTDchoice eq 'PowerLaw' AND Wind ne 1) THEN  legend, ['PL'], box = 0,
linestyle = [-99], textcolor = [fsc_color("black")], charsize=1.5,
charthick=4., /top, /right
IF (DTDchoice eq 'Gaussian' AND Wind ne 1) THEN  legend, ['NG'], box = 0,
linestyle = [-99], textcolor = [fsc_color("black")], charsize=1.5,
charthick=4., /top, /right
IF (DTDchoice eq 'BiModal' AND Wind eq 1) THEN  legend, ['BM f!Dwind!N'],
box = 0, linestyle = [-99], textcolor = [fsc_color("black")],
charsize=1.5, charthick=4., /top, /right
IF (DTDchoice eq 'PowerLaw' AND Wind eq 1) THEN  legend, ['PL f!Dwind!N'],
box = 0, linestyle = [-99], textcolor = [fsc_color("black")],
charsize=1.5, charthick=4., /top, /right
IF (DTDchoice eq 'Gaussian' AND Wind eq 1) THEN  legend, ['NG f!Dwind!N'],
box = 0, linestyle = [-99], textcolor = [fsc_color("black")],
charsize=1.5, charthick=4., /top, /right

axis, xaxis=0, xrange = xlimits, xstyle = 1, XTHICK=5., CHARSIZE=1.5
axis, yaxis=0, yrange = ylimits, ystyle = 1, YTHICK=5., CHARSIZE=1.5
axis, 0, ylimits[1], xaxis=1, xrange = xlimits, xstyle = 1, XTHICK=5.,
CHARSIZE=1.5, XTICKNAME=[' ',' ',' ',' ',' ',' ',' ',' ',' ']
axis, xlimits[1], yaxis=1, yrange = ylimits, ystyle = 1, YTHICK=5.,
CHARSIZE=1.5, YTICKNAME=[' ',' ',' ',' ',' ',' ',' ',' ']






   device,/close_file

   stop

end
