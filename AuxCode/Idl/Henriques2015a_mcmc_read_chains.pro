@procedures

;.r mcmc_read_chains.pro

;scp dph6bmh@charon.dur.ac.uk:/data/rw16/dph6bmh/allz/current/132/output/senna_g2_"*" ./mcmc_chains/


pro read_chains

  !Path = '/galformod/scratch/bmh20/plots/' + ':' + !Path  
  !Path = Expand_Path( '+' + !PATH,  /All_Dirs ) 
 
  ;mcmc_dir = '/galformod/scratch/bmh20/plots/mcmc_chains/new_postprocess_code/'
  mcmc_dir = '/galformod/scratch/bmh20/plots/mcmc_chains/'
  ;mcmc_dir = '/galformod/scratch/bmh20/SAM/Henriques2014/mcmc_chains/'

;in case multiple runs are to be combined in 1D plots
  N_MCMC=1
  chains_file=StrArr(N_MCMC)
  ;chains_file[0]='senna_gh'
  chains_file[0]='senna_g6'
  ;chains_file[0]='senna_g0'
  ;chains_file[0]='mysenna6'

  oneD_plots=1
  parameter_scalings=0
  print_stats=1
  test_convergence=1

  fix_max_like=0
  max_log_like=24.834517

  plots_in_log=0
  do_2sigma_regions=1
  plot_bestfit_value=1

  do_plot=[1,1,1,1, 1,1,1, $
           1,1,1, 1,1,1, $
           1, 1, 1,1,1, 1,1]

  Omega_M=0.315 
  Omega_Lambda=0.683

  plot_previous_model_1=1
  par_previous_model_1=[0.01,0.38,0.56,0.70,7.0e-4,0.03,280,4.0,80.0,3.2, $
                        0.2,90.0,3.2,0.3,0.03,2.0,0.]
  linestyle_previous_model_1=1
  plot_previous_model_2=1 
  ;par_previous_model_2=[0.0299, 0.12, 0.687, 1.86, 0.01, 0.0383, 755, 2.46, 440, 0.632, $
  ;                      0.639, 115, 0.74, 4.11e+10, 0.0411, 0.456, 2.46, 1.23e+04]
  par_previous_model_2=[0.055,0.38,0.56,0.70,3.2e-3,0.015,280,2.1,405,0.92, $
                         0.65,336,0.46,1.8e10,0.047,2.0,0.]
  linestyle_previous_model_2=2

;radio mode feedback parameter has new normalization
  par_previous_model_1[4]*=(0.7*0.7)
  par_previous_model_2[4]*=(0.7*0.7)

  ;par_previous_model_2=[0.0299, 0.125, 0.687, 1.86, 0.0113, 0.0383, 755, 2.46, 440, 0.632, $
  ;                      0.639, 115, 0.74, 4.11e+10, 0.0311, 0.456, 2.46, 1.23e+04] 

  ;linestyle_previous_model_2=2

;x=[0.06,0.09,0.649,0.7,7.0e-03,0.0264,432,1.53,365,1.02,0.83,208,0.44,3.79e+10,0.042,0.41,3.0,1.97e+04,6.5,6.0]
;0.06 0.09 0.649 1.90  0.012  7.0e-03 432 1.53 365 1.02  0.83 208  0.44  3.79e+10  0.042 0.41 3.0  1.97e+04  6.5 6.0          

  n_chains=120.
  ;chain_length=45000.  ;wc -l ./mcmc_chains/senna_g2_0.txt
  chain_length=900.
  N_par=17
  Nbins=40
  BurnIn=0
  nrows=n_chains*(chain_length-BurnIn)


  if(plots_in_log eq 1) then begin
     par_previous_model_1=Alog10(par_previous_model_1)
     par_previous_model_2=Alog10(par_previous_model_2)
  endif

 
  !p.charsize = 1.4
  !p.charthick = 4
  !p.thick = 4
  !x.thick = 4
  !y.thick = 4
  set_plot, 'PS'
  
  ;device, filename = 'read_chains.ps', xsize = 35, ysize = 22, /color, xoffset=1, yoffset=5
  ;device, filename = './fig/mcmc_1d.ps', xsize = 35, ysize = 8, /color, xoffset=1, yoffset=5


;Parameter Limits
   par_min=FltArr(N_par)
   par_max=FltArr(N_par)
   get_par_limits, N_par, par_min, par_max, plots_in_log


;Parameter Names
   par_name=StrArr(N_par)
   get_par_names, N_par, par_name, plots_in_log


 ;read chains
   data=DblArr(N_MCMC, N_par+2,nrows) 
   for j=0, N_MCMC-1 do begin   
      print,''
      print,'reading chain '+chains_file[j]  
      for i=0, n_chains-1 do begin
         print,''
         print,'reading chain chain number:',i  
         close,1
         file=mcmc_dir+chains_file[j]+'_'+STRTRIM(number_formatter(i, DECIMAL=0),2)+'.txt'
         openr,1,file 
         aux=fltarr(N_par+2,chain_length)
         readf,1,aux    
         data(j,*,i*(chain_length-BurnIn):(i+1)*(chain_length-BurnIn)-1)=aux(*,BurnIn:chain_length-1)
         close,1         
      endfor
      ;coldgas crit parameter has new normalization
      ;;data(j,3,*)=Alog10(10^data(j,3,*)*2.)
      ;radio mode feedback parameter has new normalization
      ;1.25 meaningless, just to get the right value
      ;;data(j,6,*)=Alog10(10^data(j,6,*)*0.673*0.673*1.25)
      ;adujst yield
      ;;data(j,16,*)=Alog10(10^data(j,16,*)*1.28)

      ;print stats
      print,''
      print,FORMAT='(%"Max Log10(Like)=%0.8g")',min(data(j,1,*))     
      ;select all the best fit steps 
      aux=data(j,*,[where(data(j,1,*) eq min(data(j,1,*)))])   
      ;print parameters from one of them 
      if(j eq 0) then $        
         for ii=0, N_par-1 do print,FORMAT='(%"%0.2g ", $)',10^aux(0,ii+2,0)  
      print,''
      print,''
        
      

   endfor
                    

   print,'marginalized stats'       
   for ii=0, N_par-1 do begin
      xmin=par_min[ii]
      xmax=par_max[ii]
      
      bin=(Alog10(xmax)-Alog10(xmin))/float(Nbins)
      hist=histogram(data[0,ii+2,*],locations=c,min=alog10(xmin)-bin*2.,max=alog10(xmax)+bin*2.,binsize=bin)
      c=10^(c+bin/2.) 
   
      sel=where(hist eq max(hist))
      
      print,FORMAT='(%"%0.2g ", $)',c[sel] 

   endfor
   print,''
   print,''

   if(plots_in_log eq 0) then begin 
      data=10^data
      ;do not correct chain weight and likelihood
      data(*,0:1,*)=Alog10(data(*,0:1,*))     
   endif

   ;data(*,15,*)=10^data(*,15,*)
   ;data(*,18,*)=10^data(*,18,*)


;test_convergence=
   loadct,6

 if(test_convergence eq 1) then begin
    pro_test_convergence, N_par, data, par_min, par_max, par_name, plots_in_log, $
       chain_length, BurnIn, n_chains
 endif

;1D plots 
   if(oneD_plots eq 1) then begin
      device, filename = './fig/mcmc_1d.ps', xsize = 35, ysize = 8, /color, xoffset=1, yoffset=5
      plot_1D, N_MCMC, N_par, data, par_min, par_max, par_name, Nbins, do_plot, $
               plots_in_log, nrows, do_2sigma_regions, plot_bestfit_value, $
               plot_previous_model_1, par_previous_model_1, linestyle_previous_model_1, $
               plot_previous_model_2, par_previous_model_2, linestyle_previous_model_2, $
               fix_max_like, max_log_like     
      device,/close_file
   endif


;statistics for MCMC nr 0
   if(print_stats eq 1) then begin
      for ii = 0,N_par-1 do begin      
         ;compute regions
         upper=ConfidVal(ii+2,(1-0.955)/2.,1, nrows, data)   
         lower=ConfidVal(ii+2,(1-0.955)/2.,0, nrows, data) 
       
         ;best fit
         if(fix_max_like eq 0) then $
            aux=data(0,ii+2,[where(data(0,1,*) eq min(data(0,1,*)))]) $
         else $
            aux=data(0,ii+2,[where(data(0,1,*) eq max_log_like)]) 
         if(ii eq 0) then print,'      par   ','     best fit','         lower','          upper'    
         if(plots_in_log eq 1) then $     
            print, ii, 10^aux(0,0,0), 10^lower, 10^upper, aux(0,0,0), lower, upper $
         else $
            print, ii, aux(0,0,0), lower, upper, Alog10(aux(0,0,0)), Alog10(lower), Alog10(upper)
      endfor     
   endif






;are these the 2D plots?
   for ix = 0,N_par-1 do begin 
      ;print,''
      ;print,'parameter',ix  
      bincounts = 0
      binlikes = 0
      binmaxlikes = 0
      binsraw = 0
      
      force_twotail=1

      num_contours=2
      contours=fltarr(num_contours)
      contours[0] = 0.683
      contours[1] = 0.955
      
      cont_lines=fltarr(N_par,2,num_contours)

      for ix1 = 0, num_contours-1 do begin
         limfrac = (1-contours(ix1))
         if (force_twotail eq 1) then $
            limfrac = limfrac/2 ;two tail 
         cont_lines(ix,1,ix1) = ConfidVal(ix+2,limfrac,1, nrows, data)        
         cont_lines(ix,0,ix1) = ConfidVal(ix+2,limfrac,0, nrows, data)        
      endfor                    ;contour lines
 
      if(plots_in_log eq 1)then begin
         ;print,'upper 2sigma=',10^cont_lines(ix,1,1)
         ;print,'upper 1sigma=',10^cont_lines(ix,1,0)      
         ;print,'lower 1sigma=',10^cont_lines(ix,0,0)
         ;print,'lower 2sigma=',10^cont_lines(ix,0,1)
      endif else begin
         ;print,'upper 2sigma=',cont_lines(ix,1,1)
         ;print,'upper 1sigma=',cont_lines(ix,1,0)      
         ;print,'lower 1sigma=',cont_lines(ix,0,0)
         ;print,'lower 2sigma=',cont_lines(ix,0,1)
      endelse

   endfor







;do plots with parameter scalings
   if(parameter_scalings eq 1) then begin
   
      ;parameters_model1=[0.011, 0.19, 0.56, 0.7, 0.0007, 0.03, 280., 4.0,80.,3.2,0.18,90.,3.2,0.3]
      H_z_1=[100.0, 170.41, 284.11, 425.86]
      ;parameters_model2=[0.055, 0.19, 0.56, 0.7, 0.0032, 0.015, 280., 2.12,405.,0.921,0.645,336.,0.46,1.81e10]
      H_z_2=[100.0, 170.41, 284.11, 425.86]

  linestyle_previous_model_1=1
  plot_previous_model_2=1 
  ;par_previous_model_2=[0.0299, 0.12, 0.687, 1.86, 0.01, 0.0383, 755, 2.46, 440, 0.632, $
  ;                      0.639, 115, 0.74, 4.11e+10, 0.0411, 0.456, 2.46, 1.23e+04]
 

      parameters_model1=par_previous_model_1
      parameters_model2=par_previous_model_2

      ;parameters_bestfit=[2.46,436.,0.63,0.64,115.,0.74,4.0e10]
      ;parameters_bestfit=[2.46,436.,0.63,10.64,115.,0.74,4.0e10]
      ;bestfit_min=[1.44,237.,0.39,0.45,91,0.40,2.6e10]
      ;bestfit_max=[3.8,727.,1.06,1.28,316.,1.3,6.5e10]

      parameters_bestfit=fltarr(14)
      bestfit_min=fltarr(14)
      bestfit_max=fltarr(14)
      ;H_z_bestfit=[100.0, 179.19, 303.41, 456.89]
      
      for i=0,13 do begin

         ;best fit
         if(fix_max_like eq 0) then $
            aux=data(0,i+2,[where(data(0,1,*) eq min(data(0,1,*)))]) $
         else $
            aux=data(0,i+2,[where(data(0,1,*) eq max_log_like)])
            
         ;2sigma regions
         upper=ConfidVal(i+2,(1-0.955)/2.,1, nrows, data)   
         lower=ConfidVal(i+2,(1-0.955)/2.,0, nrows, data)

         if(plots_in_log eq 1) then begin
            aux=10^aux
            lower=10^lower
            upper=10^upper
         endif 

         parameters_bestfit[i]=aux(0,0,0)
         
         bestfit_min[i]=lower
         bestfit_max[i]=upper
      endfor

      plot_feedback_scalings, N_two_sigma_regions=1000, nvmax=1000., parameters_model1, H_z_1, do_previous_model1,$
                              parameters_model2, H_z_2, do_previous_model2, $
                              parameters_bestfit,bestfit_min,bestfit_max,Omega_M, Omega_Lambda
    
      
   endif





stop


end



pro pro_test_convergence, N_par, data, par_min, par_max, par_name, plots_in_log, $
                          chain_length, BurnIn, n_chains

  set_plot, 'PS'
  device, filename = './fig/mcmc_test_convergence.ps', xsize = 25, ysize = 20, /color, xoffset=1, yoffset=5
 
;only the first "MCMC" in case of multiple runs/models being compared
  i=0

  Par_j=0
 
  xmin=0.
  xmax=chain_length-BurnIn
  ymin=0.01
  ymax=0.1
   
  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, $
        xtitle = 'step', ytitle = par_name[Par_j]
  
  ;for Chain_k=0, n_chains-1 do begin
  for Chain_k=0, 99 do begin
     data_Par_j_Chain_k=fltarr(chain_length-BurnIn)
     x_axis=findgen(n_elements(data_Par_j_Chain_k))
     
     data_Par_j_Chain_k[*]=data[i,Par_j+2,Chain_k*(chain_length-BurnIn):(Chain_k+1)*(chain_length-BurnIn)-1]
     
     oplot,x_axis,data_Par_j_Chain_k,linestyle=0, color=0+Chain_k*2.

  endfor

  device,/close_file
end

pro plot_1D, N_MCMC, N_Par, data, par_min, par_max, par_name, Nbins, do_plot, $
             plots_in_log, nrows, do_2sigma_regions, plot_bestfit_value, $
             plot_previous_model_1, par_previous_model_1, linestyle_previous_model_1, $
             plot_previous_model_2, par_previous_model_2, linestyle_previous_model_2, $
             fix_max_like, max_log_like

  ;!p.multi=[0,5,4]
  multiplot, /reset
  ;erase & multiplot, [7,1]
  loadct,6
  color=[200,120,70]  

  for i=0, N_par-1 do begin
     xmin=par_min[i]
     xmax=par_max[i]
     bin=(xmax-xmin)/float(Nbins)
     ymin=0.0
     ymax=1.0
     
     if(plots_in_log eq 1) then begin
        xTICKINTERVAL=0.5
        yTICKINTERVAL=0.2
        y_ntick=(ymax-ymin)/yTICKINTERVAL+1
     endif else begin
        ;xTICKINTERVAL=(xmax-xmin)/2.
        ;x_ntick=(xmax-xmin)/xTICKINTERVAL+1


  
    ;  limits[14,*] = [0.02,0.06] ;yield
    ;  limits[15,*] = [0.2, 0.7] ;Rmerger
    ;  limits[16,*] = [ 1.5, 4.0] ;Mergertime
     ; limits[17,*] = [ 7.e3, 2.e4] ;M_RamPressure



       ; xTICKV=[[0.02,0.04],[0.075,0.15],[0.5,1.0],[1.0,2.0],[0.0075,0.015],[0.02,0.06],[500.,1000.], $
       ;         [1.5,3.0],[300.,600.],[0.4,0.8],[0.3,0.6],[50.,100.],[0.6,1.2],[1.5e10,3.e10], $
       ;         [0.025,0.05],[0.3,0.6],[2.0,3.0],[7.5e3,1.5e4]]
       ; xtickname=[['0.02','0.04'],['0.075','0.15'],['0.5','1.0'],['1.0','2.0'],['0.0075','0.015'],['0.02','0.06'],['500.','1000.'], $
        ;           ['1.5','3.0'],['300.','600.'],['0.4','0.8'],['0.3','0.6'],['50','100'],['0.6','1.2'],['1.5e10','3e10'], $
        ;           ['0.025','0.05'],['0.3','0.6'],['2.0','3.0'],['7.5e3','1.5e4']]
 
       

   xTICKV=[[0.008,0.009,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08], $
           ;[0.05,0.06,0.07,0.08,0.09,0.10,0.20,0.30,0.40,0.50], $
           [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0], $
           [0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,2.00], $
           [0.40, 0.50,0.60,0.70,0.80,0.90,1.00,2.00,3.00,4.00], $

           [0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02], $
           [0.009,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09], $
           [200.,300.,400.,500.,600.,700.,800.,900,1000.,2000.], $

           [1.00,2.00,3.00,4.00,5.00,6.00,7.00,8.00,9.00,10.0], $
           [70.0,80.0,90.0,100.,200.,300.,400.,500.,600.,700.], $
           [0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,2.00], $

           [0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00], $
           [20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.,200.], $
           [0.40,0.50,0.60,0.70,0.80,0.90,1.00,2.00,3.00,4.00], $
           [1.e10,2.e10,3.e10,4.e10,5.e10,6.e10,7.e10,8.e10,9.e10,1.e11], $

           [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10], $
           ;[0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,2.00], $
           [1.00,2.00,3.00,4.00,5.00,6.00,7.00,8.00,9.00,10.0], $
           [3.e3,4.e3, 5.e3,6.e3,7.e3,8.e3,9.e3,1.e4,2.e4,3.e4]]

        xtickname=[[' ',' ','0.01',' ',' ',' ','0.05',' ',' ',' '], $
                   ['0.1',' ',' ',' ','0.5',' ',' ',' ',' ','1.0'], $
                   [' ',' ',' ','0.5',' ',' ',' ',' ','1.0',' '], $
                   [' ','0.5',' ',' ',' ',' ','1.0',' ','3.0',' '], $

                   [' ',' ',' ','0.005',' ',' ',' ',' ','0.01',' '], $
                   [' ','0.01',' ',' ',' ','0.05',' ',' ',' ',' '], $
                   [' ',' ',' ','500.',' ',' ',' ',' ','1000.',' '], $

                   ['1.0',' ',' ',' ','5.0',' ',' ',' ',' ','10.0'], $
                   [' ',' ',' ','100.',' ',' ',' ','500.',' ',' '], $
                   [' ',' ',' ','0.5',' ',' ',' ',' ','1.0',' '], $

                   ['0.1',' ',' ',' ','0.5',' ',' ',' ',' ','1.0'], $
                   [' ',' ',' ','50.',' ',' ',' ',' ','100.',' '], $
                   [' ','0.5',' ',' ',' ',' ','1.0','2.0',' ',' '], $
                   ['1.e10',' ',' ',' ','5.e10',' ',' ',' ',' ','1.e11'], $

                   ['0.01',' ',' ',' ','0.05',' ',' ',' ',' ','0.1'], $
                   ;[' ',' ',' ','0.5',' ',' ',' ',' ','1.0',' '], $
                   ['1.0',' ',' ',' ','5.0',' ',' ',' ',' ','10.0'], $
                   [' ',' ','5.e3',' ',' ',' ',' ','1.e4',' ',' ']]
 

       TICKS=9
        ;xtickname=Replicate('',x_ntick)

        ;if(i eq 1) then xtickname=[' ','','']
        ;if(i eq 5) then xtickname=[' ','','']
        ;if(i eq 6) then xtickname=[' ','','']
        ;if(i eq 15) then xtickname=[' ','','']
        yTICKINTERVAL=0.2
        y_ntick=(ymax-ymin)/yTICKINTERVAL+1
        
     endelse

     ;only do chosen plots
     if(do_plot[i] eq 1) then begin
      ;print,xTICKv[*,i]
;print,xtickname[*,i]
        if(i eq 0 or i eq 7 or i eq 14) then begin
           ;if(i eq 7 or i eq 14) then multiplot
           plot_i=0
           ytitle = 'fraction'
           if(plots_in_log eq 0) then begin
              plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
                    xstyle = 1, ystyle = 1, /xlog, $
                    xTICKS=TICKS, xTICKv=xTICKv[*,i], xtickname=xtickname[*,i], yTICKINTERVAL=yTICKINTERVAL,  $
                    xtitle = par_name[i], ytitle = ytitle, $
                    Position=[0.05+0.135*(plot_i), 0.2, 0.04+0.135*(plot_i+1.), 0.90],charsize=1.1 
           endif else begin
              plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
                    xstyle = 1, xTICKINTERVAL=xTICKINTERVAL, $
                    yTICKINTERVAL=yTICKINTERVAL, ystyle = 1, $
                    xtitle = par_name[i], ytitle = ytitle, $
                    Position=[0.05+0.135*(plot_i), 0.2, 0.04+0.135*(plot_i+1.), 0.90] 
           endelse
        endif else begin          
           plot_i+=1
           if(plots_in_log eq 0) then begin
              plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
                    xstyle = 1,  ystyle = 1, /xlog, $
                    xTICKS=TICKS, xTICKv=xTICKV[*,i], xtickname=xtickname[*,i], yTICKINTERVAL=yTICKINTERVAL,  $
                    xtitle = par_name[i], /NoErase , ytickname=Replicate(' ',y_ntick), $
                    Position=[0.05+0.135*(plot_i), 0.2, 0.04+0.135*(plot_i+1), 0.90],charsize=1.1 
           endif else begin
              plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
                    xstyle = 1, xTICKINTERVAL=xTICKINTERVAL, xtickname=xtickname,$
                    yTICKINTERVAL=yTICKINTERVAL, ystyle = 1, $
                    xtitle = par_name[i], /NoErase , ytickname=Replicate(' ',y_ntick), $
                    Position=[0.05+0.135*(plot_i), 0.2, 0.04+0.135*(plot_i+1), 0.90]
           endelse
        endelse
      
        ;plots for different MCMCs
        for j=0, N_MCMC-1 do begin
           if(plots_in_log eq 0) then begin          
              bin=(Alog10(xmax)-Alog10(xmin))/float(Nbins)
              hist=histogram(Alog10(data[j,i+2,*]),locations=c,min=alog10(xmin)-bin*2.,max=alog10(xmax)+bin*2.,binsize=bin)
              c=10^(c+bin/2.)            
              ;smooth the bining  
              for jj=0,1 do begin
                 for ii=1, n_elements(hist)-2 do hist[ii]=(hist[ii-1]+hist[ii]+hist[ii+1])/3.
              endfor
              print,max(hist)
              oplot,c,hist/(total(hist))/max(hist/(total(hist))), color=color[j], thick=6 
           endif else begin
              hist=histogram(data[j,i+2,*],locations=c,min=xmin-bin*2.,max=xmax+bin*2.,binsize=bin)
              ;smooth the bining 
              ;for jj=0,0 do begin
              ;   for ii=1, n_elements(hist)-2 do hist[ii]=(hist[ii-1]+hist[ii]+hist[ii+1])/3.
              ;endfor
              c+=bin/2.0             
              oplot,c,hist/(total(hist))/max(hist/(total(hist))), color=color[j], thick=6    
           endelse
         
        endfor
 
        ;FILL in 2sigma region - ONLY FOR FIRST SET OF CHAINS
        if(do_2sigma_regions eq 1) then begin
           ;compute regions
           upper=ConfidVal(i+2,(1-0.955)/2.,1, nrows, data)   
           lower=ConfidVal(i+2,(1-0.955)/2.,0, nrows, data)            
       ;0.955
            ;best fit
           if(fix_max_like eq 0) then $
              aux=data(0,i+2,[where(data(0,1,*) eq min(data(0,1,*)))]) $
           else $
              aux=data(0,i+2,[where(data(0,1,*) eq max_log_like)]) 
           ;if(i eq 0) then print,'      par   ','     best fit','         lower','          upper'    
           ;if(plots_in_log eq 1) then $     
           ;   print, i, 10^aux(0,0,0), 10^lower, 10^upper, aux(0,0,0), lower, upper $
           ;else $
           ;   print, i, aux(0,0,0), lower, upper, Alog10(aux(0,0,0)), Alog10(lower), Alog10(upper)
                     

           ;fill in
           ;for j=0, N_MCMC-1 do begin             
              sel=where(c gt lower and c lt upper)          
              Polyfill, [c[sel],Reverse(c[sel])], $
                        [hist[sel]/(total(hist[sel]))/max(hist[sel]/(total(hist[sel]))), fltarr(n_elements(c[sel]))], $
                        Color=fsc_color("skyblue") ;, /LINE_FILL, Spacing=0.2, Orientation=45                      
               
              if(plots_in_log eq 1) then sigma=(upper-lower)/2. else sigma=Alog10(upper/lower)/2.

              sigma=STRTRIM(number_formatter(sigma, DECIMAL=2),2)
              plot_label, xlog=1, ylog=0, type='label', label='2!7r!X!Dlogx!N='+sigma , xmin, xmax, ymin, ymax, $
                          x_percentage=2.5, x2_percentage=0., y_percentage=10.5, charsize=1.  
           ;endfor
        endif
   

        AXIS, xTICKS=TICKS, xTICKv=xTICKV[*,i], xtickname=xtickname[*,i], charsize=1.1
 
        ;plot best fit value
        for j=0, N_MCMC-1 do begin
           if(plot_bestfit_value eq 1) then begin
              x=FltArr(100)
              y=findgen(100)/100.
              if(fix_max_like eq 0) then $
                 aux=data(j,i+2,[where(data(j,1,*) eq min(data(j,1,*)))]) $
              else $
                 aux=data(j,i+2,[where(data(j,1,*) eq max_log_like)])

              if(plots_in_log eq 0) then $
                 oplot, x+aux(0,0,0), y, linestyle=0,color=color[j],thick=6 $
              else $
                 oplot, x+aux(0,0,0), y, linestyle=0,color=color[j],thick=6 
           endif 
        endfor

        ;previous models 
        x=FltArr(100)
        y=findgen(100)/100.        
        if (plot_previous_model_2 eq 1 and i ne 4) then $
           oplot, x+par_previous_model_2[i], y, linestyle=linestyle_previous_model_2 ,color=200,thick=6  
        if (plot_previous_model_1 eq 1) then $       
           oplot, x+par_previous_model_1[i], y, linestyle=linestyle_previous_model_1,color=70,thick=6     
     endif

  endfor
    


;plot labels to the right of last plot
  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 4, ystyle = 4, /NoErase , $
        Position=[0.05+0.135*(plot_i+1), 0.2, 0.04+0.135*(plot_i+2), 0.90]
  
  
  ;if(plots_in_log eq 1) then xlog=0 else xlog=0
  xlog=0
 
;THIS WORK
        plot_label, xlog=xlog, ylog=0, type='label', label='This Work-Max Like' , xmin, xmax, ymin, ymax, $
                    x_percentage=2.1, x2_percentage=0., y_percentage=8.0, charsize=1.         
        plot_label, xlog=xlog, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                    x_percentage=0.0, x2_percentage=1.6, y_percentage=8.3, $
                    linestyle=0, color=200, linethick=8

        plot_label, xlog=xlog, ylog=0, type='label', label='This Work-2!7r!X' , xmin, xmax, ymin, ymax, $
                    x_percentage=2.1, x2_percentage=0., y_percentage=6.8 , charsize=1.        
        plot_label, xlog=xlog, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                    x_percentage=0.0, x2_percentage=1.6, y_percentage=7.0, $
                    linestyle=0, color=fsc_color('sky blue'), linethick=20

;Henriques2013a, WMAP7
        plot_label, xlog=xlog, ylog=0, type='label', label='Henriques2013' , xmin, xmax, ymin, ymax, $
                    x_percentage=2.1, x2_percentage=0., y_percentage=5.5, charsize=1.
        plot_label, xlog=xlog, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                    x_percentage=0.0, x2_percentage=1.6, y_percentage=5.7, $
                    linestyle=linestyle_previous_model_2, color=200, linethick=8

;G11, WMAP7
        plot_label, xlog=xlog, ylog=0, type='label', label='Guo2013' , xmin, xmax, ymin, ymax, $
                    x_percentage=2.1, x2_percentage=0., y_percentage=4.2, charsize=1.
        plot_label, xlog=xlog, ylog=0, type='line', xmin, xmax, ymin, ymax, $
                    x_percentage=0.0, x2_percentage=1.6, y_percentage=4.4, $
                    linestyle=linestyle_previous_model_1, color=70, linethick=8
     





end


function ConfidVal, ix,limfrac, upper, nrows, data
      
  l=0
  t=nrows-1
 
  try_b = min(data(0,ix,l:t))
  try_t = max(data(0,ix,l:t))

  samps = total(data(0,0,l:t))

  lasttry = -1

  if (upper eq 1) then begin

     while(1) do begin        
        aux=data(0,ix,l:t)       
        sel=where(aux gt (try_b + try_t)/2, try)       
        if (try gt samps*limfrac) then $
           try_b = (try_b+try_t)/2 $
        else $
           try_t = (try_b+try_t)/2
        
        if (try eq lasttry) then break
        lasttry = try
      
     endwhile

  endif else begin
   
     while(1) do begin
        aux=data(0,ix,l:t)       
        sel=where(aux lt (try_b + try_t)/2, try)
    
        if (try gt samps*limfrac) then $
           try_t = (try_b+try_t)/2 $
        else $
           try_b = (try_b+try_t)/2
    
        if (try eq lasttry) then break
        lasttry = try
     endwhile

  endelse

  return, try_t
end  


pro get_par_limits, N_par, par_min, par_max, plots_in_log

   limits=fltarr(N_par, 2)

   clustering=0

   if(plots_in_log eq 1) then begin    
      limits[0,*]  = [-2.0,-1.0] ;SFE
      limits[1,*]  = [-1.6,-0.6] ;SFE_TRESH
      limits[2,*]  = [-0.9, 0.1] ;SF_burst1
      limits[3,*]  = [-0.4, 0.6] ;SF_burst2  

      ;limits[4,*]  = [-5.0,-4.0] ;KAGN
      limits[4,*]  = [-2.8,-1.8] ;KAGN
      limits[5,*]  = [-1.6,-0.6] ;FBH
      limits[6,*]  = [ 2.3, 3.3] ;VAGN
      
      limits[7,*]  = [ 0.0, 1.0] ;Edisk
      limits[8,*]  = [ 2.4, 3.4] ;Vdisk
      limits[9,*]  = [-0.8, 0.2] ;Beta1
      limits[10,*] = [-0.8, 0.2] ;Ehalo
      limits[11,*] = [ 1.9, 2.9] ;Vhalo
      limits[12,*] = [-0.8, 0.2] ;beta2
      limits[13,*] = [10.,11.] ;gama
     
      
      limits[14,*] = [-2.0,-1.0] ;yield
      ;limits[15,*] = [-1.2, -0.2] ;Rmerger
      limits[15,*] = [ -0.3, 0.7] ;Mergertime
      limits[16,*] = [ 3.5, 4.5] ;M_RamPressure
      ;limits[16,*] = [ 3.0, 6.0] ;M_RamPressure
      ;limits[18,*] = [ 0.5, 1.5] ;z0
      ;limits[19,*] = [ 0.5, 1.5] ;zi


      if(clustering eq 1) then begin
         limits[0,*]  = [-3.0, 0.0] ;SFE
         limits[1,*]  = [-3.0, 0.0] ;SFE_TRESH
         limits[2,*]  = [-2.0, 1.0] ;SF_burst1
         limits[3,*]  = [-2.0, 1.0] ;SF_burst2 
                      
         limits[4,*]  = [-5.0,-1.0] ;KAGN
         limits[5,*]  = [-3.0, 0.0] ;FBH
         limits[6,*]  = [ 1.0, 4.0] ;VAGN
      
         limits[7,*]  = [-0.5, 1.5] ;Edisk
         limits[8,*]  = [ 0.5, 3.5] ;Vdisk
         limits[9,*]  = [-0.5, 1.5] ;Beta1
         limits[10,*] = [-2.0, 1.0] ;Ehalo
         limits[11,*] = [ 0.5, 3.0] ;Vhalo
         limits[12,*] = [-1.0, 1.0] ;beta2
         limits[13,*] = [-3.0, 1.0] ;gama
         
         limits[14,*] = [-3.0, 0.0] ;yield
         ;limits[15,*] = [-2.0, 0.5] ;Rmerger
         limits[15,*] = [ -1.0, 1.0] ;Mergertime
         limits[16,*] = [ -7.0, -3.0] ;M_RamPressure
      endif

   endif else begin
      ;aaa=0
      ;if(aaa eq 1) then begin
    

      limits[0,*]  = [0.008,0.08] ;SFE
      ;limits[1,*]  = [0.05,0.5] ;SFE_TRESH
      limits[1,*]  = [0.1,1.0] ;SFE_TRESH
      limits[2,*]  = [0.2, 2.0] ;SF_burst1
      limits[3,*]  = [0.4, 4.0] ;beta_burst

      limits[4,*]  = [2.e-3,2.e-2] ;KAGN
      limits[5,*]  = [0.009,0.09] ;FBH
      limits[6,*]  = [200., 2000.] ;VAGN
      
      limits[7,*]  = [ 1., 10.0] ;Edisk
      limits[8,*]  = [ 70., 700.] ;Vdisk
      limits[9,*]  = [0.2, 2.] ;Beta1
      limits[10,*] = [0.1, 1.] ;Ehalo
      limits[11,*] = [ 20., 200.] ;Vhalo
      limits[12,*] = [0.4, 4.] ;beta2
      limits[13,*] = [1.e10,10.e10] ;gama
      
      limits[14,*] = [0.01,0.1] ;yield
      ;limits[15,*] = [0.2, 2.] ;Rmerger
      limits[15,*] = [ 1., 10.] ;Mergertime
      limits[16,*] = [ 3.e3, 3.e4] ;M_RamPressure

      ;limits[18,*] = [ 4, 14.] ;z0
      ;limits[19,*] = [ 4, 14] ;zi
    
      if(clustering eq 1) then begin
   
         limits[0,*]  = [0.01,0.1] ;SFE
         limits[1,*]  = [0.01,0.3] ;SFE_TRESH
         limits[2,*]  = [0.1, 3.0] ;SF_burst1
         limits[3,*]  = [0.2, 7.0] ;SF_burst2  
         
         limits[4,*]  = [1.e-3,5.e-2] ;KAGN
         limits[5,*]  = [0.005,0.2]   ;FBH
         limits[6,*]  = [100., 2000.] ;VAGN
         
         limits[7,*]  = [ 0.3, 10.0] ;Edisk
         limits[8,*]  = [ 30., 2000.] ;Vdisk
         limits[9,*]  = [0.05, 3.0]   ;Beta1
         limits[10,*] = [0.03, 10.]   ;Ehalo
         limits[11,*] = [ 10., 1000.] ;Vhalo
         limits[12,*] = [0.1, 10.0]   ;beta2
         limits[13,*] = [5.e9,2.e11]  ;gama
         
         limits[14,*] = [0.01,0.1] ;yield
         ;limits[15,*] = [0.05, 1.0] ;Rmerger
         limits[15,*] = [ 0.5, 10.0] ;Mergertime
         limits[16,*] = [ 1.e3, 1.e5] ;M_RamPressure
      endif
   ;endif
   endelse

   for i=0, N_par-1 do begin  
      par_min[i]=limits[i,0]
      par_max[i]=limits[i,1]   
   endfor

end


pro get_par_names, N_par, par_name, plots_in_log

  if(plots_in_log eq 0) then begin
     par_name[0]='!7a!X!DSF!N'
     par_name[1]='SF!D_threshold'
     par_name[2]='!7a!X!DSF_burst'
     par_name[3]='!7b!X!Dburst'
     par_name[4]='K!DAGN!N'
     par_name[5]='f!DBH!N'
     par_name[6]='V!DAGN!N'
     par_name[7]='!7e!X'
     par_name[8]='V!Dreheat!N'
     par_name[9]='!7b!X!D1'
     par_name[10]='!7g!X'
     par_name[11]='V!Deject!N'
     par_name[12]='!7b!X!D2'
     par_name[13]='!7c!X'
     par_name[14]='Z!Dyield!N'
     ;par_name[15]='R!Dmerger!N'
     par_name[15]='MergerTime'
     par_name[16]='M!DRamPressure!N'
     ;par_name[18]='z!D0!N'
     ;par_name[19]='z!Dr!N'
  endif else begin
      par_name[0]='log!D10!N(!7a!X!DSF!N)'
     par_name[1]='log!D10!N(SF!D_threshold!N)'
     par_name[2]='log!D10!N(!7a!X!DSF_burst)'
     par_name[3]='log!D10!N(!7b!X!Dburst)'
     par_name[4]='log!D10!N(K!DAGN!N)'
     par_name[5]='log!D10!N(f!DBH!N)'
     par_name[6]='log!D10!N(V!DAGN!N)'
     par_name[7]='log!D10!N(!7e!X)'
     par_name[8]='log!D10!N(V!Dreheat!N)'
     par_name[9]='log!D10!N(!7b!X!D!N)1'
     par_name[10]='log!D10!N(!7g!X)'
     par_name[11]='log!D10!N(V!Deject!N)'
     par_name[12]='log!D10!N(!7b!X!D2!N)'
     par_name[13]='log!D10!N(!7c!X)'
     par_name[14]='log!D10!N(Z!Dyield!N)'
     ;par_name[15]='log!D10!N(R!Dmerger!N)'
     par_name[15]='log!D10!N(MergerTime)'
     par_name[16]='log!D10!N(M!DRamPressure!N)'
     ;par_name[18]='log!D10!N(z!D0!N)'
     ;par_name[19]='log!D10!N(z!Dr!N)'
  endelse

   ;print,textoidl('\alpha')
   ;print,textoidl('\beta')
   ;print,textoidl('\gamma')
   ;print,textoidl('\epsilon')
   ;print,textoidl('\eta')

;print,textoidl('\Delta')
end




pro plot_feedback_scalings, N_two_sigma_regions=N_two_sigma_regions,nvmax=nvmax, parameters_model1, H_z_1, do_previous_model1,$
                              parameters_model2, H_z_2, do_previous_model2, $
                              parameters_bestfit,bestfit_min,bestfit_max, Omega_M, Omega_Lambda  

 

      ;!p.charsize = 2.5
      ;!p.charthick = 4
      ;!p.thick = 4
      ;!x.thick = 4
      ;!y.thick = 4
      set_plot, 'PS'
      device, filename = './fig/parameter_scalings_sn.ps', xsize = 25, ysize = 20, /color, xoffset=1, yoffset=5
 

   multiplot, /reset
   !p.multi=[0,2,2]

   loadct,6

 
  
   xmin=1.3
   xmax=3.0
    

   make_float_array, v_max, xmin, xmax, (xmax-xmin)/(nvmax*1.0)
   make_float_array, mass_ratio, 0.0, 10.0, (xmax-xmin)/(nvmax*1.0)
   make_float_array, redshift, 0.0, 10.0, (10.0)/(nvmax*1.0)

   v_max=10^v_max
  

   ;mass_ratio=findgen(nvmax)/1000. 
   ;redshift=findgen(nvmax)/100. 

   G=4.3*1.e-9

   seed=5l

  ;print,v_max


;MODEL 1
   sfe_model1=parameters_model1[0]
   sfe_thres_model1=parameters_model1[1]
   sf_burst1_model1=parameters_model1[2]
   sf_burst2_model1=parameters_model1[3]
   kagn_model1=parameters_model1[4]
   fbh_model1=parameters_model1[5]
   Vbh_model1=parameters_model1[6]
   epsilon_model1=parameters_model1[7]
   V_disk_model1=parameters_model1[8]
   beta_1_model1=parameters_model1[9]
   eta_model1=parameters_model1[10]
   V_halo_model1=parameters_model1[11]
   beta_2_model1=parameters_model1[12]
   gamma_model1=parameters_model1[13]

   epsilon_halo_model1=eta_model1*(0.5+(V_max/V_halo_model1)^(-beta_2_model1))
   epsilon_halo_model1[where(epsilon_halo_model1 gt 1.0)]=1.0

   V_SN=630.
   epsilon_disk_model1=epsilon_model1*(0.5+(V_max/V_disk_model1)^(-beta_1_model1))
   sel=where(epsilon_disk_model1 gt epsilon_halo_model1*V_SN^2/V_max^2,n)  
   epsilon_disk_model1[sel]=epsilon_halo_model1[sel]*V_SN^2/V_max[sel]^2

   eject_eff_model1=(epsilon_halo_model1*V_SN^2)/V_max^2-epsilon_disk_model1

   ;reincorporation 
   ;H_z=fltarr(4)
   H_z=H_z_1
   tdyn_0=1/(10.*H_z[0])
   tdyn_1=1/(10.*H_z[1])
   tdyn_2=1/(10.*H_z[2])
   tdyn_3=1/(10.*H_z[3])
   reinc_model1_0=1/gamma_model1*tdyn_0/(V_max/220)*1.e12
   reinc_model1_1=1/gamma_model1*tdyn_1/(V_max/220)*1.e12
   reinc_model1_2=1/gamma_model1*tdyn_2/(V_max/220)*1.e12
   reinc_model1_3=1/gamma_model1*tdyn_3/(V_max/220)*1.e12

   ;AGN feedback
   AGN_feedback_model1=kagn_model1*(10.*G*hubble_of_z(redshift, Omega_M, Omega_Lambda)*1.e11)/(200.*200.*200.*0.1)
   quasar_model1=fbh_model1/(1.+Vbh_model1/V_max)^2
 
   ;sf burst
   sf_burst_model1=sf_burst1_model1*mass_ratio^sf_burst2_model1 

;MODEL 2
   sfe_model2=parameters_model2[0]
   sfe_thres_model2=parameters_model2[1]
   sf_burst1_model2=parameters_model2[2]
   sf_burst2_model2=parameters_model2[3]
   kagn_model2=parameters_model2[4]
   fbh_model2=parameters_model2[5]
   Vbh_model2=parameters_model2[6]
   epsilon_model2=parameters_model2[7]
   V_disk_model2=parameters_model2[8]
   beta_1_model2=parameters_model2[9]
   eta_model2=parameters_model2[10]
   V_halo_model2=parameters_model2[11]
   beta_2_model2=parameters_model2[12]
   gamma_model2=parameters_model2[13]

   epsilon_halo_model2=eta_model2*(0.5+(V_max/V_halo_model2)^(-beta_2_model2))
   epsilon_halo_model2[where(epsilon_halo_model2 gt 1.0)]=1.0

   V_SN=630.
   epsilon_disk_model2=epsilon_model2*(0.5+(V_max/V_disk_model2)^(-beta_1_model2))
   sel=where(epsilon_disk_model2 gt epsilon_halo_model2*V_SN^2/V_max^2,n)  
   epsilon_disk_model2[sel]=epsilon_halo_model2[sel]*V_SN^2/V_max[sel]^2

   eject_eff_model2=(epsilon_halo_model2*V_SN^2)/V_max^2-epsilon_disk_model2

   ;reincorporation 
   ;H_z=fltarr(4)
   H_z=H_z_2
   ;G=6.672e-8
   ;G=4.3*1.e-9
   tdyn_0=1/(10.*H_z[0])
   tdyn_1=1/(10.*H_z[1])
   tdyn_2=1/(10.*H_z[2])
   tdyn_3=1/(10.*H_z[3]) 

   reinc_model2_0=gamma_model2*1.e10/(V_max^3/(10.*G*H_z[0]))
   reinc_model2_1=gamma_model2*1.e10/(V_max^3/(10.*G*H_z[1]))
   reinc_model2_2=gamma_model2*1.e10/(V_max^3/(10.*G*H_z[2]))
   reinc_model2_3=gamma_model2*1.e10/(V_max^3/(10.*G*H_z[3]))

   ;AGN feedback
   AGN_feedback_model2=kagn_model2*(10.*G*hubble_of_z(redshift, Omega_M, Omega_Lambda)*1.e11)/(200.*200.*200.*0.1)
   quasar_model2=fbh_model2/(1.+Vbh_model2/V_max)^2


   ;sf burst
   sf_burst_model2=sf_burst1_model2*mass_ratio^sf_burst2_model2




;best fit 
   sfe       =parameters_bestfit[0]
   sfe_thres = parameters_bestfit[1]
   sf_burst1 = parameters_bestfit[2]
   sf_burst2 = parameters_bestfit[3] 
   kagn      = parameters_bestfit[4]
   fbh       = parameters_bestfit[5]
   Vbh       = parameters_bestfit[6]
   epsilon   = parameters_bestfit[7]
   V_disk    = parameters_bestfit[8]
   beta_1    = parameters_bestfit[9]
   eta       = parameters_bestfit[10]
   V_halo    = parameters_bestfit[11]
   beta_2    = parameters_bestfit[12]
   gamma     = parameters_bestfit[13]
 
   epsilon_halo_best=eta*(0.5+(V_max/V_halo)^(-beta_2))
   epsilon_halo_best[where(epsilon_halo_best gt 1.0)]=1.0


   V_SN=630.
   epsilon_disk_best=epsilon*(0.5+(V_max/V_disk)^(-beta_1))
  
   sel=where(epsilon_disk_best gt epsilon_halo_best*V_SN^2/V_max^2,n)  
   epsilon_disk_best[sel]=epsilon_halo_best[sel]*V_SN^2/V_max[sel]^2

   eject_eff_best=(epsilon_halo_best*V_SN^2)/V_max^2-epsilon_disk_best

   ;AGN feedback    
   AGN_feedback_bestfit=kagn*redshift/redshift

   quasar_bestfit=fbh/(1.+Vbh/V_max)^2

   ;sf burst
   sf_burst_bestfit=sf_burst1*mass_ratio^sf_burst2

   sfe_max       = bestfit_max[0]
   sfe_thres_max = bestfit_max[1]
   sf_burst1_max = bestfit_max[2]
   sf_burst2_max = bestfit_max[3]
   kagn_max      = bestfit_max[4]
   fbh_max       = bestfit_max[5]
   Vbh_max       = bestfit_max[6]
   epsilon_max   = bestfit_max[7]
   V_disk_max    = bestfit_max[8]
   beta_1_max    = bestfit_max[9]
   eta_max       = bestfit_max[10]
   V_halo_max    = bestfit_max[11]
   beta_2_max    = bestfit_max[12]

   sfe_min       = bestfit_min[0]
   sfe_thres_min = bestfit_min[1]
   sf_burst1_min = bestfit_min[2]
   sf_burst2_min = bestfit_min[3]
   kagn_min      = bestfit_min[4]
   fbh_min       = bestfit_min[5]
   Vbh_min       = bestfit_min[6]
   epsilon_min   = bestfit_min[7]
   V_disk_min    = bestfit_min[8]
   beta_1_min    = bestfit_min[9]
   eta_min       = bestfit_min[10]
   V_halo_min    = bestfit_min[11]
   beta_2_min    = bestfit_min[12]



   ;reincorporation 
   ;H_z=fltarr(4)  
   ;G=6.672e-8
   ;G=4.3*1.e-9
   reinc_best_0=gamma*1.e10/(V_max^3/(10.*G*hubble_of_z(0.0, Omega_M, Omega_Lambda)))
   reinc_best_1=gamma*1.e10/(V_max^3/(10.*G*hubble_of_z(1.0, Omega_M, Omega_Lambda)))
   reinc_best_2=gamma*1.e10/(V_max^3/(10.*G*hubble_of_z(2.0, Omega_M, Omega_Lambda)))
   reinc_best_3=gamma*1.e10/(V_max^3/(10.*G*hubble_of_z(3.0, Omega_M, Omega_Lambda)))

   print,gamma * 1.e10*10.*G*H_z[0]/(100.^3)



;2 SIGAM REGIONS  
   for i=0, N_two_sigma_regions do begin

      if i mod 1000 eq 0 then print,i 
    
      N=1000.
      epsilon=findgen(N)/N*(epsilon_max-epsilon_min)+epsilon_min+(epsilon_max-epsilon_min)/N/2.    ;5.2  ;3.9-6.1
      V_disk=findgen(N)/N*(V_disk_max-V_disk_min)+V_disk_min+(V_disk_max-V_disk_min)/N/2. ;166.   ;99-351
      beta_1=findgen(N)/N*(beta_1_max-beta_1_min)+beta_1_min+(beta_1_max-beta_1_min)/N/2. ;0.45   ;0.31-0.66 
      eta=findgen(N)/N*(eta_max-eta_min)+eta_min+(eta_max-eta_min)/N/2.    ;0.89      ;0.46-1.61
      V_halo=findgen(N)/N*(V_halo_max-V_halo_min)+V_halo_min+(V_halo_max-V_halo_min)/N/2. ;499.   ;283-960
      beta_2=findgen(N)/N*(beta_2_max-beta_2_min)+beta_2_min+(beta_2_max-beta_2_min)/N/2.   ;0.34   ;0.32-1.1     
      ;AGN
      kagn=findgen(N)/N*(kagn_max-kagn_min)+kagn_min+(kagn_max-kagn_min)/N/2.
      fbh=findgen(N)/N*(fbh_max-fbh_min)+fbh_min+(fbh_max-fbh_min)/N/2.
      Vbh=findgen(N)/N*(Vbh_max-Vbh_min)+Vbh_min+(Vbh_max-Vbh_min)/N/2.
      ;sf burst
      sf_burst1=findgen(N)/N*(sf_burst1_max-sf_burst1_min)+sf_burst1_min+(sf_burst1_max-sf_burst1_min)/N/2.
      sf_burst2=findgen(N)/N*(sf_burst2_max-sf_burst2_min)+sf_burst2_min+(sf_burst2_max-sf_burst2_min)/N/2.

      epsilon_halo=eta[fix(randomu(seed)*1000.)]*(0.5+(V_max/V_halo[fix(randomu(seed)*1000.)])^(-beta_2[fix(randomu(seed)*1000.)]))
      epsilon_halo[where(epsilon_halo gt 1.0)]=1.0
  
      epsilon_disk=epsilon[fix(randomu(seed)*1000.)]*(0.5+(V_max/V_disk[fix(randomu(seed)*1000.)])^(-beta_1[fix(randomu(seed)*1000.)]))

      sel=where(epsilon_disk gt epsilon_halo*V_SN^2/V_max^2,n)  
      epsilon_disk[sel]=epsilon_halo[sel]*V_SN^2/V_max[sel]^2

      eject_eff=(epsilon_halo*V_SN^2)/V_max^2-epsilon_disk

      AGN_feedback=kagn[fix(randomu(seed)*1000.)]*redshift/redshift
      quasar=fbh[fix(randomu(seed)*1000.)]/(1.+Vbh[fix(randomu(seed)*1000.)]/V_max)^2
      sf_burst=sf_burst1[fix(randomu(seed)*1000.)]*mass_ratio^sf_burst2[fix(randomu(seed)*1000.)]


      if(i eq 0) then begin
         epsilon_disk_max=epsilon_disk
         epsilon_disk_min=epsilon_disk
         epsilon_halo_max=epsilon_halo
         epsilon_halo_min=epsilon_halo
         eject_eff_max=eject_eff
         eject_eff_min=eject_eff
         AGN_feedback_max=AGN_feedback
         AGN_feedback_min=AGN_feedback
         quasar_max=quasar
         quasar_min=quasar
         sf_burst_max=sf_burst
         sf_burst_min=sf_burst     
      endif else begin
         for j=0, nvmax-1 do begin
            if(epsilon_disk[j] gt epsilon_disk_max[j]) then epsilon_disk_max[j]=epsilon_disk[j]
            if(epsilon_disk[j] lt epsilon_disk_min[j]) then epsilon_disk_min[j]=epsilon_disk[j]
            if(epsilon_halo[j] gt epsilon_halo_max[j]) then epsilon_halo_max[j]=epsilon_halo[j]
            if(epsilon_halo[j] lt epsilon_halo_min[j]) then epsilon_halo_min[j]=epsilon_halo[j]
            if(eject_eff[j] gt eject_eff_max[j]) then eject_eff_max[j]=eject_eff[j]
            if(eject_eff[j] lt eject_eff_min[j]) then eject_eff_min[j]=eject_eff[j]
            if(AGN_feedback[j] gt AGN_feedback_max[j]) then AGN_feedback_max[j]=AGN_feedback[j]
            if(AGN_feedback[j] lt AGN_feedback_min[j]) then AGN_feedback_min[j]=AGN_feedback[j]
            if(quasar[j] gt quasar_max[j]) then quasar_max[j]=quasar[j]
            if(quasar[j] lt quasar_min[j]) then quasar_min[j]=quasar[j]
            if(sf_burst[j] gt sf_burst_max[j]) then sf_burst_max[j]=sf_burst[j]
            if(sf_burst[j] lt sf_burst_min[j]) then sf_burst_min[j]=sf_burst[j]           
         endfor
      endelse

      ;oplot,Alog10(v_max), epsilon_disk, linestyle=0, color=120, thick=6
  

;print,'max=',epsilon_disk_max
;print,'min=',epsilon_disk_min
   endfor


;*******************
;** PLOTS      *****
;*******************

;DISK
   ymin=-1.0
   ymax=4.5 
 
   plot, findgen(10),  /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
         xstyle = 1, ystyle = 1, xtitle = 'log!D10!N(V!Dmax!N[kms!U-1!N])', ytitle = 'log!D10!N(Mass Loading (!7e!X!Ddisk!N))'  
   Polyfill, [Alog10(v_max),Reverse(Alog10(v_max))], $
             [Alog10(epsilon_disk_min),Reverse(Alog10(epsilon_disk_max))], Color=fsc_color("sky blue"), NOCLIP = 0

   sel=where(Alog10(v_max) lt 3.675)
   oplot,Alog10(v_max[sel]), Alog10(epsilon_disk_best[sel]), linestyle=0, color=200, thick=6 
   oplot,Alog10(v_max), Alog10(epsilon_disk_model1), linestyle=1, color=200, thick=6
   oplot,Alog10(v_max), Alog10(epsilon_disk_model2), linestyle=2, color=200, thick=6
  

;DELUCIA2007
   oplot,Alog10(v_max), alog10(v_max/v_max+2.5), linestyle=2, color=120, thick=6

;Henriques2009
   ;oplot,Alog10(v_max), alog10(v_max/v_max+9.28), linestyle=2, color=70, thick=6

;BOUCHE OBSERVATIONS

   ;bouche_Vout     = [175., 500., 300., 175., 225.]
   ;bouche_Vout_err = [25.,  100., 25.,   80., 50.]
   bouche_Vout     = [131., 231., 162., 116., 240.]

   bouche_sfr      = [1.27, 0.26, 3.75, 1.96, 1.36]
   bouche_Mout     = [2.20, 6.80, 1.00, 6.00, 2.20]
   bouche_Mout_err = [1.10, 3.40, 0.50, 3.00, 1.10]
   
   bouche_loading = bouche_Mout/bouche_sfr
   bouche_loading_err = bouche_Mout_err/bouche_sfr

   symbols,2,1.0
   ;oploterror, Alog10(bouche_Vout), Alog10(bouche_Mout/bouche_sfr),Alog10(bouche_Mout_err/bouche_sfr), color = 70, errcolor = 70, psym = 8
     
   oploterror, Alog10(bouche_Vout), Alog10(bouche_loading), $
               Alog10((bouche_loading+bouche_loading_err)/bouche_loading), color = 70, $
               errcolor = 70, psym = 8, /hiba,HATLENGTH = 80.0
   oploterror, Alog10(bouche_Vout), Alog10(bouche_Mout/bouche_sfr), $
               Alog10(bouche_loading/(bouche_loading-bouche_loading_err)), color = 70, $
               errcolor = 70, psym = 8, /loba,HATLENGTH = 80.0


   plot_label, xlog=0, ylog=0, type='label', label='Bouche (2012)' , xmin, xmax, ymin, ymax, $
              x_percentage=0.8, x2_percentage=0., y_percentage=6.2, charsize=1.
   plot_label, xlog=0, ylog=0, type='symbol', label=label, xmin, xmax, ymin, ymax, $
               x_percentage=0.5, x2_percentage=0., y_percentage=6.35, $
               color=70, errcolor=70, sym_num=2, sym_size=1.0, HATLENGTH = 100, err_size=0.15
;LABELS
  
  plot_label, xlog=0, ylog=0, type='label', label='This Work - Best Fit' , xmin, xmax, ymin, ymax, $
              x_percentage=1.2, x2_percentage=0., y_percentage=9.2, charsize=1.
  plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.3, x2_percentage=1.0, y_percentage=9.35, $
              linestyle=0, color=200, linethick=8
      
  plot_label, xlog=0, ylog=0, type='label', label='This Work!M+!X2!7r!X' , xmin, xmax, ymin, ymax, $
              x_percentage=1.2, x2_percentage=0., y_percentage=8.6, charsize=1.
  plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.3, x2_percentage=1.0, y_percentage=8.75, $
              linestyle=0, color=fsc_color("sky blue"), linethick=20
      

  plot_label, xlog=0, ylog=0, type='label', label='Henriques2013' , xmin, xmax, ymin, ymax, $
              x_percentage=1.2, x2_percentage=0., y_percentage=8.0, charsize=1.
  plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.3, x2_percentage=1.0, y_percentage=8.15, $
              linestyle=2, color=200, linethick=8

  plot_label, xlog=0, ylog=0, type='label', label='Guo2013' , xmin, xmax, ymin, ymax, $
              x_percentage=1.2, x2_percentage=0., y_percentage=7.4, charsize=1.
  plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.3, x2_percentage=1.0, y_percentage=7.55, $
              linestyle=1, color=200, linethick=8

  plot_label, xlog=0, ylog=0, type='label', label='DeLucia07' , xmin, xmax, ymin, ymax, $
              x_percentage=1.2, x2_percentage=0., y_percentage=6.8, charsize=1.
  plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.3, x2_percentage=1.0, y_percentage=6.95, $
              linestyle=2, color=120, linethick=8

  plot_axis, xmin, xmax, ymin, ymax

print,textoidl('\sigma')
print,textoidl('\epsilon')
print,textoidl('\pm')

;HALO
   ymin=0
   ymax=1.1
   plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
         xstyle = 1, ystyle = 1, xtitle = 'log!D10!N(V!Dmax!N[kms!U-1!N])', ytitle = 'SN ejection efficiency (!7e!X!Dhalo!N)'  
   sel=where(Alog10(v_max) gt xmin )   
   Polyfill, [Alog10(v_max[sel]),Reverse(Alog10(v_max[sel]))], $
             [epsilon_halo_min[sel],Reverse(epsilon_halo_max[sel])], Color=fsc_color("sky blue"), NOCLIP = 0
   oplot,Alog10(v_max[sel]), epsilon_halo_best[sel], linestyle=0, color=200, thick=6 
   oplot,Alog10(v_max), epsilon_halo_model1, linestyle=1, color=200, thick=6
   oplot,Alog10(v_max), epsilon_halo_model2, linestyle=2, color=200, thick=6

;DELUCIA2007
   oplot,Alog10(v_max), v_max/v_max-1.+0.35, linestyle=2, color=120, thick=6

;Henriques2009
   ;oplot,Alog10(v_max), v_max/v_max-1.+0.53, linestyle=2, color=70, thick=6

   plot_axis, xmin, xmax, ymin, ymax

;meject/mstar (eject efficiency)
   ;xmin=0.0

   ymin=-1.0
   ymax=4.0
   plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], NOCLIP = 0, $
         xstyle = 1, ystyle = 1, xtitle = 'log!D10!N(V!Dmax!N[kms!U-1!N])', ytitle = 'log!D10!N(!7D!XM!Deject!N/!7D!XM!D*!N)' 
   sel=where(Alog10(v_max) gt xmin and eject_eff_max gt ymin and eject_eff_min gt ymin)  
 
   sel=where(Alog10(eject_eff_min) lt ymin and Alog10(eject_eff_max) gt ymin)
   eject_eff_min[sel]=10^ymin
  
   sel=where(Alog10(v_max) lt 2.41)
   Polyfill, [Alog10(v_max[sel]),Reverse(Alog10(v_max[sel]))], $
             [Alog10(eject_eff_min[sel]),Reverse(Alog10(eject_eff_max[sel]))], Color=fsc_color("sky blue"), NOCLIP = 0


   oplot,Alog10(v_max), Alog10(eject_eff_best), linestyle=0, color=200, thick=6 
   oplot,Alog10(v_max), Alog10(eject_eff_model1), linestyle=1, color=200, thick=6
   oplot,Alog10(v_max), Alog10(eject_eff_model2), linestyle=2, color=200, thick=6
 

;DELUCIA2007
   eject_eff_delucia=(0.35*V_SN^2)/V_max^2-3.5
   oplot,Alog10(v_max), Alog10(eject_eff_delucia), linestyle=2, color=120, thick=6

;Henriques2009
   ;eject_eff_hen09=(0.53*V_SN^2)/V_max^2-10.28
   ;oplot,Alog10(v_max), Alog10(eject_eff_hen09), linestyle=2, color=70, thick=6


   plot_axis, xmin, xmax, ymin, ymax

;Reinc 
   ymin=8.
   ymax=11.
   plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
         xstyle = 1, ystyle = 1, xtitle = 'log!D10!N(V!D200c!N[kms!U-1!N])', ytitle = 'Reincorporation time [log(yr)]'
   ;sel=where(Alog10(v_max) gt xmin)   
   ;Polyfill, [Alog10(v_max[sel]),Reverse(Alog10(v_max[sel]))], $
   ;          [epsilon_halo_min[sel],Reverse(epsilon_halo_max[sel])], Color=fsc_color("sky blue") 
   ;oplot,Alog10(v_max), epsilon_halo_best, linestyle=0, color=200, thick=6 

   oplot,Alog10(v_max), Alog10(reinc_model1_0), linestyle=1, color=70, thick=6
   oplot,Alog10(v_max), Alog10(reinc_model1_1), linestyle=1, color=100, thick=6
   oplot,Alog10(v_max), Alog10(reinc_model1_2), linestyle=1, color=140, thick=6
   oplot,Alog10(v_max), Alog10(reinc_model1_3), linestyle=1, color=200, thick=6

   ;oplot,Alog10(v_max), Alog10(reinc_model2_0), linestyle=2, color=70, thick=6
   ;oplot,Alog10(v_max), Alog10(reinc_model2_1), linestyle=2, color=100, thick=6
   ;oplot,Alog10(v_max), Alog10(reinc_model2_2), linestyle=2, color=140, thick=6
   ;oplot,Alog10(v_max), Alog10(reinc_model2_3), linestyle=2, color=200, thick=6

   oplot,Alog10(v_max), Alog10(reinc_best_0), linestyle=0, color=70, thick=6
   oplot,Alog10(v_max), Alog10(reinc_best_1), linestyle=0, color=100, thick=6
   oplot,Alog10(v_max), Alog10(reinc_best_2), linestyle=0, color=140, thick=6
   oplot,Alog10(v_max), Alog10(reinc_best_3), linestyle=0, color=200, thick=6

;print,Alog10(v_max),Alog10(reinc_model2_3)
;print,Alog10(v_max),Alog10(reinc_best_3)

;sel=where(Alog10(reinc_model2_0) gt Alog10(reinc_best_0))
;print,Alog10(v_max[sel])
  

;z=0  
  y=ymin+(ymax-ymin)*(9.55-0/2.5)/10.
  oplot, [xmin+(xmax-xmin)*2.9/10.,xmin+(xmax-xmin)*3.5/10.], [y,y],Color=70, thick=6
  xyouts, xmin+(xmax-xmin)*3.7/10., ymin+(ymax-ymin)*(9.4-0/2.5)/10., 'This Work /', charsize=1., charthick=4

  y=ymin+(ymax-ymin)*(9.55-0/2.5)/10.
  oplot, [xmin+(xmax-xmin)*6.5/10.,xmin+(xmax-xmin)*7.3/10.], [y,y],Color=70, thick=6, linestyle=1
  xyouts, xmin+(xmax-xmin)*7.4/10., ymin+(ymax-ymin)*(9.4-0/2.5)/10., 'G11 (z=0)', charsize=1., charthick=4


;z=1 
  y=ymin+(ymax-ymin)*(9.55-1/2.)/10.
  oplot, [xmin+(xmax-xmin)*2.9/10.,xmin+(xmax-xmin)*3.5/10.], [y,y],Color=100, thick=6
  xyouts, xmin+(xmax-xmin)*3.7/10., ymin+(ymax-ymin)*(9.4-1/2.)/10., 'This Work /', charsize=1., charthick=4

  y=ymin+(ymax-ymin)*(9.55-1/2.)/10.
  oplot, [xmin+(xmax-xmin)*6.5/10.,xmin+(xmax-xmin)*7.3/10.], [y,y],Color=100, thick=6, linestyle=1
  xyouts, xmin+(xmax-xmin)*7.4/10., ymin+(ymax-ymin)*(9.4-1/2.)/10., 'G11 (z=1)', charsize=1., charthick=4


;z=2 
  y=ymin+(ymax-ymin)*(9.55-2/2.)/10.
  oplot, [xmin+(xmax-xmin)*2.9/10.,xmin+(xmax-xmin)*3.5/10.], [y,y],Color=140, thick=6
  xyouts, xmin+(xmax-xmin)*3.7/10., ymin+(ymax-ymin)*(9.4-2/2.)/10., 'This Work /', charsize=1., charthick=4

  y=ymin+(ymax-ymin)*(9.55-2/2.)/10.
  oplot, [xmin+(xmax-xmin)*6.5/10.,xmin+(xmax-xmin)*7.3/10.], [y,y],Color=140, thick=6, linestyle=1
  xyouts, xmin+(xmax-xmin)*7.4/10., ymin+(ymax-ymin)*(9.4-2/2.)/10., 'G11 (z=2)', charsize=1., charthick=4

;z=3
  y=ymin+(ymax-ymin)*(9.55-3/2.)/10.
  oplot, [xmin+(xmax-xmin)*2.9/10.,xmin+(xmax-xmin)*3.5/10.], [y,y],Color=200, thick=6
  xyouts, xmin+(xmax-xmin)*3.7/10., ymin+(ymax-ymin)*(9.4-3/2.)/10., 'This Work /', charsize=1., charthick=4

  y=ymin+(ymax-ymin)*(9.55-3/2.)/10.
  oplot, [xmin+(xmax-xmin)*6.5/10.,xmin+(xmax-xmin)*7.3/10.], [y,y],Color=200, thick=6, linestyle=1
  xyouts, xmin+(xmax-xmin)*7.4/10., ymin+(ymax-ymin)*(9.4-3/2.)/10., 'G11 (z=3)', charsize=1., charthick=4


  plot_axis, xmin, xmax, ymin, ymax

device,/close_file

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;;       AGN feedback
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


   !p.charsize = 1.4
      !p.charthick = 4
      !p.thick = 4
      !x.thick = 4
      !y.thick = 4
      set_plot, 'PS'
      device, filename = './fig/parameter_scalings_agn.ps', xsize = 25, ysize = 10, /color, xoffset=1, yoffset=5


   multiplot, /reset
   !p.multi=[0,2,1]


;;;;;;;;;;;;;;;;;;;;;;;
;;QUASAR
;;;;;;;;;;;;;;;;;;;;;;;
  xmin=1.3
  xmax=3.0
  ymin=-5.5
  ymax=-0.5
  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, xtitle = 'log!D10!N(V!D200c!N[kms!U-1!N])', ytitle = 'log!D10!N(!7D!XM!DBH,Q!N/M!Dcold!N)'

  ;BESTFIT
  sel=where(Alog10(v_max) gt xmin)
  Polyfill, [Alog10(v_max[sel]),Reverse(Alog10(v_max[sel]))], $
             [Alog10(quasar_min[sel]),Reverse(Alog10(quasar_max[sel]))], Color=fsc_color("sky blue")
  oplot,Alog10(v_max),  Alog10(quasar_bestfit), linestyle=0, color=200, thick=6 
  ;OLD MODELS
  oplot,Alog10(v_max),  Alog10(quasar_model1), linestyle=1, color=200, thick=6
  oplot,Alog10(v_max), Alog10(quasar_model2), linestyle=2, color=200, thick=6


;LABELS
  
  plot_label, xlog=0, ylog=0, type='label', label='This Work - Best Fit' , xmin, xmax, ymin, ymax, $
              x_percentage=1.2, x2_percentage=0., y_percentage=9.2, charsize=1.3
  plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.3, x2_percentage=1.0, y_percentage=9.45, $
              linestyle=0, color=200, linethick=8
      
  plot_label, xlog=0, ylog=0, type='label', label='This Work!M+!X2!7r!X' , xmin, xmax, ymin, ymax, $
              x_percentage=1.2, x2_percentage=0., y_percentage=8.5, charsize=1.3
  plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.3, x2_percentage=1.0, y_percentage=8.725, $
              linestyle=0, color=fsc_color("sky blue"), linethick=20
      

  plot_label, xlog=0, ylog=0, type='label', label='Henriques2013' , xmin, xmax, ymin, ymax, $
              x_percentage=1.2, x2_percentage=0., y_percentage=7.7, charsize=1.3
  plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.3, x2_percentage=1.0, y_percentage=7.95, $
              linestyle=2, color=200, linethick=8

  plot_label, xlog=0, ylog=0, type='label', label='Guo2013' , xmin, xmax, ymin, ymax, $
              x_percentage=1.2, x2_percentage=0., y_percentage=7.0, charsize=1.3
  plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.3, x2_percentage=1.0, y_percentage=7.25, $
              linestyle=1, color=200, linethick=8


  plot_label, xlog=0, ylog=0, type='label', label='f!DBH!N/(1+(V!DBH!N/V!D200c!N)!U2!N)' , xmin, xmax, ymin, ymax, $
              x_percentage=4.2, x2_percentage=0., y_percentage=1.0, charsize=1.3


  plot_axis, xmin, xmax, ymin, ymax
 

;;;;;;;;;;;;;;;;;;;;;;;
;;RADIO
;;;;;;;;;;;;;;;;;;;;;;;
  
  xmin=0.0
  xmax=10.0
  ymin=-4.0
  ymax=-1.0

  ;plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
  ;      xstyle = 1, ystyle = 1, xtitle = 'redshift', ytitle = 'log!D10!N(k!DAGN!N)'
  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, xtitle = 'redshift', ytitle = 'log!D10!N(k!DAGN!N[M!D!9ng!3!Nyr!U-1!N])'
  ;BESTFIT 
  sel=where(redshift gt xmin)
  Polyfill, [redshift[sel],Reverse(redshift[sel])], $
             [Alog10(AGN_feedback_min[sel]),Reverse(Alog10(AGN_feedback_max[sel]))], Color=fsc_color("sky blue")
  oplot,redshift,  Alog10(AGN_feedback_bestfit), linestyle=0, color=200, thick=6 

  ;OLD MODELS
  oplot,redshift,  Alog10(AGN_feedback_model1), linestyle=1, color=200, thick=6
  oplot,redshift,  Alog10(AGN_feedback_model2), linestyle=2, color=200, thick=6
  

  ;plot_label, xlog=0, ylog=0, type='label', label='k!DAGN!Nxf!DBH!N/(1+(V!DBH!N/V!DVir!N)!U2!N)' , xmin, xmax, ymin, ymax, $
   ;           x_percentage=3.2, x2_percentage=0., y_percentage=1.0, charsize=1.3
  ;plot_label, xlog=0, ylog=0, type='label', label='k!DAGN!N' , xmin, xmax, ymin, ymax, $
  ;            x_percentage=8.2, x2_percentage=0., y_percentage=1.0, charsize=1.3


  ;plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], xstyle = 1, ystyle = 1


  plot_axis, xmin, xmax, ymin, ymax

device,/close_file






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;;       SF
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


   !p.charsize = 1.4
      !p.charthick = 4
      !p.thick = 4
      !x.thick = 4
      !y.thick = 4
      set_plot, 'PS'
      device, filename = './fig/parameter_scalings_sf.ps', xsize = 14, ysize = 12, /color, xoffset=1, yoffset=5


   multiplot, /reset
   !p.multi=0




;;;;;;;;;;;;;;;;;;;;;;;
;;BURST
;;;;;;;;;;;;;;;;;;;;;;;
  
 print,'star'
;print,textoidl('\star')

;print,textoidl('\alpha')
print,textoidl('\Delta')


  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  plot, findgen(10), /nodata, xrange = [xmin,xmax], yrange = [ymin, ymax], $
        xstyle = 1, ystyle = 1, xtitle = 'm!Dsat!N/m!Dcentral!N', ytitle = 'F!DMcold!N'
  ;BESTFIT
  ;sel=where(Alog10(v_max) gt xmin)

  Polyfill, [mass_ratio,Reverse(mass_ratio)], $
             [sf_burst_min,Reverse(sf_burst_max)], Color=fsc_color("sky blue"), NOCLIP = 0
  oplot,mass_ratio,  sf_burst_bestfit, linestyle=0, color=200, thick=6 

  ;OLD MODELS
  oplot,mass_ratio, sf_burst_model1, linestyle=1, color=200, thick=6
  oplot,mass_ratio, sf_burst_model2, linestyle=2, color=200, thick=6
  


;LABELS
  
  plot_label, xlog=0, ylog=0, type='label', label='This Work - Best Fit' , xmin, xmax, ymin, ymax, $
              x_percentage=1.2, x2_percentage=0., y_percentage=9.2, charsize=1.3
  plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.3, x2_percentage=1.0, y_percentage=9.4, $
              linestyle=0, color=200, linethick=8
      
  plot_label, xlog=0, ylog=0, type='label', label='This Work!M+!X2!7r!X' , xmin, xmax, ymin, ymax, $
              x_percentage=1.2, x2_percentage=0., y_percentage=8.5, charsize=1.3
  plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.3, x2_percentage=1.0, y_percentage=8.7, $
              linestyle=0, color=fsc_color("sky blue"), linethick=20





  plot_label, xlog=0, ylog=0, type='label', label='Henriques2013' , xmin, xmax, ymin, ymax, $
              x_percentage=1.2, x2_percentage=0., y_percentage=7.7, charsize=1.3
  plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.3, x2_percentage=1.0, y_percentage=7.85, $
              linestyle=2, color=200, linethick=8

  plot_label, xlog=0, ylog=0, type='label', label='Guo2013' , xmin, xmax, ymin, ymax, $
              x_percentage=1.2, x2_percentage=0., y_percentage=7.0, charsize=1.3
  plot_label, xlog=0, ylog=0, type='line', xmin, xmax, ymin, ymax, $
              x_percentage=0.3, x2_percentage=1.0, y_percentage=7.15, $
              linestyle=1, color=200, linethick=8


  plot_label, xlog=0, ylog=0, type='label', label='!7a!X!DSF,brust!N(m!Dsat!N/m!Dcentral!N)^!7b!X!DSF,brust' , xmin, xmax, ymin, ymax, $
              x_percentage=0.3, x2_percentage=0., y_percentage=5.5, charsize=1.3

  plot_axis, xmin, xmax, ymin, ymax
 

device,/close_file











end



