;@'/afs/mpa/data/robyates/Eclipse64/workspace/robyates2/awk/idl/LGalaxy'
;@'/afs/mpa/data/robyates/Eclipse64/workspace/robyates2/awk/idl/LGalaxy_classic' ;w/o YIELDS LGalaxy structure
;@'/afs/mpa/data/robyates/Eclipse64/workspace/robyates2/awk/idl/LGalaxy_original' ;w/o YIELDS, SFH and METALS LGalaxy structure

;@'/afs/mpa/data/robyates/Eclipse64/workspace/robyates2/awk/idl/LGalaxy_allElements' ;New structure, with MAINELEMENTS off
;@'/afs/mpa/data/robyates/Eclipse64/workspace/robyates2/awk/idl/LGalaxy_new' ;New structure (incorporating Bruno's changes (eg:dists from type 2 to central), with MAINELEMENTS on
;@'/afs/mpa/data/robyates/Eclipse64/workspace/peter_metals_volker/awk/idl/LGalaxy_new' ;New structure for all changes up to 15-05-12, with MAINELEMENTS on
@'/afs/mpa/data/robyates/Eclipse64/workspace/peter_metals_volker/awk/idl/LGalaxy_allElements' ;New structure for all changes up to 15-05-12, with MAINELEMENTS off

pro init_hubbletype

common  share1, ma_table, T_table, d2y


   T_table=indgen(100)/99.0 * 14 - 5

   ma_table= 0.324*(T_table+5) - 0.054*(T_table+5)^2 +0.0047*(T_table+5)^3

   d2y=spl_init(ma_table, T_table)

end

function get_hubbletype, madiff

common  share1, ma_table, T_table, d2y


   T=spl_interp(ma_table, T_table , d2y, madiff)

   ind=where(t gt 9)
   
   if ind(0) ne -1 then begin
       T(ind)=9
   endif

   return, T

end

;---------------------------------------------------------------------
function uberround,input,place,string=string
  compile_opt idl2

  if n_params() lt 2 then begin
     print,' '
     print,'  ERROR: MUST INPUT SCALAR OR ARRAY'
     print,'  AND THE DECIMAL PLACE THAT THE ARRAY'
     print,'  SHOULD BE ROUNDED TO.'
     print,' '
     res=-99
  endif else begin
     if place le 0 then begin
        print,' '
        print,'  ERROR: THE DECIMAL PLACE THAT THE'
        print,'  INPUT ARRAY SHOULD BE ROUNDED TO'
        print,'  MUST BE AN INTEGER GREATER THAN 0.'
        print,' '
        res=-99
     endif
     arr=strtrim(string(double(input)),2)
     res=dblarr(n_elements(arr))
     for i=0, n_elements(arr)-1 do begin
        dec=strpos(arr,'.')
        tmp=double(strmid(arr[i],dec+place,1))
        if tmp ge 0 and tmp le 4 then begin
           res[i]=double(strmid(strtrim(string(arr[i]),2),0,dec+place+1))+10.^(-place)
        endif else begin
           res[i]=double(strmid(strtrim(string(arr[i]),2),0,dec+place+1))
        endelse
     endfor
  endelse
  if keyword_set(string) then begin
     res=strtrim(string(res),2)
     for i=0,n_elements(res)-1 do begin
        dec=strpos(res[i],'.')
        res[i]=strmid(res[i],0,dec+place+1)
     endfor
  endif
  return,res

end
;-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
;-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

pro robs_plots_mk2

;A COPY OF: plot_data_idlde_mainelements_mk15

  ;BaseDir = '/afs/mpa/data/robyates/Eclipse64/workspace/robyates2/'
  ;BaseDir = '/afs/mpa/data/robyates/Eclipse64/workspace/robyates2_backup/'
  BaseDir = '/afs/mpa/data/robyates/Eclipse64/workspace/peter_metals_volker/'
  FileName  = 'SA'
  Location  = 'output/'
  Extn    = 'z0.00' ;'z0.00'
  
  NN = 1
  
  FirstFile = 5 ;102
  LastFile  = FirstFile + NN - 1
  
  FirstFileNr   = strcompress(FirstFile, /remove_all)  ;'10'
  LastFileNr    = strcompress(LastFile, /remove_all)   ;'150'

  MilliMil     = 0
  MSI          = 1
  MSII         = 0
  SDSSmags     = 1
  ignoremain   = 0
  loadselected = 0
  fromhome     = 0

  ;writegalaxies = 0
  ;writefilename = BaseDir+Location+'/'+FileName+'.'+Location+'.1.'+Extn+'.xyz' 
             ;ie: /afs/mpa/data/robyates/L-Galaxies_Qi/Results/SA.Results.1.z0.00.xyz

;NEW GALAXY STRUCTURE (0) OR OLD GALAXY STRUCTURE (1):
OLDGALSTRUCT = 0

;MAINELEMENTS on (0) or MAINELEMENTS off (1):
MAINELEMENTS = 1

;NEWSESTSTRUCT on (1): When using the output structure suitable after Bruno's changes.
NEWSESTSTRUCT = 1

;Define the greek letter 'alpha':
alphasymbol = '!4' + String("141B) + '!X'

;Read in Holmberg et al. (2009) MW [Fe/H] data:
Holmberg = MRDFITS('/afs/mpa/data/robyates/Chemical_Enrichment/'+'Holmberg_FeH_MWdisk_data.fit', 1)

;; secondary  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  duston = 1                    ; dust model, 0=none
  scale = 1.0                   ; to scale masses etc
  type_pivot = 0.8              ; split type by B-V colour

  LFbinwidth = 0.25             ; for luminosity function
  BARbinwidth = 0.1             ; for baryonic mass function
  
  Hubble_h = 0.73

  if MSI eq 1 then begin
    BoxSize  = 500.0
    MaxTreeFiles = 512
    MaxSnaps = 64
    PartMass = 0.0860657
  endif else if MSII eq 1 then begin
    BoxSize  = 100.0
    MaxTreeFiles = 512
    MaxSnaps = 67
    PartMass = 0.000688526
  endif else if MilliMil eq 1 then begin
    BoxSize  = 62.5
    MaxTreeFiles = 8
    MaxSnaps = 64
    PartMass = 0.0860657
  endif

  
  if fromhome eq 1 then begin   ; for faster plotting of points on slow machines
    dilute = 0.005
  endif else begin
    dilute = 0.1
  endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  DirName  = BaseDir+Location
  ModelName = DirName+FileName+'_'+Extn ;that is: /afs/mpa/data/robyates/Eclipse64/workspace/robyates2/output/SA_z0.00
  theta = findgen(30)/29.*360.
  usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill
  loadct, 39

  scalemass = scale^3.
  scaleenergy = scale^5.
  volume = (BoxSize^3.0) * 1.0*(Lastfile-Firstfile+1) / MaxTreeFiles 

  col0 = 255 & col1 = 60 & col2 = 230 & col3 = 60
  tck = 2 & tck2 = 3
  hat = 3
  ls1 = 1 & ls2 = 4

    dilute = 0.1*(16.0/(LastFile-FirstFile));0.03;
    !p.multi = 0
    !x.thick = 2
    !y.thick = 2
    !p.charthick = 2
    tck = 5
    col0 = 0 & hat = 100
    csiz = 2.0 & csiz2 = 1.5
    ssiz = 0.5 & ssiz2 = 0.5 & ssiz3 = 1.3 ;ssiz=0.5 for plotting is best
    set_plot, 'PS'
    device, color = 1, filename = BaseDir+'/figures/'+FileName+'.'+Extn+'.'+FirstFileNr+'.'+LastFileNr+'.ps'$ ;'_DiffYield'+
          , xsize = 20, ysize = 18,  yoffset = 5.0,  xoffset = 0.7, /portrait ;ysize = 18,  yoffset = 5.0
          
;Get LGalaxy structure from LGalaxy.pro
IF (OLDGALSTRUCT eq 1) THEN BEGIN
Gstruct={LGalaxy_classic}
ENDIF ELSE BEGIN
  IF (MAINELEMENTS eq 0) THEN BEGIN
    IF (NEWSESTSTRUCT eq 1) THEN BEGIN
    Gstruct={LGalaxy_new}
    ENDIF ELSE BEGIN
    Gstruct={LGalaxy}
    ENDELSE
  ENDIF ELSE BEGIN
  Gstruct={LGalaxy_allElements} ;Needs to be updated for any Merged-with-Bruno files
  ENDELSE
ENDELSE        

 ; Read the SA output files
   print, ' Now reading the input files...'
   TotNTrees = 0L & TotNGals = 0L
   close, 1
   for fnr = FirstFile, LastFile do begin
   IF (OLDGALSTRUCT eq 1) THEN BEGIN
       fname = strcompress(ModelName + '_' + string(fnr) + '_withoutYIELDS', /remove_all) ;   + '_withoutYIELDS'
   ENDIF ELSE BEGIN
       fname = strcompress(ModelName + '_' + string(fnr) + '_GaussianDTD_SNeII_7-40_Chieffi_A0_09', /remove_all) ;   + '_ManDTD_SNeII_7-40_Chieffi_A0_05'     + '_GaussianDTD_SNeII_7-40_Chieffi_A0_09'
   ENDELSE
       openr, 1, fname;, /swap_endian
       if  (FileName eq 'SA_galtree') then begin
          one = 0L & readu,1,one
          print, one
          nbytes = 0L & readu, 1, nbytes
          ngals=0L & readu, 1, ngals
          TotNGals =  TotNGals  + ngals
       endif else begin
          Ntrees = 0L     & readu, 1, Ntrees
          NtotGals = 0L   & readu, 1, NtotGals
          TotNTrees = TotNTrees + Ntrees
          TotNGals =  TotNGals  + NtotGals
       endelse
       close, 1
   endfor
   print, " Total number of groups = ", ToTNTrees
   print, " Total number of galaxies = ", TotNGals
   Ngals = lonarr(TotNTrees)
   G = replicate(Gstruct, TotNGals)
   offset = 0L
   off = 0L
   for fnr = FirstFile, LastFile do begin
       print, ' file:', fnr
   IF (OLDGALSTRUCT eq 1) THEN BEGIN
       fname = strcompress(ModelName + '_' + string(fnr) + '_withoutYIELDS', /remove_all) ;  + '_withoutYIELDS'
   ENDIF ELSE BEGIN
       fname = strcompress(ModelName + '_' + string(fnr) + '_GaussianDTD_SNeII_7-40_Chieffi_A0_09', /remove_all) ;   + '_ManDTD_SNeII_7-40_Chieffi_A0_05'     + '_GaussianDTD_SNeII_7-40_Chieffi_A0_09'
       
   ENDELSE
       openr, 1, fname;, /swap_endian
       if  (FileName eq 'SA_galtree') then begin;TREE OUTPUT
          one = 0L & readu,1,one
          nbytes = 0L & readu, 1, nbytes
          ngals=0L & readu, 1, ngals & print, ' NGals:', ngals
          nskip=nbytes/4-3
          ib=fltarr(nskip)
          readu,1,ib
          GG = replicate(Gstruct, ngals)
          readu, 1, GG
          G(offset:offset+ngals-1) = GG(*)
          offset = offset + ngals
       endif else begin;SNAP OUTPUT
          Ntrees = 0L     & readu, 1, Ntrees   & print, ' Ntrees:', Ntrees
          NtotGals = 0L   & readu, 1, NtotGals & print, ' NtotGals:', NtotGals
          GalsPerTree = lonarr(Ntrees)
          readu, 1, GalsPerTree
          Ngals(off:off+Ntrees-1) = GalsPerTree(*)
          GG = replicate(Gstruct, NtotGals)
          readu, 1, GG
          G(offset:offset+NtotGals-1) = GG(*)
          offset = offset + NtotGals
          off = off + Ntrees
       endelse
       close, 1
   endfor

  print
  print, "total galaxies considered:", TotNGals
  print  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SCALE DATA:

  G.CentralMvir = G.CentralMvir*scalemass/Hubble_h

  G.Pos  = G.Pos*scale
  G.Vel  = G.Vel*scale
  G.Mvir = G.Mvir*scalemass/Hubble_h
  G.Rvir = G.Rvir*scale
  G.Vvir = G.Vvir*scale
  G.Vmax = G.Vmax*scale

  G.ColdGas       = G.ColdGas*scalemass/Hubble_h
  ;IF (OLDGALSTRUCT eq 1) THEN BEGIN G.StellarMass   = G.StellarMass*scalemass/Hubble_h
  ;ENDIF ELSE
  G.DiskMass   = G.DiskMass*scalemass/Hubble_h
  G.BulgeMass     = G.BulgeMass*scalemass/Hubble_h
  G.HotGas        = G.HotGas*scalemass/Hubble_h
  G.EjectedMass   = G.EjectedMass*scalemass/Hubble_h
  G.BlackHoleMass = G.BlackHoleMass*scalemass/Hubble_h

  G.MetalsColdGas       = G.MetalsColdGas*scalemass/Hubble_h
  ;IF (OLDGALSTRUCT eq 1) THEN BEGIN  G.MetalsStellarMass   = G.MetalsStellarMass*scalemass/Hubble_h
  ;ENDIF ELSE
  G.MetalsDiskMass   = G.MetalsDiskMass*scalemass/Hubble_h
  G.MetalsBulgeMass     = G.MetalsBulgeMass*scalemass/Hubble_h
  G.MetalsHotGas        = G.MetalsHotGas*scalemass/Hubble_h
  G.MetalsEjectedMass   = G.MetalsEjectedMass*scalemass/Hubble_h

  G.Sfr                = G.Sfr*scalemass
  G.SfrBulge           = G.SfrBulge*scalemass
  ;G.XrayLum            = G.XrayLum + alog10(scaleenergy)
  ;G.DiskRadius         = G.DiskRadius*scale
  ;G.CoolingRadius      = G.CoolingRadius*scale
  
  G.CentralMvir = G.CentralMvir*scalemass/Hubble_h*1e10

  if (duston eq 0) then begin
    G.Mag  = G.Mag - 2.5*alog10(scaleenergy)
  endif else begin
    G.Mag  = G.MagDust - 2.5*alog10(scaleenergy)
  endelse
  G.MagBulge = G.MagBulge - 2.5*alog10(scaleenergy)

  Ratio = G.MagBulge - G.Mag
  init_hubbletype
  T = get_hubbletype(Ratio)

  zeroSfr = G.Sfr


if SDSSmags ne 1 then begin
  Mag_B =  G.Mag[0]
  Mag_V =  G.Mag[1]
  Mag_R =  G.Mag[2]
  Mag_I =  G.Mag[3]
  Mag_K =  G.Mag[4]
  Mag_Bb =  G.MagBulge[0]
  Mag_Vb =  G.MagBulge[1]
  Mag_Rb =  G.MagBulge[2]
  Mag_Ib =  G.MagBulge[3]
  Mag_Kb =  G.MagBulge[4]
endif

if SDSSmags eq 1 then begin   ; strange order - fix this!
  Mag_B =  G.Mag[1]  ; g
  Mag_V =  G.Mag[2]  ; r
  Mag_R =  G.Mag[0]  ; u
  Mag_I =  G.Mag[3]  ; i
  Mag_K =  G.Mag[4]  ; z
  Mag_Bb =  G.MagBulge[1]  ; g
  Mag_Vb =  G.MagBulge[2]  ; r
  Mag_Rb =  G.MagBulge[0]  ; u
  Mag_Ib =  G.MagBulge[3]  ; i
  Mag_Kb =  G.MagBulge[4]  ; z
endif


  blu = where (Mag_B lt -19.0 and Mag_B gt -20.0 and Mag_B-Mag_V lt type_pivot)
  red = where (Mag_B lt -19.0 and Mag_B gt -20.0 and Mag_B-Mag_V gt type_pivot)
  print, 'Ratio of red to blue at M*:', 1.0*n_elements(red)/n_elements(blu)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;BEGIN PLOTTING:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if SDSSmags eq 1 then begin
IF (OLDGALSTRUCT eq 0) THEN BEGIN 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;; plot SDSS LFs ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ymax = 0.1 & ymin = 3.0*ymax/1000000.0

;      plot, findgen(10), /nodata, charsize = csiz, $
;            xrange = [-26.0, -15.0], yrange = [ymin, ymax], /ylog, xstyle = 1, ystyle = 1, $ 
;            xtitle = 'M - 5log!D10!Nh', ytitle = '!4U!3 (h!U3!NMpc!U-3!Nmag!U-1!N)'


; g
      M = -(findgen(3000)/100)
      Mstar = -20.04
      alpha = -1.26
      phistar = 0.0206
      xx = 10.0^(0.4*(Mstar-M))
      SF = (2./5.)*alog(10.) * phistar * (xx^(alpha+1.)) * (exp(-xx))
      ;oplot, M, SF, linestyle = ls1, thick = tck, color = col3

      mag = Mag_B
      mi = floor(min(mag))-2 & ma = floor(max(mag))+2
      ;mi = -26 & ma = -15
      bins = LFbinwidth
      NB = (ma-mi)/bins
      counts = histogram(mag, min = mi, max = ma, binsize = bins)
      xaxeshisto = (indgen(n_elements(counts))+0.5)/float(NB)*(ma-mi)+mi
      ;oplot, xaxeshisto-5*alog10(Hubble_h), counts/volume/bins, thick = tck;, psym = 10


; r
      M = -(findgen(3000)/100)
      Mstar = -20.83
      alpha = -1.20
      phistar = 0.0146
      xx = 10.0^(0.4*(Mstar-M))
      SF = (2./5.)*alog(10.) * phistar * (xx^(alpha+1.)) * (exp(-xx))
      ;oplot, M, SF, linestyle = ls1, thick = tck, color = col3

      mag = Mag_V
      mi = floor(min(mag))-2 & ma = floor(max(mag))+2
      ;mi = -26 & ma = -15
      bins = LFbinwidth
      NB = (ma-mi)/bins
      counts = histogram(mag, min = mi, max = ma, binsize = bins)
      xaxeshisto = (indgen(n_elements(counts))+0.5)/float(NB)*(ma-mi)+mi
      ;oplot, xaxeshisto-5*alog10(Hubble_h), counts/volume/bins, thick = tck;, psym = 10


; i
      M = -(findgen(3000)/100)
      Mstar = -21.26
      alpha = -1.25
      phistar = 0.0128
      xx = 10.0^(0.4*(Mstar-M))
      SF = (2./5.)*alog(10.) * phistar * (xx^(alpha+1.)) * (exp(-xx))
      ;oplot, M, SF, linestyle = ls1, thick = tck, color = col3

      mag = Mag_I
      mi = floor(min(mag))-2 & ma = floor(max(mag))+2
      ;mi = -26 & ma = -15
      bins = LFbinwidth
      NB = (ma-mi)/bins
      counts = histogram(mag, min = mi, max = ma, binsize = bins)
      xaxeshisto = (indgen(n_elements(counts))+0.5)/float(NB)*(ma-mi)+mi
      ;oplot, xaxeshisto-5*alog10(Hubble_h), counts/volume/bins, thick = tck;, psym = 10


; z
      M = -(findgen(3000)/100)
      Mstar = -21.55
      alpha = -1.24
      phistar = 0.0127
      xx = 10.0^(0.4*(Mstar-M))
      SF = (2./5.)*alog(10.) * phistar * (xx^(alpha+1.)) * (exp(-xx))
      ;oplot, M, SF, linestyle = ls1, thick = tck, color = col3

      mag = Mag_K
      mi = floor(min(mag))-2 & ma = floor(max(mag))+2
      ;mi = -26 & ma = -15
      bins = LFbinwidth
      NB = (ma-mi)/bins
      counts = histogram(mag, min = mi, max = ma, binsize = bins)
      xaxeshisto = (indgen(n_elements(counts))+0.5)/float(NB)*(ma-mi)+mi
      ;oplot, xaxeshisto-5*alog10(Hubble_h), counts/volume/bins, thick = tck;, psym = 10

ENDIF ;IF (OLDGALSTRUCT eq 0) THEN BEGIN 
endif

if SDSSmags ne 1 then begin
IF (OLDGALSTRUCT eq 0) THEN BEGIN 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;; plot LF k band ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      mag = Mag_K

      ;mi = floor(min(mag))-2 & ma = floor(max(mag))+2
      mi = -30.0 & ma = -5.0
      bins = LFbinwidth
      NB = (ma-mi)/bins

      counts = histogram(mag, min = mi, max = ma, binsize = bins)
      xaxeshisto = (indgen(n_elements(counts))+0.5)/float(NB)*(ma-mi)+mi
      ymax = 0.1 & ymin = 3.0*ymax/1000000.0


;      plot, findgen(10), /nodata, charsize = csiz, $
;            xrange = [-27.5, -18.0], yrange = [ymin, ymax], /ylog, xstyle = 1, ystyle = 1, $ 
;            xtitle = 'M!DK!N - 5log!D10!Nh', ytitle = '!4U!3 (h!U3!NMpc!U-3!Nmag!U-1!N)'

      legend, ['Cole et al. (2001)'], box = 0, charsize = csiz2,  psym = [8], symsize = [ssiz3], color = [col3], /top ;, textcolor = [col1]
      legend, ['Kochanek et al. (2001)', 'Huang et al. (2003)'], box = 0, charsize = csiz2,  linestyle = [ls1, ls2], color = [col3, col3], thick = [tck, tck], /bottom, /right


; output the LF points

;       print
;       print, "k band"
;       for i = 0, n_elements(counts)-1 do print, xaxeshisto(i)-5*alog10(h), counts(i)/volume/bins
;       print


; overplot various Schechter functions from the literature

;       M = -(findgen(3000)/100)
;       Mstar = -23.44            ; Cole et al. 2001
;       alpha = -0.96
;       phistar = 0.0108
;       xx = 10.0^(0.4*(Mstar-M))
;       SF = (2./5.)*alog(10.) * phistar * (xx^(alpha+1.)) * (exp(-xx))
;       oplot, M, SF, linestyle = 1, thick = tck, color = col3

      M = -(findgen(3000)/100)
      Mstar = -23.39            ; Kochanek et al 2001
      alpha = -1.09
      phistar = 0.0116
      xx = 10.0^(0.4*(Mstar-M))
      SF = (2./5.)*alog(10.) * phistar * (xx^(alpha+1.)) * (exp(-xx))
      ;oplot, M, SF, linestyle = ls1, thick = tck, color = col3

      M = -(findgen(3000)/100)
      Mstar = -23.70            ; Huang et al 2003
      alpha = -1.37
      phistar = 0.013
      xx = 10.0^(0.4*(Mstar-M))
      SF = (2./5.)*alog(10.) * phistar * (xx^(alpha+1.)) * (exp(-xx))
      ;oplot, M, SF, linestyle = ls2, thick = tck, color = col3


; overplot the Cole et al. K band 2dFGRS LF

      fullK = [3.1315561E-03, 8.2625253E-03, 0.0000000E+00, 4.6483092E-03, 5.7576019E-03, 9.1649834E-03, 1.1232893E-02, $
               1.0536440E-02, 8.5763102E-03, 8.8181989E-03, 6.9448259E-03, 6.0896124E-03, 9.2596142E-03, 6.9631678E-03, $
               7.2867479E-03, 6.9923755E-03, 5.9844730E-03, 5.9305103E-03, 5.3865365E-03, 5.8525647E-03, 5.2373926E-03, $
               4.9635037E-03, 4.1801766E-03, 2.7171015E-03, 1.8800517E-03, 1.2136410E-03, 6.5419916E-04, 3.4594961E-04, $
               1.4771589E-04, 5.5521199E-05, 2.1283222E-05, 9.4211919E-06, 1.0871951E-06, 2.7923562E-07]

      fullErrK = [3.6377162E-03, 6.6833422E-03, 1.0000000E-10, 4.0996978E-03, 4.3155681E-03, 5.6722397E-03, 6.4211683E-03, $
                  5.7120644E-03, 4.6346937E-03, 3.8633577E-03, 2.4383855E-03, 1.6279118E-03, 1.6941463E-03, 1.1781409E-03, $
                  9.7785855E-04, 7.9027453E-04, 6.0649612E-04, 5.1598746E-04, 4.2267537E-04, 3.7395130E-04, 2.8177485E-04, $
                  2.1805518E-04, 1.6829016E-04, 1.1366483E-04, 8.1871600E-05, 5.7472309E-05, 3.6554517E-05, 2.3141622E-05, $
                  1.2801432E-05, 6.5092854E-06, 3.3540452E-06, 1.9559407E-06, 5.6035748E-07, 2.8150106E-07]

      fullMagsK = [-18.00000, -18.25000, -18.50000, -18.75000, -19.00000, -19.25000, -19.50000, -19.75000, -20.00000, $
                   -20.25000, -20.50000, -20.75000, -21.00000, -21.25000, -21.50000, -21.75000, -22.00000, -22.25000, $
                   -22.50000, -22.75000, -23.00000, -23.25000, -23.50000, -23.75000, -24.00000, -24.25000, -24.50000, $
                   -24.75000, -25.00000, -25.25000, -25.50000, -25.75000, -26.00000, -26.25000]

      oploterror, fullMagsK, fullK, fullErrK, hatlength = hat, color = col3, errcolor = col3, psym = 8, symsize = ssiz3, errthick = 2


; now plot the SA model histogram

      ;oplot, xaxeshisto-5*alog10(Hubble_h), counts/volume/bins, thick = tck;, psym = 10




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;; plot LF bJ band ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      mag = Mag_B - 0.28*(Mag_B-Mag_V) 

      ;mi = floor(min(mag))-2 & ma = floor(max(mag))+2
      mi = -30.0 & ma = -5.0     
      bins = LFbinwidth
      NB = (ma-mi)/bins


;      plot, findgen(10), /nodata, charsize = csiz, $
;            xrange = [-24.5, -15.0], yrange = [ymin, ymax], /ylog, xstyle = 1, ystyle = 1, $ 
;            xtitle = 'M!DbJ!N - 5log!D10!Nh', ytitle = '!4U!3 (h!U3!NMpc!U-3!Nmag!U-1!N)'
      
      legend, ['Norberg et al. (2002)'], box = 0, charsize = csiz2,  psym = [8], symsize = [ssiz3], color = [col3], /top ;, textcolor = [col1]



; overplot various Schechter functions from the literature
      M = -(findgen(3000)/100)

; all  - norberg et al. 2001
;       Mstar = -19.74
;       alpha = -1.16
;       phistar = 0.0155
;       xx = 10.0^(0.4*(Mstar-M))
;       SF = (2./5.)*alog(10.) * phistar * (xx^(alpha+1.)) * (exp(-xx))
;       oplot, M, SF, linestyle = 0, thick = 1, color = col3



; overplot the Norberg et al. (2002) bJ band 2dFGRS LF

      fullbJ = [0.0000000E+00, 4.9295418E-02, 4.8414569E-02, 6.1575811E-02, 4.0561352E-02, 4.4041574E-02, 3.4943096E-02, $
                3.4614738E-02, 3.2133184E-02, 2.9635418E-02, 2.3947729E-02, 2.4405226E-02, 2.3334850E-02, 2.1415474E-02, $
                1.9061793E-02, 1.8999882E-02, 1.6375385E-02, 1.5369940E-02, 1.4364687E-02, 1.3468309E-02, 1.2573769E-02, $
                1.1145771E-02, 9.4299037E-03, 8.0003440E-03, 6.1532548E-03, 4.5839567E-03, 3.0883602E-03, 1.9342416E-03, $
                1.0532817E-03, 5.1049568E-04, 2.1174403E-04, 7.1619586E-05, 2.1268283E-05, 5.9436675E-06, 1.8153961E-06, $
                5.1464536E-07, 3.0028801E-07, 2.3413641E-07]

      fullErrbJ = [1.0000000E-10, 2.4690846E-02, 1.2388939E-02, 9.7498214E-03, 5.7149841E-03, 4.5108106E-03, 2.9827626E-03, $
                   2.3375382E-03, 1.8226152E-03, 1.4552969E-03, 1.0482093E-03, 8.9600182E-04, 7.4631697E-04, 6.2729424E-04, $
                   4.9186783E-04, 4.0296704E-04, 2.9608142E-04, 2.3880773E-04, 1.9196754E-04, 1.5552285E-04, 1.2866473E-04, $
                   9.9682424E-05, 7.6998745E-05, 5.9886646E-05, 4.4532266E-05, 3.3105953E-05, 2.3582392E-05, 1.6317146E-05, $
                   1.0632872E-05, 6.6587650E-06, 4.0021373E-06, 2.2571021E-06, 1.2226549E-06, 6.4833682E-07, 3.6088997E-07, $
                   1.9372916E-07, 1.4984870E-07, 1.3478319E-07]

      fullMagsbJ = [-13.00000, -13.27500, -13.55000, -13.82500, -14.10000, -14.37500, -14.65000, -14.92500, -15.20000, -15.47500, $
                    -15.75000, -16.02500, -16.30000, -16.57500, -16.85000, -17.12500, -17.40000, -17.67500, -17.95000, -18.22500, $
                    -18.50000, -18.77500, -19.05000, -19.32500, -19.60000, -19.87500, -20.15000, -20.42500, -20.70000, -20.97500, $
                    -21.25000, -21.52500, -21.80000, -22.07500, -22.35000, -22.62500, -22.90000] 

      oploterror, fullMagsbJ, fullbJ, fullErrbJ, hatlength = hat, color = col3, errcolor = col3, psym = 8, symsize = ssiz3, errthick = 2


; calculate and plot the SA model histogram
      counts = histogram(mag, min = mi, max = ma, binsize = bins)
      xaxeshisto = (indgen(n_elements(counts))+0.5)/float(NB)*(ma-mi)+mi
      ;oplot, xaxeshisto-5*alog10(Hubble_h), counts/volume/bins, thick = tck;, psym = 10


; output the LF points

;       print
;       print, "bJ band"
;       for i = 0, n_elements(counts)-1 do print, xaxeshisto(i)-5*alog10(h), counts(i)/volume/bins
;       print




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;; plot type dependent LF bJ band ;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      mag = Mag_B - 0.28*(Mag_B-Mag_V) 

      ;mi = floor(min(mag))-2 & ma = floor(max(mag))+2
      mi = -30.0 & ma = -5.0     
      bins = LFbinwidth
      NB = (ma-mi)/bins


;      plot, findgen(10), /nodata, charsize = csiz, $
;            xrange = [-24.5, -15.0], yrange = [ymin, ymax], /ylog, xstyle = 1, ystyle = 1, $ 
;            xtitle = 'M!DbJ!N - 5log!D10!Nh', ytitle = '!4U!3 (h!U3!NMpc!U-3!Nmag!U-1!N)'
      
      legend, ['Madgwick et al. (2002)', 'early type', 'late type'], box = 0, charsize = csiz2,  $
              linestyle = [-99, 0, 0], psym = [0, 5, 6], color = [0, col2, col3], symsize = [0.0, ssiz3*0.4, ssiz3*0.4], thick = [0, 6, 6], /top



; overplot various Schechter functions from the literature
      M = -(findgen(3000)/100)

;early type - madgwick et al. 2002
;       Mstar = -19.58
;       alpha = -0.54
;       phistar = 0.0099
;       xx = 10.0^(0.4*(Mstar-M))
;       SF = (2./5.)*alog(10.) * phistar * (xx^(alpha+1.)) * (exp(-xx))
;       oplot, M, SF, linestyle = ls1, thick = tck, color = col2   

      TypeMags = [-22.7760, -22.5510, -22.3270, -22.1020, -21.8780, -21.6530, -21.4290, -21.2040, -20.9800, -20.7550, -20.5310, $
                  -20.3060, -20.0820, -19.8570, -19.6330, -19.4080, -19.1840, -18.9590, -18.7350, -18.5100, -18.2860, -18.0610, $
                  -17.8370, -17.6120, -17.3880, -17.1630, -16.9390, -16.7140, -16.4900, -16.2650, -16.0410, -15.8160, -15.5920, $
                  -15.3670, -15.1430, -14.9180, -14.6940, -14.4690, -14.2450, -14.0200, -13.7960, -13.5710, -13.3470]

      EarlyType = [1.66134e-09, 1.62423e-06, 7.97733e-06, 1.10532e-05, 2.43034e-05, 4.91651e-05, 0.000119044, 0.000269121, $
                   0.000436383, 0.000691079, 0.00107770, 0.00155109, 0.00203192, 0.00252114, 0.00313738, 0.00366140, 0.00368101, $
                   0.00416359, 0.00430455, 0.00417049, 0.00403680, 0.00369775, 0.00338008, 0.00321586, 0.00320057, 0.00286733, $
                   0.00257553, 0.00261365, 0.00254611, 0.00352320, 0.00277701, 0.00300647, 0.00273563, 0.00327098, 0.00165774, $
                   0.00431527, 0.00373609, 0.000896536, 0.00188694, 0.00375827, 0.00687318, 0.0109546, 0.0319976]

      EarlyTypeErr = [1.17474e-09,  1.14850e-06,  2.52267e-06,  2.95421e-06,  4.36525e-06,  6.19394e-06, $
                      9.62377e-06,  1.44471e-05,  1.83723e-05, 2.31007e-05,  2.88160e-05,  3.45584e-05, $
                      3.95399e-05,  4.42484e-05,  5.09783e-05,  5.91496e-05,  6.54389e-05,  7.80904e-05, $
                      9.13496e-05,  0.000104027,  0.000119408,  0.000131783,  0.000143361,  0.000159622, $
                      0.000180360,  0.000202806,  0.000232274,  0.000275526,  0.000318254,  0.000418172, $
                      0.000413994,  0.000494256,  0.000536518,  0.000682030,  0.000586125,   0.00111420, $
                      0.00132087,  0.000896521,   0.00188668,   0.00375453,   0.00683054,    0.0106302, 0.00865398]


      oploterror, TypeMags, EarlyType, EarlyTypeErr, hatlength = hat, color = col2, errcolor = col2, psym = 5, symsize = ssiz3*0.4, errthick = 2, thick = 6



;late type - madgwick et al. 2002
 ;      Mstar = -19.53
;       alpha = -0.99
;       phistar = 0.0072
;       xx = 10.0^(0.4*(Mstar-M))
;       SF = (2./5.)*alog(10.) * phistar * (xx^(alpha+1.)) * (exp(-xx))
;       Mstar = -19.17
;       alpha = -1.24
;       phistar = 0.0050
;       xx = 10.0^(0.4*(Mstar-M))
;       SF = SF + (2./5.)*alog(10.) * phistar * (xx^(alpha+1.)) * (exp(-xx))
;       Mstar = -19.15
;       alpha = -1.50
;       phistar = 0.0024
;       xx = 10.0^(0.4*(Mstar-M))
;       SF = SF + (2./5.)*alog(10.) * phistar * (xx^(alpha+1.)) * (exp(-xx))
;       oplot, M, SF, linestyle = ls1, thick = tck, color = col1

      LateType = [3.33333, 3.33333, 2.22222, 2.22222, 4.40561e-06, 1.58720e-05, 2.96734e-05, 8.90940e-05, 0.000188970, 0.000405951, $
                  0.000700606, 0.00119498, 0.00178015, 0.00255784, 0.00342355, 0.00467989, 0.00556459, 0.00676165, 0.00802202, $
                  0.00912797, 0.0101309, 0.0110775, 0.0118513, 0.0130438, 0.0146573, 0.0156096, 0.0157717, 0.0183020, 0.0184675, $
                  0.0221623, 0.0222427, 0.0211590, 0.0226055, 0.0270246, 0.0300475, 0.0317564, 0.0309249, 0.0349773, 0.0351621, $
                  0.0240171, 0.0492151, 0.0266096, 1.14692]

      LateTypeErr = [0.00000,  0.00000,  8.62398e-07,  1.11561e-06,  2.99235e-06,  5.13310e-06, $
                     6.62882e-06,  1.11406e-05,  1.76993e-05,  2.53900e-05,  3.40969e-05,  4.50930e-05, $
                     5.56902e-05,  6.80393e-05,  8.08717e-05,  0.000101552,  0.000123584,  0.000159331, $
                     0.000200303,  0.000248703,  0.000300701,  0.000356013,  0.000426286,  0.000520177, $
                     0.000634494,  0.000775759,  0.000933532,   0.00118531,   0.00137662,   0.00169484, $
                     0.00191843,   0.00216371,   0.00257518,   0.00330062,   0.00420505,   0.00516747, $
                     0.00624116,   0.00796132,    0.0100763,    0.0100451,    0.0188534,    0.0165111, 0.0114633]

      oploterror, TypeMags, LateType, LateTypeErr, hatlength = hat, color = col1, errcolor = col1, psym = 6, symsize = ssiz3*0.4, errthick = 2, thick = 6



;overplot the early and late type SA model results
      w1 = where(Mag_B-Mag_V lt type_pivot)
      w2 = where(Mag_B-Mag_V gt type_pivot)
; late types
      counts = histogram(mag(w1), min = mi, max = ma, binsize = bins)
      xaxeshisto = (indgen(n_elements(counts))+0.5)/float(NB)*(ma-mi)+mi
      ;oplot, xaxeshisto-5*alog10(Hubble_h), counts/volume/bins, thick = tck, color = col1;, psym = 10
; early types
      counts = histogram(mag(w2), min = mi, max = ma, binsize = bins)
      xaxeshisto = (indgen(n_elements(counts))+0.5)/float(NB)*(ma-mi)+mi
      ;oplot, xaxeshisto-5*alog10(Hubble_h), counts/volume/bins, thick = tck, color = col2;, psym = 10


endif
ENDIF ;IF (OLDGALSTRUCT eq 0) THEN BEGIN 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;; plot colour-mass diagram (contours) ;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;IF (OLDGALSTRUCT eq 1) THEN BEGIN
;stars = G.stellarMass
;ENDIF ELSE
stars = G.BulgeMass + G.DiskMass
g_r = G.mag[1]-G.mag[2]

hist2d, alog10(stars*1e10), g_r, hist, [8.0,12.0], [0.2,1.0], 40, 40

;print, g_r
;print, "------"
;print, hist

;print, n_elements(stars1), n_elements(g_r1), n_elements(hist)
contour, hist, findgen(40)/(40.0/(12.0-8.0)) + 8.0, findgen(40)/(40.0/(1.0-0.2)) + 0.2, nlevels = 10, $
         xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = 'g-r', thick = 6, c_color = fsc_color("black")
         
         legend, 'New model', box = 0, charsize = 1.25, linestyle = 0, color = fsc_color("black"), thick = 6, /bottom, /right

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

u_r = G.mag[0]-G.mag[2]

hist2d, alog10(stars*1e10), u_r, hist, [8.0,12.0], [0.6,3.0], 40, 40

;print, n_elements(stars1), n_elements(g_r1), n_elements(hist)
;contour, hist, findgen(40)/(40.0/(12.0-8.0)) + 8.0, findgen(40)/(40.0/(3.0+0.6)) - 0.6, nlevels = 10, $
;         xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = 'u-r'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

u_i = G.mag[0]-G.mag[3]

hist2d, alog10(stars*1e10), u_i, hist, [8.0,12.0], [0.5,3.5], 40, 40

;print, n_elements(stars1), n_elements(g_r1), n_elements(hist)
;contour, hist, findgen(40)/(40.0/(12.0-8.0)) + 8.0, findgen(40)/(40.0/(3.5+0.5)) - 0.5, nlevels = 10, $
;         xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = 'u-i'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

g_r_bulge = Mag_Bb-Mag_Vb

hist2d, alog10(stars*1e10), g_r_bulge, hist, [8.0,12.0], [0.0,1.0], 40, 40

;print, n_elements(stars1), n_elements(g_r1), n_elements(hist)
;contour, hist, findgen(40)/(40.0/(12.0-8.0)) + 8.0, findgen(40)/(40.0/(1.0+0.0)) - 0.0, nlevels = 8, $
;         xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = 'g-r (bulge)'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

u_r_bulge = Mag_Rb-Mag_Vb

hist2d, alog10(stars*1e10), u_r_bulge, hist, [8.0,12.0], [0.6,3.0], 40, 40

;print, n_elements(stars1), n_elements(g_r1), n_elements(hist)
;contour, hist, findgen(40)/(40.0/(12.0-8.0)) + 8.0, findgen(40)/(40.0/(3.0+0.6)) - 0.6, nlevels = 8, $
;         xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = 'u-r (bulge)'

;;;;;;;;;;;Gas mass vs. Stellar Mass:;;;;;;;;;;;;;;;;;;;;;;;;;
plot, findgen(10), /nodata, xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = 'log M!Dcold!N / M!D!9'+string(110B)+'!3!N', $
            charsize = csiz, xrange = [8.5, 11.5], yrange = [8.5, 11.5], xstyle = 1, ystyle = 1
            
            ;oplot, alog10(G.DiskMass*(1.0e10)), alog10(G.ColdGas*(1.0e10)), psym = 8, symsize = ssiz, color = fsc_color("grey")
            
            oplot, alog10(stars*1e10), alog10(G.ColdGas*(1.0e10)), psym = 8, symsize = ssiz, color = fsc_color("grey")
            ;w2check = WHERE(G.type eq 2)
            ;oplot, alog10(stars(w2check)*1e10), alog10(G(w2check).ColdGas*(1.0e10)), psym = 8, symsize = ssiz, color = fsc_color("red")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;; plot metallicity (coloured by SFR);;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      galtype = G.Type
      gas = G.ColdGas
 IF (OLDGALSTRUCT eq 1) THEN BEGIN  metals = G.MetalsColdGas[1]
 ENDIF ELSE BEGIN
      metals = G.MetalsColdGas[0]+G.MetalsColdGas[1]+G.MetalsColdGas[2]
      ;metalsScaleSNIa = (0.01*G.MetalsColdGas[0])+G.MetalsColdGas[1]+G.MetalsColdGas[2]
ENDELSE
      ;metals = G.MetalsColdGas[1]+G.MetalsColdGas[2] //Not including SNIa
      ;stars = G.StellarMass
      ;stars = G.BulgeMass + G.DiskMass

      plot, findgen(10), /nodata, xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = 'Z!Dcold!N', $
            charsize = csiz, xrange = [8.5, 11.5], yrange = [8.0, 9.7], xstyle = 1, ystyle = 1

      ;legend, ['Sb/c galaxies'], box = 0, charsize = csiz2, linestyle = [-99], /top

      ;ww = where(galtype eq 0 and gas/(stars+gas) gt 0.1); and randomu(2121, n_elements(galtype)) lt dilute/10.0)
      ww = where(galtype eq 0 OR galtype eq 1)
      if ww[0] ne -1 then begin      
        ;oplot, alog10(stars(ww)*1e10), alog10((metalsScaleSNIa(ww)/gas(ww))/0.0127)+9., psym = 8, symsize = ssiz      
      
;        w1 = WHERE(alog10(G(ww).sfr) lt -1.0)
;        oplot, alog10(stars(ww(w1))*1e10), alog10((metals(ww(w1))/gas(ww(w1)))/0.02)+9., psym = 8, symsize = ssiz, color = fsc_color("red") ;Wiersma's solar metallicity (0.0127)
;        w2 = WHERE(alog10(G(ww).sfr) ge -1.0 AND alog10(G(ww).sfr) le 0.5)
;        oplot, alog10(stars(ww(w2))*1e10), alog10((metals(ww(w2))/gas(ww(w2)))/0.02)+9., psym = 8, symsize = ssiz, color = fsc_color("green")
;        w3 = WHERE(alog10(G(ww).sfr) gt 0.5)
;        oplot, alog10(stars(ww(w3))*1e10), alog10((metals(ww(w3))/gas(ww(w3)))/0.02)+9., psym = 8, symsize = ssiz, color = fsc_color("blue")    
        w1 = WHERE(alog10(G(ww).sfr) lt -1.0)
        oplot, alog10(stars(ww(w1))*1e10), alog10((metals(ww(w1))/gas(ww(w1)))/0.0134)+9.0, psym = 8, symsize = ssiz, color = fsc_color("red") ;color = fsc_color("red") ;Wiersma's solar metallicity = 0.0127 (from Cloudy). Asplund's = 0.0134
        w2 = WHERE(alog10(G(ww).sfr) ge -1.0 AND alog10(G(ww).sfr) le 0.5)
        oplot, alog10(stars(ww(w2))*1e10), alog10((metals(ww(w2))/gas(ww(w2)))/0.0134)+9.0, psym = 8, symsize = ssiz, color = fsc_color("green") ;color = fsc_color("green")
        w3 = WHERE(alog10(G(ww).sfr) gt 0.5)
        oplot, alog10(stars(ww(w3))*1e10), alog10((metals(ww(w3))/gas(ww(w3)))/0.0134)+9.0, psym = 8, symsize = ssiz, color = fsc_color("blue") ;color = fsc_color("blue")   
      endif

      ;legend, ['Tremonti et al. (2003)'], box = 0, charsize = 1.25, linestyle = [0], color = [fsc_color("blue")], thick = [6], /top, /left
      legend, ['log(SFR) < -1.0','-1.0 < log(SFR) < 0.5','log(SFR) > 0.5'], box = 0, charsize = 1.25, psym = [8,8,8], color = [fsc_color("red"),fsc_color("green"),fsc_color("blue")], /top, /left

      ;Plot fit to all galaxies:
      w4 = WHERE(stars(ww) ne 0.0 AND metals(ww) ne 0.0 AND gas(ww) ne 0.0)
      ;print, "n_elements(w4) = ", n_elements(w4)
      ;print, "metals(ww(w4)) = ", metals(ww(w4))
      yyy = findgen(100)/(100.0/(12.0-8.5)) + 8.5
      FitGMZR = POLY_FIT(alog10(stars(ww(w4))*1e10), alog10((metals(ww(w4))/gas(ww(w4)))/0.0134)+8.69, 3, SIGMA=FitGasMZsigma)
      zzz = FitGMZR[0] + FitGMZR[1]*(yyy) + FitGMZR[2]*(yyy)*(yyy) + FitGMZR[3]*(yyy)*(yyy)*(yyy); + FitGMZR[4]*(yyy)*(yyy)*(yyy)*(yyy)
      www = WHERE(yyy ge 8.5 and yyy le 11.5)
      ;oplot, yyy(www), zzz(www), thick=6, linestyle = 0, color=fsc_color("black")
      
      ;Zg = alog10((metals(ww(w4))/gas(ww(w4)))/0.018)+8.69
      ;Zg = alog10((metals(ww(w4))/gas(ww(w4)))/0.0134)+8.69  
      ;Zg = alog10((metals(ww(w4))/gas(ww(w4)))/0.02)+9.0
      Zg = alog10((metals(ww(w4))/gas(ww(w4)))/0.0134)+9.0 ;Justifiable, as we're looking at a total metallicity here (not just oxygen abundance), and Asplund et al. (2009) gives Z=0.0134 (i.e. mass fraction in metals of the Sun). Also, oxygen abundance in Sun is 8.69, so full full metal abundance would be more ,like 9.0?!
 
      ;Plot 1 sigma spread in all-gal distribution:
      projmin=MIN(alog10(stars(ww(w4))*1e10))
      projmax=MAX(alog10(stars(ww(w4))*1e10))
      dim=uberround((projmax-projmin)/0.1,0)
      xx = findgen(dim)/10. + projmin
      lower1sigmag = fltarr(dim)
      upper1sigmag = fltarr(dim)     
      meanZg = fltarr(dim)
      ;medianZg = fltarr(dim)
      Zstddev1g = fltarr(dim)    
      j=0
      FOR i=projmin,projmax,0.1 DO BEGIN
        w1s = WHERE(alog10(stars(ww(w4))*1e10) ge i-0.05 AND alog10(stars(ww(w4))*1e10) lt i+0.05)
        IF(ww[0] ne -1) THEN BEGIN
          Zstddev1g[j] = stddev(Zg(w1s))
          lower1sigmag[j] = AVG(Zg(w1s)) - Zstddev1g[j]
          upper1sigmag[j] = AVG(Zg(w1s)) + Zstddev1g[j]
          meanZg[j] = MEAN(Zg(w1s))
          ;medianZg[j] = MEDIAN(Zg(w1s))
        ENDIF
        j=j+1
      ENDFOR
      www = WHERE(xx(www) ge 7.9 and xx(www) le 11.25)
      oplot, xx(www), meanZg(www), thick=6., linestyle = 0, color=fsc_color("black")
      ;oplot, xx(www), medianZg(www), thick=6., linestyle = 1, color=fsc_color("black")
      oplot, xx(www), lower1sigmag(www), thick=6., linestyle = 2, color=fsc_color("black")
      oplot, xx(www), upper1sigmag(www), thick=6., linestyle = 2, color=fsc_color("black")
      
;      print, meanZg(www)
;      print, "-----------"
;      print, lower1sigmag(www)
;      print, "-----------"
;      print, upper1sigmag(www) 

;; Old relation:
;      oplot, xx(www), [7.91404,7.96074,7.99553,8.03723,8.07996,8.11999,8.16015,8.20185,8.25475,8.28896,8.33810,8.38283,8.42819,8.47206,8.50745,8.55723,8.61917,8.65131,8.71168,8.73064,8.77592,8.82385,8.86413,8.89642,8.92927,8.95584,8.98615,8.99526,9.00173,9.01059,8.99664,9.00276,8.99802] $
;           , thick=6., linestyle = 0, color=fsc_color("red") ;color=fsc_color("dark green")
;      oplot, xx(www), [7.76847,7.81066,7.84580,7.88127,7.93088,7.97073,8.01075,8.05502,8.11491,8.14720,8.20559,8.26015,8.30106,8.35094,8.39021,8.45054,8.51545,8.54640,8.61616,8.61199,8.67647,8.71319,8.74918,8.79366,8.81178,8.84932,8.86315,8.85183,8.87039,8.86971,8.85767,8.85995,8.86571] $
;           , thick=6., linestyle = 2, color=fsc_color("red") ;color=fsc_color("dark green")
;      oplot, xx(www), [8.05961,8.11082,8.14525,8.19318,8.22905,8.26925,8.30956,8.34866,8.39460,8.43071,8.47061,8.50550,8.55531,8.59317,8.62469,8.66392,8.72288,8.75623,8.80720,8.84930,8.87538,8.93452,8.97907,8.99917,9.04675,9.06236,9.10916,9.13869,9.13307,9.15148,9.13560,9.14557,9.13033] $
;           , thick=6., linestyle = 2, color=fsc_color("red") ;color=fsc_color("dark green") 

      ;legend, ['Guo et al. (2011)', 'New model'], box = 0, charsize = 1.25, linestyle = [0,0], color = [fsc_color("red"),fsc_color("black")], thick = [6,6], /bottom, /right
      ;legend, ['Yates et al. (2012)'], box = 0, charsize = 1.25, linestyle = [0], color = [fsc_color("blue")], thick = [6], /top, /left
      legend, ['SDSS-DR7'], box = 0, charsize = 1.25, linestyle = [0], color = [fsc_color("blue")], thick = [6], /bottom, /right
      
; overplot observational results (ask Gabriella - Christy?)

;      xxx = findgen(100)/(100.0/(11.5-8.5)) + 8.5
;      yyy = -1.492 + 1.847*xxx - 0.08026*xxx*xxx
;      oplot, xxx, yyy, thick = 6, linestyle = 0, color = fsc_color("blue")
      
      xxx = findgen(100)/(100.0/(11.25-8.5)) + 8.5
      yyy = 26.6864 - 6.6399*xxx + 0.768653*xxx*xxx - 0.0282147*xxx*xxx*xxx
      ;yyy = 26.686 - 6.640*xxx + 0.769*xxx*xxx - 0.028*xxx*xxx*xxx
      oplot, xxx, yyy, thick = 6, linestyle = 0, color = fsc_color("blue")


;;---------------------------------------------------------------------------------------
;;Plot only SN-II ejecta MZR:
;      plot, findgen(10), /nodata, xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = 'Z!Dcold!N', $
;            charsize = csiz, xrange = [8.5, 12.0], yrange = [8.0, 9.7], xstyle = 1, ystyle = 1
;
;      metals = G.MetalsColdGas[1]
;
;      ;ww = where(galtype eq 0 and gas/(stars+gas) gt 0.1); and randomu(2121, n_elements(galtype)) lt dilute/10.0)
;      ww = where(galtype eq 0 OR galtype eq 1)
;      if ww[0] ne -1 then begin      
;        ;oplot, alog10(stars(ww)*1e10), alog10((metalsScaleSNIa(ww)/gas(ww))/0.0127)+9., psym = 8, symsize = ssiz      
;      
;        w1 = WHERE(alog10(G(ww).sfr) lt -1.0)
;        oplot, alog10(stars(ww(w1))*1e10), alog10((metals(ww(w1))/gas(ww(w1)))/0.02)+9., psym = 8, symsize = ssiz, color = fsc_color("red") ;Wiersma's solar metallicity (0.0127)
;        w2 = WHERE(alog10(G(ww).sfr) ge -1.0 AND alog10(G(ww).sfr) le 0.5)
;        oplot, alog10(stars(ww(w2))*1e10), alog10((metals(ww(w2))/gas(ww(w2)))/0.02)+9., psym = 8, symsize = ssiz, color = fsc_color("green")
;        w3 = WHERE(alog10(G(ww).sfr) gt 0.5)
;        oplot, alog10(stars(ww(w3))*1e10), alog10((metals(ww(w3))/gas(ww(w3)))/0.02)+9., psym = 8, symsize = ssiz, color = fsc_color("blue")       
;      endif
;
;      legend, ['Tremonti et al. (2003)'], box = 0, charsize = 1.25, linestyle = [0], color = [col1], thick = [tck], /bottom, /right
;      legend, ['log(SFR) < -1.0','-0.1 < log(SFR) < 0.5','log(SFR) > 0.5'], box = 0, charsize = 1.25, psym = [8,8,8], color = [fsc_color("red"),fsc_color("green"),fsc_color("blue")], /top, /left
;
;; overplot observational results (ask Gabriella - Christy?)
;
;      xxx = findgen(100)/20.0 + 8.0
;      yyy = -1.492 + 1.847*xxx - 0.08026*xxx*xxx
;      oplot, xxx, yyy, thick = tck, linestyle = 0, color = col1
;      
;      metals = G.MetalsColdGas[0]+G.MetalsColdGas[1]+G.MetalsColdGas[2]

;---------------------------------------------------------------------------------------
;;Real 12 + log[O/H] plot:
;
;      plot, findgen(10), /nodata, xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = '12 + log[O/H]', $
;            charsize = csiz, xrange = [8.5, 12.0], yrange = [8.0, 9.5], xstyle = 1, ystyle = 1
;
;      ;oxy = G.ColdGas_elements[2]
;      ;hyd = G.ColdGas_elements[0]
;      OH = alog10(G.ColdGas_elements[2]/G.ColdGas_elements[0]) - alog10(5.49E-3/0.7065)
;
;      legend, ['Tremonti et al. (2003)'], box = 0, charsize = csiz2, linestyle = [0], color = [col1], thick = [tck], /bottom, /right
;      legend, ['Sb/c galaxies'], box = 0, charsize = csiz2, linestyle = [-99], /top
;
;
;      w = where(galtype eq 0 and gas/(stars+gas) gt 0.1); and randomu(2121, n_elements(galtype)) lt dilute/10.0)
;      if w[0] ne -1 then begin      
;        w2 = WHERE(alog10(G(w).sfr) ge -1.0 AND alog10(G(w).sfr) le 0.5)
;        ;oplot, alog10(stars(w(w2))*1e10), alog10((oxy(w(w2))/hyd(w(w2)))/0.0127)+9., psym = 8, symsize = ssiz, color = fsc_color("green")
;        ;oplot, alog10(stars(w(w2))*1e10), alog10(OH(w2))+12., psym = 8, symsize = ssiz, color = fsc_color("green")
;        oplot, alog10(stars(w(w2))*1e10), alog10(OH(w2)/0.0127)+9., psym = 8, symsize = ssiz, color = fsc_color("green")
;        w3 = WHERE(alog10(G(w).sfr) gt 0.5)
;        ;oplot, alog10(stars(w(w3))*1e10), alog10((oxy(w(w3))/hyd(w(w3)))/0.0127)+9., psym = 8, symsize = ssiz, color = fsc_color("blue")
;        ;oplot, alog10(stars(w(w3))*1e10), alog10(OH(w3))+12., psym = 8, symsize = ssiz, color = fsc_color("blue")
;        oplot, alog10(stars(w(w3))*1e10), alog10(OH(w3)/0.0127)+9., psym = 8, symsize = ssiz, color = fsc_color("blue")
;        w1 = WHERE(alog10(G(w).sfr) lt -1.0)
;        ;oplot, alog10(stars(w(w1))*1e10), alog10((oxy(w(w1))/hyd(w(w1)))/0.0127)+9., psym = 8, symsize = ssiz, color = fsc_color("red") ;Wiersma's solar metallicity (0.0127)
;        ;oplot, alog10(stars(w(w1))*1e10), alog10(OH(w1))+12., psym = 8, symsize = ssiz, color = fsc_color("red")
;        oplot, alog10(stars(w(w1))*1e10), alog10(OH(w1)/0.0127)+9., psym = 8, symsize = ssiz, color = fsc_color("red")
;      endif
;
;
;; overplot observational results (ask Gabriella - Christy?)
;
;      xxx = findgen(100)/20.0 + 8.0
;      yyy = -1.492 + 1.847*xxx - 0.08026*xxx*xxx
;      oplot, xxx, yyy, thick = tck, linestyle = 0, color = col1

;---------------------------------------------------------------------------------------
;Plot stellar MZR:

;IF (OLDGALSTRUCT eq 1) THEN BEGIN
      ;uncorstars = G.StellarMass
      ;smetals = G.MetalsStellarMass
 ;ENDIF ELSE BEGIN
      ;uncorstars = G.DiskMass
      ;metals = G.MetalsDiskMass[0]+G.MetalsDiskMass[1]+G.MetalsDiskMass[2]
      ;uncorstars = G.DiskMass+G.BulgeMass
      smetals = G.MetalsDiskMass[0]+G.MetalsDiskMass[1]+G.MetalsDiskMass[2]+G.MetalsBulgeMass[0]+G.MetalsBulgeMass[1]+G.MetalsBulgeMass[2]
;ENDELSE

      plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = 'Z!D*!N/Z!D!9'+string(110B)+'!3!N', $
            charsize = csiz, xrange = [8.0, 11.5], yrange = [-1.5, 0.5], xstyle = 1, ystyle = 1 ;[7.5, 9.5]

      ;legend, ['Sb/c galaxies'], box = 0, charsize = csiz2, linestyle = [-99], /top

      ;ww = where(galtype eq 0 and gas/(stars+gas) gt 0.1); and randomu(2121, n_elements(galtype)) lt dilute/10.0)
      ww = where(galtype eq 0 OR galtype eq 1)
      if ww[0] ne -1 then begin      
        ;oplot, alog10(stars(ww)*1e10), alog10((metalsScaleSNIa(ww)/gas(ww))/0.0127)+9., psym = 8, symsize = ssiz      
      
        w1 = WHERE(alog10(G(ww).sfr) lt -1.0)
        oplot, alog10(stars(ww(w1))*1e10), alog10((smetals(ww(w1))/stars(ww(w1)))/0.02), psym = 8, symsize = ssiz, color = fsc_color("grey") ;Wiersma's solar metallicity (0.0127)
        w2 = WHERE(alog10(G(ww).sfr) ge -1.0 AND alog10(G(ww).sfr) le 0.5)
        oplot, alog10(stars(ww(w2))*1e10), alog10((smetals(ww(w2))/stars(ww(w2)))/0.02), psym = 8, symsize = ssiz, color = fsc_color("grey")
        w3 = WHERE(alog10(G(ww).sfr) gt 0.5)
        oplot, alog10(stars(ww(w3))*1e10), alog10((smetals(ww(w3))/stars(ww(w3)))/0.02), psym = 8, symsize = ssiz, color = fsc_color("grey")      
;        w1 = WHERE(alog10(G(ww).sfr) lt -1.0)
;        oplot, alog10(stars(ww(w1))*1e10), alog10((smetals(ww(w1))/stars(ww(w1)))/0.0127), psym = 8, symsize = ssiz, color = fsc_color("red") ;Wiersma's solar metallicity (0.0127)
;        w2 = WHERE(alog10(G(ww).sfr) ge -1.0 AND alog10(G(ww).sfr) le 0.5)
;        oplot, alog10(stars(ww(w2))*1e10), alog10((smetals(ww(w2))/stars(ww(w2)))/0.0127), psym = 8, symsize = ssiz, color = fsc_color("green")
;        w3 = WHERE(alog10(G(ww).sfr) gt 0.5)
;        oplot, alog10(stars(ww(w3))*1e10), alog10((smetals(ww(w3))/stars(ww(w3)))/0.0127), psym = 8, symsize = ssiz, color = fsc_color("blue")    
      endif
      
      ;Plot fit to all galaxies:
      w4 = WHERE(stars(ww) ne 0.0 AND smetals(ww) ne 0.0)
      yyy = findgen(100)/(100.0/(12.0-8.5)) + 8.5
      FitSMZR = POLY_FIT(alog10(stars(ww(w4))*1e10), alog10((smetals(ww(w4))/stars(ww(w4)))/0.02), 3, SIGMA=FitStellarMZsigma)
      zzz = FitSMZR[0] + FitSMZR[1]*(yyy) + FitSMZR[2]*(yyy)*(yyy) + FitSMZR[3]*(yyy)*(yyy)*(yyy); + FitSMZR[4]*(yyy)*(yyy)*(yyy)*(yyy)
      www = WHERE(yyy ge 8.5 and yyy le 11.5)
      ;oplot, yyy(www), zzz(www), thick=6, linestyle = 0, color=fsc_color("black")
      ;print, zzz
      FitSMZR = POLY_FIT(alog10(stars(ww(w4))*1e10), alog10((smetals(ww(w4))/stars(ww(w4)))/0.0127), 3, SIGMA=FitStellarMZsigma)
      zzz = FitSMZR[0] + FitSMZR[1]*(yyy) + FitSMZR[2]*(yyy)*(yyy) + FitSMZR[3]*(yyy)*(yyy)*(yyy); + FitSMZR[4]*(yyy)*(yyy)*(yyy)*(yyy)
      www = WHERE(yyy ge 8.5 and yyy le 11.5)
      ;oplot, yyy(www), zzz(www), thick=6, linestyle = 0, color=fsc_color("grey")
          
      Za = alog10((smetals(ww(w4))/stars(ww(w4)))/0.0127)  
      Zb = alog10((smetals(ww(w4))/stars(ww(w4)))/0.02)   
 
      ;Plot 1 sigma spread in all-gal distribution:
      projmin=MIN(alog10(stars(ww(w4))*1e10))
      projmax=MAX(alog10(stars(ww(w4))*1e10))
      dim=uberround((projmax-projmin)/0.1,0)
      xx = findgen(dim)/10. + projmin
      lower1sigmaa = fltarr(dim)
      lower1sigmab = fltarr(dim)
      upper1sigmaa = fltarr(dim)
      upper1sigmab = fltarr(dim)      
      meanZa = fltarr(dim)
      Zstddev1a = fltarr(dim)
      meanZb = fltarr(dim)
      Zstddev1b = fltarr(dim)      
      j=0
      FOR i=projmin,projmax,0.1 DO BEGIN
        w1s = WHERE(alog10(stars(ww(w4))*1e10) ge i-0.05 AND alog10(stars(ww(w4))*1e10) lt i+0.05)
        IF(ww[0] ne -1) THEN BEGIN
          Zstddev1a[j] = stddev(Za(w1s))
          lower1sigmaa[j] = AVG(Za(w1s)) - Zstddev1a[j]
          upper1sigmaa[j] = AVG(Za(w1s)) + Zstddev1a[j]
          meanZa[j] = MEAN(Za(w1s))
          Zstddev1b[j] = stddev(Zb(w1s))
          lower1sigmab[j] = AVG(Zb(w1s)) - Zstddev1b[j]
          upper1sigmab[j] = AVG(Zb(w1s)) + Zstddev1b[j]
          meanZb[j] = MEAN(Zb(w1s))
        ENDIF
        j=j+1
      ENDFOR
      www = WHERE(xx(www) ge 7.9 and xx(www) le 11.0)
      ;oplot, xx(www), meanZa(www), thick=6., linestyle = 0, color=fsc_color("dark grey")
      ;oplot, xx(www), lower1sigmaa(www), thick=6., linestyle = 2, color=fsc_color("dark grey")
      ;oplot, xx(www), upper1sigmaa(www), thick=6., linestyle = 2, color=fsc_color("dark grey")
      oplot, xx(www), meanZb(www), thick=6., linestyle = 0, color=fsc_color("black")
      oplot, xx(www), lower1sigmab(www), thick=6., linestyle = 2, color=fsc_color("black")
      oplot, xx(www), upper1sigmab(www), thick=6., linestyle = 2, color=fsc_color("black")
      
; Old relation:
      oplot, xx(www), [-1.12936,-1.07662,-1.02837,-0.978746,-0.931164,-0.890318,-0.850701,-0.816740,-0.764695,-0.723592,-0.690699,-0.650172,-0.610592,-0.573974,-0.532097,-0.495485,-0.444405,-0.417140,-0.379951,-0.336680,-0.308339,-0.283303,-0.254925,-0.226566,-0.210959,-0.194501,-0.167092,-0.146201,-0.125018,-0.130870,-0.114629] $
           , thick=6., linestyle = 0, color=fsc_color("red") ;color=fsc_color("dark green")
      oplot, xx(www), [-1.24992,-1.20135,-1.15437,-1.10629,-1.04805,-1.01710,-0.978461,-0.941839,-0.884364,-0.840403,-0.805147,-0.760309,-0.714601,-0.675101,-0.630491,-0.591299,-0.531214,-0.510068,-0.467629,-0.422577,-0.395979,-0.370581,-0.345908,-0.318115,-0.308524,-0.288772,-0.268637,-0.249771,-0.232671,-0.236298,-0.225512] $
           , thick=6., linestyle = 2, color=fsc_color("red") ;color=fsc_color("dark green")
      oplot, xx(www), [-1.00880,-0.951888,-0.902362,-0.851204,-0.814274,-0.763535,-0.722942,-0.691643,-0.645025,-0.606781,-0.576250,-0.540035,-0.506583,-0.472847,-0.433703,-0.399671,-0.357596,-0.324213,-0.292272,-0.250784,-0.220699,-0.196025,-0.163942,-0.135018,-0.113394,-0.100230,-0.0655471,-0.0426299,-0.0173646,-0.0254427,-0.00374631] $
           , thick=6., linestyle = 2, color=fsc_color("red") ;color=fsc_color("dark green")  

; Plot Gallazzi relation:      
      oplot, [8.91,9.11,9.31,9.51,9.72,9.91,10.11,10.31,10.51,10.72,10.91,11.11,11.31,11.51,11.72,11.91], [-0.6,-0.61,-0.65,-0.61,-0.52,-0.41,-0.23,-0.11,-0.01,0.04,0.07,0.1,0.12,0.13,0.14,0.15] $
           , linestyle = 0, thick = 6, color = fsc_color("blue")
      oplot, [8.91,9.11,9.31,9.51,9.72,9.91,10.11,10.31,10.51,10.72,10.91,11.11,11.31,11.51,11.72,11.91], [-1.11,-1.07,-1.1,-1.03,-0.97,-0.9,-0.8,-0.65,-0.41,-0.24,-0.14,-0.09,-0.06,-0.04,-0.03,-0.03] $
           , linestyle = 2, thick = 6, color = fsc_color("blue")           
      oplot, [8.91,9.11,9.31,9.51,9.72,9.91,10.11,10.31,10.51,10.72,10.91,11.11,11.31,11.51,11.72,11.91], [-0.00,-0.00,-0.05,-0.01,0.05,0.09,0.14,0.17,0.2,0.22,0.24,0.25,0.26,0.28,0.29,0.3] $
           , linestyle = 2, thick = 6, color = fsc_color("blue")                      
      
;      medianGallazzi = [-0.6,-0.61,-0.65,-0.61,-0.52,-0.41,-0.23,-0.11,-0.01,0.04,0.07,0.1,0.12,0.13,0.14,0.15]
;      lowerrors = [-1.11,-1.07,-1.1,-1.03,-0.97,-0.9,-0.8,-0.65,-0.41,-0.24,-0.14,-0.09,-0.06,-0.04,-0.03,-0.03]
;      higherrors = [-0.00,-0.00,-0.05,-0.01,0.05,0.09,0.14,0.17,0.2,0.22,0.24,0.25,0.26,0.28,0.29,0.3]
;      
;      oploterror, [8.91,9.11,9.31,9.51,9.72,9.91,10.11,10.31,10.51,10.72,10.91,11.11,11.31,11.51,11.72,11.91], [-0.6,-0.61,-0.65,-0.61,-0.52,-0.41,-0.23,-0.11,-0.01,0.04,0.07,0.1,0.12,0.13,0.14,0.15] $
;                , medianGallazzi-lowerrors, psym = 8, symsize = 4.*ssiz, color = fsc_color("blue"), errthick = 6., errcolor = fsc_color("blue"), /LOBAR 
;      oploterror, [8.91,9.11,9.31,9.51,9.72,9.91,10.11,10.31,10.51,10.72,10.91,11.11,11.31,11.51,11.72,11.91], [-0.6,-0.61,-0.65,-0.61,-0.52,-0.41,-0.23,-0.11,-0.01,0.04,0.07,0.1,0.12,0.13,0.14,0.15] $
;                , higherrors-medianGallazzi, psym = 8, symsize = 4.*ssiz, color = fsc_color("blue"), errthick = 6., errcolor = fsc_color("blue"), /HIBAR

      ;legend, ['Guo et al. (2011)', 'Z!D!9'+string(110B)+'!3!N = 0.02', 'Z!D!9'+string(110B)+'!3!N = 0.0127'], box = 0, charsize = 1.25, linestyle = [0,0,0], color = [fsc_color("dark green"),fsc_color("black"),fsc_color("dark grey")], thick = [6,6,6], /bottom, /right
      legend, ['Guo et al. (2011)', 'New model'], box = 0, charsize = 1.25, linestyle = [0,0], color = [fsc_color("red"),fsc_color("black")], thick = [6,6], /bottom, /right
      ;legend, ['log(SFR) < -1.0','-0.1 < log(SFR) < 0.5','log(SFR) > 0.5', 'Gallazzi et al. (2005)'], box = 0, charsize = 1.25, psym = [8,8,8,8], color = [fsc_color("red"),fsc_color("green"),fsc_color("blue"),fsc_color("black")], symsize=[ssiz,ssiz,ssiz,4*ssiz], /top, /left
      legend, ['Gallazzi et al. (2005)'], box = 0, charsize = 1.25, linestyle = [0], color = [fsc_color("blue")], thick = 6, /top, /left

;---------------------------------------------------------------------------------------
;---------------------------------------------------------------------------------------
IF (OLDGALSTRUCT eq 0) THEN BEGIN 
;ENHANCEMENT PLOTS:
;---------------------------------------------------------------------------------------
;---------------------------------------------------------------------------------------
;-------------------------------
;;-------------------------------
;Filter for Milky Way type haloes:
w0 = WHERE(alog10(G.CentralMvir) ge 11.75 AND alog10(G.CentralMvir) le 12.25)
ww0 = WHERE((alog10(G.CentralMvir) ge 11.75) AND (alog10(G.CentralMvir) le 12.25) AND (G.BulgeMass/(G.BulgeMass + G.DiskMass) lt 0.7) AND (G.type eq 0))
www0 = WHERE((alog10(G.CentralMvir) ge 11.75) AND (alog10(G.CentralMvir) le 12.25) AND (G.sfr ge 1.0) AND (G.sfr le 5.0) AND (G.BulgeMass/(G.BulgeMass + G.DiskMass) lt 0.7) AND (G.type eq 0))

;Filter for ellipticals:
IF (MAINELEMENTS eq 0) THEN BEGIN
w1 = WHERE((G.BulgeMass/(G.BulgeMass + G.DiskMass) ge 0.7) $
          AND (FINITE(G.DiskMass_elements[*], /NAN, SIGN=0) ne 1) $
          AND (FINITE(G.BulgeMass_elements[*], /NAN, SIGN=0) ne 1) $
          AND (G.DiskMass_elements[0] ne 0.0) $
          AND (G.DiskMass_elements[1] ne 0.0) $
          AND (G.DiskMass_elements[2] ne 0.0) $
          AND (G.DiskMass_elements[3] ne 0.0) $
          AND (G.DiskMass_elements[4] ne 0.0) $
          AND (G.BulgeMass_elements[0] ne 0.0) $
          AND (G.BulgeMass_elements[1] ne 0.0) $
          AND (G.BulgeMass_elements[2] ne 0.0) $
          AND (G.BulgeMass_elements[3] ne 0.0) $
          AND (G.BulgeMass_elements[4] ne 0.0))
          
ENDIF ELSE IF (MAINELEMENTS eq 1) THEN BEGIN   
w1 = WHERE((G.BulgeMass/(G.BulgeMass + G.DiskMass) ge 0.7) $
          AND (FINITE(G.DiskMass_elements[*], /NAN, SIGN=0) ne 1) $
          AND (FINITE(G.BulgeMass_elements[*], /NAN, SIGN=0) ne 1) $
          AND (G.DiskMass_elements[0] ne 0.0) $
          AND (G.DiskMass_elements[1] ne 0.0) $
          AND (G.DiskMass_elements[2] ne 0.0) $
          AND (G.DiskMass_elements[3] ne 0.0) $
          AND (G.DiskMass_elements[4] ne 0.0) $
          AND (G.DiskMass_elements[5] ne 0.0) $
          AND (G.DiskMass_elements[6] ne 0.0) $
          AND (G.DiskMass_elements[7] ne 0.0) $
          AND (G.DiskMass_elements[8] ne 0.0) $
          AND (G.DiskMass_elements[9] ne 0.0) $
          AND (G.DiskMass_elements[10] ne 0.0) $
          AND (G.BulgeMass_elements[0] ne 0.0) $
          AND (G.BulgeMass_elements[1] ne 0.0) $
          AND (G.BulgeMass_elements[2] ne 0.0) $
          AND (G.BulgeMass_elements[3] ne 0.0) $
          AND (G.BulgeMass_elements[4] ne 0.0) $
          AND (G.BulgeMass_elements[5] ne 0.0) $
          AND (G.BulgeMass_elements[6] ne 0.0) $
          AND (G.BulgeMass_elements[7] ne 0.0) $
          AND (G.BulgeMass_elements[8] ne 0.0) $
          AND (G.BulgeMass_elements[9] ne 0.0) $
          AND (G.BulgeMass_elements[10] ne 0.0))
ENDIF  

print, "n_elements(w1) = ", n_elements(w1)

;----------------------------------------
;PRINT DATA:
;----------------------------------------
nanfiltD = WHERE(FINITE(G(w0).DiskMass_elements[*], /NAN, SIGN=0) ne 1 $
             AND FINITE(G(w0).DiskMass, /NAN, SIGN=0) ne 1 $
             AND G(w0).DiskMass_elements[*] ne 0.0 $
             AND G(w0).DiskMass ne 0.0 $
             AND G(w0).DiskMass_elements[*] ge 0.0)
print, n_elements(nanfiltD)
nanfiltB = WHERE(FINITE(G(w0).BulgeMass_elements[*], /NAN, SIGN=0) ne 1 $
             AND FINITE(G(w0).BulgeMass, /NAN, SIGN=0) ne 1 $
             AND G(w0).BulgeMass_elements[*] ne 0.0 $
             AND G(w0).BulgeMass ne 0.0 $
             AND G(w0).BulgeMass_elements[*] ge 0.0)
print, n_elements(nanfiltB)

IF (MAINELEMENTS eq 0) THEN BEGIN
print, "DISK"
print, "M_He/M_H (code) = ", MEAN(G(www0).DiskMass_elements[1]/G(www0).DiskMass_elements[0])
print, "M_He/M_H (Asplund 2009) = ", (4.003/1.008)*10.0^(10.93-12)
print, "M_O/M_H (code) = ", MEAN(G(www0).DiskMass_elements[2]/G(www0).DiskMass_elements[0])
print, "M_O/M_H (Asplund 2009) = ", (16.00/1.008)*10.0^(8.69-12)
;print, "M_Mg/M_H (code) = ", MEAN(G[w0(nanfiltD)].DiskMass_elements[3]/G[w0(nanfiltD)].DiskMass_elements[0])
print, "M_Mg/M_H (code) = ", MEAN(G(www0).DiskMass_elements[3]/G(www0).DiskMass_elements[0])
print, "M_Mg/M_H (Asplund 2009) = ", (24.31/1.008)*10.0^(7.6-12)
;print, "M_Fe/M_H (code) = ", MEAN(G(w0(nanfiltD)).DiskMass_elements[4]/G(w0(nanfiltD)).DiskMass_elements[0])
print, "M_Fe/M_H (code) = ", MEAN(G(www0).DiskMass_elements[4]/G(www0).DiskMass_elements[0])
print, "M_Fe/M_H (Asplund 2009) = ", (55.84/1.008)*10.0^(7.5-12)
print, "H mass fraction = ", MEAN(G(www0).DiskMass_elements[0]/(G(www0).DiskMass*1.0e10/Hubble_h))
print, "He mass fraction = ", MEAN(G(www0).DiskMass_elements[1]/(G(www0).DiskMass*1.0e10/Hubble_h))
print, "O mass fraction = ", MEAN(G(www0).DiskMass_elements[2]/(G(www0).DiskMass*1.0e10/Hubble_h))
print, "Mg mass fraction = ", MEAN(G(www0).DiskMass_elements[3]/(G(www0).DiskMass*1.0e10/Hubble_h))
print, "Fe mass fraction = ", MEAN(G(www0).DiskMass_elements[4]/(G(www0).DiskMass*1.0e10/Hubble_h))

print, " "
print, "BULGE"
print, "M_He/M_H (Bulge)= ", MEAN(G(www0).BulgeMass_elements[1]/G(www0).BulgeMass_elements[0]), G(www0[0]).BulgeMass_elements[1], G(www0[0]).BulgeMass_elements[0]
print, "M_He/M_H (Disk) = ", MEAN(G(www0).DiskMass_elements[1]/G(www0).DiskMass_elements[0]), G(www0[0]).DiskMass_elements[1], G(www0[0]).DiskMass_elements[0]
print, "M_He/M_H (Cold) = ", MEAN(G(www0).ColdGas_elements[1]/G(www0).ColdGas_elements[0]), G(www0[0]).ColdGas_elements[1], G(www0[0]).ColdGas_elements[0]
print, "M_He/M_H (Hot)  = ", MEAN(G(www0).HotGas_elements[1]/G(www0).HotGas_elements[0]), G(www0[0]).HotGas_elements[1], G(www0[0]).HotGas_elements[0]
;print, "H mass fraction for one Gal = ", G(w0(nanfiltB[0])).BulgeMass_elements[0]/(G(w0(nanfiltB[0])).BulgeMass*1.0e10), " ( ", G(w0(nanfiltB[0])).BulgeMass_elements[0], " / ", (G(w0(nanfiltB[0])).BulgeMass*1.0e10), " )"
print, "H mass fraction = ", MEAN(G(w0(nanfiltB)).BulgeMass_elements[0]/(G(w0(nanfiltB)).BulgeMass*1.0e10))
;print, "H mass fraction = ", MEAN(G(w0(0)).BulgeMass_elements[0]/G(w0(0)).BulgeMass)
print, "O mass fraction = ", MEAN(G(w0(nanfiltB)).BulgeMass_elements[2]/(G(w0(nanfiltB)).BulgeMass*1.0e10))
print, "Mg mass fraction = ", MEAN(G(w0(nanfiltB)).BulgeMass_elements[3]/(G(w0(nanfiltB)).BulgeMass*1.0e10))
print, "Fe mass fraction = ", MEAN(G(w0(nanfiltB)).BulgeMass_elements[4]/(G(w0(nanfiltB)).BulgeMass*1.0e10))
print, "Mean Bulge mass = ", MEAN((G(w0(nanfiltB)).BulgeMass*1.0e10))
print, "Mean H mass = ", MEAN(G(w0(nanfiltB)).BulgeMass_elements[0])
print, "Mean Mg mass = ", MEAN(G(w0(nanfiltB)).BulgeMass_elements[3])
print, " "
print, "DISK + BULGE"
print, "H mass fraction = ", MEAN((G(www0).DiskMass_elements[0]+G(www0).BulgeMass_elements[0])/((G(www0).DiskMass+G(www0).BulgeMass)*1.0e10/Hubble_h))
print, "He mass fraction = ", MEAN((G(www0).DiskMass_elements[1]+G(www0).BulgeMass_elements[1])/((G(www0).DiskMass+G(www0).BulgeMass)*1.0e10/Hubble_h))
print, "O mass fraction = ", MEAN((G(www0).DiskMass_elements[2]+G(www0).BulgeMass_elements[2])/((G(www0).DiskMass+G(www0).BulgeMass)*1.0e10/Hubble_h))
print, "Mg mass fraction = ", MEAN((G(www0).DiskMass_elements[3]+G(www0).BulgeMass_elements[3])/((G(www0).DiskMass+G(www0).BulgeMass)*1.0e10/Hubble_h))
print, "Fe mass fraction = ", MEAN((G(www0).DiskMass_elements[4]+G(www0).BulgeMass_elements[4])/((G(www0).DiskMass+G(www0).BulgeMass)*1.0e10/Hubble_h))
;print, "Total Bulge Element Mass = ", TOTAL(G(w0(nanfiltB)).BulgeMass_elements), "Bulge Mass = ", G(w0(nanfiltB)).BulgeMass/Hubble_h*1.0e10
TotalElements = G(w0(nanfiltB)).BulgeMass_elements[*]
EleToBulge = MEAN(TotalElements/(G(w0(nanfiltB)).BulgeMass/Hubble_h*1.0e10))
print, "Mean Elements/BulgeMass = ", EleToBulge ;MEAN(TOTAL(G(w0(nanfiltB)).BulgeMass_elements)/(G(w0(nanfiltB)).BulgeMass/Hubble_h*1.0e10))
print, "Gal 0: Elements/BulgeMass = ", TOTAL(G[0].BulgeMass_elements)/(G[0].BulgeMass*1.0e10)
print, "Gal 0: (Sum Metal Elements)/TotalMetals = ", (TOTAL(G[0].BulgeMass_elements)-G[0].BulgeMass_elements[0]-G[0].BulgeMass_elements[1])/((G[0].MetalsBulgeMass[0]+G[0].MetalsBulgeMass[1]+G[0].MetalsBulgeMass[2])*1.0e10)
print, " "

print, "SFH_H = ", G[10].sfh_ElementsDiskMass[0,*]
print, " "
print, "SFH_DiskMass = ", G[10].sfh_DiskMass[*]*1.0e10
print, " "
print, "SFH: ", "    H      ", "      He      ", "      O     ", "      Mg      ", "     Fe"
print, G[10].sfh_ElementsDiskMass[*,*]
print, " "
ENDIF

IF (MAINELEMENTS eq 1) THEN BEGIN
print, "DISK"
print, "M_He/M_H (code) = ", MEAN(G(www0).DiskMass_elements[1]/G(www0).DiskMass_elements[0])
print, "M_He/M_H (Asplund 2009) = ", (4.003/1.008)*10.0^(10.93-12.)
print, "M_C/M_H (code) = ", MEAN(G(www0).DiskMass_elements[2]/G(www0).DiskMass_elements[0])
print, "M_C/M_H (Asplund 2009) = ", (12.01/1.008)*10.0^(8.43-12.)
print, "M_N/M_H (code) = ", MEAN(G(www0).DiskMass_elements[3]/G(www0).DiskMass_elements[0])
print, "M_N/M_H (Asplund 2009) = ", (14.01/1.008)*10.0^(7.83-12.)
print, "M_O/M_H (code) = ", MEAN(G(www0).DiskMass_elements[4]/G(www0).DiskMass_elements[0])
print, "M_O/M_H (Asplund 2009) = ", (16.00/1.008)*10.0^(8.69-12.)
print, "M_Ne/M_H (code) = ", MEAN(G(www0).DiskMass_elements[5]/G(www0).DiskMass_elements[0])
print, "M_Ne/M_H (Asplund 2009) = ", (20.18/1.008)*10.0^(7.93-12.)
print, "M_Mg/M_H (code) = ", MEAN(G(www0).DiskMass_elements[6]/G(www0).DiskMass_elements[0])
print, "M_Mg/M_H (Asplund 2009) = ", (24.31/1.008)*10.0^(7.6-12.)
print, "M_Si/M_H (code) = ", MEAN(G(www0).DiskMass_elements[7]/G(www0).DiskMass_elements[0])
print, "M_Si/M_H (Asplund 2009) = ", (28.09/1.008)*10.0^(7.51-12.)
print, "M_S/M_H (code) = ", MEAN(G(www0).DiskMass_elements[8]/G(www0).DiskMass_elements[0])
print, "M_S/M_H (Asplund 2009) = ", (32.07/1.008)*10.0^(7.12-12.)
print, "M_Ca/M_H (code) = ", MEAN(G(www0).DiskMass_elements[9]/G(www0).DiskMass_elements[0])
print, "M_Ca/M_H (Asplund 2009) = ", (40.08/1.008)*10.0^(6.34-12.)
print, "M_Fe/M_H (code) = ", MEAN(G(www0).DiskMass_elements[10]/G(www0).DiskMass_elements[0])
print, "M_Fe/M_H (Asplund 2009) = ", (55.84/1.008)*10.0^(7.5-12.)
print, "H mass fraction = ", MEAN(G(www0).DiskMass_elements[0]/(G(www0).DiskMass*1.0e10/Hubble_h))
print, "He mass fraction = ", MEAN(G(www0).DiskMass_elements[1]/(G(www0).DiskMass*1.0e10/Hubble_h))
print, "O mass fraction = ", MEAN(G(www0).DiskMass_elements[4]/(G(www0).DiskMass*1.0e10/Hubble_h))
print, "Mg mass fraction = ", MEAN(G(www0).DiskMass_elements[6]/(G(www0).DiskMass*1.0e10/Hubble_h))
print, "Fe mass fraction = ", MEAN(G(www0).DiskMass_elements[10]/(G(www0).DiskMass*1.0e10/Hubble_h))
MFCounta = 0.0
FOR d=0, 10 DO BEGIN
  ;MFCounta += MEAN(G(www0).DiskMass_elements[d]/(G(www0).DiskMass*1.0e10/Hubble_h))
  MFCounta += G(www0[0]).DiskMass_elements[d]/(G(www0[0]).DiskMass*1.0e10/Hubble_h)
ENDFOR
print, "SUM(Mass fractions) = ", MFCounta
print, "SUM(Elements) = ", TOTAL(G(www0[0]).DiskMass_elements[*])
print, "DiskMass = ", G(www0[0]).DiskMass*1.0e10/Hubble_h

print, " "
print, "BULGE"
print, "H mass fraction = ", MEAN(G(w0(nanfiltB)).BulgeMass_elements[0]/(G(w0(nanfiltB)).BulgeMass*1.0e10))
;print, "H mass fraction = ", MEAN(G(w0(0)).BulgeMass_elements[0]/G(w0(0)).BulgeMass)
print, "O mass fraction = ", MEAN(G(w0(nanfiltB)).BulgeMass_elements[4]/(G(w0(nanfiltB)).BulgeMass*1.0e10))
print, "Mg mass fraction = ", MEAN(G(w0(nanfiltB)).BulgeMass_elements[6]/(G(w0(nanfiltB)).BulgeMass*1.0e10))
print, "Fe mass fraction = ", MEAN(G(w0(nanfiltB)).BulgeMass_elements[10]/(G(w0(nanfiltB)).BulgeMass*1.0e10))
print, "Mean Bulge mass = ", MEAN((G(w0(nanfiltB)).BulgeMass*1.0e10))
print, "Mean H mass = ", MEAN(G(w0(nanfiltB)).BulgeMass_elements[0])
print, "Mean Mg mass = ", MEAN(G(w0(nanfiltB)).BulgeMass_elements[6])
print, " "
print, "DISK + BULGE"
print, "H mass fraction = ", MEAN((G(www0).DiskMass_elements[0]+G(www0).BulgeMass_elements[0])/((G(www0).DiskMass+G(www0).BulgeMass)*1.0e10/Hubble_h))
print, "He mass fraction = ", MEAN((G(www0).DiskMass_elements[1]+G(www0).BulgeMass_elements[1])/((G(www0).DiskMass+G(www0).BulgeMass)*1.0e10/Hubble_h))
print, "O mass fraction = ", MEAN((G(www0).DiskMass_elements[4]+G(www0).BulgeMass_elements[4])/((G(www0).DiskMass+G(www0).BulgeMass)*1.0e10/Hubble_h))
print, "Mg mass fraction = ", MEAN((G(www0).DiskMass_elements[6]+G(www0).BulgeMass_elements[6])/((G(www0).DiskMass+G(www0).BulgeMass)*1.0e10/Hubble_h))
print, "Fe mass fraction = ", MEAN((G(www0).DiskMass_elements[10]+G(www0).BulgeMass_elements[10])/((G(www0).DiskMass+G(www0).BulgeMass)*1.0e10/Hubble_h))
;print, "Total Bulge Element Mass = ", TOTAL(G(w0(nanfiltB)).BulgeMass_elements), "Bulge Mass = ", G(w0(nanfiltB)).BulgeMass/Hubble_h*1.0e10
TotalElements = G(w0(nanfiltB)).BulgeMass_elements[*]
EleToBulge = MEAN(TotalElements/(G(w0(nanfiltB)).BulgeMass/Hubble_h*1.0e10))
print, "Mean Elements/BulgeMass = ", EleToBulge ;MEAN(TOTAL(G(w0(nanfiltB)).BulgeMass_elements)/(G(w0(nanfiltB)).BulgeMass/Hubble_h*1.0e10))
print, "Gal 0: Elements/BulgeMass = ", TOTAL(G[0].BulgeMass_elements)/(G[0].BulgeMass*1.0e10)
print, "Gal 0: (Sum Metal Elements)/TotalMetals = ", (TOTAL(G[0].BulgeMass_elements)-G[0].BulgeMass_elements[0]-G[0].BulgeMass_elements[1])/((G[0].MetalsBulgeMass[0]+G[0].MetalsBulgeMass[1]+G[0].MetalsBulgeMass[2])*1.0e10)
print, " "
ENDIF

;print, "MassWeightAge = ", G[0:9].MassWeightAge
print, "MW-type galaxy discs:"
print, "SN-II: ", MEAN(G(www0).MetalsDiskMass[1]/(G(www0).MetalsDiskMass[0]+G(www0).MetalsDiskMass[1]+G(www0).MetalsDiskMass[2]))
print, "SN-Ia: ", MEAN(G(www0).MetalsDiskMass[0]/(G(www0).MetalsDiskMass[0]+G(www0).MetalsDiskMass[1]+G(www0).MetalsDiskMass[2]))
print, "AGB:   ", MEAN(G(www0).MetalsDiskMass[2]/(G(www0).MetalsDiskMass[0]+G(www0).MetalsDiskMass[1]+G(www0).MetalsDiskMass[2]))
print, "SN-II[0]: ", G(www0[0]).MetalsDiskMass[1]/(G(www0[0]).MetalsDiskMass[0]+G(www0[0]).MetalsDiskMass[1]+G(www0[0]).MetalsDiskMass[2])
print, "SN-Ia[0]: ", G(www0[0]).MetalsDiskMass[0]/(G(www0[0]).MetalsDiskMass[0]+G(www0[0]).MetalsDiskMass[1]+G(www0[0]).MetalsDiskMass[2])
print, "AGB[0]:   ", G(www0[0]).MetalsDiskMass[2]/(G(www0[0]).MetalsDiskMass[0]+G(www0[0]).MetalsDiskMass[1]+G(www0[0]).MetalsDiskMass[2])
print, "  "
print, "Elliptical galaxy bulges:"
print, "SN-II: ", MEAN(G(w1).MetalsBulgeMass[1]/(G(w1).MetalsBulgeMass[0]+G(w1).MetalsBulgeMass[1]+G(w1).MetalsBulgeMass[2]))
print, "SN-Ia: ", MEAN(G(w1).MetalsBulgeMass[0]/(G(w1).MetalsBulgeMass[0]+G(w1).MetalsBulgeMass[1]+G(w1).MetalsBulgeMass[2]))
print, "AGB:   ", MEAN(G(w1).MetalsBulgeMass[2]/(G(w1).MetalsBulgeMass[0]+G(w1).MetalsBulgeMass[1]+G(w1).MetalsBulgeMass[2]))
print, "SN-II[0]: ", G(w1[0]).MetalsBulgeMass[1]/(G(w1[0]).MetalsBulgeMass[0]+G(w1[0]).MetalsBulgeMass[1]+G(w1[0]).MetalsBulgeMass[2])
print, "SN-Ia[0]: ", G(w1[0]).MetalsBulgeMass[0]/(G(w1[0]).MetalsBulgeMass[0]+G(w1[0]).MetalsBulgeMass[1]+G(w1[0]).MetalsBulgeMass[2])
print, "AGB[0]:   ", G(w1[0]).MetalsBulgeMass[2]/(G(w1[0]).MetalsBulgeMass[0]+G(w1[0]).MetalsBulgeMass[1]+G(w1[0]).MetalsBulgeMass[2])

;----------------------------------------

;Solar mass fractions (from Wiersma et al. 2009a):
H_mf = 0.7065
He_mf = 0.2806
C_mf = 2.07E-3
N_mf = 8.36E-4
O_mf = 5.49E-3
Ne_mf = 1.41E-3
Mg_mf = 5.91E-4
Si_mf = 6.83E-4
S_mf = 4.09E-4
Ca_mf = 6.44E-5
Fe_mf = 1.1E-3

;Solar n_i/n_H ratios (from Wiersma et al. 2009a):
H_nr = 1
He_nr = 0.1
C_nr = 2.46E-4
N_nr = 8.51E-5
O_nr = 4.9E-4
Ne_nr = 1.0E-4
Mg_nr = 3.47E-5
Si_nr = 3.47E-5 ;Same as Mg
S_nr = 1.86E-5
Ca_nr = 2.29E-6
Fe_nr = 2.82E-5

;Atomic weights (from Wikipedia!):
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

IF (MAINELEMENTS eq 0) THEN BEGIN
;Disk abundance ratios ("MASS FRACTION" VERSION):
HFe_mf_disk = alog10(G.DiskMass_elements[0]/G.DiskMass_elements[4]) - alog10(H_mf/Fe_mf)
HeFe_mf_disk = alog10(G.DiskMass_elements[1]/G.DiskMass_elements[4]) - alog10(He_mf/Fe_mf)
OFe_mf_disk = alog10(G.DiskMass_elements[2]/G.DiskMass_elements[4]) - alog10(O_mf/Fe_mf)
MgFe_mf_disk = alog10(G.DiskMass_elements[3]/G.DiskMass_elements[4]) - alog10(Mg_mf/Fe_mf)
print, "OFe_mf_disk = ", MEAN(G(w1).DiskMass_elements[2]), MEAN(G(w1).DiskMass_elements[4]), MIN(G(w1).DiskMass_elements[2]), MIN(G(w1).DiskMass_elements[4])

HeH_mf_disk = alog10(G.DiskMass_elements[1]/G.DiskMass_elements[0]) - alog10(He_mf/H_mf)
OH_mf_disk = alog10(G.DiskMass_elements[2]/G.DiskMass_elements[0]) - alog10(O_mf/H_mf)
MgH_mf_disk = alog10(G.DiskMass_elements[3]/G.DiskMass_elements[0]) - alog10(Mg_mf/H_mf)
FeH_mf_disk = alog10(G.DiskMass_elements[4]/G.DiskMass_elements[0]) - alog10(Fe_mf/H_mf)

;Bulge abundance ratios ("MASS FRACTION" VERSION):
HFe_mf_bulge = alog10(G.BulgeMass_elements[0]/G.BulgeMass_elements[4]) - alog10(H_mf/Fe_mf)
HeFe_mf_bulge = alog10(G.BulgeMass_elements[1]/G.BulgeMass_elements[4]) - alog10(He_mf/Fe_mf)
OFe_mf_bulge = alog10(G.BulgeMass_elements[2]/G.BulgeMass_elements[4]) - alog10(O_mf/Fe_mf)
MgFe_mf_bulge = alog10(G.BulgeMass_elements[3]/G.BulgeMass_elements[4]) - alog10(Mg_mf/Fe_mf)

HeH_mf_bulge = alog10(G.BulgeMass_elements[1]/G.BulgeMass_elements[0]) - alog10(He_mf/H_mf)
OH_mf_bulge = alog10(G.BulgeMass_elements[2]/G.BulgeMass_elements[0]) - alog10(O_mf/H_mf)
MgH_mf_bulge = alog10(G.BulgeMass_elements[3]/G.BulgeMass_elements[0]) - alog10(Mg_mf/H_mf)
FeH_mf_bulge = alog10(G.BulgeMass_elements[4]/G.BulgeMass_elements[0]) - alog10(Fe_mf/H_mf)

;Combined disk + bulge abundance ratios ("MASS FRACTION" VERSION):
HFe_mf_comb = alog10((G.DiskMass_elements[0]+G.BulgeMass_elements[0])/(G.DiskMass_elements[4]+G.BulgeMass_elements[4])) - alog10(H_mf/Fe_mf)
HeFe_mf_comb = alog10((G.DiskMass_elements[1]+G.BulgeMass_elements[1])/(G.DiskMass_elements[4]+G.BulgeMass_elements[4])) - alog10(He_mf/Fe_mf)
OFe_mf_comb = alog10((G.DiskMass_elements[2]+G.BulgeMass_elements[2])/(G.DiskMass_elements[4]+G.BulgeMass_elements[4])) - alog10(O_mf/Fe_mf)
MgFe_mf_comb = alog10((G.DiskMass_elements[3]+G.BulgeMass_elements[3])/(G.DiskMass_elements[4]+G.BulgeMass_elements[4])) - alog10(Mg_mf/Fe_mf)

HeH_mf_comb = alog10((G.DiskMass_elements[1]+G.BulgeMass_elements[1])/(G.DiskMass_elements[0]+G.BulgeMass_elements[0])) - alog10(He_mf/H_mf)
OH_mf_comb = alog10((G.DiskMass_elements[2]+G.BulgeMass_elements[2])/(G.DiskMass_elements[0]+G.BulgeMass_elements[0])) - alog10(O_mf/H_mf)
MgH_mf_comb = alog10((G.DiskMass_elements[3]+G.BulgeMass_elements[3])/(G.DiskMass_elements[0]+G.BulgeMass_elements[0])) - alog10(Mg_mf/H_mf)
FeH_mf_comb = alog10((G.DiskMass_elements[4]+G.BulgeMass_elements[4])/(G.DiskMass_elements[0]+G.BulgeMass_elements[0])) - alog10(Fe_mf/H_mf)

;Cold Gas abundance ratios ("MASS FRACTION" VERSION):
HFe_mf_cold = alog10(G.ColdGas_elements[0]/G.ColdGas_elements[4]) - alog10(H_mf/Fe_mf)
HeFe_mf_cold = alog10(G.ColdGas_elements[1]/G.ColdGas_elements[4]) - alog10(He_mf/Fe_mf)
OFe_mf_cold = alog10(G.ColdGas_elements[2]/G.ColdGas_elements[4]) - alog10(O_mf/Fe_mf)
MgFe_mf_cold = alog10(G.ColdGas_elements[3]/G.ColdGas_elements[4]) - alog10(Mg_mf/Fe_mf)

HeH_mf_cold = alog10(G.ColdGas_elements[1]/G.ColdGas_elements[0]) - alog10(He_mf/H_mf)
OH_mf_cold = alog10(G.ColdGas_elements[2]/G.ColdGas_elements[0]) - alog10(O_mf/H_mf)
MgH_mf_cold = alog10(G.ColdGas_elements[3]/G.ColdGas_elements[0]) - alog10(Mg_mf/H_mf)
FeH_mf_cold = alog10(G.ColdGas_elements[4]/G.ColdGas_elements[0]) - alog10(Fe_mf/H_mf)

;Hot Gas abundance ratios ("MASS FRACTION" VERSION):
HFe_mf_hot = alog10(G.HotGas_elements[0]/G.HotGas_elements[4]) - alog10(H_mf/Fe_mf)
HeFe_mf_hot = alog10(G.HotGas_elements[1]/G.HotGas_elements[4]) - alog10(He_mf/Fe_mf)
OFe_mf_hot = alog10(G.HotGas_elements[2]/G.HotGas_elements[4]) - alog10(O_mf/Fe_mf)
MgFe_mf_hot = alog10(G.HotGas_elements[3]/G.HotGas_elements[4]) - alog10(Mg_mf/Fe_mf)

HeH_mf_hot = alog10(G.HotGas_elements[1]/G.HotGas_elements[0]) - alog10(He_mf/H_mf)
OH_mf_hot = alog10(G.HotGas_elements[2]/G.HotGas_elements[0]) - alog10(O_mf/H_mf)
MgH_mf_hot = alog10(G.HotGas_elements[3]/G.HotGas_elements[0]) - alog10(Mg_mf/H_mf)
FeH_mf_hot = alog10(G.HotGas_elements[4]/G.HotGas_elements[0]) - alog10(Fe_mf/H_mf)
ENDIF

IF (MAINELEMENTS eq 1) THEN BEGIN
;Disk abundance ratios ("MASS FRACTION" VERSION):
HFe_mf_disk = alog10(G.DiskMass_elements[0]/G.DiskMass_elements[10]) - alog10(H_mf/Fe_mf)
HeFe_mf_disk = alog10(G.DiskMass_elements[1]/G.DiskMass_elements[10]) - alog10(He_mf/Fe_mf)
CFe_mf_disk = alog10(G.DiskMass_elements[2]/G.DiskMass_elements[10]) - alog10(C_mf/Fe_mf)
NFe_mf_disk = alog10(G.DiskMass_elements[3]/G.DiskMass_elements[10]) - alog10(N_mf/Fe_mf)
OFe_mf_disk = alog10(G.DiskMass_elements[4]/G.DiskMass_elements[10]) - alog10(O_mf/Fe_mf)
NeFe_mf_disk = alog10(G.DiskMass_elements[5]/G.DiskMass_elements[10]) - alog10(Ne_mf/Fe_mf)
MgFe_mf_disk = alog10(G.DiskMass_elements[6]/G.DiskMass_elements[10]) - alog10(Mg_mf/Fe_mf)
SiFe_mf_disk = alog10(G.DiskMass_elements[7]/G.DiskMass_elements[10]) - alog10(Si_mf/Fe_mf)
SFe_mf_disk = alog10(G.DiskMass_elements[8]/G.DiskMass_elements[10]) - alog10(S_mf/Fe_mf)
CaFe_mf_disk = alog10(G.DiskMass_elements[9]/G.DiskMass_elements[10]) - alog10(Ca_mf/Fe_mf)
aFe_mf_disk = alog10((G.DiskMass_elements[2] + G.DiskMass_elements[4] + G.DiskMass_elements[5] + G.DiskMass_elements[7] + G.DiskMass_elements[8] + G.DiskMass_elements[9])/G.DiskMass_elements[10]) $
         - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/Fe_mf)                                                                                                                                                                                                                                                                                                                 

HeH_mf_disk = alog10(G.DiskMass_elements[1]/G.DiskMass_elements[0]) - alog10(He_mf/H_mf)
CH_mf_disk = alog10(G.DiskMass_elements[2]/G.DiskMass_elements[0]) - alog10(C_mf/H_mf)
NH_mf_disk = alog10(G.DiskMass_elements[3]/G.DiskMass_elements[0]) - alog10(N_mf/H_mf)
OH_mf_disk = alog10(G.DiskMass_elements[4]/G.DiskMass_elements[0]) - alog10(O_mf/H_mf)
NeH_mf_disk = alog10(G.DiskMass_elements[5]/G.DiskMass_elements[0]) - alog10(Ne_mf/H_mf)
MgH_mf_disk = alog10(G.DiskMass_elements[6]/G.DiskMass_elements[0]) - alog10(Mg_mf/H_mf)
SiH_mf_disk = alog10(G.DiskMass_elements[7]/G.DiskMass_elements[0]) - alog10(Si_mf/H_mf)
SH_mf_disk = alog10(G.DiskMass_elements[8]/G.DiskMass_elements[0]) - alog10(S_mf/H_mf)
CaH_mf_disk = alog10(G.DiskMass_elements[9]/G.DiskMass_elements[0]) - alog10(Ca_mf/H_mf)
FeH_mf_disk = alog10(G.DiskMass_elements[10]/G.DiskMass_elements[0]) - alog10(Fe_mf/H_mf)
aH_mf_disk = alog10((G.DiskMass_elements[2] + G.DiskMass_elements[4] + G.DiskMass_elements[5] + G.DiskMass_elements[7] + G.DiskMass_elements[8] + G.DiskMass_elements[9])/G.DiskMass_elements[0]) $
        - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/H_mf)  

;Bulge abundance ratios ("MASS FRACTION" VERSION):
HFe_mf_bulge = alog10(G.BulgeMass_elements[0]/G.BulgeMass_elements[10]) - alog10(H_mf/Fe_mf)
HeFe_mf_bulge = alog10(G.BulgeMass_elements[1]/G.BulgeMass_elements[10]) - alog10(He_mf/Fe_mf)
CFe_mf_bulge = alog10(G.BulgeMass_elements[2]/G.BulgeMass_elements[10]) - alog10(C_mf/Fe_mf)
NFe_mf_bulge = alog10(G.BulgeMass_elements[3]/G.BulgeMass_elements[10]) - alog10(N_mf/Fe_mf)
OFe_mf_bulge = alog10(G.BulgeMass_elements[4]/G.BulgeMass_elements[10]) - alog10(O_mf/Fe_mf)
NeFe_mf_bulge = alog10(G.BulgeMass_elements[5]/G.BulgeMass_elements[10]) - alog10(Ne_mf/Fe_mf)
MgFe_mf_bulge = alog10(G.BulgeMass_elements[6]/G.BulgeMass_elements[10]) - alog10(Mg_mf/Fe_mf)
SiFe_mf_bulge = alog10(G.BulgeMass_elements[7]/G.BulgeMass_elements[10]) - alog10(Si_mf/Fe_mf)
SFe_mf_bulge = alog10(G.BulgeMass_elements[8]/G.BulgeMass_elements[10]) - alog10(S_mf/Fe_mf)
CaFe_mf_bulge = alog10(G.BulgeMass_elements[9]/G.BulgeMass_elements[10]) - alog10(Ca_mf/Fe_mf)
aFe_mf_bulge = alog10((G.BulgeMass_elements[2] + G.BulgeMass_elements[4] + G.BulgeMass_elements[5] + G.BulgeMass_elements[7] + G.BulgeMass_elements[8] + G.BulgeMass_elements[9])/G.BulgeMass_elements[10]) $
         - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/Fe_mf)                                                                                                                                                                                                                                                                                                                 

HeH_mf_bulge = alog10(G.BulgeMass_elements[1]/G.BulgeMass_elements[0]) - alog10(He_mf/H_mf)
CH_mf_bulge = alog10(G.BulgeMass_elements[2]/G.BulgeMass_elements[0]) - alog10(C_mf/H_mf)
NH_mf_bulge = alog10(G.BulgeMass_elements[3]/G.BulgeMass_elements[0]) - alog10(N_mf/H_mf)
OH_mf_bulge = alog10(G.BulgeMass_elements[4]/G.BulgeMass_elements[0]) - alog10(O_mf/H_mf)
NeH_mf_bulge = alog10(G.BulgeMass_elements[5]/G.BulgeMass_elements[0]) - alog10(Ne_mf/H_mf)
MgH_mf_bulge = alog10(G.BulgeMass_elements[6]/G.BulgeMass_elements[0]) - alog10(Mg_mf/H_mf)
SiH_mf_bulge = alog10(G.BulgeMass_elements[7]/G.BulgeMass_elements[0]) - alog10(Si_mf/H_mf)
SH_mf_bulge = alog10(G.BulgeMass_elements[8]/G.BulgeMass_elements[0]) - alog10(S_mf/H_mf)
CaH_mf_bulge = alog10(G.BulgeMass_elements[9]/G.BulgeMass_elements[0]) - alog10(Ca_mf/H_mf)
FeH_mf_bulge = alog10(G.BulgeMass_elements[10]/G.BulgeMass_elements[0]) - alog10(Fe_mf/H_mf)
aH_mf_bulge = alog10((G.BulgeMass_elements[2] + G.BulgeMass_elements[4] + G.BulgeMass_elements[5] + G.BulgeMass_elements[7] + G.BulgeMass_elements[8] + G.BulgeMass_elements[9])/G.BulgeMass_elements[0]) $
        - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/H_mf)  

;Combined disk + bulge abundance ratios ("MASS FRACTION" VERSION):
HFe_mf_comb = alog10((G.DiskMass_elements[0]+G.BulgeMass_elements[0])/(G.DiskMass_elements[10]+G.BulgeMass_elements[10])) - alog10(H_mf/Fe_mf)
HeFe_mf_comb = alog10((G.DiskMass_elements[1]+G.BulgeMass_elements[1])/(G.DiskMass_elements[10]+G.BulgeMass_elements[10])) - alog10(He_mf/Fe_mf)
CFe_mf_comb = alog10((G.DiskMass_elements[2]+G.BulgeMass_elements[2])/(G.DiskMass_elements[10]+G.BulgeMass_elements[10])) - alog10(C_mf/Fe_mf)
NFe_mf_comb = alog10((G.DiskMass_elements[3]+G.BulgeMass_elements[3])/(G.DiskMass_elements[10]+G.BulgeMass_elements[10])) - alog10(N_mf/Fe_mf)
OFe_mf_comb = alog10((G.DiskMass_elements[4]+G.BulgeMass_elements[4])/(G.DiskMass_elements[10]+G.BulgeMass_elements[10])) - alog10(O_mf/Fe_mf)
NeFe_mf_comb = alog10((G.DiskMass_elements[5]+G.BulgeMass_elements[5])/(G.DiskMass_elements[10]+G.BulgeMass_elements[10])) - alog10(Ne_mf/Fe_mf)
MgFe_mf_comb = alog10((G.DiskMass_elements[6]+G.BulgeMass_elements[6])/(G.DiskMass_elements[10]+G.BulgeMass_elements[10])) - alog10(Mg_mf/Fe_mf)
SiFe_mf_comb = alog10((G.DiskMass_elements[7]+G.BulgeMass_elements[7])/(G.DiskMass_elements[10]+G.BulgeMass_elements[10])) - alog10(Si_mf/Fe_mf)
SFe_mf_comb = alog10((G.DiskMass_elements[8]+G.BulgeMass_elements[8])/(G.DiskMass_elements[10]+G.BulgeMass_elements[10])) - alog10(S_mf/Fe_mf)
CaFe_mf_comb = alog10((G.DiskMass_elements[9]+G.BulgeMass_elements[9])/(G.DiskMass_elements[10]+G.BulgeMass_elements[10])) - alog10(Ca_mf/Fe_mf)                                                                                                                                                                                                                                                                                                            
aFe_mf_comb = alog10((G.BulgeMass_elements[2] + G.BulgeMass_elements[4] + G.BulgeMass_elements[5] + G.BulgeMass_elements[7] + G.BulgeMass_elements[8] + G.BulgeMass_elements[9] $
                      + G.DiskMass_elements[2] + G.DiskMass_elements[4] + G.DiskMass_elements[5] + G.DiskMass_elements[7] + G.DiskMass_elements[8] + G.DiskMass_elements[9]) $
                      /(G.BulgeMass_elements[10]+G.DiskMass_elements[10])) $
            - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/Fe_mf) 

HeH_mf_comb = alog10((G.DiskMass_elements[1]+G.BulgeMass_elements[1])/(G.DiskMass_elements[0]+G.BulgeMass_elements[0])) - alog10(He_mf/H_mf)
CH_mf_comb = alog10((G.DiskMass_elements[2]+G.BulgeMass_elements[2])/(G.DiskMass_elements[0]+G.BulgeMass_elements[0])) - alog10(C_mf/H_mf)
NH_mf_comb = alog10((G.DiskMass_elements[3]+G.BulgeMass_elements[3])/(G.DiskMass_elements[0]+G.BulgeMass_elements[0])) - alog10(N_mf/H_mf)
OH_mf_comb = alog10((G.DiskMass_elements[4]+G.BulgeMass_elements[4])/(G.DiskMass_elements[0]+G.BulgeMass_elements[0])) - alog10(O_mf/H_mf)
NeH_mf_comb = alog10((G.DiskMass_elements[5]+G.BulgeMass_elements[5])/(G.DiskMass_elements[0]+G.BulgeMass_elements[0])) - alog10(Ne_mf/H_mf)
MgH_mf_comb = alog10((G.DiskMass_elements[6]+G.BulgeMass_elements[6])/(G.DiskMass_elements[0]+G.BulgeMass_elements[0])) - alog10(Mg_mf/H_mf)
SiH_mf_comb = alog10((G.DiskMass_elements[7]+G.BulgeMass_elements[7])/(G.DiskMass_elements[0]+G.BulgeMass_elements[0])) - alog10(Si_mf/H_mf)
SH_mf_comb = alog10((G.DiskMass_elements[8]+G.BulgeMass_elements[8])/(G.DiskMass_elements[0]+G.BulgeMass_elements[0])) - alog10(S_mf/H_mf)
CaH_mf_comb = alog10((G.DiskMass_elements[9]+G.BulgeMass_elements[9])/(G.DiskMass_elements[0]+G.BulgeMass_elements[0])) - alog10(Ca_mf/H_mf)
FeH_mf_comb = alog10((G.DiskMass_elements[10]+G.BulgeMass_elements[10])/(G.DiskMass_elements[0]+G.BulgeMass_elements[0])) - alog10(Fe_mf/H_mf)

;Cold Gas abundance ratios ("MASS FRACTION" VERSION):
HFe_mf_cold = alog10(G.ColdGas_elements[0]/G.ColdGas_elements[10]) - alog10(H_mf/Fe_mf)
HeFe_mf_cold = alog10(G.ColdGas_elements[1]/G.ColdGas_elements[10]) - alog10(He_mf/Fe_mf)
CFe_mf_cold = alog10(G.ColdGas_elements[2]/G.ColdGas_elements[10]) - alog10(C_mf/Fe_mf)
NFe_mf_cold = alog10(G.ColdGas_elements[3]/G.ColdGas_elements[10]) - alog10(N_mf/Fe_mf)
OFe_mf_cold = alog10(G.ColdGas_elements[4]/G.ColdGas_elements[10]) - alog10(O_mf/Fe_mf)
NeFe_mf_cold = alog10(G.ColdGas_elements[5]/G.ColdGas_elements[10]) - alog10(Ne_mf/Fe_mf)
MgFe_mf_cold = alog10(G.ColdGas_elements[6]/G.ColdGas_elements[10]) - alog10(Mg_mf/Fe_mf)
SiFe_mf_cold = alog10(G.ColdGas_elements[7]/G.ColdGas_elements[10]) - alog10(Si_mf/Fe_mf)
SFe_mf_cold = alog10(G.ColdGas_elements[8]/G.ColdGas_elements[10]) - alog10(S_mf/Fe_mf)
CaFe_mf_cold = alog10(G.ColdGas_elements[9]/G.ColdGas_elements[10]) - alog10(Ca_mf/Fe_mf)
aFe_mf_cold = alog10((G.ColdGas_elements[2] + G.ColdGas_elements[4] + G.ColdGas_elements[5] + G.ColdGas_elements[7] + G.ColdGas_elements[8] + G.ColdGas_elements[9])/G.ColdGas_elements[10]) $
         - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/Fe_mf)                                                                                                                                                                                                                                                                                                                 

HeH_mf_cold = alog10(G.ColdGas_elements[1]/G.ColdGas_elements[0]) - alog10(He_mf/H_mf)
CH_mf_cold = alog10(G.ColdGas_elements[2]/G.ColdGas_elements[0]) - alog10(C_mf/H_mf)
NH_mf_cold = alog10(G.ColdGas_elements[3]/G.ColdGas_elements[0]) - alog10(N_mf/H_mf)
OH_mf_cold = alog10(G.ColdGas_elements[4]/G.ColdGas_elements[0]) - alog10(O_mf/H_mf)
;NeH_mf_cold = alog10(G.ColdGas_elements[5]/G.ColdGas_elements[0]) - alog10(Ne_mf/H_mf)
MgH_mf_cold = alog10(G.ColdGas_elements[6]/G.ColdGas_elements[0]) - alog10(Mg_mf/H_mf)
SiH_mf_cold = alog10(G.ColdGas_elements[7]/G.ColdGas_elements[0]) - alog10(Si_mf/H_mf)
SH_mf_cold = alog10(G.ColdGas_elements[8]/G.ColdGas_elements[0]) - alog10(S_mf/H_mf)
CaH_mf_cold = alog10(G.ColdGas_elements[9]/G.ColdGas_elements[0]) - alog10(Ca_mf/H_mf)
FeH_mf_cold = alog10(G.ColdGas_elements[10]/G.ColdGas_elements[0]) - alog10(Fe_mf/H_mf)
aH_mf_cold = alog10((G.ColdGas_elements[2] + G.ColdGas_elements[4] + G.ColdGas_elements[5] + G.ColdGas_elements[7] + G.ColdGas_elements[8] + G.ColdGas_elements[9])/G.ColdGas_elements[0]) $
        - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/H_mf)  

;Hot Gas abundance ratios ("MASS FRACTION" VERSION):
HFe_mf_hot = alog10(G.HotGas_elements[0]/G.HotGas_elements[10]) - alog10(H_mf/Fe_mf)
HeFe_mf_hot = alog10(G.HotGas_elements[1]/G.HotGas_elements[10]) - alog10(He_mf/Fe_mf)
CFe_mf_hot = alog10(G.HotGas_elements[2]/G.HotGas_elements[10]) - alog10(C_mf/Fe_mf)
NFe_mf_hot = alog10(G.HotGas_elements[3]/G.HotGas_elements[10]) - alog10(N_mf/Fe_mf)
OFe_mf_hot = alog10(G.HotGas_elements[4]/G.HotGas_elements[10]) - alog10(O_mf/Fe_mf)
NeFe_mf_hot = alog10(G.HotGas_elements[5]/G.HotGas_elements[10]) - alog10(Ne_mf/Fe_mf)
MgFe_mf_hot = alog10(G.HotGas_elements[6]/G.HotGas_elements[10]) - alog10(Mg_mf/Fe_mf)
SiFe_mf_hot = alog10(G.HotGas_elements[7]/G.HotGas_elements[10]) - alog10(Si_mf/Fe_mf)
SFe_mf_hot = alog10(G.HotGas_elements[8]/G.HotGas_elements[10]) - alog10(S_mf/Fe_mf)
CaFe_mf_hot = alog10(G.HotGas_elements[9]/G.HotGas_elements[10]) - alog10(Ca_mf/Fe_mf)
aFe_mf_hot = alog10((G.HotGas_elements[2] + G.HotGas_elements[4] + G.HotGas_elements[5] + G.HotGas_elements[7] + G.HotGas_elements[8] + G.HotGas_elements[9])/G.HotGas_elements[10]) $
         - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/Fe_mf)                                                                                                                                                                                                                                                                                                                 

HeH_mf_hot = alog10(G.HotGas_elements[1]/G.HotGas_elements[0]) - alog10(He_mf/H_mf)
CH_mf_hot = alog10(G.HotGas_elements[2]/G.HotGas_elements[0]) - alog10(C_mf/H_mf)
NH_mf_hot = alog10(G.HotGas_elements[3]/G.HotGas_elements[0]) - alog10(N_mf/H_mf)
OH_mf_hot = alog10(G.HotGas_elements[4]/G.HotGas_elements[0]) - alog10(O_mf/H_mf)
NeH_mf_hot = alog10(G.HotGas_elements[5]/G.HotGas_elements[0]) - alog10(Ne_mf/H_mf)
MgH_mf_hot = alog10(G.HotGas_elements[6]/G.HotGas_elements[0]) - alog10(Mg_mf/H_mf)
SiH_mf_hot = alog10(G.HotGas_elements[7]/G.HotGas_elements[0]) - alog10(Si_mf/H_mf)
SH_mf_hot = alog10(G.HotGas_elements[8]/G.HotGas_elements[0]) - alog10(S_mf/H_mf)
CaH_mf_hot = alog10(G.HotGas_elements[9]/G.HotGas_elements[0]) - alog10(Ca_mf/H_mf)
FeH_mf_hot = alog10(G.HotGas_elements[10]/G.HotGas_elements[0]) - alog10(Fe_mf/H_mf)
aH_mf_hot = alog10((G.HotGas_elements[2] + G.HotGas_elements[4] + G.HotGas_elements[5] + G.HotGas_elements[7] + G.HotGas_elements[8] + G.HotGas_elements[9])/G.HotGas_elements[0]) $
        - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/H_mf) 
ENDIF


;;Bulge abundance ratios ("CONVERTING FROM n_i/n_H TO MASS FRACTION" VERSION):
;HFe_nr = alog10(G.BulgeMass_elements[0]/G.BulgeMass_elements[4]) - alog10((H_nr/Fe_nr)*(H_aw/Fe_aw))
;HeFe_nr = alog10(G.BulgeMass_elements[1]/G.BulgeMass_elements[4]) - alog10((He_nr/Fe_nr)*(He_aw/Fe_aw))
;CFe_nr = alog10(G.BulgeMass_elements[2]/G.BulgeMass_elements[4]) - alog10((C_nr/Fe_nr)*(C_aw/Fe_aw))
;NFe_nr = alog10(G.BulgeMass_elements[3]/G.BulgeMass_elements[4]) - alog10((N_nr/Fe_nr)*(N_aw/Fe_aw))
;OFe_nr = alog10(G.BulgeMass_elements[2]/G.BulgeMass_elements[4]) - alog10((O_nr/Fe_nr)*(O_aw/Fe_aw))
;NeFe_nr = alog10(G.BulgeMass_elements[4]/G.BulgeMass_elements[4]) - alog10((Ne_nr/Fe_nr)*(Ne_aw/Fe_aw))
;MgFe_nr = alog10(G.BulgeMass_elements[3]/G.BulgeMass_elements[4]) - alog10((Mg_nr/Fe_nr)*(Mg_aw/Fe_aw))
;SiFe_nr = alog10(G.BulgeMass_elements[7]/G.BulgeMass_elements[4]) - alog10((Si_nr/Fe_nr)*(Si_aw/Fe_aw))
;SFe_nr = alog10(G.BulgeMass_elements[8]/G.BulgeMass_elements[4]) - alog10((S_nr/Fe_nr)*(S_aw/Fe_aw))
;CaFe_nr = alog10(G.BulgeMass_elements[9]/G.BulgeMass_elements[4]) - alog10((Ca_nr/Fe_nr)*(Ca_aw/Fe_aw))
;aFe_nr = alog10((G.BulgeMass_elements[2] $
;                 + G.BulgeMass_elements[4] $
;                 + G.BulgeMass_elements[3] $
;                 + G.BulgeMass_elements[7] $
;                 + G.BulgeMass_elements[8] $                                                   
;                 + G.BulgeMass_elements[9])/G.BulgeMass_elements[4]) - alog10((O_nr/Fe_nr)*(O_aw/Fe_aw) $
;                                                                               + (Ne_mf/Fe_nr)*(O_aw/Fe_aw) $
;                                                                               + (Mg_mf/Fe_nr)*(Mg_aw/Fe_aw) $
;                                                                               + (Si_mf/Fe_nr)*(Si_aw/Fe_aw) $
;                                                                               + (S_mf/Fe_nr)*(S_aw/Fe_aw) $   
;                                                                               + (Ca_mf/Fe_nr)*(Ca_aw/Fe_aw))  
;
;HeH_nr = alog10(G.BulgeMass_elements[1]/G.BulgeMass_elements[0]) - alog10((He_nr/H_nr)*(He_aw/H_aw))
;CH_nr = alog10(G.BulgeMass_elements[2]/G.BulgeMass_elements[0]) - alog10((C_nr/H_nr)*(C_aw/H_aw))
;NH_nr = alog10(G.BulgeMass_elements[3]/G.BulgeMass_elements[0]) - alog10((N_nr/H_nr)*(N_aw/H_aw))
;OH_nr = alog10(G.BulgeMass_elements[2]/G.BulgeMass_elements[0]) - alog10((O_nr/H_nr)*(O_aw/H_aw))
;NeH_nr = alog10(G.BulgeMass_elements[4]/G.BulgeMass_elements[0]) - alog10((Ne_nr/H_nr)*(Ne_aw/H_aw))
;MgH_nr = alog10(G.BulgeMass_elements[3]/G.BulgeMass_elements[0]) - alog10((Mg_nr/H_nr)*(Mg_aw/H_aw))
;SiH_nr = alog10(G.BulgeMass_elements[7]/G.BulgeMass_elements[0]) - alog10((Si_nr/H_nr)*(Si_aw/H_aw))
;SH_nr = alog10(G.BulgeMass_elements[8]/G.BulgeMass_elements[0]) - alog10((S_nr/H_nr)*(S_aw/H_aw))
;CaH_nr = alog10(G.BulgeMass_elements[9]/G.BulgeMass_elements[0]) - alog10((Ca_nr/H_nr)*(Ca_aw/H_aw))
;FeH_nr = alog10(G.BulgeMass_elements[4]/G.BulgeMass_elements[0]) - alog10((Fe_nr/H_nr)*(Fe_aw/H_aw))
;aH_nr = alog10((G.BulgeMass_elements[2] $
;                 + G.BulgeMass_elements[4] $
;                 + G.BulgeMass_elements[3] $
;                 + G.BulgeMass_elements[7] $
;                 + G.BulgeMass_elements[8] $                                                   
;                 + G.BulgeMass_elements[9])/G.BulgeMass_elements[0]) - alog10((O_nr/H_nr)*(O_aw/H_aw) $
;                                                                               + (Ne_mf/H_nr)*(O_aw/H_aw) $
;                                                                               + (Mg_mf/H_nr)*(Mg_aw/H_aw) $
;                                                                               + (Si_mf/H_nr)*(Si_aw/H_aw) $
;                                                                               + (S_mf/H_nr)*(S_aw/H_aw) $   
;                                                                               + (Ca_mf/H_nr)*(Ca_aw/H_aw))  

;!P.MULTI = [0,2,2,0,0]

;------------------------------------------------------------------
;------------------------------------------------------------------
;DISK (MILKY WAY TYPE HALOES):
;------------------------------------------------------------------
;------------------------------------------------------------------

;-------------------------------
;PLOT [O/Fe]_DiskMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[O/Fe]!Ddisk!N', $
            charsize = csiz, xrange = [MIN(alog10(stars(w1)*1e10)), 11.5], yrange = [MIN(OFe_mf_disk(w0)), 0.5], xstyle = 1, ystyle = 1 ;xrange = [MIN(alog10(stars(w1)*1e10)), 11.5], yrange = [MIN(OFe_mf_disk(w0)), 0.5]

oplot, alog10(stars(w0)*1e10), OFe_mf_disk(w0), psym = 8, symsize = ssiz, color = fsc_color("black")
oplot, alog10(stars(ww0)*1e10), OFe_mf_disk(ww0), psym = 8, symsize = ssiz, color = fsc_color("blue")
oplot, alog10(stars(www0)*1e10), OFe_mf_disk(www0), psym = 8, symsize = ssiz, color = fsc_color("red")
legend, ['In MW haloes', '...and are central discs', '...and have 1.0<SFR<5.0'], box = 0, charsize = 1.5, linestyle = [-99,-99,-99], textcolor=[fsc_color("black"), fsc_color("blue"), fsc_color("red")], /top, /left

;-------------------------------
;PLOT [Mg/Fe]_DiskMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[Mg/Fe]!Ddisk!N', $
            charsize = csiz, xrange = [MIN(alog10(stars(w1)*1e10)), 11.5], yrange = [MIN(MgFe_mf_disk(w0)), 0.5], xstyle = 1, ystyle = 1 ;yrange = [MIN(MgFe_mf_disk(w0)), 0.5]
            
oplot, alog10(stars(w0)*1e10), MgFe_mf_disk(w0), psym = 8, symsize = ssiz, color = fsc_color("black")
oplot, alog10(stars(ww0)*1e10), MgFe_mf_disk(ww0), psym = 8, symsize = ssiz, color = fsc_color("blue")
oplot, alog10(stars(www0)*1e10), MgFe_mf_disk(www0), psym = 8, symsize = ssiz, color = fsc_color("red")
;legend, ['MW haloes'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left
legend, ['In MW haloes', '...and are central discs', '...and have 1.0<SFR<5.0'], box = 0, charsize = 1.5, linestyle = [-99,-99,-99], textcolor=[fsc_color("black"), fsc_color("blue"), fsc_color("red")], /top, /left

;;-------------------------------
;;PLOT [alpha/Fe]_DiskMass vs. M*:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = '[alpha/Fe]!Ddisk!N', $
;            charsize = 1.0, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [0.0,1.0], xstyle = 1, ystyle = 1 ;[MIN(aFe_mf_disk(w0)), MAX(aFe_mf_disk(w0))]
;            
;;oplot, alog10(stars(w0)*1e10), aFe_mf_disk(w0), psym = 8, symsize = ssiz, color = fsc_color("black")
;legend, ['No alpha with MAINELEMENTS'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left

;;-------------------------------
;;PLOT [O/Fe]_DiskMass vs. [Fe/H]_DiskMass:
;;-------------------------------
;www0 = WHERE(FeH_mf_disk(w0) ge 0.0)
;plot, findgen(10), /nodata, xtitle = '[Fe/H]!Ddisk!N', ytitle = '[O/Fe]!Ddisk!N', $
;            charsize = 1.0, xrange = [MIN(FeH_mf_disk(w0(www0))), MAX(FeH_mf_disk(w0(www0)))], yrange = [MIN(OFe_mf_disk(w0(www0))), MAX(OFe_mf_disk(w0(www0)))], xstyle = 1, ystyle = 1 ;xrange = [-1.5, 0.5], yrange = [-0.1, 0.65], 
;;print, "MIN(FeH_mf_disk(w0(ww0))) = ", MIN(FeH_mf_disk(w0(ww0))), "MAX(FeH_mf_disk(w0(ww0))) = ", MAX(FeH_mf_disk(w0(ww0)))      
;oplot, FeH_mf_disk(w0(www0)), OFe_mf_disk(w0(www0)), psym = 8, symsize = ssiz, color = fsc_color("black")
;oplot, FeH_mf_disk(ww0(www0)), OFe_mf_disk(ww0(www0)), psym = 8, symsize = ssiz, color = fsc_color("blue")
;legend, ['MW haloes'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left

;;-------------------------------
;;PLOT [O/Fe]_DiskMass vs. [Fe/H]_DiskMass:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = '[Fe/H]!Ddisk!N', ytitle = '[O/Fe]!Ddisk!N', $
;            charsize = csiz, xrange = [-3.0, 1.0], yrange = [-0.25, 0.6], xstyle = 1, ystyle = 1
;;w0b = w0(WHERE(alog10(G(w0[i]).sfh_ElementsDiskMass[10,*]/G(w0[i]).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf) gt -3.0 AND alog10(G(w0[i]).sfh_ElementsDiskMass[4,*]/G(w0[i]).sfh_ElementsDiskMass[10,*]) - alog10(O_mf/Fe_mf) lt 0.6)) 
;w0b = w0
;FOR i=0, n_elements(w0b)-1 DO BEGIN
;oplot, alog10(G(w0b[i]).sfh_ElementsDiskMass[10,*]/G(w0b[i]).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf), alog10(G(w0b[i]).sfh_ElementsDiskMass[4,*]/G(w0b[i]).sfh_ElementsDiskMass[10,*]) - alog10(O_mf/Fe_mf), psym = 8, symsize = ssiz, color = fsc_color("black")
;ENDFOR
;FOR i=0, n_elements(ww0)-1 DO BEGIN
;;plots, alog10(G(w0[i]).sfh_ElementsDiskMass[10,*]/G(w0[i]).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf), alog10(G(w0[i]).sfh_ElementsDiskMass[4,*]/G(w0[i]).sfh_ElementsDiskMass[10,*]) - alog10(O_mf/Fe_mf), psym = 8, symsize = ssiz, color = fsc_color("black")
;oplot, alog10(G(ww0[i]).sfh_ElementsDiskMass[10,*]/G(ww0[i]).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf), alog10(G(ww0[i]).sfh_ElementsDiskMass[4,*]/G(ww0[i]).sfh_ElementsDiskMass[10,*]) - alog10(O_mf/Fe_mf), psym = 8, symsize = ssiz, color = fsc_color("blue")
;;print, i, n_elements(ww0)
;;plots, alog10(G(www0[i]).sfh_ElementsDiskMass[10,*]/G(www0[i]).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf), alog10(G(www0[i]).sfh_ElementsDiskMass[4,*]/G(www0[i]).sfh_ElementsDiskMass[10,*]) - alog10(O_mf/Fe_mf), psym = 8, symsize = ssiz, color = fsc_color("red")
;ENDFOR
;FOR i=0, n_elements(www0)-1 DO BEGIN
;oplot, alog10(G(www0[i]).sfh_ElementsDiskMass[10,*]/G(www0[i]).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf), alog10(G(www0[i]).sfh_ElementsDiskMass[4,*]/G(www0[i]).sfh_ElementsDiskMass[10,*]) - alog10(O_mf/Fe_mf), psym = 8, symsize = ssiz, color = fsc_color("red")
;ENDFOR
;legend, ['At z=0: In MW haloes', '...and are central discs', '...and have 1.0<SFR<5.0'], box = 0, charsize = 1.5, linestyle = [-99,-99,-99], textcolor=[fsc_color("black"), fsc_color("blue"), fsc_color("red")], /top, /right
;;print, alog10(G(www0[0]).sfh_ElementsDiskMass[10,*]/G(www0[0]).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf)
;print, "............... ", n_elements(alog10(G(ww0[0]).sfh_ElementsDiskMass[10,*]/G(ww0[0]).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf))
;print, "........ ", n_elements(ww0), n_elements(www0)

;-------------------------------
;PLOT [O/Fe]_DiskMass vs. [Fe/H]_DiskMass (LOW RES):
;-------------------------------
plot, findgen(10), /nodata, xtitle = '[Fe/H]!Ddisc!N', ytitle = '[O/Fe]!Ddisc!N', $
      xrange = [-2.0, 0.5], yrange = [-0.1, 0.5], xstyle = 1, ystyle = 1, CHARSIZE=1.75, CHARTHICK=4., XTHICK=5., YTHICK=5. ;charsize = csiz

;Plottables:
;TheFeH = fltarr(n_elements(ww0),20)
;TheOFe = fltarr(n_elements(ww0),20)
TheFeH = fltarr(n_elements(ww0)*20)
TheOFe = fltarr(n_elements(ww0)*20)
k=0
FOR i=0, n_elements(ww0)-1 DO BEGIN
  FOR j=0, 20-1 DO BEGIN
    ;TheFeH[i,j] = alog10(G(ww0[i]).sfh_ElementsDiskMass[10,j]/G(ww0[i]).sfh_ElementsDiskMass[0,j]) - alog10(Fe_mf/H_mf)
    ;TheOFe[i,j] = alog10(G(ww0[i]).sfh_ElementsDiskMass[4,j]/G(ww0[i]).sfh_ElementsDiskMass[10,j]) - alog10(O_mf/Fe_mf)
    TheFeH[k] = alog10(G(ww0[i]).sfh_ElementsDiskMass[10,j]/G(ww0[i]).sfh_ElementsDiskMass[0,j]) - alog10(Fe_mf/H_mf)
    TheOFe[k] = alog10(G(ww0[i]).sfh_ElementsDiskMass[4,j]/G(ww0[i]).sfh_ElementsDiskMass[10,j]) - alog10(O_mf/Fe_mf)
    k = k+1
  ENDFOR
ENDFOR

loadct, 13

;Parameters:
FeHmin = -3.0 
FeHmax = 1.0
;FeHbinwidth = 4./20.
;FeHbinno = long((FeHmax-FeHmin)*FeHbinwidth)
FeHbinno = 75.
FeHbinwidth = (FeHmax-FeHmin) / FeHbinno
OFemin = -0.3
OFemax = 0.7
;OFebinwidth = 4./20.
;OFebinno = long((OFemax-OFemin)*OFebinwidth)
OFebinno = 75.
OFebinwidth = (OFemax-OFemin) / OFebinno

;Bin galaxies by [Fe/H] and [O/Fe]:
AveFeH = fltarr(FeHbinno*OFebinno)
AveOFe = fltarr(FeHbinno*OFebinno)
count = intarr(FeHbinno*OFebinno)
deleg = intarr(n_elements(ww0)*20)
k=0
FOR i=FeHmin, FeHmax-FeHbinwidth, FeHbinwidth DO BEGIN
    FOR j=OFemin, OFemax-OFebinwidth, OFebinwidth DO BEGIN
        w99 = WHERE(TheFeH ge i and TheFeH lt i+FeHbinwidth and TheOFe ge j and TheOFe lt j+OFebinwidth)
        IF(w99[0] ne -1) THEN BEGIN
            IF(n_elements(w99) ge 5) THEN deleg(w99) = 1
            ;deleg(w99) = 1
            AveFeH[k] = i + (0.5*FeHbinwidth)
            AveOFe[k] = j + (0.5*OFebinwidth)
            count[k] = n_elements(w99)
        ENDIF
        k+=1
    ENDFOR
ENDFOR

ww2 = WHERE(count ge 5)

;Make plotting symbol:
usersym, [-1.1,1.1,1.1,-1.1], [1.1,1.1,-1.1,-1.1], /fill
;usersym, [-0.85,0.85,0.85,-0.85], [0.85,0.85,-0.85,-0.85], /fill
;usersym, [-0.5*FeHbinwidth,0.5*FeHbinwidth,0.5*FeHbinwidth,-0.5*FeHbinwidth], [0.5*OFebinwidth,0.5*OFebinwidth,-0.5*OFebinwidth,-0.5*OFebinwidth], /fill
ctload, 0, /REVERSE, CLIP=[50,200] ;CLIP=[50,225]

;print, "n_elements(deleg(WHERE (deleg ne 1))) = ", n_elements(deleg(WHERE (deleg ne 1)))
;print, "n_elements(deleg(WHERE (deleg eq 0))) = ", n_elements(deleg(WHERE (deleg eq 0)))
plots, AveFeH(ww2), AveOFe(ww2), psym=8, color=bytscl(count(ww2))

;goodcount = count(WHERE(count gt 0))
;print, goodcount[0:99] 

loadct, 13
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill

ww99 = WHERE(deleg eq 0)
;print, "n_elements(ww99) = ", n_elements(ww99)
oplot, TheFeH[ww99], TheOFe[ww99], psym = 8, symsize = ssiz, color=fsc_color("grey")

;Plot solar lines:
oplot, [-4.0, 2.0], [0.0, 0.0], thick=6, linestyle = 2, color=fsc_color("black")
oplot, [0.0, 0.0], [-0.5, 0.8], thick=6, linestyle = 2, color=fsc_color("black")

;Title:
xyouts, -1.8, 0.45, "MW-type galaxies", charsize=2., charthick=4.

;Plot colour table
     xran = [-1.8,-0.5] ;[-1.7,-1.8]
     yran = [0.05,0.1] ;[-0.07,0.0]
     delta_yran = yran[1] - yran[0]
     delta_xran = xran[1] - xran[0]
     xarr = findgen(100)/100. * (delta_xran+0.01) + xran[0] ;* 0.50 + delta_xran * 0.25 + xran[0]
     yarr = findgen(2)/2. * delta_yran + yran[0] ;* 0.12 + delta_yran * 0.805 + yran[0]
     zarr = FLTARR([(SIZE(xarr))[1], (SIZE(yarr))[1]])
     zarr[*,0] = findgen(100)/100.
     zarr[*,1] = findgen(100)/100.  
     loadct, 0, /silent
     contour, zarr, xarr, yarr, /fill, /over, levels=findgen(175) / 175., c_colors=reverse(indgen(175)) ;levels=findgen(256) / 256., c_colors=reverse(indgen(256))
     plots, [xran(0), xran(0), xran(1), xran(1), xran(0)], [yran(0), (delta_yran/2.)+yran(0), (delta_yran/2.)+yran(0), yran(0), yran(0)], thick=2.
print, "MIN(count(ww2)) = ", MIN(count(ww2))
print, "MEAN(count(ww2)) = ", MEAN(count(ww2))
print, "MEDIAN(count(ww2)) = ", MEDIAN(count(ww2))
print, "MID(count(ww2)) = ", ((MAX(count(ww2))-MIN(count(ww2)))/2.0)+MIN(count(ww2))
print, "MAX(count(ww2)) = ", MAX(count(ww2))
print, "  "
     xyouts, (xran[1]-xran[0])/2. + xran[0]-0.12, 0.08, 'Number', charsize=1., charthick=4.  
     xyouts, xran[0]-0.03, 0.03, strcompress(STRING(FIX(MIN(count(ww2)))), /remove_all), charsize=1., charthick=3.
     xyouts, (xran[1]-xran[0])/2. + xran[0] - 0.06, 0.03, strcompress(STRING(FIX(((MAX(count(ww2))-MIN(count(ww2)))/2.0)+MIN(count(ww2)))), /remove_all), charsize=1., charthick=3.
     xyouts, xran[1]-0.07, 0.03, strcompress(STRING(FIX(MAX(count(ww2)))), /remove_all), charsize=1., charthick=3.
     plots, [xran[0],xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [(xran[1]-xran[0])/2. + xran[0],(xran[1]-xran[0])/2. + xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [xran[1],xran[1]], [yran[0], yran[0]-0.006], thick=2.

;-------------------------------
;PLOT [Mg/Fe]_DiskMass vs. [Fe/H]_DiskMass (LOW RES):
;-------------------------------
plot, findgen(10), /nodata, xtitle = '[Fe/H]!Ddisc!N', ytitle = '[Mg/Fe]!Ddisc!N', $
      xrange = [-2.0, 0.5], yrange = [-0.1, 0.5], xstyle = 1, ystyle = 1, CHARSIZE=1.75, CHARTHICK=4., XTHICK=5., YTHICK=5. ;charsize = csiz

;Plottables:
TheFeH = fltarr(n_elements(ww0)*20)
TheMgFe = fltarr(n_elements(ww0)*20)
k=0
FOR i=0, n_elements(ww0)-1 DO BEGIN
  FOR j=0, 20-1 DO BEGIN
    TheFeH[k] = alog10(G(ww0[i]).sfh_ElementsDiskMass[10,j]/G(ww0[i]).sfh_ElementsDiskMass[0,j]) - alog10(Fe_mf/H_mf)
    TheMgFe[k] = alog10(G(ww0[i]).sfh_ElementsDiskMass[6,j]/G(ww0[i]).sfh_ElementsDiskMass[10,j]) - alog10(Mg_mf/Fe_mf)
    k = k+1
  ENDFOR
ENDFOR

loadct, 13

;Parameters:
FeHmin = -3.0 
FeHmax = 1.0
FeHbinno = 75.
FeHbinwidth = (FeHmax-FeHmin) / FeHbinno
MgFemin = -0.3
MgFemax = 0.7
MgFebinno = 75.
MgFebinwidth = (MgFemax-MgFemin) / MgFebinno

;Bin galaxies by [Fe/H] and [Mg/Fe]:
AveFeH = fltarr(FeHbinno*MgFebinno)
AveMgFe = fltarr(FeHbinno*MgFebinno)
count = intarr(FeHbinno*MgFebinno)
deleg = intarr(n_elements(ww0)*20)
k=0
FOR i=FeHmin, FeHmax-FeHbinwidth, FeHbinwidth DO BEGIN
    FOR j=MgFemin, MgFemax-MgFebinwidth, MgFebinwidth DO BEGIN
        w99 = WHERE(TheFeH ge i and TheFeH lt i+FeHbinwidth and TheMgFe ge j and TheMgFe lt j+MgFebinwidth)
        IF(w99[0] ne -1) THEN BEGIN
            IF(n_elements(w99) ge 5) THEN deleg(w99) = 1
            ;deleg(w99) = 1
            AveFeH[k] = i + (0.5*FeHbinwidth)
            AveMgFe[k] = j + (0.5*MgFebinwidth)
            count[k] = n_elements(w99)
        ENDIF
        k+=1
    ENDFOR
ENDFOR

ww2 = WHERE(count ge 5)

;Make plotting symbol:
usersym, [-1.1,1.1,1.1,-1.1], [1.1,1.1,-1.1,-1.1], /fill
ctload, 0, /REVERSE, CLIP=[50,200] ;CLIP=[50,225]

plots, AveFeH(ww2), AveMgFe(ww2), psym=8, color=bytscl(count(ww2))

loadct, 13
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill

ww99 = WHERE(deleg eq 0)
oplot, TheFeH[ww99], TheMgFe[ww99], psym = 8, symsize = ssiz, color=fsc_color("grey")

;Plot solar lines:
oplot, [-4.0, 2.0], [0.0, 0.0], thick=6, linestyle = 2, color=fsc_color("black")
oplot, [0.0, 0.0], [-0.5, 0.8], thick=6, linestyle = 2, color=fsc_color("black")

;Title:
xyouts, -1.8, 0.45, "MW-type galaxies", charsize=2., charthick=4.

;Plot colour table
     xran = [-1.8,-0.5] ;[-1.7,-1.8]
     yran = [0.05,0.1] ;[-0.07,0.0]
     delta_yran = yran[1] - yran[0]
     delta_xran = xran[1] - xran[0]
     xarr = findgen(100)/100. * (delta_xran+0.01) + xran[0] ;* 0.50 + delta_xran * 0.25 + xran[0]
     yarr = findgen(2)/2. * delta_yran + yran[0] ;* 0.12 + delta_yran * 0.805 + yran[0]
     zarr = FLTARR([(SIZE(xarr))[1], (SIZE(yarr))[1]])
     zarr[*,0] = findgen(100)/100.
     zarr[*,1] = findgen(100)/100.  
     loadct, 0, /silent
     contour, zarr, xarr, yarr, /fill, /over, levels=findgen(175) / 175., c_colors=reverse(indgen(175)) ;levels=findgen(256) / 256., c_colors=reverse(indgen(256))
     plots, [xran(0), xran(0), xran(1), xran(1), xran(0)], [yran(0), (delta_yran/2.)+yran(0), (delta_yran/2.)+yran(0), yran(0), yran(0)], thick=2.
print, "MIN(count(ww2)) = ", MIN(count(ww2))
print, "MEAN(count(ww2)) = ", MEAN(count(ww2))
print, "MEDIAN(count(ww2)) = ", MEDIAN(count(ww2))
print, "MID(count(ww2)) = ", ((MAX(count(ww2))-MIN(count(ww2)))/2.0)+MIN(count(ww2))
print, "MAX(count(ww2)) = ", MAX(count(ww2))
print, "  "
     xyouts, (xran[1]-xran[0])/2. + xran[0]-0.12, 0.08, 'Number', charsize=1., charthick=4.  
     xyouts, xran[0]-0.03, 0.03, strcompress(STRING(FIX(MIN(count(ww2)))), /remove_all), charsize=1., charthick=3.
     xyouts, (xran[1]-xran[0])/2. + xran[0] - 0.06, 0.03, strcompress(STRING(FIX(((MAX(count(ww2))-MIN(count(ww2)))/2.0)+MIN(count(ww2)))), /remove_all), charsize=1., charthick=3.
     xyouts, xran[1]-0.07, 0.03, strcompress(STRING(FIX(MAX(count(ww2)))), /remove_all), charsize=1., charthick=3.
     plots, [xran[0],xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [(xran[1]-xran[0])/2. + xran[0],(xran[1]-xran[0])/2. + xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [xran[1],xran[1]], [yran[0], yran[0]-0.006], thick=2.

;-------------------------------
;PLOT [Si/Fe]_DiskMass vs. [Fe/H]_DiskMass (LOW RES):
;-------------------------------
plot, findgen(10), /nodata, xtitle = '[Fe/H]!Ddisc!N', ytitle = '[Si/Fe]!Ddisc!N', $
      xrange = [-2.0, 0.5], yrange = [-0.1, 0.5], xstyle = 1, ystyle = 1, CHARSIZE=1.75, CHARTHICK=4., XTHICK=5., YTHICK=5. ;charsize = csiz

;Plottables:
TheFeH = fltarr(n_elements(ww0)*20)
TheSiFe = fltarr(n_elements(ww0)*20)
k=0
FOR i=0, n_elements(ww0)-1 DO BEGIN
  FOR j=0, 20-1 DO BEGIN
    TheFeH[k] = alog10(G(ww0[i]).sfh_ElementsDiskMass[10,j]/G(ww0[i]).sfh_ElementsDiskMass[0,j]) - alog10(Fe_mf/H_mf)
    TheSiFe[k] = alog10(G(ww0[i]).sfh_ElementsDiskMass[7,j]/G(ww0[i]).sfh_ElementsDiskMass[10,j]) - alog10(Si_mf/Fe_mf)
    k = k+1
  ENDFOR
ENDFOR

loadct, 13

;Parameters:
FeHmin = -3.0 
FeHmax = 1.0
FeHbinno = 75.
FeHbinwidth = (FeHmax-FeHmin) / FeHbinno
SiFemin = -0.3
SiFemax = 0.7
SiFebinno = 75.
SiFebinwidth = (SiFemax-SiFemin) / SiFebinno

;Bin galaxies by [Fe/H] and [Si/Fe]:
AveFeH = fltarr(FeHbinno*SiFebinno)
AveSiFe = fltarr(FeHbinno*SiFebinno)
count = intarr(FeHbinno*SiFebinno)
deleg = intarr(n_elements(ww0)*20)
k=0
FOR i=FeHmin, FeHmax-FeHbinwidth, FeHbinwidth DO BEGIN
    FOR j=SiFemin, SiFemax-SiFebinwidth, SiFebinwidth DO BEGIN
        w99 = WHERE(TheFeH ge i and TheFeH lt i+FeHbinwidth and TheSiFe ge j and TheSiFe lt j+SiFebinwidth)
        IF(w99[0] ne -1) THEN BEGIN
            IF(n_elements(w99) ge 5) THEN deleg(w99) = 1
            ;deleg(w99) = 1
            AveFeH[k] = i + (0.5*FeHbinwidth)
            AveSiFe[k] = j + (0.5*SiFebinwidth)
            count[k] = n_elements(w99)
        ENDIF
        k+=1
    ENDFOR
ENDFOR

ww2 = WHERE(count ge 5)

;Make plotting symbol:
usersym, [-1.1,1.1,1.1,-1.1], [1.1,1.1,-1.1,-1.1], /fill
ctload, 0, /REVERSE, CLIP=[50,200] ;CLIP=[50,225]

plots, AveFeH(ww2), AveSiFe(ww2), psym=8, color=bytscl(count(ww2))

loadct, 13
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill

ww99 = WHERE(deleg eq 0)
oplot, TheFeH[ww99], TheSiFe[ww99], psym = 8, symsize = ssiz, color=fsc_color("grey")

;Plot solar lines:
oplot, [-4.0, 2.0], [0.0, 0.0], thick=6, linestyle = 2, color=fsc_color("black")
oplot, [0.0, 0.0], [-0.5, 0.8], thick=6, linestyle = 2, color=fsc_color("black")

;Title:
xyouts, -1.8, 0.45, "MW-type galaxies", charsize=2., charthick=4.

;Plot colour table
     xran = [-1.8,-0.5] ;[-1.7,-1.8]
     yran = [0.05,0.1] ;[-0.07,0.0]
     delta_yran = yran[1] - yran[0]
     delta_xran = xran[1] - xran[0]
     xarr = findgen(100)/100. * (delta_xran+0.01) + xran[0] ;* 0.50 + delta_xran * 0.25 + xran[0]
     yarr = findgen(2)/2. * delta_yran + yran[0] ;* 0.12 + delta_yran * 0.805 + yran[0]
     zarr = FLTARR([(SIZE(xarr))[1], (SIZE(yarr))[1]])
     zarr[*,0] = findgen(100)/100.
     zarr[*,1] = findgen(100)/100.  
     loadct, 0, /silent
     contour, zarr, xarr, yarr, /fill, /over, levels=findgen(175) / 175., c_colors=reverse(indgen(175)) ;levels=findgen(256) / 256., c_colors=reverse(indgen(256))
     plots, [xran(0), xran(0), xran(1), xran(1), xran(0)], [yran(0), (delta_yran/2.)+yran(0), (delta_yran/2.)+yran(0), yran(0), yran(0)], thick=2.
print, "MIN(count(ww2)) = ", MIN(count(ww2))
print, "MEAN(count(ww2)) = ", MEAN(count(ww2))
print, "MEDIAN(count(ww2)) = ", MEDIAN(count(ww2))
print, "MID(count(ww2)) = ", ((MAX(count(ww2))-MIN(count(ww2)))/2.0)+MIN(count(ww2))
print, "MAX(count(ww2)) = ", MAX(count(ww2))
print, "  "
     xyouts, (xran[1]-xran[0])/2. + xran[0]-0.12, 0.08, 'Number', charsize=1., charthick=4.  
     xyouts, xran[0]-0.03, 0.03, strcompress(STRING(FIX(MIN(count(ww2)))), /remove_all), charsize=1., charthick=3.
     xyouts, (xran[1]-xran[0])/2. + xran[0] - 0.06, 0.03, strcompress(STRING(FIX(((MAX(count(ww2))-MIN(count(ww2)))/2.0)+MIN(count(ww2)))), /remove_all), charsize=1., charthick=3.
     xyouts, xran[1]-0.07, 0.03, strcompress(STRING(FIX(MAX(count(ww2)))), /remove_all), charsize=1., charthick=3.
     plots, [xran[0],xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [(xran[1]-xran[0])/2. + xran[0],(xran[1]-xran[0])/2. + xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [xran[1],xran[1]], [yran[0], yran[0]-0.006], thick=2.

;-------------------------------
;PLOT [S/Fe]_DiskMass vs. [Fe/H]_DiskMass (LOW RES):
;-------------------------------
plot, findgen(10), /nodata, xtitle = '[Fe/H]!Ddisc!N', ytitle = '[S/Fe]!Ddisc!N', $
      xrange = [-2.0, 0.5], yrange = [-0.1, 0.5], xstyle = 1, ystyle = 1, CHARSIZE=1.75, CHARTHICK=4., XTHICK=5., YTHICK=5. ;charsize = csiz

;Plottables:
TheFeH = fltarr(n_elements(ww0)*20)
TheSFe = fltarr(n_elements(ww0)*20)
k=0
FOR i=0, n_elements(ww0)-1 DO BEGIN
  FOR j=0, 20-1 DO BEGIN
    TheFeH[k] = alog10(G(ww0[i]).sfh_ElementsDiskMass[10,j]/G(ww0[i]).sfh_ElementsDiskMass[0,j]) - alog10(Fe_mf/H_mf)
    TheSFe[k] = alog10(G(ww0[i]).sfh_ElementsDiskMass[8,j]/G(ww0[i]).sfh_ElementsDiskMass[10,j]) - alog10(S_mf/Fe_mf)
    k = k+1
  ENDFOR
ENDFOR

loadct, 13

;Parameters:
FeHmin = -3.0 
FeHmax = 1.0
FeHbinno = 75.
FeHbinwidth = (FeHmax-FeHmin) / FeHbinno
SFemin = -0.3
SFemax = 0.7
SFebinno = 75.
SFebinwidth = (SFemax-SFemin) / SFebinno

;Bin galaxies by [Fe/H] and [S/Fe]:
AveFeH = fltarr(FeHbinno*SFebinno)
AveSFe = fltarr(FeHbinno*SFebinno)
count = intarr(FeHbinno*SFebinno)
deleg = intarr(n_elements(ww0)*20)
k=0
FOR i=FeHmin, FeHmax-FeHbinwidth, FeHbinwidth DO BEGIN
    FOR j=SFemin, SFemax-SFebinwidth, SFebinwidth DO BEGIN
        w99 = WHERE(TheFeH ge i and TheFeH lt i+FeHbinwidth and TheSFe ge j and TheSFe lt j+SFebinwidth)
        IF(w99[0] ne -1) THEN BEGIN
            IF(n_elements(w99) ge 5) THEN deleg(w99) = 1
            ;deleg(w99) = 1
            AveFeH[k] = i + (0.5*FeHbinwidth)
            AveSFe[k] = j + (0.5*SFebinwidth)
            count[k] = n_elements(w99)
        ENDIF
        k+=1
    ENDFOR
ENDFOR

ww2 = WHERE(count ge 5)

;Make plotting symbol:
usersym, [-1.1,1.1,1.1,-1.1], [1.1,1.1,-1.1,-1.1], /fill
ctload, 0, /REVERSE, CLIP=[50,200] ;CLIP=[50,225]

plots, AveFeH(ww2), AveSFe(ww2), psym=8, color=bytscl(count(ww2))

loadct, 13
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill

ww99 = WHERE(deleg eq 0)
oplot, TheFeH[ww99], TheSFe[ww99], psym = 8, symsize = ssiz, color=fsc_color("grey")

;Plot solar lines:
oplot, [-4.0, 2.0], [0.0, 0.0], thick=6, linestyle = 2, color=fsc_color("black")
oplot, [0.0, 0.0], [-0.5, 0.8], thick=6, linestyle = 2, color=fsc_color("black")

;Title:
xyouts, -1.8, 0.45, "MW-type galaxies", charsize=2., charthick=4.

;Plot colour table
     xran = [-1.8,-0.5] ;[-1.7,-1.8]
     yran = [0.05,0.1] ;[-0.07,0.0]
     delta_yran = yran[1] - yran[0]
     delta_xran = xran[1] - xran[0]
     xarr = findgen(100)/100. * (delta_xran+0.01) + xran[0] ;* 0.50 + delta_xran * 0.25 + xran[0]
     yarr = findgen(2)/2. * delta_yran + yran[0] ;* 0.12 + delta_yran * 0.805 + yran[0]
     zarr = FLTARR([(SIZE(xarr))[1], (SIZE(yarr))[1]])
     zarr[*,0] = findgen(100)/100.
     zarr[*,1] = findgen(100)/100.  
     loadct, 0, /silent
     contour, zarr, xarr, yarr, /fill, /over, levels=findgen(175) / 175., c_colors=reverse(indgen(175)) ;levels=findgen(256) / 256., c_colors=reverse(indgen(256))
     plots, [xran(0), xran(0), xran(1), xran(1), xran(0)], [yran(0), (delta_yran/2.)+yran(0), (delta_yran/2.)+yran(0), yran(0), yran(0)], thick=2.
print, "MIN(count(ww2)) = ", MIN(count(ww2))
print, "MEAN(count(ww2)) = ", MEAN(count(ww2))
print, "MEDIAN(count(ww2)) = ", MEDIAN(count(ww2))
print, "MID(count(ww2)) = ", ((MAX(count(ww2))-MIN(count(ww2)))/2.0)+MIN(count(ww2))
print, "MAX(count(ww2)) = ", MAX(count(ww2))
print, "  "
     xyouts, (xran[1]-xran[0])/2. + xran[0]-0.12, 0.08, 'Number', charsize=1., charthick=4.  
     xyouts, xran[0]-0.03, 0.03, strcompress(STRING(FIX(MIN(count(ww2)))), /remove_all), charsize=1., charthick=3.
     xyouts, (xran[1]-xran[0])/2. + xran[0] - 0.06, 0.03, strcompress(STRING(FIX(((MAX(count(ww2))-MIN(count(ww2)))/2.0)+MIN(count(ww2)))), /remove_all), charsize=1., charthick=3.
     xyouts, xran[1]-0.07, 0.03, strcompress(STRING(FIX(MAX(count(ww2)))), /remove_all), charsize=1., charthick=3.
     plots, [xran[0],xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [(xran[1]-xran[0])/2. + xran[0],(xran[1]-xran[0])/2. + xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [xran[1],xran[1]], [yran[0], yran[0]-0.006], thick=2.

;-------------------------------
;PLOT [Ca/Fe]_DiskMass vs. [Fe/H]_DiskMass (LOW RES):
;-------------------------------
plot, findgen(10), /nodata, xtitle = '[Fe/H]!Ddisc!N', ytitle = '[Ca/Fe]!Ddisc!N', $
      xrange = [-2.0, 0.5], yrange = [-0.1, 0.5], xstyle = 1, ystyle = 1, CHARSIZE=1.75, CHARTHICK=4., XTHICK=5., YTHICK=5. ;charsize = csiz

;Plottables:
TheFeH = fltarr(n_elements(ww0)*20)
TheCaFe = fltarr(n_elements(ww0)*20)
k=0
FOR i=0, n_elements(ww0)-1 DO BEGIN
  FOR j=0, 20-1 DO BEGIN
    TheFeH[k] = alog10(G(ww0[i]).sfh_ElementsDiskMass[10,j]/G(ww0[i]).sfh_ElementsDiskMass[0,j]) - alog10(Fe_mf/H_mf)
    TheCaFe[k] = alog10(G(ww0[i]).sfh_ElementsDiskMass[9,j]/G(ww0[i]).sfh_ElementsDiskMass[10,j]) - alog10(Ca_mf/Fe_mf)
    k = k+1
  ENDFOR
ENDFOR

loadct, 13

;Parameters:
FeHmin = -3.0 
FeHmax = 1.0
FeHbinno = 75.
FeHbinwidth = (FeHmax-FeHmin) / FeHbinno
CaFemin = -0.3
CaFemax = 0.7
CaFebinno = 75.
CaFebinwidth = (CaFemax-CaFemin) / CaFebinno

;Bin galaxies by [Fe/H] and [Ca/Fe]:
AveFeH = fltarr(FeHbinno*CaFebinno)
AveCaFe = fltarr(FeHbinno*CaFebinno)
count = intarr(FeHbinno*CaFebinno)
deleg = intarr(n_elements(ww0)*20)
k=0
FOR i=FeHmin, FeHmax-FeHbinwidth, FeHbinwidth DO BEGIN
    FOR j=CaFemin, CaFemax-CaFebinwidth, CaFebinwidth DO BEGIN
        w99 = WHERE(TheFeH ge i and TheFeH lt i+FeHbinwidth and TheCaFe ge j and TheCaFe lt j+CaFebinwidth)
        IF(w99[0] ne -1) THEN BEGIN
            IF(n_elements(w99) ge 5) THEN deleg(w99) = 1
            ;deleg(w99) = 1
            AveFeH[k] = i + (0.5*FeHbinwidth)
            AveCaFe[k] = j + (0.5*CaFebinwidth)
            count[k] = n_elements(w99)
        ENDIF
        k+=1
    ENDFOR
ENDFOR

ww2 = WHERE(count ge 5)

;Make plotting symbol:
usersym, [-1.1,1.1,1.1,-1.1], [1.1,1.1,-1.1,-1.1], /fill
ctload, 0, /REVERSE, CLIP=[50,200] ;CLIP=[50,225]

plots, AveFeH(ww2), AveCaFe(ww2), psym=8, color=bytscl(count(ww2))

loadct, 13
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill

ww99 = WHERE(deleg eq 0)
oplot, TheFeH[ww99], TheCaFe[ww99], psym = 8, symsize = ssiz, color=fsc_color("grey")

;Plot solar lines:
oplot, [-4.0, 2.0], [0.0, 0.0], thick=6, linestyle = 2, color=fsc_color("black")
oplot, [0.0, 0.0], [-0.5, 0.8], thick=6, linestyle = 2, color=fsc_color("black")

;Title:
xyouts, -1.8, 0.45, "MW-type galaxies", charsize=2., charthick=4.

;Plot colour table
     xran = [-1.8,-0.5] ;[-1.7,-1.8]
     yran = [0.05,0.1] ;[-0.07,0.0]
     delta_yran = yran[1] - yran[0]
     delta_xran = xran[1] - xran[0]
     xarr = findgen(100)/100. * (delta_xran+0.01) + xran[0] ;* 0.50 + delta_xran * 0.25 + xran[0]
     yarr = findgen(2)/2. * delta_yran + yran[0] ;* 0.12 + delta_yran * 0.805 + yran[0]
     zarr = FLTARR([(SIZE(xarr))[1], (SIZE(yarr))[1]])
     zarr[*,0] = findgen(100)/100.
     zarr[*,1] = findgen(100)/100.  
     loadct, 0, /silent
     contour, zarr, xarr, yarr, /fill, /over, levels=findgen(175) / 175., c_colors=reverse(indgen(175)) ;levels=findgen(256) / 256., c_colors=reverse(indgen(256))
     plots, [xran(0), xran(0), xran(1), xran(1), xran(0)], [yran(0), (delta_yran/2.)+yran(0), (delta_yran/2.)+yran(0), yran(0), yran(0)], thick=2.
print, "MIN(count(ww2)) = ", MIN(count(ww2))
print, "MEAN(count(ww2)) = ", MEAN(count(ww2))
print, "MEDIAN(count(ww2)) = ", MEDIAN(count(ww2))
print, "MID(count(ww2)) = ", ((MAX(count(ww2))-MIN(count(ww2)))/2.0)+MIN(count(ww2))
print, "MAX(count(ww2)) = ", MAX(count(ww2))
print, "  "
     xyouts, (xran[1]-xran[0])/2. + xran[0]-0.12, 0.08, 'Number', charsize=1., charthick=4.  
     xyouts, xran[0]-0.03, 0.03, strcompress(STRING(FIX(MIN(count(ww2)))), /remove_all), charsize=1., charthick=3.
     xyouts, (xran[1]-xran[0])/2. + xran[0] - 0.06, 0.03, strcompress(STRING(FIX(((MAX(count(ww2))-MIN(count(ww2)))/2.0)+MIN(count(ww2)))), /remove_all), charsize=1., charthick=3.
     xyouts, xran[1]-0.07, 0.03, strcompress(STRING(FIX(MAX(count(ww2)))), /remove_all), charsize=1., charthick=3.
     plots, [xran[0],xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [(xran[1]-xran[0])/2. + xran[0],(xran[1]-xran[0])/2. + xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [xran[1],xran[1]], [yran[0], yran[0]-0.006], thick=2.

;;-------------------------------
;;PLOT [a/Fe]_DiskMass vs. [Fe/H]_DiskMass:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = '[Fe/H]!Ddisk!N', ytitle = '[alpha/Fe]!Ddisk!N', $
;            charsize = csiz, xrange = [-3.0, 1.0], yrange = [-0.25, 0.6], xstyle = 1, ystyle = 1
;
;;aFew0 = alog10((G(w0).sfh_ElementsDiskMass[2,*] + G(w0).sfh_ElementsDiskMass[4,*] + G(w0).sfh_ElementsDiskMass[5,*] + G(w0).sfh_ElementsDiskMass[7,*] + G(w0).sfh_ElementsDiskMass[8,*] + G(w0).sfh_ElementsDiskMass[9,*])/G(w0).sfh_ElementsDiskMass[10,*]) $
;;         - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/Fe_mf)       
;;aFeww0 = alog10((G(ww0).sfh_ElementsDiskMass[2,*] + G(ww0).sfh_ElementsDiskMass[4,*] + G(ww0).sfh_ElementsDiskMass[5,*] + G(ww0).sfh_ElementsDiskMass[7,*] + G(ww0).sfh_ElementsDiskMass[8,*] + G(ww0).sfh_ElementsDiskMass[9,*])/G(ww0).sfh_ElementsDiskMass[10,*]) $
;;         - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/Fe_mf) 
;;aFewww0 = alog10((G(www0).sfh_ElementsDiskMass[2,*] + G(www0).sfh_ElementsDiskMass[4,*] + G(www0).sfh_ElementsDiskMass[5,*] + G(www0).sfh_ElementsDiskMass[7,*] + G(www0).sfh_ElementsDiskMass[8,*] + G(www0).sfh_ElementsDiskMass[9,*])/G(www0).sfh_ElementsDiskMass[10,*]) $
;;         - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/Fe_mf) 
;
;FOR i=0, n_elements(w0)-1 DO BEGIN
;oplot, alog10(G(w0[i]).sfh_ElementsDiskMass[10,*]/G(w0[i]).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf), alog10((G(w0[i]).sfh_ElementsDiskMass[2,*] + G(w0[i]).sfh_ElementsDiskMass[4,*] + G(w0[i]).sfh_ElementsDiskMass[5,*] + G(w0[i]).sfh_ElementsDiskMass[7,*] + G(w0[i]).sfh_ElementsDiskMass[8,*] + G(w0[i]).sfh_ElementsDiskMass[9,*])/G(w0[i]).sfh_ElementsDiskMass[10,*]) $
;         - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/Fe_mf) $
;         , psym = 8, symsize = ssiz, color = fsc_color("black")
;ENDFOR
;FOR i=0, n_elements(ww0)-1 DO BEGIN
;oplot, alog10(G(ww0[i]).sfh_ElementsDiskMass[10,*]/G(ww0[i]).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf), alog10((G(ww0[i]).sfh_ElementsDiskMass[2,*] + G(ww0[i]).sfh_ElementsDiskMass[4,*] + G(ww0[i]).sfh_ElementsDiskMass[5,*] + G(ww0[i]).sfh_ElementsDiskMass[7,*] + G(ww0[i]).sfh_ElementsDiskMass[8,*] + G(ww0[i]).sfh_ElementsDiskMass[9,*])/G(ww0[i]).sfh_ElementsDiskMass[10,*]) $
;         - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/Fe_mf) $
;         , psym = 8, symsize = ssiz, color = fsc_color("blue")
;ENDFOR
;FOR i=0, n_elements(www0)-1 DO BEGIN
;oplot, alog10(G(www0[i]).sfh_ElementsDiskMass[10,*]/G(www0[i]).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf), alog10((G(www0[i]).sfh_ElementsDiskMass[2,*] + G(www0[i]).sfh_ElementsDiskMass[4,*] + G(www0[i]).sfh_ElementsDiskMass[5,*] + G(www0[i]).sfh_ElementsDiskMass[7,*] + G(www0[i]).sfh_ElementsDiskMass[8,*] + G(www0[i]).sfh_ElementsDiskMass[9,*])/G(www0[i]).sfh_ElementsDiskMass[10,*]) $
;         - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/Fe_mf)  $
;         , psym = 8, symsize = ssiz, color = fsc_color("red")
;ENDFOR
;legend, ['At z=0: In MW haloes', '...and are central discs', '...and have 1.0<SFR<5.0'], box = 0, charsize = 1.5, linestyle = [-99,-99,-99], textcolor=[fsc_color("black"), fsc_color("blue"), fsc_color("red")], /top, /right

;-------------------------------
;PLOT [alpha/Fe]_DiskMass vs. [Fe/H]_DiskMass (LOW RES): NB: Bovy et al. 2012 use different [a/Fe] (i.e. mean of [Mg/Fe], [Si/Fe] and [Ca/Fe]):
;-------------------------------
plot, findgen(10), /nodata, xtitle = '[Fe/H]!Ddisc!N', ytitle = '['+alphasymbol+'/Fe]!Ddisc!N', $
      xrange = [-2.0, 0.5], yrange = [-0.1, 0.5], xstyle = 1, ystyle = 1, CHARSIZE=1.75, CHARTHICK=4., XTHICK=5., YTHICK=5. ;charsize = csiz

;Plottables:
TheFeH = fltarr(n_elements(ww0)*20)
TheaFe = fltarr(n_elements(ww0)*20)
k=0
FOR i=0, n_elements(ww0)-1 DO BEGIN
  FOR j=0, 20-1 DO BEGIN
    TheFeH[k] = alog10(G(ww0[i]).sfh_ElementsDiskMass[10,j]/G(ww0[i]).sfh_ElementsDiskMass[0,j]) - alog10(Fe_mf/H_mf)
    ;TheaFe[k] = MEAN([alog10(G(ww0[i]).sfh_ElementsDiskMass[6,j]/G(ww0[i]).sfh_ElementsDiskMass[10,j]) - alog10(Mg_mf/Fe_mf), $
    ;                 alog10(G(ww0[i]).sfh_ElementsDiskMass[7,j]/G(ww0[i]).sfh_ElementsDiskMass[10,j]) - alog10(Si_mf/Fe_mf), $
    ;                 alog10(G(ww0[i]).sfh_ElementsDiskMass[9,j]/G(ww0[i]).sfh_ElementsDiskMass[10,j]) - alog10(Ca_mf/Fe_mf)])
    TheaFe[k] = alog10((G(ww0[i]).sfh_ElementsDiskMass[2,j] + G(ww0[i]).sfh_ElementsDiskMass[4,j] + G(ww0[i]).sfh_ElementsDiskMass[5,j] + G(ww0[i]).sfh_ElementsDiskMass[7,j] + G(ww0[i]).sfh_ElementsDiskMass[8,j] + G(ww0[i]).sfh_ElementsDiskMass[9,j])/G(ww0[i]).sfh_ElementsDiskMass[10,j]) $
              - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/Fe_mf)    
    k = k+1
  ENDFOR
ENDFOR

loadct, 13

;Parameters:
FeHmin = -3.0 
FeHmax = 1.0
FeHbinno = 75.
FeHbinwidth = (FeHmax-FeHmin) / FeHbinno
aFemin = -0.3
aFemax = 0.7
aFebinno = 75.
aFebinwidth = (OFemax-OFemin) / OFebinno

;Bin galaxies by [Fe/H] and [a/Fe]:
AveFeH = fltarr(FeHbinno*aFebinno)
AveaFe = fltarr(FeHbinno*aFebinno)
count = intarr(FeHbinno*aFebinno)
deleg = intarr(n_elements(ww0)*20)
k=0
FOR i=FeHmin, FeHmax-FeHbinwidth, FeHbinwidth DO BEGIN
    FOR j=aFemin, aFemax-aFebinwidth, aFebinwidth DO BEGIN
        w99 = WHERE(TheFeH ge i and TheFeH lt i+FeHbinwidth and TheaFe ge j and TheaFe lt j+aFebinwidth)
        IF(w99[0] ne -1) THEN BEGIN
            IF(n_elements(w99) ge 5) THEN deleg(w99) = 1
            ;deleg(w99) = 1
            AveFeH[k] = i + (0.5*FeHbinwidth)
            AveaFe[k] = j + (0.5*aFebinwidth)
            count[k] = n_elements(w99)
        ENDIF
        k+=1
    ENDFOR
ENDFOR

ww2 = WHERE(count ge 5)

;Make plotting symbol:
usersym, [-1.1,1.1,1.1,-1.1], [1.1,1.1,-1.1,-1.1], /fill
ctload, 0, /REVERSE, CLIP=[50,200] ;CLIP=[50,225]

plots, AveFeH(ww2), AveaFe(ww2), psym=8, color=bytscl(count(ww2))

loadct, 13, /silent
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill

ww99 = WHERE(deleg eq 0)
oplot, TheFeH[ww99], TheaFe[ww99], psym = 8, symsize = ssiz, color=fsc_color("grey")

;Plot solar lines:
oplot, [-4.0, 2.0], [0.0, 0.0], thick=6, linestyle = 2, color=fsc_color("black")
oplot, [0.0, 0.0], [-0.5, 0.8], thick=6, linestyle = 2, color=fsc_color("black")

;Title:
xyouts, -1.8, 0.45, "MW-type galaxies", charsize=2., charthick=4.

;Plot colour table
     xran = [-1.8,-0.5] ;[-1.7,-1.8]
     yran = [0.05,0.1] ;[-0.07,0.0]
     delta_yran = yran[1] - yran[0]
     delta_xran = xran[1] - xran[0]
     xarr = findgen(100)/100. * (delta_xran+0.01) + xran[0] ;* 0.50 + delta_xran * 0.25 + xran[0]
     yarr = findgen(2)/2. * delta_yran + yran[0] ;* 0.12 + delta_yran * 0.805 + yran[0]
     zarr = FLTARR([(SIZE(xarr))[1], (SIZE(yarr))[1]])
     zarr[*,0] = findgen(100)/100.
     zarr[*,1] = findgen(100)/100.  
     loadct, 0, /silent
     contour, zarr, xarr, yarr, /fill, /over, levels=findgen(175) / 175., c_colors=reverse(indgen(175)) ;levels=findgen(256) / 256., c_colors=reverse(indgen(256))
     plots, [xran(0), xran(0), xran(1), xran(1), xran(0)], [yran(0), (delta_yran/2.)+yran(0), (delta_yran/2.)+yran(0), yran(0), yran(0)], thick=2.
print, "MIN(count(ww2)) = ", MIN(count(ww2))
print, "MEAN(count(ww2)) = ", MEAN(count(ww2))
print, "MEDIAN(count(ww2)) = ", MEDIAN(count(ww2))
print, "MID(count(ww2)) = ", ((MAX(count(ww2))-MIN(count(ww2)))/2.0)+MIN(count(ww2))
print, "MAX(count(ww2)) = ", MAX(count(ww2))
print, "  "
     xyouts, (xran[1]-xran[0])/2. + xran[0]-0.12, 0.08, 'Number', charsize=1., charthick=4.  
     xyouts, xran[0]-0.03, 0.03, '5', charsize=1., charthick=3.
     xyouts, (xran[1]-xran[0])/2. + xran[0] - 0.06, 0.03, '578', charsize=1., charthick=3.
     xyouts, xran[1]-0.07, 0.03, '1151', charsize=1., charthick=3.
     plots, [xran[0],xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [(xran[1]-xran[0])/2. + xran[0],(xran[1]-xran[0])/2. + xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [xran[1],xran[1]], [yran[0], yran[0]-0.006], thick=2.

;-------------------------------
;PLOT [alpha/Fe]_DiskMass vs. [Fe/H]_DiskMass (LOW RES) (COLOURED BY MASS IN STARS FORMED):
;-------------------------------
plot, findgen(10), /nodata, xtitle = '[Fe/H]!Ddisc!N', ytitle = '['+alphasymbol+'/Fe]!Ddisc!N', $
      xrange = [-2.0, 0.5], yrange = [-0.1, 0.5], xstyle = 1, ystyle = 1, CHARSIZE=1.75, CHARTHICK=4., XTHICK=5., YTHICK=5. ;charsize = csiz

;Plottables:
TheFeH = fltarr(n_elements(ww0)*20)
TheaFe = fltarr(n_elements(ww0)*20)
StarMassFormed = fltarr(n_elements(ww0)*20)
AgeFormed = fltarr(n_elements(ww0)*20) ;i.e. Lookback time (in years) from z=0 to middle of the SFH bin
k=0
FOR i=0, n_elements(ww0)-1 DO BEGIN
  FOR j=0, 20-1 DO BEGIN
    TheFeH[k] = alog10(G(ww0[i]).sfh_ElementsDiskMass[10,j]/G(ww0[i]).sfh_ElementsDiskMass[0,j]) - alog10(Fe_mf/H_mf)
    TheaFe[k] = alog10((G(ww0[i]).sfh_ElementsDiskMass[2,j] + G(ww0[i]).sfh_ElementsDiskMass[4,j] + G(ww0[i]).sfh_ElementsDiskMass[5,j] + G(ww0[i]).sfh_ElementsDiskMass[7,j] + G(ww0[i]).sfh_ElementsDiskMass[8,j] + G(ww0[i]).sfh_ElementsDiskMass[9,j])/G(ww0[i]).sfh_ElementsDiskMass[10,j]) $
              - alog10((O_mf + Ne_mf + Mg_mf + Si_mf + S_mf + Ca_mf)/Fe_mf)
    StarMassFormed[k] = G(ww0[i]).sfh_DiskMass[j]
    AgeFormed[k] = G(ww0[i]).sfh_time[j]
    k = k+1
  ENDFOR
ENDFOR

loadct, 13

;Parameters:
FeHmin = -3.0 
FeHmax = 1.0
FeHbinno = 75.
FeHbinwidth = (FeHmax-FeHmin) / FeHbinno
aFemin = -0.3
aFemax = 0.7
aFebinno = 75.
aFebinwidth = (OFemax-OFemin) / OFebinno

;Bin galaxies by [Fe/H] and [a/Fe]:
AveFeH = fltarr(FeHbinno*aFebinno)
AveaFe = fltarr(FeHbinno*aFebinno)
AveStarMassFormed = fltarr(FeHbinno*aFebinno)
;TotalStarMassFormed = fltarr(FeHbinno*aFebinno)
AveAgeFormed = fltarr(FeHbinno*aFebinno)
count = intarr(FeHbinno*aFebinno)
deleg = intarr(n_elements(ww0)*20)
k=0
FOR i=FeHmin, FeHmax-FeHbinwidth, FeHbinwidth DO BEGIN
    FOR j=aFemin, aFemax-aFebinwidth, aFebinwidth DO BEGIN
        w99 = WHERE(TheFeH ge i and TheFeH lt i+FeHbinwidth and TheaFe ge j and TheaFe lt j+aFebinwidth)
        IF(w99[0] ne -1) THEN BEGIN
            IF(n_elements(w99) ge 5) THEN deleg(w99) = 1
            ;deleg(w99) = 1
            AveFeH[k] = i + (0.5*FeHbinwidth)
            AveaFe[k] = j + (0.5*aFebinwidth)
            AveStarMassFormed[k] = MEDIAN(StarMassFormed(w99))
            AveAgeFormed[k] = MEDIAN(AgeFormed(w99))
            ;TotalStarMassFormed[k] = TOTAL(StarMassFormed(w99))
            count[k] = n_elements(w99)
        ENDIF
        k+=1
    ENDFOR
ENDFOR

ww2 = WHERE(count ge 5)

;Make plotting symbol:
usersym, [-1.1,1.1,1.1,-1.1], [1.1,1.1,-1.1,-1.1], /fill
ctload, 0, /REVERSE, CLIP=[50,200] ;CLIP=[50,225]

plots, AveFeH(ww2), AveaFe(ww2), psym=8, color=bytscl(AveStarMassFormed(ww2))

;loadct, 13
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill

ww99 = WHERE(deleg eq 0)
oplot, TheFeH[ww99], TheaFe[ww99], psym = 8, symsize = ssiz, color=fsc_color("grey") ;color=bytscl(StarMassFormed[ww99]) ;NB: Scale won't be the same for the colour here as in the low-res plot line just above

loadct, 13, /silent

oplot, [-4.0, 2.0], [0.0, 0.0], thick=6, linestyle = 2, color=fsc_color("black")
oplot, [0.0, 0.0], [-0.5, 0.8], thick=6, linestyle = 2, color=fsc_color("black")

;Title:
xyouts, -1.8, 0.45, "MW-type galaxies", charsize=2., charthick=4.

;Plot colour table
     xran = [-1.8,-0.5] ;[-1.7,-1.8]
     yran = [0.05,0.1] ;[-0.07,0.0]
     delta_yran = yran[1] - yran[0]
     delta_xran = xran[1] - xran[0]
     xarr = findgen(100)/100. * (delta_xran+0.01) + xran[0] ;* 0.50 + delta_xran * 0.25 + xran[0]
     yarr = findgen(2)/2. * delta_yran + yran[0] ;* 0.12 + delta_yran * 0.805 + yran[0]
     zarr = FLTARR([(SIZE(xarr))[1], (SIZE(yarr))[1]])
     zarr[*,0] = findgen(100)/100.
     zarr[*,1] = findgen(100)/100.  
     loadct, 0, /silent
     contour, zarr, xarr, yarr, /fill, /over, levels=findgen(175) / 175., c_colors=reverse(indgen(175)) ;levels=findgen(256) / 256., c_colors=reverse(indgen(256))
     plots, [xran(0), xran(0), xran(1), xran(1), xran(0)], [yran(0), (delta_yran/2.)+yran(0), (delta_yran/2.)+yran(0), yran(0), yran(0)], thick=2.
print, "MIN(AveStarMassFormed(ww2)) = ", MIN(AveStarMassFormed(ww2))*1.0e10/Hubble_h
print, "MEAN(AveStarMassFormed(ww2)) = ", MEAN(AveStarMassFormed(ww2))*1.0e10/Hubble_h
print, "MEDIAN(AveStarMassFormed(ww2)) = ", MEDIAN(AveStarMassFormed(ww2))*1.0e10/Hubble_h
print, "MID(AveStarMassFormed(ww2)) = ", (((MAX(AveStarMassFormed(ww2))-MIN(AveStarMassFormed(ww2)))/2.0)+MIN(AveStarMassFormed(ww2)))*1.0e10/Hubble_h
print, "MAX(AveStarMassFormed(ww2)) = ", MAX(AveStarMassFormed(ww2))*1.0e10/Hubble_h
print, "  "
     xyouts, (xran[1]-xran[0])/2. + xran[0]-0.25, 0.08, 'Mass formed [M!D!9'+string(110B)+'!3!N]', charsize=1., charthick=4.  
     xyouts, xran[0]-0.12, 0.03, '1.3x10!E5!N', charsize=1., charthick=3.
     xyouts, (xran[1]-xran[0])/2. + xran[0]-0.12, 0.03, '8.7x10!E9!N', charsize=1., charthick=3.
     xyouts, xran[1]-0.12, 0.03, '1.7x10!E10!N', charsize=1., charthick=3.
     plots, [xran[0],xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [(xran[1]-xran[0])/2. + xran[0],(xran[1]-xran[0])/2. + xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [xran[1],xran[1]], [yran[0], yran[0]-0.006], thick=2.

;-------------------------------
;PLOT [alpha/Fe]_DiskMass vs. [Fe/H]_DiskMass (LOW RES) (COLOURED BY AGE):

plot, findgen(10), /nodata, xtitle = '[Fe/H]!Ddisk!N', ytitle = '['+alphasymbol+'/Fe]!Ddisk!N', $
      xrange = [-2.0, 0.5], yrange = [-0.1, 0.5], xstyle = 1, ystyle = 1, CHARSIZE=1.75, CHARTHICK=4., XTHICK=5., YTHICK=5. ;charsize = csiz
            
;Make plotting symbol:
usersym, [-1.1,1.1,1.1,-1.1], [1.1,1.1,-1.1,-1.1], /fill
ctload, 0, /REVERSE, CLIP=[50,200] ;CLIP=[50,225]

plots, AveFeH(ww2), AveaFe(ww2), psym=8, color=bytscl(AveAgeFormed(ww2))

;loadct, 13
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill

ww99 = WHERE(deleg eq 0)
oplot, TheFeH[ww99], TheaFe[ww99], psym = 8, symsize = ssiz, color=fsc_color("grey") ;color=bytscl(StarMassFormed[ww99]) ;NB: Scale won't be the same for the colour here as in the low-res plot line just above

loadct, 13, /silent

oplot, [-4.0, 2.0], [0.0, 0.0], thick=6, linestyle = 2, color=fsc_color("black")
oplot, [0.0, 0.0], [-0.5, 0.8], thick=6, linestyle = 2, color=fsc_color("black")

;Title:
xyouts, -1.8, 0.45, "MW-type galaxies", charsize=2., charthick=4.

;Plot colour table
     xran = [-1.8,-0.5] ;[-1.7,-1.8]
     yran = [0.05,0.1] ;[-0.07,0.0]
     delta_yran = yran[1] - yran[0]
     delta_xran = xran[1] - xran[0]
     xarr = findgen(100)/100. * (delta_xran+0.01) + xran[0] ;* 0.50 + delta_xran * 0.25 + xran[0]
     yarr = findgen(2)/2. * delta_yran + yran[0] ;* 0.12 + delta_yran * 0.805 + yran[0]
     zarr = FLTARR([(SIZE(xarr))[1], (SIZE(yarr))[1]])
     zarr[*,0] = findgen(100)/100.
     zarr[*,1] = findgen(100)/100.  
     loadct, 0, /silent
     contour, zarr, xarr, yarr, /fill, /over, levels=findgen(175) / 175., c_colors=reverse(indgen(175)) ;levels=findgen(256) / 256., c_colors=reverse(indgen(256))
     plots, [xran(0), xran(0), xran(1), xran(1), xran(0)], [yran(0), (delta_yran/2.)+yran(0), (delta_yran/2.)+yran(0), yran(0), yran(0)], thick=2.
print, "MIN(AveAgeFormed(ww2)) = ", MIN(AveAgeFormed(ww2))
print, "MEAN(AveAgeFormed(ww2)) = ", MEAN(AveAgeFormed(ww2))
print, "MEDIAN(AveAgeFormed(ww2)) = ", MEDIAN(AveAgeFormed(ww2))
print, "MID(AveAgeFormed(ww2)) = ", ((MAX(AveAgeFormed(ww2))-MIN(AveAgeFormed(ww2)))/2.0)+MIN(AveAgeFormed(ww2))
print, "MAX(AveAgeFormed(ww2)) = ", MAX(AveAgeFormed(ww2))
print, "  "
     xyouts, (xran[1]-xran[0])/2. + xran[0]-0.14, 0.08, 'Age [Gyrs]', charsize=1., charthick=4.  
     xyouts, xran[0]-0.07, 0.03, '0.04', charsize=1., charthick=3.
     xyouts, (xran[1]-xran[0])/2. + xran[0] - 0.07, 0.03, '6.33', charsize=1., charthick=3.
     xyouts, xran[1]-0.07, 0.03, '12.6', charsize=1., charthick=3.
     plots, [xran[0],xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [(xran[1]-xran[0])/2. + xran[0],(xran[1]-xran[0])/2. + xran[0]], [yran[0], yran[0]-0.006], thick=2.
     plots, [xran[1],xran[1]], [yran[0], yran[0]-0.006], thick=2.

;-------------------------------
;PLOT [Fe/H]_Disk distribution histogram:
;-------------------------------
plot, findgen(10), /nodata, xtitle = '[Fe/H]!Ddisc!N', ytitle = 'Fraction', $
      xrange = [-2.0, 0.5], yrange = [0.0,0.2], xstyle = 1, ystyle = 1, CHARSIZE=1.75, CHARTHICK=4., XTHICK=5., YTHICK=5. ;charsize = csiz

;Bin galaxies by [Fe/H]
FeHHistx = fltarr(FeHbinno)
HistCount = intarr(FeHbinno)
obsHistCount = intarr(FeHbinno)
k=0
FOR i=FeHmin, FeHmax-FeHbinwidth, FeHbinwidth DO BEGIN
    w99 = WHERE(TheFeH ge i and TheFeH lt i+FeHbinwidth)
    w99o = WHERE(Holmberg.Fe_H ge i and Holmberg.Fe_H lt i+FeHbinwidth)
    IF(w99[0] ne -1) THEN BEGIN
        FeHHistx[k] = i + (0.5*FeHbinwidth)
        HistCount[k] = n_elements(w99)
        obsHistCount[k] = n_elements(w99o)
    ENDIF
    k+=1
ENDFOR

;Normalise HistCount:
NormHistCount = fltarr(FeHbinno)
NormHistCount = HistCount/TOTAL(HistCount)
NormobsHistCount = fltarr(FeHbinno)
NormobsHistCount = obsHistCount/TOTAL(obsHistCount)
;print, FeHHistx[0:9]
;print, HistCount[0:9]
;print, NormHistCount[0:9]

;Make histogram arrays:
xM = fltarr(n_elements(FeHHistx)*2.)
j=0
FOR i=0, n_elements(FeHHistx)-1 DO BEGIN
    xM[j] = FeHHistx[i]-0.5*FeHbinwidth
    xM[j+1] = FeHHistx[i]+0.5*FeHbinwidth
    j=j+2
ENDFOR
yM = fltarr(n_elements(NormHistCount)*2.)
yyM = fltarr(n_elements(NormobsHistCount)*2.)
j=0
FOR i=0, n_elements(NormHistCount)-1 DO BEGIN
    yM[j] = NormHistCount[i]
    yM[j+1] = NormHistCount[i]
    yyM[j] = NormobsHistCount[i]
    yyM[j+1] = NormobsHistCount[i]
    j=j+2
ENDFOR    

oplot, xM, yyM, thick=6, linestyle = 0, color=fsc_color("red")
oplot, xM, yM, thick=6, linestyle = 0, color=fsc_color("black")

legend, ['Model MW-type galaxies', 'Holmberg et al. (2009)'], box = 0, charsize=1.75, charthick=4., linestyle = [-99,-99], textcolor=[fsc_color("black"),fsc_color("red")], /top, /left

;-------------------------------
;PLOT [Fe/H]_Bulge distribution histogram:
;-------------------------------
plot, findgen(10), /nodata, xtitle = '[Fe/H]!Dbulge!N', ytitle = 'Fraction', $
      xrange = [-2.0, 0.5], yrange = [0.0,0.2], xstyle = 1, ystyle = 1, CHARSIZE=1.75, CHARTHICK=4., XTHICK=5., YTHICK=5. ;charsize = csiz

;Plottables:
TheFeH_bulge = fltarr(n_elements(ww0)*20)
k=0
FOR i=0, n_elements(ww0)-1 DO BEGIN
  FOR j=0, 20-1 DO BEGIN
    TheFeH_bulge[k] = alog10(G(ww0[i]).sfh_ElementsBulgeMass[10,j]/G(ww0[i]).sfh_ElementsBulgeMass[0,j]) - alog10(Fe_mf/H_mf)
    k = k+1
  ENDFOR
ENDFOR

;Bin galaxies by [Fe/H]
FeHHistx = fltarr(FeHbinno)
HistCount = intarr(FeHbinno)
k=0
FOR i=FeHmin, FeHmax-FeHbinwidth, FeHbinwidth DO BEGIN
    w99 = WHERE(TheFeH_bulge ge i and TheFeH_bulge lt i+FeHbinwidth)
    IF(w99[0] ne -1) THEN BEGIN
        FeHHistx[k] = i + (0.5*FeHbinwidth)
        HistCount[k] = n_elements(w99)
    ENDIF
    k+=1
ENDFOR

;Normalise HistCount:
NormHistCount = fltarr(FeHbinno)
NormHistCount = HistCount/TOTAL(HistCount)

print, FeHHistx[0:9]
print, HistCount[0:9]
print, NormHistCount[0:9]

;Make histogram arrays:
xM = fltarr(n_elements(FeHHistx)*2.)
j=0
FOR i=0, n_elements(FeHHistx)-1 DO BEGIN
    xM[j] = FeHHistx[i]-0.5*FeHbinwidth
    xM[j+1] = FeHHistx[i]+0.5*FeHbinwidth
    j=j+2
ENDFOR
yM = fltarr(n_elements(NormHistCount)*2.)
j=0
FOR i=0, n_elements(NormHistCount)-1 DO BEGIN
    yM[j] = NormHistCount[i]
    yM[j+1] = NormHistCount[i]
    j=j+2
ENDFOR    

oplot, xM, yM, thick=6, linestyle = 0, color=fsc_color("black")

;-------------------------------
;PLOT [O/Fe]_Disk+BulgeMass vs. [Fe/H]_Disk+BulgeMass:
;-------------------------------
;www0a = WHERE(FeH_mf_disk(w0) ge 0.0)
;www0b = WHERE(FeH_mf_bulge(w0) ge 0.0)
plot, findgen(10), /nodata, xtitle = '[Fe/H]', ytitle = '[O/Fe]', $
            charsize = csiz, xrange = [-1.5, 0.5], yrange = [-0.2, 0.4], xstyle = 1, ystyle = 1 ; xrange = [MIN(FeH_mf_bulge(w0)), MAX(FeH_mf_disk(w0))], yrange = [MIN(OFe_mf_disk(w0)), MAX(OFe_mf_bulge(w0))]
oplot, FeH_mf_disk(w0), OFe_mf_disk(w0), psym = 8, symsize = ssiz, color = fsc_color("black")
oplot, FeH_mf_disk(ww0), OFe_mf_disk(ww0), psym = 8, symsize = ssiz, color = fsc_color("blue")
oplot, FeH_mf_disk(www0), OFe_mf_disk(www0), psym = 8, symsize = ssiz, color = fsc_color("red")
oplot, FeH_mf_bulge(w0), OFe_mf_bulge(w0), psym = 8, symsize = ssiz, color = fsc_color("black")
oplot, FeH_mf_bulge(ww0), OFe_mf_bulge(ww0), psym = 8, symsize = ssiz, color = fsc_color("blue")
oplot, FeH_mf_bulge(www0), OFe_mf_bulge(www0), psym = 8, symsize = ssiz, color = fsc_color("red")
legend, ['In MW haloes', '...and are central discs', '...and have 1.0<SFR<5.0'], box = 0, charsize = 1.5, linestyle = [-99,-99,-99], textcolor=[fsc_color("black"), fsc_color("blue"), fsc_color("red")], /top, /left

;-------------------------------
;PLOT Age vs. [Fe/H]_DiskMass:
;-------------------------------
age1 = (G(ww0).sfh_time[*]+(0.5*G(ww0).sfh_dt[*]))/1.0e9

print, alog10(G(ww0[0]).sfh_time[*]/1.0e9)

plot, findgen(10), /nodata, xtitle = 'Age [Gyr]', ytitle = '[Fe/H]!Ddisk!N', $
            charsize = csiz, xrange = [-2.2, 1.3], yrange = [-3.0, 1.0], xstyle = 1, ystyle = 1
FOR i=0, n_elements(ww0)-1 DO BEGIN
;plots, alog10(G(w0[i]).sfh_ElementsDiskMass[10,*]/G(w0[i]).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf), alog10(G(w0[i]).sfh_ElementsDiskMass[4,*]/G(w0[i]).sfh_ElementsDiskMass[10,*]) - alog10(O_mf/Fe_mf), psym = 8, symsize = ssiz, color = fsc_color("black")
c = G(ww0[i])
d = c.sfh_ElementsDiskMass
;plots, alog10(age1), alog10(G(ww0[i]).sfh_ElementsDiskMass[10,*]/G(ww0[i]).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf), psym = 8, symsize = ssiz, color = fsc_color("blue")
plots, alog10(age1), alog10(d[10,*]/d[0,*]) - alog10(Fe_mf/H_mf), psym = 8, symsize = ssiz, color = fsc_color("blue")
;plots, alog10(G(www0[i]).sfh_ElementsDiskMass[10,*]/G(www0[i]).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf), alog10(G(www0[i]).sfh_ElementsDiskMass[4,*]/G(www0[i]).sfh_ElementsDiskMass[10,*]) - alog10(O_mf/Fe_mf), psym = 8, symsize = ssiz, color = fsc_color("red")
ENDFOR
;legend, ['In MW haloes', '...and are central discs', '...and have 1.0<SFR<5.0'], box = 0, charsize = 1.5, linestyle = [-99,-99,-99], textcolor=[fsc_color("black"), fsc_color("blue"), fsc_color("red")], /top, /left
;print, alog10(G(www0[0]).sfh_ElementsDiskMass[10,*]/G(www0[0]).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf)


;-------------------------------
;PLOT [Mg/Fe]_Disk+BulgeMass vs. [Fe/H]_Disk+BulgeMass:
;-------------------------------
;www0a = WHERE(FeH_mf_disk(w0) ge 0.0)
;www0b = WHERE(FeH_mf_bulge(w0) ge 0.0)
plot, findgen(10), /nodata, xtitle = '[Fe/H]', ytitle = '[Mg/Fe]', $
            charsize = csiz, xrange = [-1.5, 0.5], yrange = [-0.2, 0.4], xstyle = 1, ystyle = 1 ; xrange = [MIN(FeH_mf_bulge(w0)), MAX(FeH_mf_disk(w0))], yrange = [MIN(OFe_mf_disk(w0)), MAX(OFe_mf_bulge(w0))]
oplot, FeH_mf_disk(w0), MgFe_mf_disk(w0), psym = 8, symsize = ssiz, color = fsc_color("black")
oplot, FeH_mf_disk(ww0), MgFe_mf_disk(ww0), psym = 8, symsize = ssiz, color = fsc_color("blue")
oplot, FeH_mf_disk(www0), MgFe_mf_disk(www0), psym = 8, symsize = ssiz, color = fsc_color("red")
oplot, FeH_mf_bulge(w0), MgFe_mf_bulge(w0), psym = 8, symsize = ssiz, color = fsc_color("black")
oplot, FeH_mf_bulge(ww0), MgFe_mf_bulge(ww0), psym = 8, symsize = ssiz, color = fsc_color("blue")
oplot, FeH_mf_bulge(www0), MgFe_mf_bulge(www0), psym = 8, symsize = ssiz, color = fsc_color("red")
legend, ['In MW haloes', '...and are central discs', '...and have 1.0<SFR<5.0'], box = 0, charsize = 1.5, linestyle = [-99,-99,-99], textcolor=[fsc_color("black"), fsc_color("blue"), fsc_color("red")], /top, /left

;------------------------------------------------------------------
;------------------------------------------------------------------
;BULGE (MILKY WAY TYPE HALOES):
;------------------------------------------------------------------
;------------------------------------------------------------------

;-------------------------------
;PLOT [O/Fe]_BulgeMass vs. M*:
;-------------------------------
;plot, findgen(10), /nodata, xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = '[O/Fe]!Dbulge!N', $
;            charsize = csiz, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [MIN(OFe_mf_bulge(w0)), MAX(OFe_mf_bulge(w0))], xstyle = 1, ystyle = 1 
;
;oplot, alog10(stars(w0)*1e10), OFe_mf_bulge(w0), psym = 8, symsize = ssiz, color = fsc_color("black")
;legend, ['MW haloes'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left

plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[O/Fe]!Dbulge!N', $
            charsize = csiz, xrange = [MIN(alog10(stars(w1)*1e10)), 11.5], yrange = [MIN(OFe_mf_disk(w0)), 0.5], xstyle = 1, ystyle = 1 ;[MIN(OFe_mf_disk(w0)), MAX(OFe_mf_disk(w0))] ;[MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))] 

oplot, alog10(stars(w0)*1e10), OFe_mf_bulge(w0), psym = 8, symsize = ssiz, color = fsc_color("black")
oplot, alog10(stars(ww0)*1e10), OFe_mf_bulge(ww0), psym = 8, symsize = ssiz, color = fsc_color("blue")
oplot, alog10(stars(www0)*1e10), OFe_mf_bulge(www0), psym = 8, symsize = ssiz, color = fsc_color("red")
legend, ['In MW haloes', '...and are central discs', '...and have 1.0<SFR<5.0'], box = 0, charsize = 1.5, linestyle = [-99,-99,-99], textcolor=[fsc_color("black"), fsc_color("blue"), fsc_color("red")], /top, /left


;-------------------------------
;PLOT [Mg/Fe]_BulgeMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[Mg/Fe]!Dbulge!N', $
            charsize = csiz, xrange = [MIN(alog10(stars(w1)*1e10)), 11.5], yrange = [MIN(MgFe_mf_disk(w0)), 0.5], xstyle = 1, ystyle = 1 
            
oplot, alog10(stars(w0)*1e10), MgFe_mf_bulge(w0), psym = 8, symsize = ssiz, color = fsc_color("black")
oplot, alog10(stars(ww0)*1e10), MgFe_mf_bulge(ww0), psym = 8, symsize = ssiz, color = fsc_color("blue")
oplot, alog10(stars(www0)*1e10), MgFe_mf_bulge(www0), psym = 8, symsize = ssiz, color = fsc_color("red")
;legend, ['MW haloes'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left
legend, ['In MW haloes', '...and are central discs', '...and have 1.0<SFR<5.0'], box = 0, charsize = 1.5, linestyle = [-99,-99,-99], textcolor=[fsc_color("black"), fsc_color("blue"), fsc_color("red")], /top, /left


;;-------------------------------
;;PLOT [alpha/Fe]_BulgeMass vs. M*:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = '[alpha/Fe]!Dbulge!N', $
;            charsize = 1.0, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [0.0,1.0], xstyle = 1, ystyle = 1 ;[MIN(aFe_mf_disk(w0)), MAX(aFe_mf_disk(w0))]
;            
;;oplot, alog10(stars(w0)*1e10), aFe_mf_bulge(w0), psym = 8, symsize = ssiz, color = fsc_color("black")
;legend, ['No alpha with MAINELEMENTS'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left
;
;;-------------------------------
;;PLOT [O/Fe]_BulgeMass vs. [Fe/H]_BulgeMass:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = '[Fe/H]!DBulge!N', ytitle = '[O/Fe]!Dbulge!N', $
;            charsize = 1.0, xrange = [MIN(FeH_mf_bulge(w0)), MAX(FeH_mf_bulge(w0))], yrange = [MIN(OFe_mf_bulge(w0)), MAX(OFe_mf_bulge(w0))], xstyle = 1, ystyle = 1 
;            
;oplot, FeH_mf_bulge(w0), OFe_mf_bulge(w0), psym = 8, symsize = ssiz, color = fsc_color("black")
;legend, ['MW haloes'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left

;------------------------------------------------------------------
;------------------------------------------------------------------
;DISK + BULGE (MILKY WAY TYPE HALOES):
;------------------------------------------------------------------
;------------------------------------------------------------------
;
;;-------------------------------
;;PLOT [O/Fe]_Disk+BulgeMass vs. M*:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = '[O/Fe]!Ddisk+bulge!N', $
;            charsize = 1.0, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [MIN(OFe_mf_comb(w0)), MAX(OFe_mf_comb(w0))], xstyle = 1, ystyle = 1 
;
;oplot, alog10(stars(w0)*1e10), OFe_mf_comb(w0), psym = 8, symsize = ssiz, color = fsc_color("black")
;legend, ['MW haloes'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left
;
;;-------------------------------
;;PLOT [Mg/Fe]_Disk+BulgeMass vs. M*:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = '[Mg/Fe]!Ddisk+bulge!N', $
;            charsize = 1.0, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [MIN(MgFe_mf_comb(w0)), MAX(MgFe_mf_comb(w0))], xstyle = 1, ystyle = 1 
;            
;oplot, alog10(stars(w0)*1e10), MgFe_mf_comb(w0), psym = 8, symsize = ssiz, color = fsc_color("black")
;legend, ['MW haloes'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left
;
;;-------------------------------
;;PLOT [alpha/Fe]_Disk+BulgeMass vs. M*:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = '[alpha/Fe]!Ddisk+bulge!N', $
;            charsize = 1.0, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [0.0,1.0], xstyle = 1, ystyle = 1 ;[MIN(aFe_mf_disk(w0)), MAX(aFe_mf_disk(w0))]
;            
;;oplot, alog10(stars(w0)*1e10), aFe_mf_bulge(w0), psym = 8, symsize = ssiz, color = fsc_color("black")
;legend, ['No alpha with MAINELEMENTS'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left
;
;;-------------------------------
;;PLOT [O/Fe]_Disk+BulgeMass vs. [Fe/H]_Disk+BulgeMass:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = '[Fe/H]!DBulge!N', ytitle = '[O/Fe]!Ddisk+bulge!N', $
;            charsize = 1.0, xrange = [MIN(FeH_mf_comb(w0(ww0))), MAX(FeH_mf_comb(w0(ww0)))], yrange = [MIN(OFe_mf_comb(w0(ww0))), MAX(OFe_mf_comb(w0(ww0)))], xstyle = 1, ystyle = 1 
;            
;oplot, FeH_mf_comb(w0(ww0)), OFe_mf_comb(w0(ww0)), psym = 8, symsize = ssiz, color = fsc_color("black")
;legend, ['MW haloes'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left


;------------------------------------------------------------------
;------------------------------------------------------------------
;DISK (ELLIPTICALS):
;------------------------------------------------------------------
;------------------------------------------------------------------

wexp1 = WHERE((alog10(stars(w1)*1e10) gt 11.0) AND (MgFe_mf_disk(w1) gt 0.1) AND (G(w1).sfh_DiskMass[15] eq 0.0) AND (G(w1).sfh_DiskMass[14] eq 0.0) AND (G(w1).sfh_DiskMass[13] eq 0.0) AND (G(w1).sfh_DiskMass[12] eq 0.0) AND (G(w1).sfh_DiskMass[11] eq 0.0) AND (G(w1).sfh_DiskMass[10] eq 0.0) AND (G(w1).sfh_DiskMass[9] eq 0.0) AND (G(w1).sfh_DiskMass[8] eq 0.0) AND (G(w1).sfh_DiskMass[7] eq 0.0) AND (G(w1).sfh_DiskMass[6] eq 0.0) AND (G(w1).sfh_DiskMass[5] eq 0.0) AND (G(w1).sfh_DiskMass[4] eq 0.0) AND (G(w1).sfh_DiskMass[3] eq 0.0) AND (G(w1).sfh_DiskMass[2] eq 0.0) AND (G(w1).sfh_DiskMass[1] ne 0.0))
print, wexp1
;print, G(w1(144)).sfh_DiskMass
;print, G(w1(319)).sfh_DiskMass
;print, G(w1(61)).sfh_DiskMass
;print, G(w1(238)).sfh_DiskMass
;print, G(w1(481)).sfh_DiskMass

;-------------------------------
;PLOT Age vs. M*:
;-------------------------------
;w0 = WHERE(G.Mg_yield_DiskMass ge 0.0 && G.Fe_yield_DiskMass gt 0.0)
;print, n_elements(w0)

Age = G.MassWeightAge
;print, MIN(Age), MAX(Age)

plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = 'log(Age) [Gyrs]', $
            charsize = csiz, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [MIN(alog10(Age(w1)))-0.02,MAX(alog10(Age(w1)))+0.02], xstyle = 1, ystyle = 1 ;[8.5, 12.0] ;[0.01692, 0.01696]
            
oplot, alog10(stars(w1)*1e10), alog10(Age(w1)), psym = 8, symsize = ssiz   
;oplot, [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], [alog10(13.6), alog10(13.6)], linestyle = 2
legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /bottom, /right

;-------------------------------
;PLOT [O/Fe]_DiskMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[O/Fe]!Ddisk!N', $
            charsize = csiz, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [MIN(OFe_mf_disk(w1)), MAX(OFe_mf_disk(w1))], xstyle = 1, ystyle = 1 

oplot, alog10(stars(w1)*1e10), OFe_mf_disk(w1), psym = 8, symsize = ssiz, color = fsc_color("black")

usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
plots, alog10(stars(w1(0))*1e10), OFe_mf_disk(w1(0)), psym = 8, symsize = ssiz*3., color = fsc_color("blue")
;plots, alog10(stars(w1(61))*1e10), OFe_mf_disk(w1(61)), psym = 8, symsize = ssiz*3., color = fsc_color("red")
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill

legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right
print, OFe_mf_disk(w1[0:9])
print, "MIN(OFe_mf_disk(w1)) = ", MIN(OFe_mf_disk(w1)), "MAX(OFe_mf_disk(w1)) = ", MAX(OFe_mf_disk(w1))
;-------------------------------
;PLOT [Mg/Fe]_DiskMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[Mg/Fe]!Ddisk!N', $
            charsize = csiz, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [MIN(MgFe_mf_disk(w1)), MAX(MgFe_mf_disk(w1))], xstyle = 1, ystyle = 1 
            
oplot, alog10(stars(w1)*1e10), MgFe_mf_disk(w1), psym = 8, symsize = ssiz, color = fsc_color("black")

usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
plots, alog10(stars(w1(0))*1e10), MgFe_mf_disk(w1(0)), psym = 8, symsize = ssiz*3., color = fsc_color("blue")
;plots, alog10(stars(w1(61))*1e10), MgFe_mf_disk(w1(61)), psym = 8, symsize = ssiz*3., color = fsc_color("red")
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill

legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

;;-------------------------------
;;PLOT [alpha/Fe]_DiskMass vs. M*:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = '[alpha/Fe]!Ddisk!N', $
;            charsize = 1.0, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [0.0,1.0], xstyle = 1, ystyle = 1 ;[MIN(aFe_mf_disk(w0)), MAX(aFe_mf_disk(w0))]
;            
;;oplot, alog10(stars(w1)*1e10), aFe_mf_disk(w1), psym = 8, symsize = ssiz, color = fsc_color("black")
;
;;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
;;plots, alog10(stars(w1(0))*1e10), aFe_mf_disk(w1(0)), psym = 8, symsize = ssiz*3., color = fsc_color("blue")
;;plots, alog10(stars(2434)*1e10), aFe_mf_disk(2434), psym = 8, symsize = ssiz*3., color = fsc_color("red")
;;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill
;
;legend, ['No alpha with MAINELEMENTS'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left
;
;;-------------------------------
;;PLOT [O/Fe]_DiskMass vs. [Fe/H]_DiskMass:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = '[Fe/H]!Ddisk!N', ytitle = '[O/Fe]!Ddisk!N', $
;            charsize = 1.0, xrange = [MIN(FeH_mf_disk(w1)), MAX(FeH_mf_disk(w1))], yrange = [MIN(OFe_mf_disk(w1)), MAX(OFe_mf_disk(w1))], xstyle = 1, ystyle = 1
;            
;oplot, FeH_mf_disk(w1), OFe_mf_disk(w1), psym = 8, symsize = ssiz, color = fsc_color("black")
;
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
;plots, FeH_mf_disk(w1(0)), OFe_mf_disk(w1(0)), psym = 8, symsize = ssiz*3., color = fsc_color("blue")
;plots, FeH_mf_disk(2434), OFe_mf_disk(2434), psym = 8, symsize = ssiz*3., color = fsc_color("red")
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill
;
;legend, ['Ellipticals'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left

;------------------------------------------------------------------
;------------------------------------------------------------------
;BULGE (ELLIPTICALS):
;------------------------------------------------------------------
;------------------------------------------------------------------
;print, "n_elements(w1) = ", n_elements(w1)
;-------------------------------
;PLOT [O/Fe]_BulgeMass vs. M*:
;-------------------------------
;plot, findgen(10), /nodata, xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = '[O/Fe]!Dbulge!N', $
;            charsize = 1.0, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [MIN(OFe_mf_bulge(w1)), MAX(OFe_mf_bulge(w1))], xstyle = 1, ystyle = 1 
;
;oplot, alog10(stars(w1)*1e10), OFe_mf_bulge(w1), psym = 8, symsize = ssiz, color = fsc_color("black")
;
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
;plots, alog10(stars(w1(0))*1e10), OFe_mf_bulge(w1(0)), psym = 8, symsize = ssiz*3., color = fsc_color("blue")
;plots, alog10(stars(2434)*1e10), OFe_mf_bulge(2434), psym = 8, symsize = ssiz*3., color = fsc_color("red")
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill
;
;legend, ['Ellipticals'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left

plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[O/Fe]!Dbulge!N', $
            charsize = csiz, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [MIN(OFe_mf_disk(w1)), MAX(OFe_mf_disk(w1))], xstyle = 1, ystyle = 1 

oplot, alog10(stars(w1)*1e10), OFe_mf_bulge(w1), psym = 8, symsize = ssiz, color = fsc_color("black")

usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
plots, alog10(stars(w1(0))*1e10), OFe_mf_bulge(w1(0)), psym = 8, symsize = ssiz*3., color = fsc_color("blue")
;plots, alog10(stars(w1(61))*1e10), OFe_mf_bulge(w1(61)), psym = 8, symsize = ssiz*3., color = fsc_color("red")
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill

legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

;-------------------------------
;PLOT [Mg/Fe]_BulgeMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[Mg/Fe]!Dbulge!N', $
            charsize = csiz, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [MIN(MgFe_mf_disk(w1)), MAX(MgFe_mf_disk(w1))], xstyle = 1, ystyle = 1 
            
oplot, alog10(stars(w1)*1e10), MgFe_mf_bulge(w1), psym = 8, symsize = ssiz, color = fsc_color("black")

usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
plots, alog10(stars(w1(0))*1e10), MgFe_mf_bulge(w1(0)), psym = 8, symsize = ssiz*3., color = fsc_color("blue")
;plots, alog10(stars(w1(61))*1e10), MgFe_mf_bulge(w1(61)), psym = 8, symsize = ssiz*3., color = fsc_color("red")
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill

legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

;;-------------------------------
;;PLOT [alpha/Fe]_BulgeMass vs. M*:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = '[alpha/Fe]!Dbulge!N', $
;            charsize = 1.0, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [0.0,1.0], xstyle = 1, ystyle = 1 ;[MIN(aFe_mf_disk(w0)), MAX(aFe_mf_disk(w0))]
;            
;;oplot, alog10(stars(w1)*1e10), aFe_mf_bulge(w1), psym = 8, symsize = ssiz, color = fsc_color("black")
;
;;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
;;plots, alog10(stars(w1(0))*1e10), aFe_mf_bulge(w1(0)), psym = 8, symsize = ssiz*3., color = fsc_color("blue")
;;plots, alog10(stars(2434)*1e10), aFe_mf_bulge(2434), psym = 8, symsize = ssiz*3., color = fsc_color("red")
;;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill
;
;legend, ['No alpha with MAINELEMENTS'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left
;
;;-------------------------------
;;PLOT [O/Fe]_BulgeMass vs. [Fe/H]_BulgeMass:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = '[Fe/H]!DBulge!N', ytitle = '[O/Fe]!Dbulge!N', $
;            charsize = 1.0, xrange = [MIN(FeH_mf_bulge(w1)), MAX(FeH_mf_bulge(w1))], yrange = [MIN(OFe_mf_bulge(w1)), MAX(OFe_mf_bulge(w1))], xstyle = 1, ystyle = 1 
;            
;oplot, FeH_mf_bulge(w1), OFe_mf_bulge(w1), psym = 8, symsize = ssiz, color = fsc_color("black")
;
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
;plots, FeH_mf_bulge(w1(0)), OFe_mf_bulge(w1(0)), psym = 8, symsize = ssiz*3., color = fsc_color("blue")
;plots, FeH_mf_bulge(2434), OFe_mf_bulge(2434), psym = 8, symsize = ssiz*3., color = fsc_color("red")
;plots, FeH_mf_bulge(w1(200)), OFe_mf_bulge(w1(200)), psym = 8, symsize = ssiz*3., color = fsc_color("green")
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill
;
;legend, ['Ellipticals'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left

;------------------------------------------------------------------
;------------------------------------------------------------------
;DISC AND BULGE (ELLIPTICALS):
;------------------------------------------------------------------
;------------------------------------------------------------------

;-------------------------------
;PLOT CLOSE-UP [H/Fe]_Disc+BulgeMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[H/Fe]!Ddisk+bulge!N', $
            charsize = csiz, xrange = [9.5, 11.5], yrange = [-0.2, 0.4], xstyle = 1, ystyle = 1
            
oplot, alog10(stars(w1)*1e10), HFe_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")

;Plot fit to all galaxies:
comp = WHERE(alog10(stars(w1)*1e10) ge 10.0 AND alog10(stars(w1)*1e10) le 11.3)
xxx = findgen(100)/(100.0/(11.3-10.0)) + 10.0
FitHFe = POLY_FIT(alog10(stars(w1(comp))*1e10), HFe_mf_comb(w1(comp)), 1, SIGMA=FitHFeSigma)
yyy = FitHFe[0] + FitHFe[1]*(xxx); + FitNFe[2]*(xxx)*(xxx) + FitNFe[3]*(xxx)*(xxx)*(xxx)
oplot, xxx, yyy, thick=6, linestyle = 0, color=fsc_color("black")

legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

;-------------------------------
;PLOT CLOSE-UP [He/Fe]_Disc+BulgeMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[He/Fe]!Ddisk+bulge!N', $
            charsize = csiz, xrange = [9.5, 11.5], yrange = [-0.2, 0.4], xstyle = 1, ystyle = 1
            
oplot, alog10(stars(w1)*1e10), HeFe_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")

;Plot fit to all galaxies:
comp = WHERE(alog10(stars(w1)*1e10) ge 10.0 AND alog10(stars(w1)*1e10) le 11.3)
xxx = findgen(100)/(100.0/(11.3-10.0)) + 10.0
FitHeFe = POLY_FIT(alog10(stars(w1(comp))*1e10), HeFe_mf_comb(w1(comp)), 1, SIGMA=FitHeFeSigma)
yyy = FitHeFe[0] + FitHeFe[1]*(xxx); + FitNFe[2]*(xxx)*(xxx) + FitNFe[3]*(xxx)*(xxx)*(xxx)
oplot, xxx, yyy, thick=6, linestyle = 0, color=fsc_color("black")

legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

IF (MAINELEMENTS eq 1) THEN BEGIN
;-------------------------------
;PLOT CLOSE-UP [C/Fe]_Disc+BulgeMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[C/Fe]!Ddisk+bulge!N', $
            charsize = csiz, xrange = [9.5, 11.5], yrange = [-0.2, 0.4], xstyle = 1, ystyle = 1
            
oplot, alog10(stars(w1)*1e10), CFe_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")

;Plot Jonas' fit:
logSigma = findgen(100)/(100.0/(2.45-1.8)) + 1.8
logMstar = (2.2*logSigma) + 5.85 ;MY log(sigma) TO log(M*) CONVERTER! USING FIG. 2 FROM WAKE ET AL. (2012)
yJ = 0.35*(logSigma) - 0.55 ;Slope and intercept from Johansson et al. (2009), Fig. 12
oplot, logMstar, yJ, linestyle = 0, thick = 6., color = fsc_color("orange")

;Plot fit to all galaxies:
comp = WHERE(alog10(stars(w1)*1e10) ge 10.0 AND alog10(stars(w1)*1e10) le 11.3)
xxx = findgen(100)/(100.0/(11.3-10.0)) + 10.0
FitCFe = POLY_FIT(alog10(stars(w1(comp))*1e10), CFe_mf_comb(w1(comp)), 1, SIGMA=FitCFeSigma)
yyy = FitCFe[0] + FitCFe[1]*(xxx); + FitCFe[2]*(xxx)*(xxx) + FitCFe[3]*(xxx)*(xxx)*(xxx)
oplot, xxx, yyy, thick=6, linestyle = 0, color=fsc_color("black")

legend, 'Johansson et al. (2009)', box = 0, charsize = 1.25, linestyle = 0, color = fsc_color("orange"), thick = 6, /bottom, /left
legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

;-------------------------------
;PLOT CLOSE-UP [N/Fe]_Disc+BulgeMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[N/Fe]!Ddisk+bulge!N', $
            charsize = csiz, xrange = [9.5, 11.5], yrange = [-0.2, 0.4], xstyle = 1, ystyle = 1
            
oplot, alog10(stars(w1)*1e10), NFe_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")
;print, "[N/Fe]  ", NFe_mf_comb(w1[0:9])

;Plot Jonas' fit:
logSigma = findgen(100)/(100.0/(2.45-1.8)) + 1.8
logMstar = (2.2*logSigma) + 5.85 ;MY log(sigma) TO log(M*) CONVERTER! USING FIG. 2 FROM WAKE ET AL. (2012)
yJ = 0.48*(logSigma) - 0.99 ;Slope and intercept from Johansson et al. (2009), Fig. 12
oplot, logMstar, yJ, linestyle = 0, thick = 6., color = fsc_color("orange")

;Plot fit to all galaxies:
comp = WHERE(alog10(stars(w1)*1e10) ge 10.0 AND alog10(stars(w1)*1e10) le 11.3)
xxx = findgen(100)/(100.0/(11.3-10.0)) + 10.0
FitNFe = POLY_FIT(alog10(stars(w1(comp))*1e10), NFe_mf_comb(w1(comp)), 1, SIGMA=FitNFeSigma)
yyy = FitNFe[0] + FitNFe[1]*(xxx); + FitNFe[2]*(xxx)*(xxx) + FitNFe[3]*(xxx)*(xxx)*(xxx)
oplot, xxx, yyy, thick=6, linestyle = 0, color=fsc_color("black")

legend, 'Johansson et al. (2009)', box = 0, charsize = 1.25, linestyle = 0, color = fsc_color("orange"), thick = 6, /bottom, /left
legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right


ENDIF ;IF (MAINELEMENTS eq 1) THEN BEGIN

;;-------------------------------
;;PLOT [O/Fe]_Disc+BulgeMass vs. M*:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[O/Fe]!Ddisk+bulge!N', $
;            charsize = csiz, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [-0.44, 0.6], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk(w1)), MAX(OFe_mf_disk(w1))] 
;
;oplot, alog10(stars(w1)*1e10), OFe_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")
;
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
;plots, alog10(stars(w1(0))*1e10), OFe_mf_comb(w1(0)), psym = 8, symsize = ssiz*3., color = fsc_color("blue")
;;plots, alog10(stars(w1(61))*1e10), OFe_mf_comb(w1(61)), psym = 8, symsize = ssiz*3., color = fsc_color("red")
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill
;
;legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

;-------------------------------
;PLOT CLOSE-UP [O/Fe]_Disc+BulgeMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[O/Fe]!Ddisk+bulge!N', $
            charsize = csiz, xrange = [9.5, 11.5], yrange = [0.0, 0.6], xstyle = 1, ystyle = 1
            
oplot, alog10(stars(w1)*1e10), OFe_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")

;Plot Jonas' fit:

logSigma = findgen(100)/(100.0/(2.45-1.8)) + 1.8
logMstar = (2.2*logSigma) + 5.85 ;MY log(sigma) TO log(M*) CONVERTER! USING FIG. 2 FROM WAKE ET AL. (2012)
yJ = 0.22*(logSigma) - 0.26 ;Slope and intercept from Johansson et al. (2009), Fig. 11
oplot, logMstar, yJ, linestyle = 0, thick = 6., color = fsc_color("orange")

;Plot fit to all galaxies:
;xxx = findgen(100)/(100.0/(11.3-10.0)) + 10.0
;FitOFe = POLY_FIT(alog10(stars(w1)*1e10), OFe_mf_comb(w1), 1, SIGMA=FitOFeSigma)
;yyy = FitOFe[0] + FitOFe[1]*(xxx); + FitMgFe[2]*(xxx)*(xxx) + FitMgFe[3]*(xxx)*(xxx)*(xxx)
;oplot, xxx, yyy, thick=6, linestyle = 0, color=fsc_color("black")
comp = WHERE(alog10(stars(w1)*1e10) ge 10.0 AND alog10(stars(w1)*1e10) le 11.3)
xxx = findgen(100)/(100.0/(11.3-10.0)) + 10.0
FitOFe = POLY_FIT(alog10(stars(w1(comp))*1e10), OFe_mf_comb(w1(comp)), 1, SIGMA=FitOFeSigma)
yyy = FitOFe[0] + FitOFe[1]*(xxx); + FitMgFe[2]*(xxx)*(xxx) + FitMgFe[3]*(xxx)*(xxx)*(xxx)
oplot, xxx, yyy, thick=6, linestyle = 0, color=fsc_color("black")

legend, 'Johansson et al. (2009)', box = 0, charsize = 1.25, linestyle = 0, color = fsc_color("orange"), thick = 6, /bottom, /left
legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

IF (MAINELEMENTS eq 1) THEN BEGIN
;-------------------------------
;PLOT CLOSE-UP [Ne/Fe]_Disc+BulgeMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[Ne/Fe]!Ddisk+bulge!N', $
            charsize = csiz, xrange = [9.5, 11.5], yrange = [-0.2, 0.4], xstyle = 1, ystyle = 1
            
oplot, alog10(stars(w1)*1e10), NeFe_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")

;Plot fit to all galaxies:
comp = WHERE(alog10(stars(w1)*1e10) ge 10.0 AND alog10(stars(w1)*1e10) le 11.3)
xxx = findgen(100)/(100.0/(11.3-10.0)) + 10.0
FitNeFe = POLY_FIT(alog10(stars(w1(comp))*1e10), NeFe_mf_comb(w1(comp)), 1, SIGMA=FitNeFeSigma)
yyy = FitNeFe[0] + FitNeFe[1]*(xxx); + FitCFe[2]*(xxx)*(xxx) + FitCFe[3]*(xxx)*(xxx)*(xxx)
oplot, xxx, yyy, thick=6, linestyle = 0, color=fsc_color("black")

legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

ENDIF ;IF (MAINELEMENTS eq 1) THEN BEGIN

;;-------------------------------
;;PLOT [Mg/Fe]_Disc+BulgeMass vs. M*:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[Mg/Fe]!Ddisk+bulge!N', $
;            charsize = csiz, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [MIN(MgFe_mf_disk(w1)), MAX(MgFe_mf_disk(w1))], xstyle = 1, ystyle = 1  ;xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [MIN(MgFe_mf_disk(w1)), MAX(MgFe_mf_disk(w1))]
;            
;oplot, alog10(stars(w1)*1e10), MgFe_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")
;
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
;plots, alog10(stars(w1(0))*1e10), MgFe_mf_comb(w1(0)), psym = 8, symsize = ssiz*3., color = fsc_color("blue")
;;plots, alog10(stars(w1(61))*1e10), MgFe_mf_comb(w1(61)), psym = 8, symsize = ssiz*3., color = fsc_color("red")
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill
;
;legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

;-------------------------------
;PLOT CLOSE-UP [Mg/Fe]_Disc+BulgeMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[Mg/Fe]!Ddisk+bulge!N', $
            charsize = csiz, xrange = [9.5, 11.5], yrange = [-0.44, 1.0], xstyle = 1, ystyle = 1  ;xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [MIN(MgFe_mf_disk(w1)), MAX(MgFe_mf_disk(w1))]
            
oplot, alog10(stars(w1)*1e10), MgFe_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")

;Plot Jonas' and Genevieve's fits:
;logSigma = findgen(100)/(100.0/(2.6-1.2)) + 1.2
logSigma = findgen(100)/(100.0/(2.45-1.8)) + 1.8
logMstar = (2.2*logSigma) + 5.85 ;MY log(sigma) TO log(M*) CONVERTER! USING FIG. 2 FROM WAKE ET AL. (2012)
;logMstar = findgen(100)/(100.0/(11.5-9.0)) + 9.0
;logSigma = (logMstar-5.85)/2.2
yJ = 0.30*(logSigma) - 0.41 ;Slope and intercept from Johansson et al. (2009), Fig. 11
yG = 0.375*(logSigma) - 0.63 ;Slope and intercept derived from Graves et al. (2009), Fig. 4
oplot, logMstar, yJ, linestyle = 0, thick = 6., color = fsc_color("orange")
oplot, logMstar, yG, linestyle = 0, thick = 6., color = fsc_color("red")

;Plot fit to all galaxies:
xxx = findgen(100)/(100.0/(11.5-10.0)) + 10.0
FitMgFe = POLY_FIT(alog10(stars(w1)*1e10), MgFe_mf_comb(w1), 1, SIGMA=FitMgFeSigma)
yyy = FitMgFe[0] + FitMgFe[1]*(xxx); + FitMgFe[2]*(xxx)*(xxx) + FitMgFe[3]*(xxx)*(xxx)*(xxx)
oplot, xxx, yyy, thick=6, linestyle = 0, color=fsc_color("black")

;ShortSFt = WHERE((G(w1).sfh_DiskMass[3] + G(w1).sfh_BulgeMass[3] eq 0.0) AND (G(w1).sfh_DiskMass[4] + G(w1).sfh_BulgeMass[4] eq 0.0) AND (G(w1).sfh_DiskMass[5] + G(w1).sfh_BulgeMass[5] eq 0.0) AND (G(w1).sfh_DiskMass[6] + G(w1).sfh_BulgeMass[6] eq 0.0) AND  (G(w1).sfh_DiskMass[7] + G(w1).sfh_BulgeMass[7] eq 0.0) AND  (G(w1).sfh_DiskMass[8] + G(w1).sfh_BulgeMass[8] eq 0.0) AND  (G(w1).sfh_DiskMass[9] + G(w1).sfh_BulgeMass[9] eq 0.0) AND  (G(w1).sfh_DiskMass[10] + G(w1).sfh_BulgeMass[10] eq 0.0) AND (G(w1).sfh_DiskMass[11] + G(w1).sfh_BulgeMass[11] eq 0.0) AND (G(w1).sfh_DiskMass[12] + G(w1).sfh_BulgeMass[12] eq 0.0) AND (G(w1).sfh_DiskMass[13] + G(w1).sfh_BulgeMass[13] eq 0.0) AND (G(w1).sfh_DiskMass[14] + G(w1).sfh_BulgeMass[14] eq 0.0) AND (G(w1).sfh_DiskMass[15] + G(w1).sfh_BulgeMass[15] eq 0.0))
;print, "n_elements(ShortSFt) = ", n_elements(ShortSFt)
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
;oplot, alog10(stars(w1(ShortSFt))*1e10), MgFe_mf_comb(w1(ShortSFt)), psym = 8, symsize = ssiz*3., color = fsc_color("red")
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill

legend, ['Johansson et al. (2009)', 'Graves et al. (2009)'], box = 0, charsize = 1.25, linestyle = [0,0], color = [fsc_color("orange"),fsc_color("red")], thick = [6,6], /bottom, /left
legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

IF (MAINELEMENTS eq 1) THEN BEGIN
;-------------------------------
;PLOT CLOSE-UP [Si/Fe]_Disc+BulgeMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[Si/Fe]!Ddisk+bulge!N', $
            charsize = csiz, xrange = [9.5, 11.5], yrange = [-0.2, 0.4], xstyle = 1, ystyle = 1
            
oplot, alog10(stars(w1)*1e10), SiFe_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")

;Plot fit to all galaxies:
comp = WHERE(alog10(stars(w1)*1e10) ge 10.0 AND alog10(stars(w1)*1e10) le 11.3)
xxx = findgen(100)/(100.0/(11.3-10.0)) + 10.0
FitSiFe = POLY_FIT(alog10(stars(w1(comp))*1e10), SiFe_mf_comb(w1(comp)), 1, SIGMA=FitSiFeSigma)
yyy = FitSiFe[0] + FitSiFe[1]*(xxx); + FitCaFe[2]*(xxx)*(xxx) + FitCaFe[3]*(xxx)*(xxx)*(xxx)
oplot, xxx, yyy, thick=6, linestyle = 0, color=fsc_color("black")

legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

;-------------------------------
;PLOT CLOSE-UP [S/Fe]_Disc+BulgeMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[S/Fe]!Ddisk+bulge!N', $
            charsize = csiz, xrange = [9.5, 11.5], yrange = [-0.2, 0.4], xstyle = 1, ystyle = 1
            
oplot, alog10(stars(w1)*1e10), SFe_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")

;Plot fit to all galaxies:
comp = WHERE(alog10(stars(w1)*1e10) ge 10.0 AND alog10(stars(w1)*1e10) le 11.3)
xxx = findgen(100)/(100.0/(11.3-10.0)) + 10.0
FitSFe = POLY_FIT(alog10(stars(w1(comp))*1e10), SFe_mf_comb(w1(comp)), 1, SIGMA=FitSiFeSigma)
yyy = FitSFe[0] + FitSFe[1]*(xxx); + FitCaFe[2]*(xxx)*(xxx) + FitCaFe[3]*(xxx)*(xxx)*(xxx)
oplot, xxx, yyy, thick=6, linestyle = 0, color=fsc_color("black")

legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

;-------------------------------
;PLOT CLOSE-UP [Ca/Fe]_Disc+BulgeMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[Ca/Fe]!Ddisk+bulge!N', $
            charsize = csiz, xrange = [9.5, 11.5], yrange = [-0.2, 0.4], xstyle = 1, ystyle = 1
            
oplot, alog10(stars(w1)*1e10), CaFe_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")
;print, "[N/Fe]  ", NFe_mf_comb(w1[0:9])

;Plot Jonas' fit:
logSigma = findgen(100)/(100.0/(2.45-1.8)) + 1.8
logMstar = (2.2*logSigma) + 5.85 ;MY log(sigma) TO log(M*) CONVERTER! USING FIG. 2 FROM WAKE ET AL. (2012)
yJ = 0.15*(logSigma) - 0.23 ;Slope and intercept from Johansson et al. (2009), Fig. 12
oplot, logMstar, yJ, linestyle = 0, thick = 6., color = fsc_color("orange")

;Plot fit to all galaxies:
comp = WHERE(alog10(stars(w1)*1e10) ge 10.0 AND alog10(stars(w1)*1e10) le 11.3)
xxx = findgen(100)/(100.0/(11.3-10.0)) + 10.0
FitCaFe = POLY_FIT(alog10(stars(w1(comp))*1e10), CaFe_mf_comb(w1(comp)), 1, SIGMA=FitCaFeSigma)
yyy = FitCaFe[0] + FitCaFe[1]*(xxx); + FitCaFe[2]*(xxx)*(xxx) + FitCaFe[3]*(xxx)*(xxx)*(xxx)
oplot, xxx, yyy, thick=6, linestyle = 0, color=fsc_color("black")

legend, 'Johansson et al. (2009)', box = 0, charsize = 1.25, linestyle = 0, color = fsc_color("orange"), thick = 6, /bottom, /left
legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

;-------------------------------
;PLOT CLOSE-UP [alpha/Fe]_Disc+BulgeMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[alpha/Fe]!Ddisk+bulge!N', $
            charsize = csiz, xrange = [9.5, 11.5], yrange = [0.0, 0.6], xstyle = 1, ystyle = 1
            
oplot, alog10(stars(w1)*1e10), aFe_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")
;print, "[N/Fe]  ", NFe_mf_comb(w1[0:9])

;;Plot Jonas' fit:
;logSigma = findgen(100)/(100.0/(2.45-1.8)) + 1.8
;logMstar = (2.2*logSigma) + 5.85 ;MY log(sigma) TO log(M*) CONVERTER! USING FIG. 2 FROM WAKE ET AL. (2012)
;yJ = 0.15*(logSigma) - 0.23 ;Slope and intercept from Johansson et al. (2009), Fig. 12
;oplot, logMstar, yJ, linestyle = 0, thick = 6., color = fsc_color("orange")

;Plot fit to all galaxies:
comp = WHERE(alog10(stars(w1)*1e10) ge 10.0 AND alog10(stars(w1)*1e10) le 11.3)
xxx = findgen(100)/(100.0/(11.3-10.0)) + 10.0
FitaFe = POLY_FIT(alog10(stars(w1(comp))*1e10), aFe_mf_comb(w1(comp)), 1, SIGMA=FitCaFeSigma)
yyy = FitaFe[0] + FitaFe[1]*(xxx); + FitaFe[2]*(xxx)*(xxx) + FitaFe[3]*(xxx)*(xxx)*(xxx)
oplot, xxx, yyy, thick=6, linestyle = 0, color=fsc_color("black")

;legend, 'Johansson et al. (2009)', box = 0, charsize = 1.25, linestyle = 0, color = fsc_color("orange"), thick = 6, /bottom, /left
legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

ENDIF ;IF (MAINELEMENTS eq 1) THEN BEGIN

;-------------------------------
;PLOT CLOSE-UP [Mg/H]_Disc+BulgeMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[Mg/H]!Ddisk+bulge!N', $
            charsize = csiz, xrange = [9.5, 11.5], yrange = [-0.6, 0.3], xstyle = 1, ystyle = 1  ;xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [MIN(MgFe_mf_disk(w1)), MAX(MgFe_mf_disk(w1))]
            
oplot, alog10(stars(w1)*1e10), MgH_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")

;Plot Genevieve's fit:
logSigma = findgen(100)/(100.0/(2.4-1.9)) + 1.9
logMstar = (2.2*logSigma) + 5.85 ;MY log(sigma) TO log(M*) CONVERTER! USING FIG. 2 FROM WAKE ET AL. (2012)
yG = 0.7*(logSigma) - 1.48 ;Slope and intercept derived from Graves et al. (2009), Fig. 4
oplot, logMstar, yG, linestyle = 0, thick = 6., color = fsc_color("red")

legend, ['Graves et al. (2009)'], box = 0, charsize = 1.25, linestyle = [0], color = [fsc_color("red")], thick = 6, /bottom, /left
legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

;-------------------------------
;PLOT CLOSE-UP [Fe/H]_Disc+BulgeMass vs. M*:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'log(M!D*!N) [M!D!9'+string(110B)+'!3!N]', ytitle = '[Fe/H]!Ddisk+bulge!N', $
            charsize = csiz, xrange = [9.5, 11.5], yrange = [-0.6, 0.3], xstyle = 1, ystyle = 1  ;xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [MIN(MgFe_mf_disk(w1)), MAX(MgFe_mf_disk(w1))]
            
oplot, alog10(stars(w1)*1e10), FeH_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")

;Plot Genevieve's fit:
logSigma = findgen(100)/(100.0/(2.4-1.9)) + 1.9
logMstar = (2.2*logSigma) + 5.85 ;MY log(sigma) TO log(M*) CONVERTER! USING FIG. 2 FROM WAKE ET AL. (2012)
yG = 0.35*(logSigma) - 0.51 ;Slope and intercept derived from Graves et al. (2009), Fig. 4
oplot, logMstar, yG, linestyle = 0, thick = 6., color = fsc_color("red")

;Plot Jonas' fit:
logSigma = findgen(100)/(100.0/(2.45-1.8)) + 1.8
logMstar = (2.2*logSigma) + 5.85 ;MY log(sigma) TO log(M*) CONVERTER! USING FIG. 2 FROM WAKE ET AL. (2012)
yJ = -0.03*(logSigma) - 0.02 ;Slope and intercept from Johansson et al. (2009), Fig. 12
oplot, logMstar, yJ, linestyle = 0, thick = 6., color = fsc_color("orange")

legend, ['Johansson et al. (2009)', 'Graves et al. (2009)'], box = 0, charsize = 1.25, linestyle = [0,0], color = [fsc_color("orange"),fsc_color("red")], thick = [6,6], /bottom, /left
legend, ['Ellipticals'], box = 0, charsize = 1.75, linestyle = [-99], /top, /right

;-------------------------------
;PLOT deltaM* vs. Lookback time:
;-------------------------------
sfhtime = G(w1).sfh_time/1.0e9 ;Lookback time from z=0 to middle of SFH bin, in Gyrs. It's the same for all gals.
;print, " 1 = ", G(w1[0]).sfh_time/1.0e9
;print, " 2 = ", G(w1[10]).sfh_time/1.0e9
;
SFH_diskE = ((G(w1).sfh_DiskMass+G(w1).sfh_BulgeMass)*1.0e10/Hubble_h)
;print, "SFH_diskE[*,0] = ", SFH_diskE[*,0]
;whistE = WHERE(G(w1).sfh_DiskMass eq 0.0)
;SFH_diskE(whistE) = -99.0

;print, "G(w1[0]).sfh_DiskMass = ", G(w1[0]).sfh_DiskMass*1.0e10/Hubble_h
;print, "G(w1[0]).sfh_DiskMass/sfhdt = ", (G(w1[0]).sfh_DiskMass*1.0e10/Hubble_h)/sfhdt

;wcleanE = WHERE(SFH_diskE ne -99.0)
;wdirtyE = WHERE(SFH_diskE eq -99.0)

plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = 'Delta M!D*!N [M!D!9'+string(110B)+'!3!N]', $ ;ytitle = 'log(Delta M!D*!N) [M!D!9'+string(110B)+'!3!N]'
            charsize = csiz, xrange = [MAX(sfhtime), MIN(sfhtime)], yrange = [0.0, 12.0], xstyle = 1, ystyle = 1 ;yrange = [MIN(TotalSFR), MAX(TotalSFR)]

;oplot, sfhtime, alog10(SFH_diskE), psym = 8, symsize = ssiz*3.0, color = fsc_color("black")
;oplot, sfhtime, alog10(SFH_diskE), linestyle = 0, thick = 6., color = fsc_color("black")
  
;------------------  
  
;;lMass = WHERE(alog10(((G(w1(wcleanE)).DiskMass+G(w1(wcleanE)).BulgeMass)*1.0e10/Hubble_h)) lt 10.0)
;;hMass = WHERE(alog10(((G(w1(wcleanE)).DiskMass+G(w1(wcleanE)).BulgeMass)*1.0e10/Hubble_h)) ge 10.0)  
lMass = WHERE(alog10(stars(w1)*1e10) lt 10.0)
hMass = WHERE(alog10(stars(w1)*1e10) ge 10.0)
print, "n_elements(lMass) = ", n_elements(lMass), "n_elements(hMass) = ", n_elements(hMass)

sfhtimelM = G(w1(lMass)).sfh_time/1.0e9
SFH_diskElM = alog10((G(w1(lMass)).sfh_DiskMass+G(w1(lMass)).sfh_BulgeMass)*1.0e10/Hubble_h)
sfhtimehM = G(w1(hMass)).sfh_time/1.0e9
SFH_diskEhM = alog10((G(w1(hMass)).sfh_DiskMass+G(w1(hMass)).sfh_BulgeMass)*1.0e10/Hubble_h)
;  
;;plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = 'log(Delta M!D*!N) [M!D!9'+string(110B)+'!3!N]', $ ;'log(Stars formed!Ddisk!N) [10!E10!N Msun]'
;;            charsize = csiz, xrange = [MAX(sfhtime), MIN(sfhtime)], yrange = [0.0, 12.0], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean)), MAX(OFe_mf_disk_hist)]
;
oplot, sfhtimehM, SFH_diskEhM, psym = 8, symsize = ssiz*3.0, color = fsc_color("red")
oplot, sfhtimehM, SFH_diskEhM, linestyle = 0, thick = 6., color = fsc_color("red")
oplot, sfhtimelM, SFH_diskElM, psym = 8, symsize = ssiz*3.0, color = fsc_color("blue")
oplot, sfhtimelM, SFH_diskElM, linestyle = 0, thick = 6., color = fsc_color("blue")
;print, "SFH_diskE(hMass[0]) = ", SFH_diskE(hMass[0])

;-------------------------------
;PLOT [Mg/Fe] vs. Lookback time:
;-------------------------------
;MgFe_mf_disk_hist = alog10(G(w1).sfh_ElementsDiskMass[3,*]/G(w1).sfh_ElementsDiskMass[4,*]) - alog10(Mg_mf/Fe_mf)
;print, "MgFe_mf_disk_hist = ", MgFe_mf_disk_hist[0,*]
if (MAINELEMENTS eq 0) THEN BEGIN
MgFe_mf_disk_histlM = alog10(G(w1(lMass)).sfh_ElementsDiskMass[3,*]/G(w1(lMass)).sfh_ElementsDiskMass[4,*]) - alog10(Mg_mf/Fe_mf)
MgFe_mf_disk_histhM = alog10(G(w1(hMass)).sfh_ElementsDiskMass[3,*]/G(w1(hMass)).sfh_ElementsDiskMass[4,*]) - alog10(Mg_mf/Fe_mf)
;print, "MgFe_mf_disk_histlM[0,*] = ", MgFe_mf_disk_histlM[0,*]
;print, "MgFe_mf_disk_histhM[0,*] = ", MgFe_mf_disk_histhM[0,*]
ENDIF ELSE IF (MAINELEMENTS eq 1) THEN BEGIN
MgFe_mf_disk_histlM = alog10(G(w1(lMass)).sfh_ElementsDiskMass[6,*]/G(w1(lMass)).sfh_ElementsDiskMass[10,*]) - alog10(Mg_mf/Fe_mf)
MgFe_mf_disk_histhM = alog10(G(w1(hMass)).sfh_ElementsDiskMass[6,*]/G(w1(hMass)).sfh_ElementsDiskMass[10,*]) - alog10(Mg_mf/Fe_mf)
ENDIF
print, n_elements(MgFe_mf_disk_histlM)

plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = '[Mg/Fe]!Ddisc+bulge!N', $ ;ytitle = 'log(Delta M!D*!N) [M!D!9'+string(110B)+'!3!N]'
            charsize = csiz, xrange = [MAX(sfhtime), MIN(sfhtime)], yrange = [-1.0, 2.0], xstyle = 1, ystyle = 1

oplot, sfhtimelM, MgFe_mf_disk_histlM, psym = 8, symsize = ssiz*3.0, color = fsc_color("blue")
oplot, sfhtimelM, MgFe_mf_disk_histlM, linestyle = 0, thick = 6., color = fsc_color("blue")
oplot, sfhtimehM, MgFe_mf_disk_histhM, psym = 8, symsize = ssiz*3.0, color = fsc_color("red")
oplot, sfhtimehM, MgFe_mf_disk_histhM, linestyle = 0, thick = 6., color = fsc_color("red")
;oplot, G(w1(lMass[0])).sfh_time/1.0e9, MgFe_mf_disk_histlM[0,*], psym = 8, symsize = ssiz*3.0, color = fsc_color("sky blue")
;oplot, G(w1(lMass[0])).sfh_time/1.0e9, MgFe_mf_disk_histlM[0,*], linestyle = 0, thick = 6., color = fsc_color("sky blue")
;oplot, G(w1(hMass[0])).sfh_time/1.0e9, MgFe_mf_disk_histhM[0,*], psym = 8, symsize = ssiz*3.0, color = fsc_color("pink")
;oplot, G(w1(hMass[0])).sfh_time/1.0e9, MgFe_mf_disk_histhM[0,*], linestyle = 0, thick = 6., color = fsc_color("pink")

;highMgFe = WHERE(MgFe_mf_comb(w1(hMass)) gt 0.7)
;oplot, G(w1(hMass[highMgFe])).sfh_time/1.0e9, MgFe_mf_disk_histhM[highMgFe,*], psym = 8, symsize = ssiz*3.0, color = fsc_color("pink")
;oplot, G(w1(hMass[highMgFe])).sfh_time/1.0e9, MgFe_mf_disk_histhM[highMgFe,*], linestyle = 0, thick = 6., color = fsc_color("pink")


;;-------------------------------
;;PLOT SFH vs. Lookback time:
;;-------------------------------
;;sfhtime = G(w1).sfh_time/1.0e9 ;Lookback time from z=0 to middle of SFH bin, in Gyrs.
;sfhdt = [1.905159,3.218720,2.268308,2.395445,1.147430,1.070113,0.498352,0.469735,0.223084,0.110342,0.105373,0.052686,0.052686,0.026343,0.013172,0.013172,0.000000,0.000000,0.000000,0.000000]*1.0e9
;;sfhdt = [1.905159,3.218720,2.268308,2.395445,1.147430,1.070113,0.498352,0.469735,0.223084,0.110342,0.105373,0.052686,0.052686,0.026343,0.013172,0.013172,1.000000,1.000000,1.000000,1.000000]*1.0e9
;DiskSFR = fltarr(n_elements(G(w1)),20)
;BulgeSFR = fltarr(n_elements(G(w1)),20)
;sfhtime = fltarr(n_elements(G(w1)),20)
;FOR i=0,i<n_elements(G(w1)) DO BEGIN
;  DiskSFR[i,*] = ((G(w1[i]).sfh_DiskMass[*])*1.0e10/Hubble_h)/sfhdt[*]
;  BulgeSFR[i,*] = ((G(w1[i]).sfh_BulgeMass[*])*1.0e10/Hubble_h)/sfhdt[*]
;  sfhtime[i,*] = G(w1[i]).sfh_time/1.0e9
;ENDFOR
;TotalSFR = DiskSFR+BulgeSFR
;;;print, sfhtime[0,*]
;;print, "G(w1[0]).sfh_time = ", G(w1[0]).sfh_time/1.0e9
;;print, "sfhtime[*,0] = ", sfhtime[*,0]
;print, "G(w1[0]).sfh_DiskMass = ", G(w1[0]).sfh_DiskMass*1.0e10/Hubble_h
;print, "G(w1[0]).sfh_DiskMass/sfhdt = ", (G(w1[0]).sfh_DiskMass*1.0e10/Hubble_h)/sfhdt
;print,  "DiskSFR[0,*] = ", DiskSFR[0,*]
;print,  "BulgeSFR[0,*] = ", BulgeSFR[0,*]
;print,  "TotalSFR[0,*] = ", TotalSFR[0,*]
;;
;;SFH_diskE = ((G(w1).sfh_DiskMass+G(w1).sfh_BulgeMass)*1.0e10/Hubble_h)
;;SFH_diskE = ((G(w1).sfh_DiskMass+G(w1).sfh_BulgeMass)*1.0e10/Hubble_h)
;;print, "SFH_diskE[*,0] = ", SFH_diskE[*,0]
;;whistE = WHERE(G(w1).sfh_DiskMass eq 0.0)
;;SFH_diskE(whistE) = -99.0
;
;;print, "G(w1[0]).sfh_DiskMass = ", G(w1[0]).sfh_DiskMass*1.0e10/Hubble_h
;;print, "G(w1[0]).sfh_DiskMass/sfhdt = ", (G(w1[0]).sfh_DiskMass*1.0e10/Hubble_h)/sfhdt
;
;;wcleanE = WHERE(SFH_diskE ne -99.0)
;;wdirtyE = WHERE(SFH_diskE eq -99.0)
;
;plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = 'SFR [M!D!9'+string(110B)+'!3!N/yr]', $ ;ytitle = 'log(Delta M!D*!N) [M!D!9'+string(110B)+'!3!N]'
;            charsize = csiz, xrange = [MAX(sfhtime), MIN(sfhtime)], yrange = [0.0, 1.0], xstyle = 1, ystyle = 1 ;yrange = [MIN(TotalSFR), MAX(TotalSFR)]
;
;FOR i=0,i<n_elements(G(w1)) DO BEGIN
;oplot, sfhtime[i,*], TotalSFR[i,*], psym = 8, symsize = ssiz*3.0, color = fsc_color("black")
;oplot, sfhtime[i,*], TotalSFR[i,*], linestyle = 0, thick = 6., color = fsc_color("black")
;ENDFOR
;;oplot, sfhtime, SFH_diskE, psym = 8, symsize = ssiz*3.0, color = fsc_color("black")
;;oplot, sfhtime, SFH_diskE, linestyle = 0, thick = 6., color = fsc_color("black")
;;oplot, sfhtime[*,0], SFH_diskE[*,0], psym = 8, symsize = ssiz*3.0, color = fsc_color("black")
;;oplot, sfhtime[*,0], SFH_diskE[*,0], linestyle = 0, thick = 6., color = fsc_color("black")
;
;
;  
;------------------  
  
;;lMass = WHERE(alog10(((G(w1(wcleanE)).DiskMass+G(w1(wcleanE)).BulgeMass)*1.0e10/Hubble_h)) lt 10.0)
;;hMass = WHERE(alog10(((G(w1(wcleanE)).DiskMass+G(w1(wcleanE)).BulgeMass)*1.0e10/Hubble_h)) ge 10.0)  
;lMass = WHERE(alog10(stars(w1)*1e10) lt 10.0)
;hMass = WHERE(alog10(stars(w1)*1e10) ge 10.0)
;print, "n_elements(hMass) = ", n_elements(hMass)
;
;sfhtimelM = G(w1(lMass)).sfh_time/1.0e9
;SFH_diskElM = alog10(((G(w1(lMass)).sfh_DiskMass+G(w1(lMass)).sfh_BulgeMass)*1.0e10/Hubble_h)/sfhdt)
;sfhtimehM = G(w1(hMass)).sfh_time/1.0e9
;SFH_diskEhM = alog10(((G(w1(hMass)).sfh_DiskMass+G(w1(hMass)).sfh_BulgeMass)*1.0e10/Hubble_h)/sfhdt)
;  
;;plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = 'log(Delta M!D*!N) [M!D!9'+string(110B)+'!3!N]', $ ;'log(Stars formed!Ddisk!N) [10!E10!N Msun]'
;;            charsize = csiz, xrange = [MAX(sfhtime), MIN(sfhtime)], yrange = [0.0, 12.0], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean)), MAX(OFe_mf_disk_hist)]
;
;oplot, sfhtimehM, SFH_diskEhM, psym = 8, symsize = ssiz*3.0, color = fsc_color("red")
;oplot, sfhtimehM, SFH_diskEhM, linestyle = 0, thick = 6., color = fsc_color("red")
;oplot, sfhtimelM, SFH_diskElM, psym = 8, symsize = ssiz*3.0, color = fsc_color("blue")
;oplot, sfhtimelM, SFH_diskElM, linestyle = 0, thick = 6., color = fsc_color("blue")
;print, "SFH_diskE(hMass[0]) = ", SFH_diskE(hMass[0])

;;-------------------------------
;;PLOT [alpha/Fe]_Disc+BulgeMass vs. M*:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = 'log M!Dstellar!N / M!D!9'+string(110B)+'!3!N', ytitle = '[alpha/Fe]!Ddisk+bulge!N', $
;            charsize = 1.0, xrange = [MIN(alog10(stars(w1)*1e10)), MAX(alog10(stars(w1)*1e10))], yrange = [0.0,1.0], xstyle = 1, ystyle = 1 ;[MIN(aFe_mf_disk(w0)), MAX(aFe_mf_disk(w0))]
;            
;;oplot, alog10(stars(w1)*1e10), aFe_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")
;
;;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
;;plots, alog10(stars(w1(0))*1e10), aFe_mf_comb(w1(0)), psym = 8, symsize = ssiz*3., color = fsc_color("blue")
;;plots, alog10(stars(2434)*1e10), aFe_mf_comb(2434), psym = 8, symsize = ssiz*3., color = fsc_color("red")
;;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill
;
;legend, ['No alpha with MAINELEMENTS'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left
;
;;-------------------------------
;;PLOT [O/Fe]_Disc+BulgeMass vs. [Fe/H]_combMass:
;;-------------------------------
;plot, findgen(10), /nodata, xtitle = '[Fe/H]!DBulge!N', ytitle = '[O/Fe]!Ddisk+bulge!N', $
;            charsize = 1.0, xrange = [MIN(FeH_mf_comb(w1)), MAX(FeH_mf_comb(w1))], yrange = [MIN(OFe_mf_comb(w1)), MAX(OFe_mf_comb(w1))], xstyle = 1, ystyle = 1 
;            
;oplot, FeH_mf_comb(w1), OFe_mf_comb(w1), psym = 8, symsize = ssiz, color = fsc_color("black")
;
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
;plots, FeH_mf_comb(w1(0)), OFe_mf_comb(w1(0)), psym = 8, symsize = ssiz*3., color = fsc_color("blue")
;plots, FeH_mf_comb(2434), OFe_mf_comb(2434), psym = 8, symsize = ssiz*3., color = fsc_color("red")
;plots, FeH_mf_comb(w1(200)), OFe_mf_comb(w1(200)), psym = 8, symsize = ssiz*3., color = fsc_color("green")
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill
;
;legend, ['Ellipticals'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left

;------------------------------------------------------------------
;------------------------------------------------------------------
;DISK AND BULGE HISTORIES (ELLIPTICALS: Galaxy 0):
;------------------------------------------------------------------
;------------------------------------------------------------------
IF (MAINELEMENTS eq 2) THEN BEGIN

sfh_dt = [5.12,2.56,2.56,1.28,0.64,0.64,0.32,0.16,0.08,0.08,0.04,0.04,0.02,0.01,0.01,0.01,0.0,0.0,0.0,0.0]*1.0e9 ;IN YEARS

;sfh_time = alog10(G(w1(0)).sfh_time[*]) 
;IF (NEWSESTSTRUCT eq 1) THEN BEGIN
;;sfh_time1 = (G(w1(0)).sfh_time[*])/1.0e9
;sfh_time1 = (1.36e10 - ((G(w1(0)).sfh_time[*] * 0.73 / (10.0^(11.978))) - (0.5*sfh_dt[*])) - (0.5*sfh_dt[*]))/1.0e9
;ENDIF ELSE BEGIN
sfh_time1 = (G(w1(0)).sfh_time[*]+(0.5*sfh_dt[*]))/1.0e9
;ENDELSE
;
;sfh_time1 = (G(w1(0)).sfh_time[*]+(0.5*sfh_dt[*]))/1.0e9
whist1 = WHERE(G(w1(0)).sfh_DiskMass eq 0.0)
bhist1 = WHERE(G(w1(0)).sfh_BulgeMass eq 0.0)
;-------
OFe_mf_disk_hist1 = alog10(G(w1(0)).sfh_ElementsDiskMass[2,*]/G(w1(0)).sfh_ElementsDiskMass[4,*]) - alog10(O_mf/Fe_mf)
;print, G(w1(0)).sfh_ElementsDiskMass[2,*]/G(w1(0)).sfh_ElementsDiskMass[4,*]
;FOR i=0, 19 DO BEGIN
;    print, "O: ", G(w1(0)).sfh_ElementsDiskMass[2,i]
;    print, "Fe:", G(w1(0)).sfh_ElementsDiskMass[4,i]
;    print, "--"
;ENDFOR  
OFe_mf_disk_hist1(whist1) = -99.0

MgFe_mf_disk_hist1 = alog10(G(w1(0)).sfh_ElementsDiskMass[3,*]/G(w1(0)).sfh_ElementsDiskMass[4,*]) - alog10(Mg_mf/Fe_mf)
MgFe_mf_disk_hist1(whist1) = -99.0

FeH_mf_disk_hist1 = alog10(G(w1(0)).sfh_ElementsDiskMass[4,*]/G(w1(0)).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf)
FeH_mf_disk_hist1(whist1) = -99.0

;SNII_disk_hist1 = G(w1(0)).sfh_MetalsDiskMass[1,*]/G(w1(0)).sfh_DiskMass[*]
SNII_disk_hist1 = G(w1(0)).sfh_MetalsDiskMass[1,*]/(G(w1(0)).sfh_MetalsDiskMass[0,*]+G(w1(0)).sfh_MetalsDiskMass[1,*]+G(w1(0)).sfh_MetalsDiskMass[2,*])
SNII_disk_hist1(whist1) = -99.0

SNIa_disk_hist1 = G(w1(0)).sfh_MetalsDiskMass[0,*]/(G(w1(0)).sfh_MetalsDiskMass[0,*]+G(w1(0)).sfh_MetalsDiskMass[1,*]+G(w1(0)).sfh_MetalsDiskMass[2,*])
SNIa_disk_hist1(whist1) = -99.0

AGB_disk_hist1 = G(w1(0)).sfh_MetalsDiskMass[2,*]/(G(w1(0)).sfh_MetalsDiskMass[0,*]+G(w1(0)).sfh_MetalsDiskMass[1,*]+G(w1(0)).sfh_MetalsDiskMass[2,*])
AGB_disk_hist1(whist1) = -99.0

SFH_disk1 = alog10(((G(w1(0)).sfh_DiskMass[*]+G(w1(0)).sfh_BulgeMass[*])*1.0e10/Hubble_h)/sfh_dt)
SFH_disk1(whist1) = -99.0
;print, SFH_disk1
;print, G(w1(0)).sfr, alog10(G(w1(0)).sfr)
;;--------
;OFe_mf_bulge_hist = alog10(G(w1(0)).sfh_ElementsBulgeMass[2,*]/G(w1(0)).sfh_ElementsBulgeMass[4,*]) - alog10(O_mf/Fe_mf)
;OFe_mf_bulge_hist(bhist) = 0.0
;
;MgFe_mf_bulge_hist = alog10(G(w1(0)).sfh_ElementsBulgeMass[3,*]/G(w1(0)).sfh_ElementsBulgeMass[4,*]) - alog10(Mg_mf/Fe_mf)
;MgFe_mf_bulge_hist(bhist) = 0.0
;
;FeH_mf_bulge_hist = alog10(G(w1(0)).sfh_ElementsBulgeMass[4,*]/G(w1(0)).sfh_ElementsBulgeMass[0,*]) - alog10(Fe_mf/H_mf)
;FeH_mf_bulge_hist(bhist) = 0.0
;
;SNII_bulge_hist = G(w1(0)).sfh_MetalsBulgeMass[1,*]/G(w1(0)).sfh_BulgeMass[*]
;SNII_bulge_hist(bhist) = -99.0
;
;SNIa_bulge_hist = G(w1(0)).sfh_MetalsBulgeMass[0,*]/G(w1(0)).sfh_BulgeMass[*]
;SNIa_bulge_hist(bhist) = -99.0
;
;AGB_bulge_hist = G(w1(0)).sfh_MetalsBulgeMass[2,*]/G(w1(0)).sfh_BulgeMass[*]
;AGB_bulge_hist(bhist) = -99.0
;
;SFH_bulge = alog10(G(w1(0)).sfh_BulgeMass[*])
;SFH_bulge(bhist) = 0.0
;

;IF (NEWSESTSTRUCT eq 1) THEN BEGIN
;;sfh_time2 = (G(w1(61)).sfh_time[*])/1.0e9
;sfh_time2 = (1.36e10 - ((G(w1(61)).sfh_time[*] * 0.73 / (10.0^(11.978))) - (0.5*sfh_dt[*])) - (0.5*sfh_dt[*]))/1.0e9
;print, sfh_time2
;ENDIF ELSE BEGIN
sfh_time2 = (G(w1(61)).sfh_time[*]+(0.5*sfh_dt[*]))/1.0e9
;ENDELSE

whist2 = WHERE(G(w1(61)).sfh_DiskMass eq 0.0)
bhist2 = WHERE(G(w1(61)).sfh_BulgeMass eq 0.0)
;-------
;Ohist = G(w1(61)).sfh_ElementsDiskMass[2,*]
;Fehist = G(w1(61)).sfh_ElementsDiskMass[4,*]

OFe_mf_disk_hist2 = alog10(G(w1(61)).sfh_ElementsDiskMass[2,*]/G(w1(61)).sfh_ElementsDiskMass[4,*]) - alog10(O_mf/Fe_mf)
;print, OFe_mf_disk_hist2, alog10(Ohist/Fehist) - alog10(O_mf/Fe_mf)
OFe_mf_disk_hist2(whist2) = -99.0

MgFe_mf_disk_hist2 = alog10(G(w1(61)).sfh_ElementsDiskMass[3,*]/G(w1(61)).sfh_ElementsDiskMass[4,*]) - alog10(Mg_mf/Fe_mf)
MgFe_mf_disk_hist2(whist2) = -99.0

FeH_mf_disk_hist2 = alog10(G(w1(61)).sfh_ElementsDiskMass[4,*]/G(w1(61)).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf)
FeH_mf_disk_hist2(whist2) = -99.0

;SNII_disk_hist2 = G(w1(61)).sfh_MetalsDiskMass[1,*]/G(w1(61)).sfh_DiskMass[*]
SNII_disk_hist2 = G(w1(61)).sfh_MetalsDiskMass[1,*]/(G(w1(61)).sfh_MetalsDiskMass[0,*]+G(w1(61)).sfh_MetalsDiskMass[1,*]+G(w1(61)).sfh_MetalsDiskMass[2,*])
SNII_disk_hist2(whist2) = -99.0

SNIa_disk_hist2 = G(w1(61)).sfh_MetalsDiskMass[0,*]/(G(w1(61)).sfh_MetalsDiskMass[0,*]+G(w1(61)).sfh_MetalsDiskMass[1,*]+G(w1(61)).sfh_MetalsDiskMass[2,*])
SNIa_disk_hist2(whist2) = -99.0

AGB_disk_hist2 = G(w1(61)).sfh_MetalsDiskMass[2,*]/(G(w1(61)).sfh_MetalsDiskMass[0,*]+G(w1(61)).sfh_MetalsDiskMass[1,*]+G(w1(61)).sfh_MetalsDiskMass[2,*])
AGB_disk_hist2(whist2) = -99.0

SFH_disk2 = alog10(((G(w1(61)).sfh_DiskMass[*]+G(w1(61)).sfh_BulgeMass[*])*1.0e10/Hubble_h)/sfh_dt)
SFH_disk2(whist2) = -99.0
SFH_disk2(bhist2) = -99.0
;print, SFH_disk1
;print, G(w1(61)).sfr, alog10(G(w1(61)).sfr)
;print, SFH_disk2
;print, G(w1(61)).sfh_DiskMass
;print, G(w1(61)).sfh_BulgeMass

;FOR i=0, 20 DO BEGIN
;    print, TOTAL(G(w1(0)).sfh_ElementsDiskMass[*,i])
;    print, "--", G(w1(0)).sfh_DiskMass[i]*1.0e10/Hubble_h
;ENDFOR    

;-------------------------------
;PLOT [O/Fe]_DiskMass vs. Lookback time:
;-------------------------------
wclean1 = WHERE(OFe_mf_disk_hist1 ne -99.0)
wdirty1 = WHERE(OFe_mf_disk_hist1 eq -99.0)
wclean2 = WHERE(OFe_mf_disk_hist2 ne -99.0)
wdirty2 = WHERE(OFe_mf_disk_hist2 eq -99.0)
;bclean = WHERE(OFe_mf_bulge_hist ne 0.0)
;bdirty = WHERE(OFe_mf_bulge_hist eq 0.0)

;plot, findgen(10), /nodata, xtitle = 'log(Lookback Time) [yrs]', ytitle = '[O/Fe]!Ddisk!N', $
;            charsize = 1.0, xrange = [MAX(sfh_time), MIN(sfh_time)], yrange = [MIN(OFe_mf_disk_hist(wclean))-0.001, MAX(OFe_mf_disk_hist(wclean))+0.001], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean)), MAX(OFe_mf_disk_hist)]
plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = '[O/Fe]!Ddisk!N', $ ;xtitle = 'log(Lookback Time) [yrs]'
            charsize = csiz, xrange = [MAX(sfh_time1), MIN(sfh_time1)], yrange = [MIN(OFe_mf_disk_hist1(wclean1))-0.01, MAX(OFe_mf_disk_hist1(wclean1))+0.01], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean))-0.001, MAX(OFe_mf_disk_hist(wclean))+0.001]

oplot, sfh_time1(wclean1), OFe_mf_disk_hist1(wclean1), linestyle = 0, thick = 6., color = fsc_color("blue")
oplot, sfh_time1, OFe_mf_disk_hist1, psym = 8, symsize = ssiz*3.0, color = fsc_color("blue")
;oplot, sfh_time2(wclean2), OFe_mf_disk_hist2(wclean2), linestyle = 0, thick = 6., color = fsc_color("red")
oplot, sfh_time2, OFe_mf_disk_hist2, psym = 8, symsize = ssiz*3.0, color = fsc_color("red")

;legend, ['Galaxy 0'], box = 0, charsize = 1.75, textcolor = fsc_color("blue"), linestyle = [-99], /top, /right

;-------------------------------
;PLOT O_DiskMass vs. Lookback time:
;-------------------------------
plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = 'X!Dmass fraction!N', $ ;xtitle = 'log(Lookback Time) [yrs]'
            charsize = csiz, xrange = [MAX(sfh_time1), MIN(sfh_time1)], yrange = [0.0, 0.012], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean))-0.001, MAX(OFe_mf_disk_hist(wclean))+0.001]
oplot, sfh_time1, (G(w1(0)).sfh_ElementsDiskMass[2,*]/(G(w1(0)).sfh_DiskMass[*]*1.0e10/Hubble_h)), psym = 8, symsize = ssiz*3.0, color = fsc_color("blue")
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
oplot, sfh_time1, (G(w1(0)).sfh_ElementsDiskMass[4,*]/(G(w1(0)).sfh_DiskMass[*]*1.0e10/Hubble_h)), psym = 8, symsize = ssiz*3.0, color = fsc_color("blue")
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill
;print, (G(w1(0)).sfh_ElementsDiskMass[2,*]/(G(w1(0)).sfh_DiskMass[*]*1.0e10/Hubble_h))
;print, (G(w1(0)).sfh_ElementsDiskMass[2,*]/(G(w1(0)).sfh_DiskMass[*]*1.0e10/Hubble_h))/(G(w1(0)).sfh_ElementsDiskMass[4,*]/(G(w1(0)).sfh_DiskMass[*]*1.0e10/Hubble_h))
oplot, sfh_time2, (G(w1(61)).sfh_ElementsDiskMass[2,*]/(G(w1(61)).sfh_DiskMass[*]*1.0e10/Hubble_h)), psym = 8, symsize = ssiz*3.0, color = fsc_color("red")
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
oplot, sfh_time2, (G(w1(61)).sfh_ElementsDiskMass[4,*]/(G(w1(61)).sfh_DiskMass[*]*1.0e10/Hubble_h)), psym = 8, symsize = ssiz*3.0, color = fsc_color("red")
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill

;print, (G(w1(1347)).sfh_ElementsDiskMass[2,*]/(G(w1(1347)).sfh_DiskMass[*]*1.0e10/Hubble_h))
legend, ['Oxygen'], box = 0, charsize = 1.75, psym = 8, symsize = ssiz*3., color = fsc_color("black"), /top, /left
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
legend, ['Iron'], box = 0, charsize = 1.75, psym = 8, symsize = ssiz*3., color = fsc_color("black"), /bottom, /left
usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill
;legend, ['Oxygen','Iron'], box = 0, charsize = 1.75, psym = [8,88], symsize = [ssiz*3.,ssiz*3.], color = [fsc_color("black"),fsc_color("black")], /top, /left
;-------------------------------
;PLOT [Mg/Fe]_DiskMass vs. Lookback time:
;-------------------------------
wclean1 = WHERE(MgFe_mf_disk_hist1 ne -99.0)
wdirty1 = WHERE(MgFe_mf_disk_hist1 eq -99.0)
wclean2 = WHERE(MgFe_mf_disk_hist2 ne -99.0)
wdirty2 = WHERE(MgFe_mf_disk_hist2 eq -99.0)

plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = '[Mg/Fe]!Ddisk!N', $ ;xtitle = 'log(Lookback Time) [yrs]'
            charsize = csiz, xrange = [MAX(sfh_time1), MIN(sfh_time1)], yrange = [MIN(MgFe_mf_disk_hist1(wclean1))-0.01, MAX(MgFe_mf_disk_hist1(wclean1))+0.01], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean))-0.001, MAX(OFe_mf_disk_hist(wclean))+0.001]

oplot, sfh_time1(wclean1), MgFe_mf_disk_hist1(wclean1), linestyle = 0, thick = 6., color = fsc_color("blue")
oplot, sfh_time1, MgFe_mf_disk_hist1, psym = 8, symsize = ssiz*3.0, color = fsc_color("blue")
;oplot, sfh_time2(wclean2), MgFe_mf_disk_hist2(wclean2), linestyle = 0, thick = 6., color = fsc_color("red")
oplot, sfh_time2, MgFe_mf_disk_hist2, psym = 8, symsize = ssiz*3.0, color = fsc_color("red")

;;-------------------------------
;;PLOT [Fe/H]_DiskMass vs. Lookback time:
;;-------------------------------
;wclean = WHERE(FeH_mf_disk_hist ne 0.0)
;wdirty = WHERE(FeH_mf_disk_hist eq 0.0)
;
;plot, findgen(10), /nodata, xtitle = 'log(Lookback Time) [yrs]', ytitle = '[Fe/H]!Ddisk!N', $
;            charsize = 1.0, xrange = [MAX(sfh_time), MIN(sfh_time)], yrange = [MIN(FeH_mf_disk_hist(wclean))-0.001, MAX(FeH_mf_disk_hist(wclean))+0.001], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean)), MAX(OFe_mf_disk_hist)]
;
;oplot, sfh_time, FeH_mf_disk_hist, psym = 8, symsize = ssiz*2.0, color = fsc_color("black")
;oplot, sfh_time(wdirty), FeH_mf_disk_hist(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;legend, ['Galaxy 0'], box = 0, charsize = 1.5, linestyle = [-99], /top, /right

;-------------------------------
;PLOT SFH_DiskMass and SFH_BulgeMass vs. Lookback time:
;-------------------------------
wclean1 = WHERE(SFH_disk1 ne -99.0)
wdirty1 = WHERE(SFH_disk1 eq -99.0)
wclean2 = WHERE(SFH_disk2 ne -99.0)
wdirty2 = WHERE(SFH_disk2 eq -99.0)

plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = 'log(SFR) [M!D!9'+string(110B)+'!3!N/yr]', $ ;'log(Stars formed!Ddisk!N) [10!E10!N Msun]'
            charsize = csiz, xrange = [MAX(sfh_time1), MIN(sfh_time1)], yrange = [MIN(SFH_disk1(wclean1))-0.1, MAX(SFH_disk2(wclean2))+0.1], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean)), MAX(OFe_mf_disk_hist)]

oplot, sfh_time1, SFH_disk1, psym = 8, symsize = ssiz*3.0, color = fsc_color("blue")
oplot, sfh_time1(wclean1), SFH_disk1(wclean1), linestyle = 0, thick = 6., color = fsc_color("blue")
oplot, sfh_time2, SFH_disk2, psym = 8, symsize = ssiz*3.0, color = fsc_color("red")
;oplot, sfh_time2(wclean2), SFH_disk2(wclean2), linestyle = 0, thick = 6., color = fsc_color("red")
;print, SFH_disk2(wclean2)
;legend, ['Galaxy 0'], box = 0, charsize = 1.5, linestyle = [-99], /top, /right

;-------------------------------
;PLOT Mode_ejecta_DiskMass vs. Lookback time:
;-------------------------------
wclean1 = WHERE(SNII_disk_hist1 ne -99.0)
wdirty1 = WHERE(SNII_disk_hist1 eq -99.0)

plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = 'Metal Fraction!Ddisk!N', $
            charsize = csiz, xrange = [MAX(sfh_time1), MIN(sfh_time1)], yrange = [0.0, 0.1], xstyle = 1, ystyle = 1 ;0.001 ;MAX(SNII_disk_hist(wclean))+0.01

oplot, sfh_time1(wclean1), SNII_disk_hist1(wclean1)-0.9, psym = 8, symsize = ssiz*3.0, color = fsc_color("blue")
;oplot, sfh_time1(wdirty1), SNII_disk_hist1(wdirty1), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, thick = 4.0
oplot, sfh_time1(wclean1), SNIa_disk_hist1(wclean1), psym = 5, symsize = ssiz*2.5, thick = 4.0, color = fsc_color("blue")
;oplot, sfh_time1(wdirty1), SNIa_disk_hist1(wdirty1), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;usersym, cos(theta * !dtor)*0.6, sin(theta * !dtor)*0.6, /fill
oplot, sfh_time1(wclean1), AGB_disk_hist1(wclean1), psym = 6, symsize = ssiz*2.5, thick = 4.0, color = fsc_color("blue")
;oplot, sfh_time1(wdirty1), AGB_disk_hist1(wdirty1), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")

;oplot, sfh_time2(wclean2), SNII_disk_hist2(wclean2)-0.9, psym = 8, symsize = ssiz*3.0, color = fsc_color("red")
;oplot, sfh_time2(wclean2), SNIa_disk_hist2(wclean2), psym = 5, symsize = ssiz*2.5, thick = 4.0, color = fsc_color("red")
;oplot, sfh_time2(wclean2), AGB_disk_hist2(wclean2), psym = 6, symsize = ssiz*2.5, thick = 4.0, color = fsc_color("red")
;legend, ['Galaxy 0'], box = 0, charsize = 1.5, linestyle = [-99], /top, /left
legend, ['(SNe-II) - 0.9', 'SNe-Ia', 'AGB'], box = 0, charsize = 1.75, psym = [8,5,6], symsize = [ssiz*3.,ssiz*2.5,ssiz*2.5], color = [fsc_color("black"),fsc_color("black"),fsc_color("black")], /top, /left
;
;;------------------------------------------------------------------
;;------------------------------------------------------------------
;;DISK AND BULGE HISTORIES (ELLIPTICALS: Galaxy 200):
;;------------------------------------------------------------------
;;------------------------------------------------------------------
;
;;sfh_time = alog10(G(w1(200)).sfh_time[*]) 
;sfh_time = G(w1(200)).sfh_time[*]/1.0e9
;whist = WHERE(G(w1(200)).sfh_DiskMass eq 0.0)
;bhist = WHERE(G(w1(200)).sfh_BulgeMass eq 0.0)
;;-------
;OFe_mf_disk_hist = alog10(G(w1(200)).sfh_ElementsDiskMass[2,*]/G(w1(200)).sfh_ElementsDiskMass[4,*]) - alog10(O_mf/Fe_mf)
;OFe_mf_disk_hist(whist) = 0.0
;
;MgFe_mf_disk_hist = alog10(G(w1(200)).sfh_ElementsDiskMass[3,*]/G(w1(200)).sfh_ElementsDiskMass[4,*]) - alog10(Mg_mf/Fe_mf)
;MgFe_mf_disk_hist(whist) = 0.0
;
;FeH_mf_disk_hist = alog10(G(w1(200)).sfh_ElementsDiskMass[4,*]/G(w1(200)).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf)
;FeH_mf_disk_hist(whist) = 0.0
;
;SNII_disk_hist = G(w1(200)).sfh_MetalsDiskMass[1,*]/G(w1(200)).sfh_DiskMass[*]
;SNII_disk_hist(whist) = -99.0
;
;SNIa_disk_hist = G(w1(200)).sfh_MetalsDiskMass[0,*]/G(w1(200)).sfh_DiskMass[*]
;SNIa_disk_hist(whist) = -99.0
;
;AGB_disk_hist = G(w1(200)).sfh_MetalsDiskMass[2,*]/G(w1(200)).sfh_DiskMass[*]
;AGB_disk_hist(whist) = -99.0
;
;SFH_disk = alog10(G(w1(200)).sfh_DiskMass[*])
;SFH_disk(whist) = 0.0
;
;;--------
;OFe_mf_bulge_hist = alog10(G(w1(200)).sfh_ElementsBulgeMass[2,*]/G(w1(200)).sfh_ElementsBulgeMass[4,*]) - alog10(O_mf/Fe_mf)
;OFe_mf_bulge_hist(bhist) = 0.0
;
;MgFe_mf_bulge_hist = alog10(G(w1(200)).sfh_ElementsBulgeMass[3,*]/G(w1(200)).sfh_ElementsBulgeMass[4,*]) - alog10(Mg_mf/Fe_mf)
;MgFe_mf_bulge_hist(bhist) = 0.0
;
;FeH_mf_bulge_hist = alog10(G(w1(200)).sfh_ElementsBulgeMass[4,*]/G(w1(200)).sfh_ElementsBulgeMass[0,*]) - alog10(Fe_mf/H_mf)
;FeH_mf_bulge_hist(bhist) = 0.0
;
;SNII_bulge_hist = G(w1(200)).sfh_MetalsBulgeMass[1,*]/G(w1(200)).sfh_BulgeMass[*]
;SNII_bulge_hist(bhist) = -99.0
;
;SNIa_bulge_hist = G(w1(200)).sfh_MetalsBulgeMass[0,*]/G(w1(200)).sfh_BulgeMass[*]
;SNIa_bulge_hist(bhist) = -99.0
;
;AGB_bulge_hist = G(w1(200)).sfh_MetalsBulgeMass[2,*]/G(w1(200)).sfh_BulgeMass[*]
;AGB_bulge_hist(bhist) = -99.0
;
;SFH_bulge = alog10(G(w1(200)).sfh_BulgeMass[*])
;SFH_bulge(bhist) = 0.0
;;-------------------------------
;;PLOT [O/Fe]_DiskMass vs. Lookback time:
;;-------------------------------
;wclean = WHERE(OFe_mf_disk_hist ne 0.0)
;wdirty = WHERE(OFe_mf_disk_hist eq 0.0)
;bclean = WHERE(OFe_mf_bulge_hist ne 0.0)
;bdirty = WHERE(OFe_mf_bulge_hist eq 0.0)
;
;;plot, findgen(10), /nodata, xtitle = 'log(Lookback Time) [yrs]', ytitle = '[O/Fe]!Ddisk!N', $
;;            charsize = 1.0, xrange = [MAX(sfh_time), MIN(sfh_time)], yrange = [MIN(OFe_mf_disk_hist(wclean))-0.001, MAX(OFe_mf_disk_hist(wclean))+0.001], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean)), MAX(OFe_mf_disk_hist)]
;plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = '[O/Fe]!Ddisk!N', $ ;xtitle = 'log(Lookback Time) [yrs]'
;            charsize = 1.0, xrange = [MAX(sfh_time), MIN(sfh_time)], yrange = [MIN(OFe_mf_disk_hist(wclean))-0.001, MAX(OFe_mf_disk_hist(wclean))+0.001], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean))-0.001, MAX(OFe_mf_disk_hist(wclean))+0.001]
;
;
;oplot, sfh_time, OFe_mf_disk_hist, psym = 8, symsize = ssiz*2.0, color = fsc_color("black")
;oplot, sfh_time(wdirty), OFe_mf_disk_hist(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;;oplot, sfh_time, OFe_mf_bulge_hist, psym = 6, symsize = ssiz*2.0, color = fsc_color("black")
;;oplot, sfh_time(bdirty), OFe_mf_bulge_hist(bdirty), psym = 6, symsize = ssiz*2.0, color = fsc_color("white")
;
;legend, ['Galaxy 0'], box = 0, charsize = 1.5, linestyle = [-99], /top, /right
;
;;-------------------------------
;;PLOT [Mg/Fe]_DiskMass vs. Lookback time:
;;-------------------------------
;wclean = WHERE(MgFe_mf_disk_hist ne 0.0)
;wdirty = WHERE(MgFe_mf_disk_hist eq 0.0)
;
;plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = '[Mg/Fe]!Ddisk!N', $
;            charsize = 1.0, xrange = [MAX(sfh_time), MIN(sfh_time)], yrange = [MIN(MgFe_mf_disk_hist(wclean))-0.001, MAX(MgFe_mf_disk_hist(wclean))+0.001], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean)), MAX(OFe_mf_disk_hist)]
;
;oplot, sfh_time, MgFe_mf_disk_hist, psym = 8, symsize = ssiz*2.0, color = fsc_color("black")
;oplot, sfh_time(wdirty), MgFe_mf_disk_hist(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;legend, ['Galaxy 0'], box = 0, charsize = 1.5, linestyle = [-99], /top, /right
;
;;;-------------------------------
;;;PLOT [Fe/H]_DiskMass vs. Lookback time:
;;;-------------------------------
;;wclean = WHERE(FeH_mf_disk_hist ne 0.0)
;;wdirty = WHERE(FeH_mf_disk_hist eq 0.0)
;;
;;plot, findgen(10), /nodata, xtitle = 'log(Lookback Time) [yrs]', ytitle = '[Fe/H]!Ddisk!N', $
;;            charsize = 1.0, xrange = [MAX(sfh_time), MIN(sfh_time)], yrange = [MIN(FeH_mf_disk_hist(wclean))-0.001, MAX(FeH_mf_disk_hist(wclean))+0.001], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean)), MAX(OFe_mf_disk_hist)]
;;
;;oplot, sfh_time, FeH_mf_disk_hist, psym = 8, symsize = ssiz*2.0, color = fsc_color("black")
;;oplot, sfh_time(wdirty), FeH_mf_disk_hist(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;;legend, ['Galaxy 0'], box = 0, charsize = 1.5, linestyle = [-99], /top, /right
;
;;-------------------------------
;;PLOT SFH_DiskMass and SFH_BulgeMass vs. Lookback time:
;;-------------------------------
;wclean = WHERE(SFH_disk ne 0.0)
;wdirty = WHERE(SFH_disk eq 0.0)
;bclean = WHERE(SFH_bulge ne 0.0)
;bdirty = WHERE(SFH_bulge eq 0.0)
;
;plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = 'log(Stars formed!Ddisk!N) [10!E10!N Msun]', $
;            charsize = 1.0, xrange = [MAX(sfh_time), MIN(sfh_time)], yrange = [MIN(SFH_disk(wclean))-0.1, MAX(SFH_disk(wclean))+0.1], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean)), MAX(OFe_mf_disk_hist)]
;
;oplot, sfh_time, SFH_disk, psym = 8, symsize = ssiz*2.0, color = fsc_color("black")
;oplot, sfh_time(wdirty), SFH_disk(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;;oplot, sfh_time, SFH_bulge, psym = 6, symsize = ssiz, color = fsc_color("black")
;;oplot, sfh_time(bdirty), SFH_bulge(bdirty), psym = 6, symsize = ssiz, color = fsc_color("white")
;
;legend, ['Galaxy 0'], box = 0, charsize = 1.5, linestyle = [-99], /top, /right
;
;;-------------------------------
;;PLOT Mode_ejecta_DiskMass vs. Lookback time:
;;-------------------------------
;wclean = WHERE(SNII_disk_hist ne -99.0)
;wdirty = WHERE(SNII_disk_hist eq -99.0)
;
;plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = 'Metal Fraction!Ddisk!N', $
;            charsize = 1.0, xrange = [MAX(sfh_time), MIN(sfh_time)], yrange = [MIN(AGB_disk_hist(wclean))-0.0001, MAX(SNII_disk_hist(wclean))+0.0001], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean)), MAX(OFe_mf_disk_hist)]
;
;oplot, sfh_time, SNII_disk_hist, psym = 8, symsize = ssiz*2.0, color = fsc_color("black")
;oplot, sfh_time(wdirty), SNII_disk_hist(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;oplot, sfh_time, SNIa_disk_hist, psym = 8, symsize = ssiz*2.0, color = fsc_color("red")
;oplot, sfh_time(wdirty), SNIa_disk_hist(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;oplot, sfh_time, AGB_disk_hist, psym = 8, symsize = ssiz*2.0, color = fsc_color("blue")
;oplot, sfh_time(wdirty), AGB_disk_hist(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;legend, ['Galaxy 0'], box = 0, charsize = 1.5, linestyle = [-99], /top, /right
;legend, ['SNe-II', 'SNe-Ia', 'AGB'], box = 0, charsize = 1., psym = [8,8,8], symsize = [ssiz*2.5,ssiz*2.5,ssiz*2.5], color = [fsc_color("black"),fsc_color("red"),fsc_color("blue")], /bottom, /right
;
;;print, SFH_disk
;;------------------------------------------------------------------
;;------------------------------------------------------------------
;;DISK HISTORIES (ELLIPTICALS: Galaxy 242):
;;------------------------------------------------------------------
;;------------------------------------------------------------------
;
;;sfh_time = alog10(G(2434).sfh_time[*])
;sfh_time = G(2434).sfh_time[*]/1.0e9
;whist = WHERE(G(2434).sfh_DiskMass eq 0.0)
;bhist = WHERE(G(2434).sfh_BulgeMass eq 0.0)
;;--------
;OFe_mf_disk_hist = alog10(G(2434).sfh_ElementsDiskMass[2,*]/G(2434).sfh_ElementsDiskMass[4,*]) - alog10(O_mf/Fe_mf)
;OFe_mf_disk_hist(whist) = 0.0
;
;MgFe_mf_disk_hist = alog10(G(2434).sfh_ElementsDiskMass[3,*]/G(2434).sfh_ElementsDiskMass[4,*]) - alog10(Mg_mf/Fe_mf)
;MgFe_mf_disk_hist(whist) = 0.0
;
;FeH_mf_disk_hist = alog10(G(2434).sfh_ElementsDiskMass[4,*]/G(2434).sfh_ElementsDiskMass[0,*]) - alog10(Fe_mf/H_mf)
;FeH_mf_disk_hist(whist) = 0.0
;
;SNII_disk_hist = G(2434).sfh_MetalsDiskMass[1,*]/G(2434).sfh_DiskMass[*]
;SNII_disk_hist(whist) = -99.0
;
;SNIa_disk_hist = G(2434).sfh_MetalsDiskMass[0,*]/G(2434).sfh_DiskMass[*]
;SNIa_disk_hist(whist) = -99.0
;
;AGB_disk_hist = G(2434).sfh_MetalsDiskMass[2,*]/G(2434).sfh_DiskMass[*]
;AGB_disk_hist(whist) = -99.0
;
;SFH_disk = alog10(G(2434).sfh_DiskMass[*])
;SFH_disk(whist) = 0.0
;
;;a = [4.,5.,6.,7.]
;;b = [6.,7.,8.,9.]
;;print, a[*]/b[*]
;;print, a/b
;;print, G(0).sfh_MetalsDiskMass[1,*], G(0).sfh_DiskMass[*]
;;print, G(0).sfh_MetalsDiskMass[1,1]/G(0).sfh_DiskMass[1]
;
;;--------
;OFe_mf_bulge_hist = alog10(G(2434).sfh_ElementsBulgeMass[2,*]/G(2434).sfh_ElementsBulgeMass[4,*]) - alog10(O_mf/Fe_mf)
;OFe_mf_bulge_hist(bhist) = 0.0
;
;MgFe_mf_bulge_hist = alog10(G(2434).sfh_ElementsBulgeMass[3,*]/G(2434).sfh_ElementsBulgeMass[4,*]) - alog10(Mg_mf/Fe_mf)
;MgFe_mf_bulge_hist(bhist) = 0.0
;
;FeH_mf_bulge_hist = alog10(G(2434).sfh_ElementsBulgeMass[4,*]/G(2434).sfh_ElementsBulgeMass[0,*]) - alog10(Fe_mf/H_mf)
;FeH_mf_bulge_hist(bhist) = -99.0
;
;SNII_bulge_hist = G(2434).sfh_MetalsBulgeMass[1,*]/G(2434).sfh_BulgeMass[*]
;SNII_bulge_hist(bhist) = -99.0
;
;SNIa_bulge_hist =G(2434).sfh_MetalsBulgeMass[0,*]/G(2434).sfh_BulgeMass[*]
;SNIa_bulge_hist(bhist) = -99.0
;
;AGB_bulge_hist = G(2434).sfh_MetalsBulgeMass[2,*]/G(2434).sfh_BulgeMass[*]
;AGB_bulge_hist(bhist) = -99.0
;
;SFH_bulge = alog10(G(2434).sfh_BulgeMass[*])
;SFH_bulge(bhist) = 0.0
;;-------------------------------
;;PLOT [O/Fe]_DiskMass vs. Lookback time:
;;-------------------------------
;wclean = WHERE(OFe_mf_disk_hist ne 0.0)
;wdirty = WHERE(OFe_mf_disk_hist eq 0.0)
;bclean = WHERE(OFe_mf_bulge_hist ne 0.0)
;bdirty = WHERE(OFe_mf_bulge_hist eq 0.0)
;
;plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = '[O/Fe]!Ddisk!N', $ ;xtitle = 'log(Lookback Time) [yrs]'
;            charsize = 1.0, xrange = [MAX(sfh_time), MIN(sfh_time)], yrange = [MIN(OFe_mf_disk_hist(wclean))-0.001, MAX(OFe_mf_disk_hist(wclean))+0.001], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean))-0.001, MAX(OFe_mf_disk_hist(wclean))+0.001]
;
;oplot, sfh_time, OFe_mf_disk_hist, psym = 8, symsize = ssiz*2.0, color = fsc_color("black")
;oplot, sfh_time(wdirty), OFe_mf_disk_hist(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;oplot, sfh_time, OFe_mf_bulge_hist, psym = 6, symsize = ssiz, color = fsc_color("black")
;oplot, sfh_time(bdirty), OFe_mf_bulge_hist(bdirty), psym = 6, symsize = ssiz, color = fsc_color("white")
;
;legend, ['Galaxy 2434'], box = 0, charsize = 1.5, linestyle = [-99], /top, /right
;
;;-------------------------------
;;PLOT [Mg/Fe]_DiskMass vs. Lookback time:
;;-------------------------------
;wclean = WHERE(MgFe_mf_disk_hist ne 0.0)
;wdirty = WHERE(MgFe_mf_disk_hist eq 0.0)
;bclean = WHERE(MgFe_mf_bulge_hist ne 0.0)
;bdirty = WHERE(MgFe_mf_bulge_hist eq 0.0)
;
;plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = '[Mg/Fe]!Ddisk!N', $
;            charsize = 1.0, xrange = [MAX(sfh_time), MIN(sfh_time)], yrange = [MIN(MgFe_mf_disk_hist(wclean))-0.001, MAX(MgFe_mf_disk_hist(wclean))+0.001], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean)), MAX(OFe_mf_disk_hist)]
;
;oplot, sfh_time, MgFe_mf_disk_hist, psym = 8, symsize = ssiz*2.0, color = fsc_color("black")
;oplot, sfh_time(wdirty), MgFe_mf_disk_hist(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;oplot, sfh_time, MgFe_mf_bulge_hist, psym = 6, symsize = ssiz, color = fsc_color("black")
;oplot, sfh_time(bdirty), MgFe_mf_bulge_hist(bdirty), psym = 6, symsize = ssiz, color = fsc_color("white")
;
;legend, ['Galaxy 2434'], box = 0, charsize = 1.5, linestyle = [-99], /top, /right
;
;;;-------------------------------
;;;PLOT [Fe/H]_DiskMass vs. Lookback time:
;;;-------------------------------
;;wclean = WHERE(FeH_mf_disk_hist ne 0.0)
;;wdirty = WHERE(FeH_mf_disk_hist eq 0.0)
;;bclean = WHERE(FeH_mf_bulge_hist ne 0.0)
;;bdirty = WHERE(FeH_mf_bulge_hist eq 0.0)
;;
;;plot, findgen(10), /nodata, xtitle = 'log(Lookback Time) [yrs]', ytitle = '[Fe/H]!Ddisk!N', $
;;            charsize = 1.0, xrange = [MAX(sfh_time), MIN(sfh_time)], yrange = [MIN(FeH_mf_disk_hist(wclean))-0.001, MAX(FeH_mf_disk_hist(wclean))+0.001], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean)), MAX(OFe_mf_disk_hist)]
;;
;;oplot, sfh_time, FeH_mf_disk_hist, psym = 8, symsize = ssiz*2.0, color = fsc_color("black")
;;oplot, sfh_time(wdirty), FeH_mf_disk_hist(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;;oplot, sfh_time, FeH_mf_bulge_hist, psym = 6, symsize = ssiz, color = fsc_color("black")
;;;oplot, sfh_time(bdirty), FeH_mf_bulge_hist(bdirty), psym = 6, symsize = ssiz, color = fsc_color("white")
;;
;;legend, ['Galaxy 2434'], box = 0, charsize = 1.5, linestyle = [-99], /top, /right
;
;;-------------------------------
;;PLOT SFH_DiskMass and SFH_BulgeMass vs. Lookback time:
;;-------------------------------
;wclean = WHERE(FeH_mf_disk_hist ne 0.0)
;wdirty = WHERE(FeH_mf_disk_hist eq 0.0)
;bclean = WHERE(SFH_disk ne 0.0)
;bdirty = WHERE(SFH_disk eq 0.0)
;
;plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = 'log(Stars formed!Ddisk!N) [10!E10!N Msun]', $
;            charsize = 1.0, xrange = [MAX(sfh_time), MIN(sfh_time)], yrange = [MIN(SFH_disk(wclean))-0.1, MAX(SFH_disk(wclean))+0.1], xstyle = 1, ystyle = 1 ;yrange = [MIN(OFe_mf_disk_hist(wclean)), MAX(OFe_mf_disk_hist)]
;
;oplot, sfh_time, SFH_disk, psym = 8, symsize = ssiz*2.0, color = fsc_color("black")
;oplot, sfh_time(wdirty), SFH_disk(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;oplot, sfh_time, SFH_bulge, psym = 6, symsize = ssiz, color = fsc_color("black")
;oplot, sfh_time(bdirty), SFH_bulge(bdirty), psym = 6, symsize = ssiz, color = fsc_color("white")
;
;legend, ['Galaxy 2434'], box = 0, charsize = 1.5, linestyle = [-99], /top, /right
;
;;-------------------------------
;;PLOT Mode_ejecta_DiskMass vs. Lookback time:
;;-------------------------------
;wclean = WHERE(SNII_disk_hist ne -99.0)
;wdirty = WHERE(SNII_disk_hist eq -99.0)
;bclean = WHERE(SNII_bulge_hist ne -99.0)
;bdirty = WHERE(SNII_bulge_hist eq -99.0)
;
;plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = 'Metal Fraction!Ddisk!N', $
;            charsize = 1.0, xrange = [MAX(sfh_time), MIN(sfh_time)], yrange = [MIN(AGB_disk_hist(wclean))-0.001, MAX(SNII_disk_hist(wclean))+0.001], xstyle = 1, ystyle = 1 ;MAX(SNII_disk_hist(wclean))+0.001  ;, /ylog ;yrange = [0.0, 0.001]
;
;oplot, sfh_time, SNII_disk_hist, psym = 8, symsize = ssiz*2.0, color = fsc_color("black")
;oplot, sfh_time(wdirty), SNII_disk_hist(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;oplot, sfh_time, SNIa_disk_hist, psym = 8, symsize = ssiz*2.0, color = fsc_color("red")
;oplot, sfh_time(wdirty), SNIa_disk_hist(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;oplot, sfh_time, AGB_disk_hist, psym = 8, symsize = ssiz*2.0, color = fsc_color("blue")
;oplot, sfh_time(wdirty), AGB_disk_hist(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;
;;oplot, sfh_time, SNII_bulge_hist, psym = 6, symsize = ssiz, color = fsc_color("black")
;;oplot, sfh_time(bdirty), SNII_bulge_hist(bdirty), psym = 6, symsize = ssiz, color = fsc_color("white")
;;oplot, sfh_time, SNIa_bulge_hist, psym = 6, symsize = ssiz, color = fsc_color("red")
;;oplot, sfh_time(bdirty), SNIa_bulge_hist(bdirty), psym = 6, symsize = ssiz, color = fsc_color("white")
;;oplot, sfh_time, AGB_bulge_hist, psym = 6, symsize = ssiz, color = fsc_color("blue")
;;oplot, sfh_time(bdirty), AGB_bulge_hist(bdirty), psym = 6, symsize = ssiz, color = fsc_color("white")
;
;legend, ['Galaxy 2434'], box = 0, charsize = 1.5, linestyle = [-99], /top, /right
;legend, ['SNe-II', 'SNe-Ia', 'AGB'], box = 0, charsize = 1., psym = [8,8,8], symsize = [ssiz*2.5,ssiz*2.5,ssiz*2.5], color = [fsc_color("black"),fsc_color("red"),fsc_color("blue")], /bottom, /right
;
;;print, SNII_disk_hist
;;print, "----"
;;print, G[0].sfh_time
;;print, "----"
;
;;-------------------------------
;;PLOT Fractions vs. Lookback time:
;;-------------------------------
;
;;SNIIfrac = G(2434).sfh_MetalsDiskMass[1,*]/G(2434).sfh_DiskMass[*]
;SNIIfrac_disk = G(2434).sfh_MetalsDiskMass[1,*]/(G(2434).sfh_MetalsDiskMass[0,*]+G(2434).sfh_MetalsDiskMass[1,*]+G(2434).sfh_MetalsDiskMass[2,*])
;SNIafrac_disk = G(2434).sfh_MetalsDiskMass[0,*]/(G(2434).sfh_MetalsDiskMass[0,*]+G(2434).sfh_MetalsDiskMass[1,*]+G(2434).sfh_MetalsDiskMass[2,*])
;AGBfrac_disk = G(2434).sfh_MetalsDiskMass[2,*]/(G(2434).sfh_MetalsDiskMass[0,*]+G(2434).sfh_MetalsDiskMass[1,*]+G(2434).sfh_MetalsDiskMass[2,*])
;
;whist = WHERE((G(2434).sfh_MetalsDiskMass[0,*]+G(2434).sfh_MetalsDiskMass[1,*]+G(2434).sfh_MetalsDiskMass[2,*]) eq 0.0)
;SNIIfrac_disk(whist) = -99.0
;SNIafrac_disk(whist) = -99.0
;AGBfrac_disk(whist) = -99.0
;
;wclean = WHERE(SNIIfrac_disk ne -99.0)
;wdirty = WHERE(SNIIfrac_disk eq -99.0)
;;print, SNIafrac_disk, SNIIfrac_disk, AGBfrac_disk
;plot, findgen(10), /nodata, xtitle = 'Lookback Time [Gyrs]', ytitle = 'M!DZ,mode!N/M!DZ!N', $
;            charsize = 1.0, xrange = [MAX(sfh_time), MIN(sfh_time)], yrange = [0.0, 0.05], xstyle = 1, ystyle = 1 ;, MAX(SNIIfrac_disk(wclean))+0.1
;
;oplot, sfh_time, SNIIfrac_disk, psym = 8, symsize = ssiz*2.0, color = fsc_color("black")
;oplot, sfh_time(wdirty), SNIIfrac_disk(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;oplot, sfh_time, SNIafrac_disk, psym = 8, symsize = ssiz*2.0, color = fsc_color("red")
;oplot, sfh_time(wdirty), SNIafrac_disk(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;oplot, sfh_time, AGBfrac_disk, psym = 8, symsize = ssiz*2.0, color = fsc_color("blue")
;oplot, sfh_time(wdirty), AGBfrac_disk(wdirty), psym = 8, symsize = ssiz*2.5, color = fsc_color("white")
;
;
;;print, "----"
;;print, G(2434).sfh_ElementsDiskMass[0,*]
;;print, "----"
;;print, G(2434).sfh_ElementsDiskMass[2,*]
;;print, "----"
;;print, G(2434).sfh_ElementsDiskMass[4,*]
;;print, G(2434).sfh_time
ENDIF ;IF (MAINELEMENTS eq 0) THEN BEGIN

ENDIF ;  IF (OLDGALSTRUCT eq 0) THEN BEGIN 
;-------------------------------

device, /close
print, "DONE!"  
end

pro hist2d,x,y,hist,xrange,yrange,nxbins,nybins
;+
; NAME:
;    HIST2D
;
; PURPOSE:
;    like the IDL built in function hist_2d
;    but better since it can accept floats
;    output histogram is a longarray(nxbins,nybins)
;    uses only the relevent data in range
;
;  Dave Johnston
;-
if n_params() eq 0 then begin
        print,'-syntax hist2d,x,y,hist,xrange,yrange,nxbins,nybins'
        return
endif

xrange=float(xrange)
yrange=float(yrange)

xmin=xrange(0)
xmax=xrange(1)
ymin=yrange(0)
ymax=yrange(1)

w=where(x gt xmin and x lt xmax and y gt ymin and y lt ymax,count)
if count eq 0 then begin
        hist=-1
        return
endif

xind=floor((x(w)-xmin)*(nxbins/(xmax-xmin)))
yind=floor((y(w)-ymin)*(nybins/(ymax-ymin)))

ind=xind+nxbins*yind
h=histogram(ind,min=0l,max=long(nxbins*nybins)-1)

hist=lonarr(nxbins,nybins)
hist(*)=h

return
end