pro squish,dataset

; Because the full data-set is so large, it can take a long time to
; read in on many machines and can over-run memory.  This procedure
; reads in the data a bit at a time and then write out desired
; quantities into individual arrays which can later be read in by
; read_squish
; Currently this procedure must be edited by hand to specify the
; source files, the destination directory, and the quantities that
; are sought.  As it will only be run occasionally, I have not sought
; to make it any cleverer.

; Input:
;   dataset - 'HT09' or 'DLB07', defines which dataset to use

; Define which files you want to read in
if dataset eq 'HT09' then begin
    format='HT09'
    DirNameRead = '/export/virgo/SAM/snaps/HT09/Z000_63/'
    DirNameWrite = '../output/HT09/'
endif else begin
    format='DeLucia'
    DirNameRead = '/export/virgo/SAM/snaps/DLB07/Z000_63/'
    DirNameWrite = '../output/DLB07/'
endelse
FileName = 'SA_z0.00'
FirstFile = 0
LastFile = 511
BlockSize=64
 
;------------------------------------------------------------------------------
; Create output data files and set the galaxy number equal to zero
; Add in more of these, as desired (no need to repeat those already done)
;Type
fname = strcompress(DirNameWrite+FileName+'_Type', /remove_all)
if keyword_set(swap_endian) then openw, 1, fname, /swap_endian else openw, 1, fname
writeu,1,0L
close, 1
;Mvir
fname = strcompress(DirNameWrite+FileName+'_Mvir', /remove_all)
if keyword_set(swap_endian) then openw, 1, fname, /swap_endian else openw, 1, fname
writeu,1,0L
close, 1
;StellarMass
fname = strcompress(DirNameWrite+FileName+'_StellarMass', /remove_all)
if keyword_set(swap_endian) then openw, 1, fname, /swap_endian else openw, 1, fname
writeu,1,0L
close, 1
;Sfr
fname = strcompress(DirNameWrite+FileName+'_Sfr', /remove_all)
if keyword_set(swap_endian) then openw, 1, fname, /swap_endian else openw, 1, fname
writeu,1,0L
close, 1
;MetalsStellarMass
fname = strcompress(DirNameWrite+FileName+'_MetalsStellarMass', /remove_all)
if keyword_set(swap_endian) then openw, 1, fname, /swap_endian else openw, 1, fname
writeu,1,0L
close, 1
;MassweightAge
fname = strcompress(DirNameWrite+FileName+'_MassWeightAge', /remove_all)
if keyword_set(swap_endian) then openw, 1, fname, /swap_endian else openw, 1, fname
writeu,1,0L
close, 1
if dataset eq 'HT09' then begin 
;ICL
fname = strcompress(DirNameWrite+FileName+'_ICLmass', /remove_all)
if keyword_set(swap_endian) then openw, 1, fname, /swap_endian else openw, 1, fname
writeu,1,0L
close, 1
endif

;-------------------------------------------------------------------------------
;Loop over blocks of data files
print,'Processing block...'
for FirstFileBlock=FirstFile,LastFile,BlockSize do begin
    LastFileBlock=FirstFileBlock+BlockSize-1 < LastFile
    print,'   ',FirstFileBlock,' - ',LastFileBlock

;-------------------------------------------------------------------------------
; Read in data
    read_lgal,DirNameRead+FileName,GalStruct $
         ,FirstFile=FirstFileBlock,LastFile=LastFileBlock $
         ,format=format                ;,/swap_endian

;-------------------------------------------------------------------------------
; Write out simple arrays containing useful data from structure array
    NGals=n_elements(GalStruct)
;Type
    fname = strcompress(DirNameWrite+FileName+'_Type', /remove_all)
    if keyword_set(swap_endian) then openu, 1, fname, /swap_endian else openu, 1, fname
    Ngals_old=0L & readu,1,NGals_old
    point_lun,1,0 & writeu,1,NGals_old+Ngals
    point_lun,1,(fstat(1)).size & writeu,1,GalStruct.Type
    close, 1
;Mvir
    fname = strcompress(DirNameWrite+FileName+'_Mvir', /remove_all)
    if keyword_set(swap_endian) then openu, 1, fname, /swap_endian else openu, 1, fname
    Ngals_old=0L & readu,1,NGals_old
    point_lun,1,0 & writeu,1,NGals_old+Ngals
    point_lun,1,(fstat(1)).size & writeu,1,GalStruct.Mvir
    close, 1
;StellarMass
    fname = strcompress(DirNameWrite+FileName+'_StellarMass', /remove_all)
    if keyword_set(swap_endian) then openu, 1, fname, /swap_endian else openu, 1, fname
    Ngals_old=0L & readu,1,NGals_old
    point_lun,1,0 & writeu,1,NGals_old+Ngals
    point_lun,1,(fstat(1)).size & writeu,1,GalStruct.StellarMass
    close, 1
;Sfr
    fname = strcompress(DirNameWrite+FileName+'_Sfr', /remove_all)
    if keyword_set(swap_endian) then openu, 1, fname, /swap_endian else openu, 1, fname
    Ngals_old=0L & readu,1,NGals_old
    point_lun,1,0 & writeu,1,NGals_old+Ngals
    point_lun,1,(fstat(1)).size & writeu,1,GalStruct.Sfr
    close, 1
;MetalsStellarMass
    fname = strcompress(DirNameWrite+FileName+'_MetalsStellarMass', /remove_all)
    if keyword_set(swap_endian) then openu, 1, fname, /swap_endian else openu, 1, fname
    Ngals_old=0L & readu,1,NGals_old
    point_lun,1,0 & writeu,1,NGals_old+Ngals
    point_lun,1,(fstat(1)).size & writeu,1,GalStruct.MetalsStellarMass
    close, 1
;MassWeightAge
    fname = strcompress(DirNameWrite+FileName+'_MassWeightAge', /remove_all)
    if keyword_set(swap_endian) then openu, 1, fname, /swap_endian else openu, 1, fname
    Ngals_old=0L & readu,1,NGals_old
    point_lun,1,0 & writeu,1,NGals_old+Ngals
    point_lun,1,(fstat(1)).size & writeu,1,GalStruct.MassWeightAge
    close, 1
;ICL
    if dataset eq 'HT09' then begin
    fname = strcompress(DirNameWrite+FileName+'_ICLmass', /remove_all)
    if keyword_set(swap_endian) then openu, 1, fname, /swap_endian else openu, 1, fname
    Ngals_old=0L & readu,1,NGals_old
    point_lun,1,0 & writeu,1,NGals_old+Ngals
    point_lun,1,(fstat(1)).size & writeu,1,GalStruct.ICLmass
    close, 1
    endif
; Nullify galaxy structure
    GalStruct=0
    
;-------------------------------------------------------------------------------
; End loop over blocks
endfor

end
