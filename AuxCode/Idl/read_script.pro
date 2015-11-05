;=========================================================================
;
;  Script to read in L-galaxies snapshot data
;
;  To force a re-read of the data do "delvar, GalStruct
;
;-------------------------------------------------------------------------

; Parameters

; Template for L-Galaxies data
;@'/export/data/virgo/MillGas/WMAP7/data/500_sam/Guo10/snap_template.pro'
@'../output/bursttest/snap_template.pro'
; Define which files you want to read in
;
DirName = '../output/bursttest/'
FileName = 'SA_z0.99'
;FirstFile = 0
FirstFile = 511
LastFile = 511

;-------------------------------------------------------------------------------

ModelName = DirName + FileName

; read in data
if n_elements(GalStruct) eq 0 then read_lgal,ModelName,Template,G $
         ,FileNr=FileNr,TreeNr=TreeNr,NGalsPerTree=NGalsPerTree $
         ,FirstFile=FirstFile,LastFile=LastFile ;,/swap_endian

;-------------------------------------------------------------------------------
