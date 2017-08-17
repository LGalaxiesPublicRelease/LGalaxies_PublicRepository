;=========================================================================
;
;  Script to read in L-galaxies galaxytree data
;
;-------------------------------------------------------------------------

; Define which files you want to read in
DirName = '../output/'
FileName = 'SA_galtree'
ModelName = DirName + FileName
;FirstFile = 0
;LastFile = 511
FirstFile = 287
LastFile = 287
@'../output/galaxytree_template.pro'

;-------------------------------------------------------------------------------

; read in data

read_lgaltree,ModelName,Template,Gnew $
         ,FileNr=FileNr $
         ,FirstFile=FirstFile,LastFile=LastFile ;,/swap_endian

;-------------------------------------------------------------------------------
