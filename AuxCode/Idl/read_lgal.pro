pro read_lgal,FileRoot, Galstruct, GalStructArray $
             ,FileNr=Filenr, TreeNr=TreeNr, NGalsPerTree=NGalsPerTree $
             ,FirstFile=FirstFile,LastFile=LastFile $
             ,swap_endian=swap_endian

; Procedure to read in an L-galaxies output file and return as an
; array of galaxy structures, GalStructArray

; Input:
;    FileRoot - root file name (will append .* where * runs over all
;               numbers  between FirstFile and LastFile
;    Galstruct - template for galaxy structure
;    FirstFile - First file number
;    LastFile - Last file number
;    swap_endian - flag to swap the endianness of the data upon read
; Output:
;    GalStructArray - an array of Galaxy structures containing the
;                     galaxy properties
;    FileNr - the file corresponding to each galaxy
;    TreeNr - the tree corresponding to each galaxy
;    NGalsPerTree - the number of galaxies in each tree

compile_opt idl2

; Ideally change the following to detect the file numbers present in
; the output directory
if n_elements(FirstFile) eq 0 then FirstFile=0
if n_elements(LastFile) eq 0 then LastFile=511
NFiles = LastFile - FirstFile + 1

; First loop to count trees and galaxies
TotNTrees = 0L & TotNGals = 0L
for fnr = FirstFile, LastFile do begin
   fname = strcompress(FileRoot+'_'+string(fnr), /remove_all)
   if keyword_set(swap_endian) then $
      openr, 1, fname, /swap_endian $
   else $
      openr, 1, fname
   NTrees = 0L     & readu, 1, NTrees
   NGals = 0L   & readu, 1, NGals
   TotNTrees = TotNTrees + NTrees
   TotNGals =  TotNGals  + NGals
   close, 1
endfor
print, " Total number of trees = ", TotNTrees
print, " Total number of galaxies = ", TotNGals

; Define output arrays
NGalsPerTree = lonarr(TotNTrees)
FileNr = lonarr(TotNgals)  
TreeNr = lonarr(TotNgals)  
GalStructArray = replicate(GalStruct, TotNGals)

; Second loop to populate structure array
offset = 0L
off = 0L
for fnr = FirstFile, LastFile do begin
   print, ' file:', fnr
   fname = strcompress(FileRoot+'_'+string(fnr), /remove_all)
   if keyword_set(swap_endian) then $
      openr, 1, fname, /swap_endian $
   else $
      openr, 1, fname
   NTrees = 0L     & readu, 1, NTrees   & print, ' Ntrees:', Ntrees
   NGals = 0L   & readu, 1, NGals & print, ' NGals:', NGals
   ; cannot read into array sections, so need these temporary variables
   temp = lonarr(NTrees)
   readu, 1, temp
   NGalsPerTree[off:off+NTrees-1] = temp
   if Ngals gt 0 then begin
      temp = replicate(GalStruct, Ngals)
      readu, 1, temp
      GalStructArray[offset:offset+NGals-1] = temp
      skip=0L
      for i=0L, Ntrees-1 do begin
         if NGalsPerTree[off+i] gt 0 then begin
            FileNr[offset+skip:offset+skip+NGalsPerTree[off+i]-1] = fnr
            TreeNr[offset+skip:offset+skip+NGalsPerTree[off+i]-1] = i
            skip = skip + NGalsPerTree[off+i]
         endif
      endfor
   endif
   offset = offset + NGals
   off = off + Ntrees
   close, 1
endfor

end
