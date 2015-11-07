pro read_lgaltree,FileRoot, Galstruct, GalStructArray $
             ,FileNr=Filenr $
             ,StructSize=StructSize $
             ,FirstFile=FirstFile,LastFile=LastFile $
             ,swap_endian=swap_endian

; Procedure to read in an L-galaxies galaxy tree output file and return as an
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

compile_opt idl2

; Ideally change the following to detect the file numbers present in
; the output directory
if n_elements(FirstFile) eq 0 then FirstFile=0
if n_elements(LastFile) eq 0 then LastFile=511
NFiles = LastFile - FirstFile + 1

; First loop to count trees and galaxies
TotNGals = 0L
StructSize = 0L
Ngals = lonarr(Nfiles)
for fnr = FirstFile, LastFile do begin
   fname = strcompress(FileRoot+'_'+string(fnr), /remove_all)
   if keyword_set(swap_endian) then $
      openr, 1, fname, /swap_endian $
   else $
      openr, 1, fname
   One = 0L     & readu, 1, One & print,'One=',One
   readu, 1, StructSize
   if fnr eq FirstFile then print,'StructSize=',StructSize,' bytes'
   temp=0L & readu, 1,temp & NGals[fnr-FirstFile]=temp
   print,'Ngals=',Ngals[fnr-FirstFile]
   close, 1
endfor
TotNGals =  nint(total(NGals,/integer))
print, " Total number of galaxies = ", TotNGals

; Define output arrays
FileNr = lonarr(TotNgals)  
GalStructArray = replicate(GalStruct, TotNGals)

; Second loop to populate structure array
offset = 0L
for fnr = FirstFile, LastFile do begin
   print, ' file:', fnr
   fname = strcompress(FileRoot+'_'+string(fnr), /remove_all)
   if keyword_set(swap_endian) then $
      openr, 1, fname, /swap_endian $
   else $
      openr, 1, fname
   ; read number of bytes equal to one galaxy structure as an offset
   temp=replicate(GalStruct,1)
   readu,1,temp
   ; cannot read into array sections, so need these temporary variables
   Ngal=Ngals[fnr-FirstFile]
   if Ngal gt 0 then begin
      temp = replicate(GalStruct, Ngal)
      readu, 1, temp
      GalStructArray[offset:offset+NGal-1] = temp
      FileNr[offset:offset+NGal-1] = fnr
   endif
   offset = offset + NGal
   close, 1
endfor

end
