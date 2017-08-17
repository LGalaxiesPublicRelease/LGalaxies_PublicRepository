pro write_lgaltree,FileRoot, GalStructArray ,StructSize, swap_endian=swap_endian

; Procedure to write an L-galaxies output file.

; Input:
;    FileRoot - root file name (will append .* where * runs over all
;               present file numbers
;    GalStructArray - an array of Galaxy structures containing the
;                     galaxy properties

compile_opt idl2

; Find the files to write
FirstFile=min(GalstructArray.FileTreeNr)
LastFile=max(GalstructArray.FileTreeNr)
ngal=n_elements(GalStructArray)

for fnr = FirstFile, LastFile do begin
   print,'Writing file ',fnr
   index=where(GalStructArray.FileTreeNr eq fnr)
; First write the header
   fname = strcompress(FileRoot+'_'+string(fnr), /remove_all)
   if keyword_set(swap_endian) then $
      openw, 1, fname, /swap_endian $
   else $
      openw, 1, fname
   One = 0L     & writeu, 1, One
   writeu, 1, StructSize
   writeu, 1, n_elements(index)
   Zero = lonarr(StructSize/4-3) & writeu, 1, Zero
   ; Now write out all the Structures
   writeu, 1, GalStructArray[index]
   close, 1
endfor

end
