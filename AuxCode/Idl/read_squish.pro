function read_squish,FileRoot,prop_name,type,swap_endian=swap_endian

; Function to read in the data written by squish.pro

; Inputs:
;   FileRoot - Filename, excluding property
;   prop_name - property name
;   (Fileroot+'_'+prop_name should give file name)
;   type - property type
;   (The type of a variable can be found from size(variable,/type))
;   swap_endian - keyword to caused binary data to be read with
;                 oppostive endianness
; Output:
;   An array containing the required data

fname = strcompress(FileRoot+'_'+prop_name, /remove_all)
if keyword_set(swap_endian) then $
  openr, 1, fname, /swap_endian $
else $
  openr, 1, fname

NGals = 0L     & readu, 1, NGals
prop_data=make_array(Ngals,type=type) & readu, 1, prop_data

close, 1

return, prop_data

end

