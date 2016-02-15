;+
; NAME:
;   is_file_readable
; CALLING SEQUENCE:
;   good= is_file_readable(filename,/compress)
; OUTPUT:
;   good  - 1 if good, 0 if bad
; BUGS:
;   - NOT extensively tested.
; REVISION HISTORY:
;   2005-06-03  converted from Schlegel email to code - Hogg
;-
function is_file_readable, filename,compress=compress
openr, rlun,filename,/get_lun,compress=compress
qgoodfile = 0B                ; Default to saying that the file is bad
on_ioerror, bad_file_label
image1 = uintarr(300, 300)
readu, rlun, image1
qgoodfile = 1B          ; We only get to this line if the file is good
bad_file_label: free_lun, rlun
return, qgoodfile 
end
