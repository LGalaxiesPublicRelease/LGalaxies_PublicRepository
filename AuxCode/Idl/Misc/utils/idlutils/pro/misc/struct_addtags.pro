;+
; NAME:
;   struct_addtags
;
; PURPOSE:
;   Add tags from one structure (array) to another
;
; CALLING SEQUENCE:
;   outstruct = struct_addtags(astruct, bstruct)
;
; INPUTS:
;   astruct    - First structure, which can be an array
;   bstruct    - Second structure, which can be an array
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   outstruct  - Ouput structure array
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The dimensions of the output array match that of ASTRUCT.  For example,
;   if ASTRUCT has dimensions [5,10], and BSTRUCT has dimensions [2,25],
;   the output structure has dimensions [5,10].
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   copy_struct
;   copy_struct_inx
;
; REVISION HISTORY:
;   28-Jun-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function struct_addtags, astruct, bstruct

   if (N_elements(astruct) EQ 0) then $
    return, bstruct

   num1 = N_elements(astruct)
   num2 = N_elements(bstruct)
   if (num1 NE num2) then $
    message, 'Both structures must have the same number of elements'

   ;----------
   ; Create an empty structure with all the tags from both structures

   obj1 = create_struct(astruct[0], bstruct[0])
   dims = size(astruct,/dimens)
   outstruct = make_array(dimension=dims, value=obj1)

   ;----------
   ; Assign elements from ASTRUCT into the new output structure

   copy_struct, astruct, outstruct

   ;----------
   ; Assign elements from BSTRUCT into the new output structure

   copy_struct_inx, bstruct, outstruct

   return, outstruct
end
;------------------------------------------------------------------------------
