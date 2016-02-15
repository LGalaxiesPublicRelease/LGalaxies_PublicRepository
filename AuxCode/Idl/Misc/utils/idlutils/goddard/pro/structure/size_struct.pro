function size_struct, structure, struct, Max_Field_Size, PRINT=print 
;+
; NAME:
;	SIZE_STRUCT
; PURPOSE:
;	Obtain the size in bytes of an IDL structure definition.    
; EXPLANATION:
;	For most applications this function is superceded by use 
;	of the /LENGTH keyword to the intrinsic N_TAGS function 
;	introduced in IDL V2.3.0
;
; CALLING SEQUENCE:
;			bytes = size_struct( structure )
;	examples:
;			print, size_struct( "fdq_sdf" )
; INPUTS:
;		structure = a structure variable or
;				a string giving the structure name
;				as known by IDL (help,/struct,variable).
;		/PRINT = to print all sub structure sizes.
;
; inputs/outputs used recursively:
;		struct = the structure VARIABLE currently analyzed.
;		Max_Field_Size = size of the largest field found in structure.
; RESULT:
;		Function returns the total size in bytes of a structure element.
; PROCEDURE:
;		Strategy is to call size_struct recursively if
;		structure contains sub-structures.
;		Otherwise just add up the field sizes.
;
; MODIFICATION HISTORY:
;	written 1990 Frank Varosi STX @ NASA/GSFC (using align_struct).
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
	if N_tags( struct ) LE 0 then begin

		s = size( structure )

		CASE s[s[0]+1] OF
		7: BEGIN
			struct_name = strupcase( structure )
			com = "struct = replicate( {" + struct_name + "}, 1 )"
			status = execute( com )
			if NOT status then begin
				message," unable to define structure: " + $
							struct_name,/CONTIN
				return,(-1)
			   endif
		     END
		8: BEGIN
			struct_name = "?"
			struct = structure[0]
		     END
		else: BEGIN
			message,"specify STRUCTURE VARIABLE or NAME",/CONTIN
			return,(-1)
			END
		ENDCASE

		if keyword_set( print  ) then $
			print," Structure Name    Size   Fields  Max Field Size"

	  endif else begin

		s = size( structure )
		if  ( s[s[0]+1] EQ 7 )  then struct_name = structure $
					else struct_name = "?"
	   endelse

	Tags = Tag_Names( struct )
	Ntag = N_tags( struct )
	Max_Field_Size = 1
	struct_size = 0

	for t=0,Ntag-1 do begin

		s = size( struct.(t) )
		N_elem = s[ s[0]+2 ]
		Field_Type = s[ s[0]+1 ]

		if (Field_Type EQ 8) then begin

			ST_siz = size_struct( Tags[t], struct.(t), MaxF_siz, $
								   PRINT=print )
			struct_size = struct_size + ST_siz * N_elem
			Max_Field_Size = Max_Field_Size > MaxF_siz

		  endif else begin

			CASE Field_Type OF
			1: 	Field_Size = 1
			2: 	Field_Size = 2
			3: 	Field_Size = 4
			4: 	Field_Size = 4
			5: 	Field_Size = 8
			6: 	Field_Size = 8
			else: 	Field_Size = 0
			ENDCASE

			struct_size = struct_size + Field_Size * N_elem
			Max_Field_Size = Max_Field_Size > Field_Size

		   endelse
	  endfor

;write out the sizes of structure(s) :

	if keyword_set( print  ) then $
	   print, struct_name, struct_size, Ntag, Max_Field_Size, $
		FORM="(A15,I8,I8,I8)"

return, struct_size
end
