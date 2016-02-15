function tab_null,values
;+
; NAME:
;	TAB_NULL  
; PURPOSE:
;	function to locate null values within a vector of values from
;	an STSDAS table.
;
; CALLING SEQUENCE
;	result = tab_null(values)
;
; INPUTS:
;	values - data value(s)
;
; OUTPUTS:
;	a boolean variable is returned with the same length as values.
;	1 indicates that the corresponding value was null
;
; OPERATIONAL NOTES:
;	Boolean columns in an STSDAS table does not presently
;	have the capability to flag null values.
;
; HISTORY:
;	version 1   D. Lindler   April 1989
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;------------------------------------------------------------------------------
;
; determine what the null value is
;
	s=size(values) & dtype=s[s[0]+1]
	case dtype of
		7: null=''		;string
		1: null=999		;byte/boolean (not implemented)
		2: null='80000000'XL	;integer*2 (should not occur)
		4: null=1.6e38		;real*4
		3: null='80000000'XL	;longword integer
		5: null=1.6d38		;real*8
	endcase
;
	if dtype ne 1 then return,values eq null
	nvals=n_elements(values)
	ndim=s[0]
	if ndim eq 0 then return,strtrim(nulltrim(values)) eq ''
	nulls=bytarr(nvals)
	for i=0,nvals-1 do nulls[i]=strtrim(nulltrim(values[i])) eq ''
	return,nulls
end
