function stsub,unit,x1,x2,y1,y2,step
;+
; NAME:
;	STSUB
; PURPOSE:
;	Subroutine of STSUBIM to read a subset of a SDAS image file.   
; EXPLANATION:
;	User can specify a subimage range or a step size    Called by STSUBIM
;
; CALLING SEQUENCE:
;	Result =  stsub( unit, x1, x2, y1, y2, step)
;
; INPUTS:
;	UNIT      =  Unit number of file, must be from 1 to 9.
;		     Unit must have been opened with SXOPEN.
;       x1        =  lower x value
;       x2        =  upper x value
;       y1        =  lower y value
;       y2        =  upper y value
;	step      =  used to extract every nth pixel.  If step = 1, a full res.
; 	             subimage is extracted; step = 2, every other pixel is
;	             extracted, etc.  Defaults to 1.  The minimum value is 1.
; OUTPUTS:
;	Result of function = array constructed from designated record.
;
; COMMON BLOCKS:
;	Uses idl common stcommn to access parameters (see SXOPEN)
;
; MODIFICATION HISTORY:
;	Written, M. Greason, STX, July 1990.
;	Remove initialization of array for increased efficiency  January, 1991
;	Removed call to STSUBC.EXE, do it all in IDL
;                                      - K. Venkatakrishna, STX, April 1992
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
; common block containing description of file (see sxopen)
;
common stcommn,result,filename
;
; check if unit open
;
if ( unit LT 1 ) or ( unit GT 9 ) then $
	                message,'Invalid unit number'
;
if n_elements(result) EQ 0 then result=0
if ( N_elements(result NE 200) ) or (result[0,unit] NE 121147) then $
        message,'Header file has not been opened'
;
; check the step parameter.
;
if N_params() LT 6 then step = 1
step = step > 1
;
; get image descriptors.
;
desc = result[*,unit]		      ;description for unit
npix = desc[10]			      ;number of pixels/line
ysize = (y2 - y1) / step + 1          ;Y size of subarray
dtype = desc[8]                       ;Datatype
bpix = abs( desc[2] ) /8l
xsize = (x2 - x1) / step + 1          ;X size of subarray
skip = npix*step                      ;Bytes to skip between lines

; Create output subarray 

    sub = make_array( size = [2,[xsize,ysize], dtype, 0], /NOZERO)  


;  read the image every n-th line (where n=step), one line at a time 
;   and load the output buffer

 offset = y1*npix + x1
 line = make_array(TYPE = dtype, DIMEN = (xsize-1)*step+1,/NOZERO)
 if step GT 1 then  begin 
    index = indgen(xsize)*step
    for i= 0l, ysize - 1 do begin
          point_lun, unit, (offset + i*skip)*bpix
           readu,unit,line
           sub[0,i] = line[index] 
     endfor 
 endif else begin 
    for i= 0l, ysize - 1 do begin
          point_lun, unit, (offset + i*skip)*bpix
           readu,unit,line
           sub[0,i] = line
     endfor 
  endelse

return, sub
;
end
