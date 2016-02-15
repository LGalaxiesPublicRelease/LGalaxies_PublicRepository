pro CONV_UNIX_VAX, variable, SOURCE_ARCH=source
;+
; NAME:
;      CONV_UNIX_VAX
; PURPOSE:
;      To convert Unix IDL data types to Vax IDL data types. 
; EXPLANATION:
;      CONV_UNIX_VAX assumes the Unix IDL data type is IEEE standard in either
;      big-endian or little-endian format.
;
; CALLING SEQUENCE:
;      CONV_UNIX_VAX, variable, [ SOURCE_ARCH = ]
;
; PARAMETERS:
;      variable - The data variable to be converted.  This may be a scalar
;            or an array.  Valid datatypes are integer, longword,
;            floating point, and double precision. The result of the 
;            conversion is passed back in the original variable.
; OPTIONAL INPUT KEYWORD:  
;      SOURCE_ARCH = name (string) of source architecture
;            if using this function on a VAX, otherwise
;            !VERSION.ARCH is used to determine the conversion.
;            **If run on a VAX, the default is to assume the source to be
;            a little-endian machine with IEEE floating point
;            (e.g. MIPSEL or Alpha***).
; RESTRICTIONS:
;      Requires that data be from IEEE standard Unix machines
;      (e.g. SUN, MIPSEL, or Alpha).
; EXAMPLE:
;      Read a 100 by 100 matrix of floating point numbers from a data
;      file created on a Sun.  Then convert the matrix values into
;      VAX format.
;
;      IDL> openr,1,'vax_float.dat
;      IDL> data = fltarr(100,100)
;      IDL> forrd,1,data
;      IDL> CONV_UNIX_VAX,data,SOURCE_ARCH='sparc'
;
; MODIFICATION HISTORY:
;      Version 1      By John Hoegy            13-Jun-88
;      04-May-90 - WTT:  Created CONV_UNIX_VAX from VAX2SUN,
;                         reversing floating point procedure.
;       Modified  P. Keegstra             September 1994
;           Implemented MIPSEL and ALPHA architecture,
;           distinguishing VMS and OSF
;       Modified  P. Keegstra             February 1995
;           Added 386 PC based architectures
;       If since V5.1 then VMS is always little endian    June 1998
;       Convert to IDL V5.0   W. Landsman                 June 1998
;-                                   
;****************************************************************************
;
;  Check to see if VARIABLE is defined.
;
 if N_params() LT 1 then begin
      print,'Syntax - CONV_UNIX_VAX, variable, [ SOURCE_ARCH = ]
      return
 endif

 if n_elements(variable) eq 0 then begin
      print,'*** VARIABLE not defined, routine CONV_UNIX_VAX.'
      retall
 endif

if N_elements( source ) EQ 1 then arch = source  else arch = !VERSION.ARCH 
 little_endian = 0

CASE arch OF

"sparc":                    ;Assume default big-endian

; Demo version of PV-WAVE for Linux reports itself as arch="i386".
; IDL for MS-WINDOWS reports itself as arch="3.1".

'i386': little_endian = 1
'3.1':  little_endian = 1
'mipsel': little_endian = 1
'386':  little_endian = 1
'386i': little_endian = 1
'x86': little_endian = 1

 "vax": BEGIN
        message,"machine is VAX, " + $
                "will assume source has little-endian " + $
                "architecture and IEEE floating point",/CONTIN
      little_endian = 1
      END

 "alpha": BEGIN
        IF !VERSION.OS EQ 'vms' THEN BEGIN
            if !VERSION.RELEASE LT '5.1' then $
                    message,"machine is alpha running VMS, " + $
                    "will assume source has little-endian " + $
                    "architecture and IEEE floating point",/CONTIN
          little_endian = 1
        ENDIF ELSE little_endian = 1
      END

 else:                  ;default is to assume big endian architecture
 ENDCASE
;
 if little_endian then begin
      swap_ints = 0
      swap_float = 2
 endif else begin
      swap_ints = 1
      swap_float = 1
 endelse
 
var_chars = size(variable)
var_type = var_chars[var_chars[0]+1]
 
;
case var_type of
  1: return                             ; byte

  2: if (swap_ints GT 0) then byteorder,variable,/SSWAP    ;integer

  3: if (swap_ints GT 0) then byteorder,variable,/LSWAP         ;longword

  4: BEGIN                             ; floating point
        scalar = (var_chars[0] eq 0)
        var_elems = long(var_chars[var_chars[0]+2])
        byte_elems = var_elems*4L
        byte_eq = byte(variable, 0, byte_elems)
    ;
        if (swap_float GT 1) then byteorder, byte_eq, /LSWAP
    ;
        i1 = lindgen(byte_elems/4L)*4L
        i2 = i1 + 1L
        biased = byte((byte_eq[i1] AND '7F'X) * 2) OR byte(byte_eq[i2]/128L)
        i = where(biased ne 0)
        if (i[0] ne -1) then biased[i] = byte(biased[i] + 2)
        byte_eq[i1] = byte(byte_eq[i1] AND '80'X) OR byte(biased/2)
        byte_eq[i2] = byte(byte_eq[i2] AND '7F'X) OR byte(biased*128)
    ; 
    ; swap bytes
    ;
        byte_elems = byte_elems + 3L
        byteorder, byte_eq, /SSWAP
;
        if scalar then begin
           tmp = fltarr(1)
           tmp[0] = float(byte_eq, 0, var_elems)
           variable = tmp[0]
           endif else variable[0] = float(byte_eq, 0, var_elems)
        return
     END

  5: BEGIN                         ; double precision
        var_elems = long(var_chars[var_chars[0]+2])
        byte_elems = var_elems*8L
      scalar = (var_chars[0] eq 0)
        if scalar then begin
             tmp = dblarr(1)
           tmp[0] = variable
                 byte_eq = byte(tmp, 0, byte_elems)
        endif else byte_eq = byte(variable, 0, byte_elems)
    ;
    ;  Bring it up to at least a double-precision level.
    ;
        byte_elems = byte_elems + 7L
        i1 = lindgen(byte_elems/8L)*8L
        i2 = i1 + 1L
        i3 = i2 + 1L
        i4 = i3 + 1L
        i5 = i4 + 1L
        i6 = i5 + 1L
        i7 = i6 + 1L
        i8 = i7 + 1L
    ;
      if (swap_float GT 1) then begin
             byte_eq2     = bytarr(byte_elems)
           byte_eq2[i1] = byte_eq[i8]
           byte_eq2[i2] = byte_eq[i7]
           byte_eq2[i3] = byte_eq[i6]
           byte_eq2[i4] = byte_eq[i5]
           byte_eq2[i5] = byte_eq[i4]
           byte_eq2[i6] = byte_eq[i3]
           byte_eq2[i7] = byte_eq[i2]
           byte_eq2[i8] = byte_eq[i1]
             byte_eq      = byte_eq2
        endif
    ;
    ;  Bring it up to at least a double-precision level.
    ;

        exponent = fix( ((byte_eq[i1] AND '7F'X)*16) OR $
                 ((byte_eq[i2] AND 'F0'X)/16) )
        i = where(exponent ne 0)
        if (i[0] ne -1) then exponent[i] = exponent[i] + 128 - 1022
        tmp1 = byte_eq[i8]
        byte_eq[i8] = ((byte_eq[i7] and '1f'x)*8) or ((tmp1 and 'e0'x)/32)
        tmp2 = byte_eq[i7]
        byte_eq[i7] = (tmp1 and '1f'x)*8
        tmp3 = byte_eq[i6]
        byte_eq[i6] = ((byte_eq[i5] and '1f'x)*8) or ((tmp3 and 'e0'x)/32)
        tmp1 = byte_eq[i5]
        byte_eq[i5] = ((tmp3 and '1f'x)*8) or ((tmp2 and 'e0'x)/32)
        tmp2 = byte_eq[i4]
        byte_eq[i4] = ((byte_eq[i3] and '1f'x)*8) or ((tmp2 and 'e0'x)/32)
        tmp3 = byte_eq[i3]
        byte_eq[i3] = ((tmp2 and '1f'x)*8) or ((tmp1 and 'e0'x)/32)
        tmp1 = byte_eq[i2]
        byte_eq[i2] = (byte_eq[i1] and '80'x) or byte((exponent and 'fe'x)/2)
        byte_eq[i1] = byte((exponent and '1'x)*128) or ((tmp1 and 'f'x)*8) or $
             ((tmp3 and 'e0'x)/32)
;
        if scalar then begin
           tmp = dblarr(1)
           tmp[0] = double(byte_eq, 0, var_elems)
           variable = tmp[0]
           endif else variable[0] = double(byte_eq, 0, var_elems)
        return
     END

  6: begin                  ; complex
       rvalue = float(variable)
       ivalue = imaginary(variable)
       conv_unix_vax,rvalue, SOURCE_ARCH = source
       conv_unix_vax,ivalue, SOURCE_ARCH = source
       variable = complex(rvalue,ivalue)
       end
 
  7: return                     ; string

  else: begin                   ; unknown
       print,'*** Data type ' + strtrim(var_type,2) + $
                  ' unknown, routine CONV_UNIX_VAX.'
       retall
       end
  endcase
return
end
