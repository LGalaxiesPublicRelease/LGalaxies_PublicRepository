pro host_to_ieee, data, IDLTYPE = idltype
;+
; NAME:
;     HOST_TO_IEEE
; PURPOSE:
;     Translate an IDL variable from host to IEEE representation 
; EXPLANATION:
;     The variable is converted from the format used by the host architecture
;     into IEEE-754 representation (as used, for example, in FITS data ).
;
;     Duplicates most of the functionality of the SWAP_ENDIAN_INPLACE procedure
;     introduced in V5.6, with the addition of the IDLTYPE keyword.
; CALLING SEQUENCE:
;     HOST_TO_IEEE, data, [ IDLTYPE = ]
;
; INPUT-OUTPUT PARAMETERS:
;     data - any IDL variable, scalar or vector.   It will be modified by
;             HOST_TO_IEEE to convert from host to IEEE representation.  Byte 
;             and string variables are returned by HOST_TO_IEEE unchanged
;
; OPTIONAL KEYWORD INPUTS:
;     IDLTYPE - scalar integer (1-15) specifying the IDL datatype according
;               to the code given by the SIZE function.      This keyword
;               will usually be used when supplying a byte array that needs
;               to be interpreted as another data type (e.g. FLOAT).
;
; EXAMPLE:
;     Suppose FITARR is a 2880 element byte array to be converted to a FITS
;     record and interpreted a FLOAT data.
;
;       IDL> host_to_ieee, FITARR, IDLTYPE = 4
;
; METHOD:
;     The BYTEORDER procedure is called with the appropriate keywords
;
; MODIFICATION HISTORY:
;      Adapted from CONV_UNIX_VAX, W. Landsman   Hughes/STX    January, 1992
;      Version for IDL V5.0  August 1997
;      Converted to IDL V5.0   W. Landsman   September 1997
;      Added new integer datatypes  C. Markwardt/W. Landsman  July 2000
;      Use /SWAP_IF_LITTLE_ENDIAN keyword for 64bit types W. Landsman Feb 2003
;-
 On_error,2 

 if N_params() EQ 0 then begin
    print,'Syntax - HOST_TO_IEEE, data, [IDLTYPE = ]'
    return
 endif  

 npts = N_elements( data )
 if npts EQ 0 then $
     message,'ERROR - IDL data variable (first parameter) not defined'

 if not keyword_set( idltype) then idltype = size(data,/type)

 case idltype of

      1: return                             ;byte

      2: byteorder, data, /HTONS            ;integer

      3: byteorder, data, /HTONL            ;long

      4: byteorder, data, /FTOXDR           ;float

      5: byteorder,data,/DTOXDR              ;double
 
      6: byteorder, data, /FTOXDR
     
      7: return                             ;string

      8: BEGIN                              ;structure

        Ntag = N_tags( data )

        for t=0,Ntag-1 do  begin
          temp = data.(t)
          host_to_ieee, temp
          data.(t) = temp
        endfor 
       END

     9: byteorder, data, /DTOXDR
 
     12: byteorder, data, /HTONS

     13: byteorder, data, /HTONL

     14: byteorder, data, /L64swap, /SWAP_IF_LITTLE

     15: byteorder, data, /L64swap, /SWAP_IF_LITTLE

     else: message,'Unrecognized datatype ' + strtrim(idltype,2)

 ENDCASE

 return
 end 
