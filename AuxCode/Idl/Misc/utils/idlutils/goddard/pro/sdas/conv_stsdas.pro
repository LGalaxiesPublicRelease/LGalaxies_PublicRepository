pro conv_stsdas, sdas_name, FROM_IEEE = from_ieee
;+
; NAME:
;       CONV_STSDAS
; PURPOSE:
;       Convert internal format of an STSDAS image to host machine architecture
; EXPLANATION:
;       Converts the internal format of an STSDAS image (.hhh and .hhd file)
;       to the host machine architecture.     Useful for copying STSDAS files
;       between different machines.     If the host is not a VMS machine, then
;       by default CONV_STSDAS assumes the image originated on VMS.   If the
;       host is VMS, then CONV_STSDAS assumes that the image originated on
;       an IEEE machine (e.g. SparcStation).
;
; CALLING SEQUENCE:
;       CONV_STSDAS, sdas_name, [ /FROM_IEEE]
;
; INPUTS:
;       sdas_name - scalar string giving name of the STSDAS image
;               CONV_STSDAS assumes a default header extension of .hhh -- 
;               otherwise the header extension should be included in sdas_name.
;               The internal format of the file will be modified by CONV_STSDAS.
;
; OPTIONAL KEYWORD INPUT:
;       /FROM_IEEE - On little endian machines (OSF, windows) this keyword
;               indicates that the STSDAS file originated on an IEEE machine
;               (e.g SparcStation) rather than a VMS machine
;
; EXAMPLE:
;       Suppose files test.hhd and test.hhh have been copied with FTP from
;       a Vax to a Sparcstation.   Convert these files to the SparcStation
;       internal format.
;
;       IDL> conv_stsdas, 'test'
;
; METHOD:
;       CONV_STSDAS reads each group image and parameter block and uses 
;       IEEE_TO_HOST or CONV_VAX_UNIX to convert the internal format.   The
;       converted images and parameter blocks are written back to the orginal
;       file.
;
; PROCEDURE CALLS
;       sxopen, fdecomp, sxgpar(), sxpar(), ieee_to_host, conv_vax_unix()
;
; NOTES:
;       (1)  When copying STSDAS files to VMS, be sure the .hhh file is 
;               formatted as fixed block 80 byte.
;       (2)  CONV_STSDAS has no way of knowing if a file really came from
;               a different machine architecture.    If it is applied to a file
;               that already has the correct internal format, then CONV_STSDAS
;               will "convert" this file and corrupt the internal format.
;       (3)  Note that CONV_STSDAS currently does not support conversion *from*
;               a little-endian machine (OSF, windows)          
;
; REVISION HISTORY:
;       Written   W. Landsman                     January, 1993
;       Don't require .hhh extension            April, 1993
;       Increase speed by calling SXGINFO       May, 1993
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Replace DATATYPE() with size(/TNAME)  W. Landsman   November 2001
;-
 On_error,2

 if N_params() EQ 0 then begin
     print,'Syntax - CONV_STSDAS, sdas_name, [ /FROM_IEEE ]'
     return
 endif

 common stcommn, result, filename      ;SXOPEN common block

 fdecomp,sdas_name,disk,dir,name,ext       ;Decompose filename
 sxopen, 1, name + '.' + ext, h            ;Open file and find stcommn values
 close,1 
 desc = result[*, 1]
 ndimen = desc[3]                      ;Number of dimensions of each image
 dtype  =  desc[8]                     ;Datatype of each image
 dimen = desc[10:9+ndimen]             ;Vector of image dimensions

 ext = strmid(ext,0,2) + 'd'
 openu,lun, name + '.' + ext, /BLOCK, /GET_LUN
 gcount = sxpar( h, 'GCOUNT' ) > 1     ;Number of groups
 pcount = sxpar( h, 'PCOUNT' )         ;Number of group parameters

 fromieee = (!VERSION.OS eq "vms") or keyword_set( FROM_IEEE)

 for group = 0,gcount-1 do begin       ;Loop over each group
 
 message,'Processing group '+ strtrim( group , 2), /CONT

 if ( pcount GT 0 ) then begin         ;Parameter block in file?

   parrec = assoc(lun, bytarr(desc[7]),(group+1)*desc[9]-desc[7])
   par = parrec[0]
   if group EQ 0 then sxginfo, h, par, pdtype, pdsbyte, nbyte

   for i = 0, pcount-1 do begin          ;Loop over each parameter

       parval = sxgpar(h, par,i+1, pdtype[i], pdsbyte[i], nbyte[i] )
       if fromieee then ieee_to_host, parval $
                  else parval = conv_vax_unix( parval)
       if size( parval,/TNAME ) NE 'STRING' then $
          par[pdsbyte[i] ] = byte( parval, 0, nbyte[i] )
  endfor

  parrec[0] = par          ;Write back converted parameter block

 endif

  sbyte = long(group)*desc[9]
  rec  =  assoc( lun, make_array(size=[ndimen,dimen>1,dtype,0],/nozero), sbyte)
  im = rec[0]
  if fromieee then ieee_to_host, im $
             else im = conv_vax_unix( im)
  rec[0] = im
  endfor

  free_lun, lun
  return
  end
