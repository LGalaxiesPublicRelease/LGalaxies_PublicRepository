pro AFhread,HdrFile,hdr
;+
; NAME:
;      AFhread
; PURPOSE:
;      Subroutine of WFPCREAD to read a GEIS header from an HST STSDAS image.
; EXPLANATION:
;       This procedure reads a GEIS header from an HST image.   It then looks
;       if a .SHH file is present for FOC images to calculate better 
;       astrometry by getting the current PSANGLV3 from this file.   Called by
;        WFPCREAD.PRO
;
; CALLING SEQUENCE:
;       AFhread, HdrFile, hdr
;
; INPUTS:
;       HdrFile - scalar string giving name of STSDAS header for an FOC image   
;
; OUTPUTS:
;       hdr - string array, FITS header for the FOC image.    The position
;               angle of the V3 axis of HST (PSANGLV3) is added, if it could 
;               be found in the .SHH file       
; PROCEDURE CALLS:
;       STRN(), SXADDPAR, SXHREAD, SXPAR()
; REVISION HISTORY:
;       Written         Eric W. Deutsch  (U. of Washington)    June, 1994
;       Documentation update   W. Landsman  (HSTX)             July, 1994
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Removed call to EXIST() function   W. Landsman        April 1999
;-

  if (n_params(0) lt 2) then begin
    print,'Call> AFhread,Header_File,Returned_Header_Array'
    print,"e.g.> AFhread,'X0C35123T.D0H',h"
    return
    endif

  sxhread,HdrFile,hdr
  INSTRUME = sxpar(hdr,'INSTRUME')
  if (strn(INSTRUME) ne 'FOC') then return

  i=strlen(HdrFile)-1
  while (i gt 0) and (strmid(HdrFile,i,1) ne '.') and $
    (strmid(HdrFile,i,1) ne ']') and (strmid(HdrFile,i,1) ne ':') do i=i-1
  SHHFile=HdrFile
  if (strmid(HdrFile,i,1) eq '.') then SHHFile=strmid(HdrFile,0,i)
  SHHFile=SHHFile+'.shh'

  test = findfile(SHHFile, Count = count)
  if count GT 0 then sxhread,SHHFile,SHH else return
  PSANGLV3=sxpar(SHH,'PSANGLV3')
  if (!err lt 0) then begin
    print,'Unable to find keyword PSANGLV3 in ',SHHFile
    return
    endif
  sxaddpar,hdr,'PSANGLV3',PSANGLV3,' Position angle of V3 axis of HST'

  return
  end
