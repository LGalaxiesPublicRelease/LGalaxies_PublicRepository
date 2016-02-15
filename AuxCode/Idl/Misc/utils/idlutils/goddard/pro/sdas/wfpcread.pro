pro wfpcread,file,chip,wfhdr,wfimg,par
;+
; NAME:
;       WFPCREAD
; PURPOSE:
;       Read designated header and chip of a WFPC1 image
; EXPLANATION:
;       This procedure is designed to read the designated header and chip of a
;       WFPC image.  If the PAR input parameter is supplied, then the group
;       PARameter byte array is is returned.  If it is not, then the  header 
;       is modified by placing all the group parameters in the header as data 
;       cards.
;
;       Use the procedure WFPC2_READ to read WFPC2 images.
; CALLING SEQUENCE:
;       WFPCREAD, file, chip, wfhdr, wfimg, par
;
; INPUT:
;       FILE - The filename of the Header file of the image
;       CHIP - The chip number to read (usually 0-3)
;
; OUTPUT:
;       WFHDR -  Returned WF/PC header in a string array
;       WFIMG -  Returned WF/PC float image array
;
; OPTIONAL OUTPUT:
;       PAR  -  PARameter byte array (for group format header)
;
; HISTORY:
;       25-JUN-1990 Version 1 written
;       2-APR-1992 Added code to add CAM and CHIP onto the FILTNAM1    EWD
;       27-JUL-1992 Proper Header finally added  (E. Deutsch)
;       Converted to IDL V5.0   W. Landsman   September 1997
;-

  arg=n_params(0)
  if (arg lt 3) then begin
    print,'Call: IDL> WFPCREAD,header_file,chip,header,image,[params]'
    print,"e.g.: IDL> WFPCREAD,'W034234et.c0h',0,wfhdr,wfimg,wfpar"
    return
    endif

  asfix=1
  if (chip lt 0) then begin
    asfix=0 & chip=-chip-10
    endif

  AFhread,file,wfhdr
  sxopen,1,file
  wfimg=sxread(1,chip,par)
  close,1

  if (arg lt 5) then extgrp,wfhdr,par
  if (asfix eq 1) then astrmfix,wfhdr,chip+1

  pmode=sxpar(wfhdr,'PHOTMODE')
  filtr=strn(sxpar(wfhdr,'FILTNAM1'))
  if (strpos(filtr,' ') eq -1) then begin
    filtr=strn(filtr)+' '+strmid(pmode,0,2)+strmid(pmode,3,1)
    sxaddpar,wfhdr,'FILTNAM1',filtr,' First filter name and Chip No.'
    endif

  return
end
