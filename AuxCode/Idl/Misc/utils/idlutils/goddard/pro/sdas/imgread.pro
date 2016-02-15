pro IMGread,image,h,filename,group,NoAssoc=NoAssoc,silent=silent, $
            AstrmFix=AstrmFix
;+
; NAME:
;	IMGREAD
; PURPOSE:
;	Read a WFPC or FOC file into IDL image and data arrays
; EXPLANATION:
;	Open an SDAS/GEIS file and read the image into a data array of 
;	appropriate type and read the header into a string array.  This 
;	procedure was designed to be more versatile than the STRD procedure 
;	and to be specifically useful to WF/PC and FOC data, as well as all 
;	other GEIS images.  IMGread supports multiple GROUPS (i.e. in STSDAS 
;	format).
;
; CALLING SEQEUNCE:
;	IMGread,image,hdr,[filename],[groupno],[/NoAssoc,/silent,Astrmfix=]
;
; OPTIONAL INPUT:
;	FILENAME  The filename of the HEADER file (must have extention .xxh 
;		where xx may be any two alphanumerics but it usually hh.)  If 
;		there is no extention supplied, .hhh and .hhd are assumed.  If
;		this parameter is not supplied, a filename is prompted for, 
;		with the option of pressing [ENTER] to call the GETFILE() 
;		function provides a menu listing of available *.*h files.  If 
;		widgets are available, the function PICKFILE() is called 
;		instead.
;	GROUP - This parameter specifies the GROUP number image to read from a
;               file which contains multiple groups.  For example, for WF/PC
;               images where all four chips are contained in one file, one
;               specify a GROUP of 0 to read PC5, 1 for PC6, 3 for PC8,
;               0 to read WF1 for a WF image, etc.  therefore, the range of
;               GROUP is 0 to GCOUNT-1 (where GCOUNT is a header keyword.)
; OUTPUT:
;	IMAGE - The returned array which contains the pixel information.  
;		IMAGE will be of whatever datatype the header indicates (or 
;		seems to... i.e. if BITPIX=32 but there is no DATATYPE keyword,
;               IMGread assumes REAL*4 if BZERO is 0 or non-existant and
;               INTEGER*4 if BZERO is not 0.  This is usually right, but not
;               always.)
;	H - The returned string array containing the image header 
;		information as if SXHREAD were used.
;
; OPTIONAL KEYWORDS:
;	NoAssoc -  This keyword controls how IDL reads the file.    If NoAssoc
;		is set and non-zero then the READU function is used instead
;		of the ASSOC function.    The user can select the type of
;		read that gives the best performance on his particular setup.
;		In general, the ASSOC function seems to be faster, but is more
;		demanding on virtual memory.
;	SILENT - If this keyword is set and non-zero, then the "Loading..." 
;		message will not be ;               printed.
;	ASTRMFIX  Controls whether the procedure AstrmFix is run.  AstrmFix
;               calculates an astrometric solution from the HST Spacecraft
;               angle in the header.  CRPIXn and CRVALn are left alone.  Only
;               CDn_n are changed.  The Default is currently set to 1 since
;               correct astrometry still does not come with the headers.  Once
;               the astrometric fix is implemented in PODPS, the default should
;               be switched to 0.
;
; SIDE EFFECTS:
;	For an image with group parameters, all parameters are extracted from
;	the .HHD file and values are inserted into the returned header variable.
;	To get the original header, use SXHREAD for these type of image files.
;	The EXTGRP procedure takes care of this process.
;
; EXAMPLE:
;	Read the WF/PC file named 'w0hd0203t.c1h' into IDL variables, IM and H.
;
;	IDL> IMGREAD, im,h,'w0hd0203t.c1h'
;
; OTHER PROCEDURES CALLED:
;	SXPAR, SXADDPAR, SXOPEN, SXHREAD, FDECOMP, WFPCREAD, PICKFILE, EXTGRP
;
; HISTORY:
;	09-JUL-92 Header finally added to this procedure which has been in use
;	 for two or more years.  All versions and header by Eric W. Deutsch
;	01-APR-93 Made a few minor adjustments.  EWD.  (No, really)
;	July 93 Added /NoAssoc, MAKE_ARRAY, removed GET_FILE W. Landsman (HSTX)
;	Converted to IDL V5.0   W. Landsman   September 1997
;       Remove use of !ERR  W. Landsman   January 2001
;-

  On_error,2
  err='[IMGread] Error: '
  warn='[IMGread] Warning: '

  arg=n_params(0)
  if (arg lt 2) then begin
    print,'Call> IMGread,imagearray,header,[filename],[groupno],[/NoAssoc,/silent,Astrmfix=]'
    print,"e.g.> IMGread,img1,h1,'test.hhh'"
    return
    endif

  if (n_elements(silent) eq 0) then silent=0
  if (n_elements(AstrmFix) eq 0) then AstrmFix=1
  if (arg lt 3) then filename=''
  s=size(filename) & if (s[0] ne 0) and (s[0] ne 7) then filename=''

  if (arg lt 4) then group=-1
  chipno=group

GET_FILE:
   if filename EQ '' then begin
      if (!d.flags and 65536) EQ 65536 then begin 
	   message,'Select SDAS/GEIS header filename to read',/INF
           tmp = pickfile(filter='*.*h',/read)
      endif else begin
	    print,'Enter SDAS/GEIS filename to READ: ' +  $
		'(Include extension if not .hhh)'
	    print,'Hit [RETURN] to Cancel'
	    tmp='' & read,'Filename: ',tmp
      endelse
    if (tmp eq '') then return else filename = tmp
  endif

  fdecomp,filename,disk,dir,name,ext,ver
  if (ver ne '') then filename=disk+dir+name+'.'+ext
  if (ext eq '') then begin
    ext=ext+'.hhh' & filename=filename+ext
    endif

  if (strupcase(strmid(ext,2,1)) ne 'H') then begin
    print,err,'SDAS filename must have extension .xxh'
    filename='' & goto,GET_FILE
    endif

  find=findfile(filename,count=i)    ;Does file exist?
  if (i lt 1) then begin
    print,err,'Unable to find file '+filename
    filename='' & goto,GET_FILE
    endif

  on_ioerror,IOERR
  sxhread,filename,h

  SIMPLE=sxpar(h,'SIMPLE') & GROUPS=sxpar(h,'GROUPS') & PCOUNT=sxpar(h,'PCOUNT')
  if (SIMPLE eq 0) and (GROUPS eq 1) and (PCOUNT ne 0) then begin
    GCOUNT=sxpar(h,'GCOUNT')
    if (GCOUNT eq 0) then begin 
      print,err,'Unable to find keyword GCOUNT or GCOUNT=0' & return & endif
    chip=0
    if (chipno ne -1) then chip=chipno
    if (GCOUNT eq 1) then chip=-10
    if (GCOUNT gt 1) and (chipno eq -1) then begin
      print,'This image has several groups: Enter group to load: 0-',strn(GCOUNT-1)
      tmp='' & read,'Chip: ',tmp & if (tmp eq '') then return
      chip=fix(tmp)
      if (chip lt 0) or (chip ge GCOUNT) then begin
        message,'ERROR - Invalid chip number',/CON & return & endif
      endif

    NAXIS=sxpar(h,'NAXIS') & NAX=sxpar(h,'NAXIS*')
    st = strn(nax[0])
    if NAXIS GT 1 then for i=1,naxis-1 do st = st + ' x ' + strn(nax[i])
    if not silent then message,/INF,'Loading '+ filename + ' (' + st + ')...'
    if (AstrmFix eq 0) and (chip ge 0) then chip=-10-chip
    wfpcread,filename,chip,h,image
    return
    endif

  if (SIMPLE eq 1) or (PCOUNT eq 0) then begin
    NAXIS=sxpar(h,'NAXIS') & NAXIS1=sxpar(h,'NAXIS1') & NAXIS2=sxpar(h,'NAXIS2')
    BSCALE=sxpar(h,'BSCALE') & BZERO=sxpar(h,'BZERO')
    ORIGIN='?' & tmp1 = sxpar(h,'ORIGIN', Count = N_Origin) 
    if (N_Origin ge 1) then ORIGIN=tmp1

    dtype=0
    DATATYPE=strn(sxpar(h,'DATATYPE', Count = N_Datatype))
    if (N_Datatype GE 1) then begin
      case DATATYPE of        ;Convert datatype to type code
        'BYTE':                 dtype=1
        'LOGICAL*1':            dtype=1 ;Byte
        'INTEGER*1':            dtype=1
        'REAL*4':               dtype=4
        'INTEGER*2':            dtype=2
        'UNSIGNED*2':           dtype=2
        'INTEGER*4':            dtype=3
        'UNSIGNED*4':           dtype=3
        'REAL*8':               dtype=5
        'COMPLEX*8':            dtype=6
        else:                   message,'Invalid DATATYPE'
      endcase
    endif else begin
      BITPIX=sxpar(h,'BITPIX')
      case BITPIX of
          8:                    dtype=1              ;byte
         16:                    dtype=2              ;integer*2
         32:                    dtype=3              ;integer*4
        -32:                    dtype=4
         64:                    dtype=5
        -64:                    dtype=5
        else:                   message,'Unable to determine data type'
        endcase
      if (BITPIX eq 32) and (BZERO ne 0) then begin
        print,warn,'BITPIX=32 and no DATATYPE.  Assuming INTEGER*4' & dtype=3
        endif
      if (BITPIX eq 32) and (BZERO eq 0) then begin
        print,warn,'BITPIX=32 and no DATATYPE.  Assuming REAL*4' & dtype=4
        endif
      endelse

    chip=0
    if (chipno ne -1) then chip=chipno
    if (NAXIS eq 3) and (chipno eq -1) then begin
      NAXIS3=sxpar(h,'NAXIS3')
      if (NAXIS3 gt 1) then begin
        print,'This image has several planes: Enter plane to load: 0-',strn(NAXIS3-1)
        intmp='' & read,'Plane: ',intmp
        if (intmp eq '') then return
        chip=fix(intmp)
        if (chip lt 0) or (chip ge NAXIS3) then begin
          print,'Invalid plane number.' & return
          endif
      endif else chip=0
      endif

    imtyp=['NOTHING','BYTE','INTEGER','LONG INTEGER','FLOAT','DOUBLE FLOAT']
    if not silent then print,'Loading ',filename,' (',strn(NAXIS1),'x',strn(NAXIS2),' ', $
      imtyp[dtype],')...'

    size_arr = [2, NAXIS1, NAXIS2, dtype, NAXIS1*NAXIS2]

    if not keyword_set( NoAssoc ) then begin
      inpfil=8
      sxopen,inpfil,filename
      tmp = assoc( inpfil, make_array( SIZE = size_arr, /NoZero ) ) 
      image=tmp[0]
      close,inpfil

      endif else begin

; this method is used in STRD, but seems to be less memory efficient than the
;   the ASSOC method.  EWD 16-APR-1992
	image = make_array( SIZE = size_arr, /NoZero)

; this only works for the first chip... I doubt it's worth fixing... EWD
      dataname=strmid(filename,0,strlen(filename)-1)+'d'
      openr,unit,dataname,/BLOCK,/GET_LUN
      readu,unit,image
      free_lun,unit

      endelse

    if (dtype eq 3) and (BZERO ne 0) then begin
      print,'Applying BSCALE and BZERO...'
      image=TEMPORARY(image)*BSCALE+BZERO
      sxaddpar,h,'BSCALE',1 & sxaddpar,h,'BZERO',0
      Check_FITS,image,h,/sdas,/update
      endif

    ; FILTNAM1 gets modified here!  This is probably dishonest, but I'm getting
    ;  paid to do it.  EWD.
    if (strn(sxpar(h,'INSTRUME')) eq 'WFPC') then begin
      pmode=sxpar(h,'PHOTMODE')
      filtr=strn(sxpar(h,'FILTNAM1'))
      if (strpos(filtr,' ') eq -1) then begin
        filtr=strn(filtr)+' '+strmid(pmode,0,2)+strmid(pmode,3,1)
        sxaddpar,h,'FILTNAM1',filtr,' First filter name and Chip No.'
        endif
      endif

    return
    endif

  message,'Conflicting keywords in header.  Unrecognized save method.',/CON
  return

IOERR:
  print,err,'I/O Error.  Unable to read file ',filename
  return

end
