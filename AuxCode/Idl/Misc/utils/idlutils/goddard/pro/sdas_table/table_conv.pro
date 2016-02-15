 pro table_hconv,name, convtype = convtype
;
; NAME:
;	TABLE_HCONV
; PURPOSE:
;	Convert STSDAS Tables header into the host format, 
;	This is not a standalone procedure but called by TABLE_CONV
; CALLING SEQUENCE:
;     	TABLE_HCONV, name, convtype=convtype
; INPUT PARAMETERS:
;	name - name of table, scalar string
; REQUIRED INPUT KEYWORD:
;	CONVTYPE  - integer scalar specifying conversion direction
;	1 -  Host is big endian, table origin is VMS
;	2 -  Host is big endian, table origin is little endian
;	3 -  Host is VMS, table origin is little endian 
;	4 -  Host is VMS, table origin is big endian
;	5 -  Host is little endian, table origin is VMS
;	6 -  Host is little endian, table origin is big endian
; NOTE:
;		tcb(*,i) contains description of column i
;		   word 0	column number
;			1	offset for start of row in units of 2-bytes
;			2	width or column in 2-byte units
;			3	data type
;					6 = real*4
;					7 = real*8
;					4 = integer*4
;					1 = boolean*4
;					2 = character string
;			4-8	ascii column name up to 19 characters
;			9-13	column units (up to 19 characters)
;			14-15	format string
; HISTORY:
;	Written,  W. Landsman         
;	Adapted from GHRS version by Keith Feggans
;
;
;------------------------------------------------------------------------

; open file

 get_lun,unit
 fname = name
 if strpos(fname,'.') lt 1 then fname = strtrim(fname)+'.tab'
 openu,unit,fname,/block,error=err                  
 if err LT 0 then begin
	free_lun,unit
	print,'TABLE_CONV - Error opening file ' + fname
	print,strmessage(-!err)                                
	return
 end

; read header record

 lrec = assoc(unit,lonarr(12),0)
 h = lrec[0]
 case convtype of 
		1: h = conv_vax_unix(h)
		2: h = swap_endian(h)
		3: conv_unix_vax,h
		4: ieee_to_host,h
		5: h = conv_vax_unix(h)
		6: ieee_to_host,h
		else: message,' Conversion for this architecture is not avail.'
 end
 lrec[0]=h
 maxcol=h[5]

; Convert the TCB

 sbyte=12*4+h[1]*80L			;starting byte of column descriptions
 lrec = assoc(unit,lonarr(16,maxcol),sbyte)
 tcb=lrec[0]			;read col. descriptions
;		   word 0	column number
;			1	offset for start of row in units of 2-bytes
;			2	width or column in 2-byte units
;			3	data type
;					6 = real*4
;					7 = real*8
;					4 = integer*4
;					1 = boolean*4
;					2 = character string
;			4-8	ascii column name up to 19 characters
;			9-13	column units (up to 19 characters)
;			14-15	format string
 for ii=0,3 do begin
	tmp = tcb[ii,*]
	; do not convert TEXT fields(4 thru n)
	case convtype of 
		1: tmp = conv_vax_unix(tmp)
		2: tmp = swap_endian(tmp)
		3: conv_unix_vax,tmp
		4: ieee_to_host,tmp
		5: tmp = conv_vax_unix(tmp)
		6: ieee_to_host,tmp
		else: message,' Conversion for this architecture is not avail.'
	end
	tcb[ii,*]= tmp
 end
 lrec[0]=tcb

 free_lun,unit
 return
end


 pro table_conv,filespec,from_vms = from_vms, from_little = from_little
;+
; NAME:
;	TABLE_CONV
;
; PURPOSE:
;	Convert STSDAS table(s) to the host format
;
; EXPLANATION:
;	If on a BIG_ENDIAN machine (e.g. SparcStation), assumes table came
;		from a little endian machine unless /FROM_VMS keyword is set
;	If on a LITTLE_ENDIAN machine (e.g. OSF, Windows), assumes table came
;		from a big endian machine unless /FROM_VMS keyword is set
;	If on a VMS machine, assumes table came from a big endian machine
;		unless the /FROM_LITTLE keyword is set
;
; CALLING SEQUENCE:
;	TABLE_CONV, filespec, [ /FROM_VMS, /FROM_LITTLE ]
;
; INPUT PARAMETERS:
;	filespec - file specification for table(s), scalar string.
;		Can include wildcard values, e.g. '*.tab'
;
; EXAMPLE:
;	(1) An STSDAS table "calspec.tab" has been FTP'ed from a Sparcstation 
;	to a VMS machine.   Convert the table to the host VMS format.
;	(The FTP mode should be set to binary when copying STSDAS tables)
;
;	IDL> table_conv, 'calspec.tab'
;
;	(2) A set of files '*.tab' have been FTP'ed from VMS machine to a 
;	Sparcstation.   Convert all the files to the host format
;
;	IDL> table_conv, '*.tab', /FROM_VMS
;
; NOTES:
;	TABLE_CONV does not check whether byte-swapping is actually needed.
;	If this procedure is applied to a file that is already in the host
;	format, then that file will be corrupted.
;
; PROCEDURES CALLED:
;	CONV_VAX_UNIX(), CONV_UNIX_VAX, FDECOMP, IS_IEEE_BIG(), SWAP_ENDIAN(), 
;	TABLE_HCONV, TAB_PUT, TAB_READ, TAB_SIZE, TAB_WRITE
; MODIFICATION HISTORY:
;	W. Landsman,  Hughes STX/Goddad              July 1996 
;	Adapted from GHRS version by Don Lindler, Keith Feggans 
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-------------------------------------------------------------------------------

 if N_params() lt 1 then begin
	print,'Syntax: TABLE_CONV, filespec, [/FROM_VMS, /FROM_LITTLE ]'
	return
 endif

 list = findfile(filespec, count = n)
 if n EQ 0 then message, $
	'ERROR - no files found meeting specification ' + filespec
 big_end = is_ieee_big()
 if !VERSION.OS EQ 'vms' then vms = 1 else vms = 0
 little_end = (not vms) and (not big_end)
 if N_elements(from_vms) EQ 0 then from_vms = 0
 if N_elements(from_little) EQ 0 then from_little = 0

 if big_end then if from_vms then convtype = 1 else convtype = 2
 if vms then if from_little then convtype = 3 else convtype = 4
 if little_end then if from_vms then convtype = 5 else convtype = 6

 case convtype of
 1: message,/INF,'Converting table from VMS to big endian format'
 2: message,/INF,'Converting table from little endian to big endian format'
 3: message,/INF,'Converting table from little endian to VMS format'
 4: message,/INF,'Converting table from big endian to VMS format'
 5: message,/INF,'Converting table from VMS to little endian format'
 6: message,/INF,'Converting table from big endian to little endian format'
 else: message,'ERROR - unrecognized conversion number ' + strtrim(convtype)
 endcase

 for j=0,n_elements(list)-1 do begin
    if !VERSION.OS EQ 'vms' then begin    ;Remove version number
	fdecomp,list[j],disk,dir,name,ext
	list[j] =  disk + dir + name + '.' + ext
    endif
    print,list[j]
    table_hconv,list[j],convtype = convtype  ; byte swap the header.
    tab_read,list[j],tcb,tab,h
    tab_size,tcb,nrows,ncols
    for i=1,ncols do begin
        x = tab_val(tcb,tab,i)
	case convtype of 
		1: x = conv_vax_unix(x)
		2: x = swap_endian(x)
		3: conv_unix_vax,x
		4: ieee_to_host,x
		5: x = conv_vax_unix(x)
		6: ieee_to_host,x
	end
        tab_put,i,x,tcb,tab
    end
    tab_write,list[j],tcb,tab,h
 end
 return
 end
