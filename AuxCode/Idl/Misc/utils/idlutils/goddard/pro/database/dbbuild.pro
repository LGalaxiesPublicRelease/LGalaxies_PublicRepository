pro dbbuild,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18, $
    v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30, $
    NOINDEX = noindex, STATUS=STATUS, SILENT=SILENT
;+
; NAME:
;	DBBUILD
; PURPOSE:
;	Build a database by appending new values for every item.  
; EXPLANATION:
;	The database must be opened for update (with DBOPEN) before calling 
;	DBBUILD.
;
; CALLING SEQUENCE:
;	DBBUILD, [ v1, v2, v3, v4......v30, /NOINDEX, /SILENT, STATUS =  ]
;
; INPUTS:
;	v1,v2....v30 - vectors containing values for all items in the database.
;         V1 contains values for the first item, V2 for the second, etc.
;         The number of vectors supplied must equal the number of items
;         (excluding entry number) in the database.  The number of elements 
;         in each vector should be the same.   A multiple valued item
;         should be dimensioned NVALUE by NENTRY, where NVALUE is the number
;         of values, and NENTRY is the number of entries.
;
; OPTIONAL INPUT KEYWORDS:
;	NOINDEX - If this keyword is supplied and non-zero then DBBUILD will
;             *not* create an indexed file.    Useful to save time if
;             DBBUILD is to be called several times and the indexed file need
;             only be created on the last call
;
;	SILENT  - If the keyword SILENT is set and non-zero, then DBBUILD
;	      will not print a message when the index files are generated
;
; OPTIONAL OUTPUT KEYWORD:
;	STATUS - Returns a status code denoting whether the operation was
;	      successful (1) or unsuccessful (0).  Useful when DBBUILD is
;	      called from within other applications.
;
; EXAMPLE:
;	Suppose a database named STARS contains the four items NAME,RA,DEC, and 
;	FLUX.   Assume that one already has the four vectors containing the
;	values, and that the database definition (.DBD) file already exists.
;
;	IDL> !PRIV=2                  ;Writing to database requires !PRIV=2
;	IDL> dbcreate,'stars',1,1   ;Create database (.DBF) & index (.DBX) file
;	IDL> dbopen,'stars',1         ;Open database for update
;	IDL> dbbuild,name,ra,dec,flux ;Write 4 vectors into the database
;
; NOTES:
;	Do not call DBCREATE before DBBUILD if you want to append entries to
;	an existing database
;
;	DBBUILD checks that each value vector matches the idl type given in the
;	database definition (.DBD) file, and that character strings are the 
;	proper length. 
; REVISION HISTORY:
;	Written          W. Landsman           March, 1989
;	Added /NOINDEX keyword           W. Landsman        November, 1992
;	User no longer need supply all items   W. Landsman  December, 1992 
;	Added STATUS keyword, William Thompson, GSFC, 1 April 1994
;	Added /SILENT keyword, William Thompson, GSFC, October 1995
;	Allow up to 30 items, fix problem if first item was multiple value
;				  W. Landsman    GSFC, July 1996
;	Faster build of external databases on big endian machines 
;				  W. Landsman    GSFC, November 1997  
;	Converted to IDL V5.0   W. Landsman 24-Nov-1997
;       Use SIZE(/TNAME) for error mesage display  W.Landsman   July 2001
;       Fix message display error introduced July 2001  W. Landsman   Oct. 2001 
;-
  COMPILE_OPT IDL2
  On_error,2                            ;Return to caller
  npar = N_params()
  if npar LT 1 then begin
    print,'Syntax - DBBUILD, v1, [ v2, v3, v4, v5, ... v30,' 
    print,'         /NOINDEX, /SILENT, STATUS =  ]'
    return
  endif

 dtype = ['UNDEFINED','BYTE','INTEGER*2','INTEGER*4','REAL*4','REAL*8', $
        'COMPLEX','STRING','STRUCTURE','DCOMPLEX','POINTER','OBJECT', $ 
        'UNSIGNED*2', 'UNSIGNED*4', 'INTEGER*8','UNSIGNED*8']

 
;  Initialize STATUS as unsuccessful (0).  If the routine is successful, this
;  will be updated below.

  status = 0

  nitem = db_info( 'ITEMS' )
  if N_elements( nitem ) EQ 0 then return

   items = indgen(nitem)
   db_item, items, itnum, ivalnum, idltype, sbyte, numvals, nbyte
  for i = 1,npar do begin    ;Get the dimensions and type of each input vector

    ii = strtrim(i,2)
    test = execute('s=size(v' + ii +')' )

    if s[s[0] + 1] NE idltype[i] then begin
        message, 'Item ' + strtrim( db_item_info('NAME',i),2) + $
           ' - parameter '+strtrim(i,2) + ' - has an incorrect data type',/INF
        message, 'Required data type is ' + dtype[idltype[i]], /INF
        message, 'Supplied data type is ' + dtype[s[s[0]+1]], /INF
        return
     endif

  endfor
  external = db_info('external',0)
  if external then noconvert = is_ieee_big() else noconvert = 1b

  nitems = ( (npar<nitem) GE indgen(31))
  entry = make_array( DIMEN = db_info('LENGTH'),/BYTE ) ;Empty entry array
  nvalues = long( db_item_info( 'NVALUES' ) )       ;# of values per item
  nbyte = nbyte*nvalues                             ;Number of bytes per item
  Nv = N_elements(v1)/nvalues[1]                   
  for i = 0l, Nv - 1 do begin

       i1 = i*nvalues         
       i2 = i1 + nvalues -1
        dbxput,0l,entry,idltype[0],sbyte[0],nbyte[0]
        dbxput,v1[i1[1]:i2[1]],entry,idltype[1],sbyte[1],nbyte[1]
       if nitems[2] then begin
        dbxput,v2[i1[2]:i2[2]],entry,idltype[2],sbyte[2],nbyte[2]
       if nitems[3] then begin 
        dbxput,v3[i1[3]:i2[3]],entry,idltype[3],sbyte[3],nbyte[3]
       if nitems[4] then begin 
        dbxput,v4[i1[4]:i2[4]],entry,idltype[4],sbyte[4],nbyte[4]
       if nitems[5] then begin 
        dbxput,v5[i1[5]:i2[5]],entry,idltype[5],sbyte[5],nbyte[5]
       if nitems[6] then begin 
        dbxput,v6[i1[6]:i2[6]],entry,idltype[6],sbyte[6],nbyte[6]
       if nitems[7] then begin 
        dbxput,v7[i1[7]:i2[7]],entry,idltype[7],sbyte[7],nbyte[7]
       if nitems[8] then begin 
        dbxput,v8[i1[8]:i2[8]],entry,idltype[8],sbyte[8],nbyte[8]
       if nitems[9] then begin 
        dbxput,v9[i1[9]:i2[9]],entry,idltype[9],sbyte[9],nbyte[9]
       if nitems[10] then begin 
        dbxput,v10[i1[10]:i2[10]],entry,idltype[10],sbyte[10],nbyte[10]
       if nitems[11] then begin 
        dbxput,v11[i1[11]:i2[11]],entry,idltype[11],sbyte[11],nbyte[11]
       if nitems[12] then begin 
        dbxput,v12[i1[12]:i2[12]],entry,idltype[12],sbyte[12],nbyte[12]
       if nitems[13] then begin 
        dbxput,v13[i1[13]:i2[13]],entry,idltype[13],sbyte[13],nbyte[13]
       if nitems[14] then begin
        dbxput,v14[i1[14]:i2[14]],entry,idltype[14],sbyte[14],nbyte[14]
       if nitems[15] then begin
        dbxput,v15[i1[15]:i2[15]],entry,idltype[15],sbyte[15],nbyte[15]
       if nitems[16] then begin
        dbxput,v16[i1[16]:i2[16]],entry,idltype[16],sbyte[16],nbyte[16]
       if nitems[17] then begin
        dbxput,v17[i1[17]:i2[17]],entry,idltype[17],sbyte[17],nbyte[17]
       if nitems[18] then begin
        dbxput,v18[i1[18]:i2[18]],entry,idltype[18],sbyte[18],nbyte[18]
       if nitems[19] then begin
          dbxput,v19[i1[19]:i2[19]],entry,idltype[19],sbyte[19], nbyte[19]
       if nitems[20] then begin
          dbxput, v20[ i1[20]:i2[20] ],entry,idltype[20],sbyte[20],nbyte[20]
       if nitems[21] then begin
          dbxput,v21[i1[21]:i2[21]],entry,idltype[21],sbyte[21],nbyte[21]
       if nitems[22] then begin
          dbxput,v22[i1[22]:i2[22]],entry,idltype[22],sbyte[22],nbyte[22]
       if nitems[23] then begin
          dbxput,v23[i1[23]:i2[23]],entry,idltype[23],sbyte[23],nbyte[23]
       if nitems[24] then begin
          dbxput,v24[i1[24]:i2[24]],entry,idltype[24],sbyte[24],nbyte[24]
       if nitems[25] then $
          dbxput,v25[i1[25]:i2[25]],entry,idltype[25],sbyte[25],nbyte[25]
       if nitems[26] then $
          dbxput,v26[i1[26]:i2[26]],entry,idltype[26],sbyte[26],nbyte[26]
       if nitems[27] then $
          dbxput,v27[i1[27]:i2[27]],entry,idltype[27],sbyte[27],nbyte[27]
       if nitems[28] then $
          dbxput,v28[i1[28]:i2[28]],entry,idltype[28],sbyte[28],nbyte[28]
       if nitems[29] then $
          dbxput,v29[i1[29]:i2[29]],entry,idltype[29],sbyte[29],nbyte[29]
       if nitems[30] then $
          dbxput,v30[i1[30]:i2[30]],entry,idltype[30],sbyte[30],nbyte[30]

     endif & endif & endif & endif & endif & endif & endif & endif & endif
     endif & endif & endif & endif & endif & endif & endif & endif & endif
     endif & endif & endif & endif & endif

     dbwrt,entry,noconvert=noconvert        ;Write the entry into the database

  endfor

  if not keyword_set( NOINDEX ) then begin

      indexed = db_item_info( 'INDEX' )      ;Need to create an indexed file?
      if total(indexed) GE 1 then begin
	   if not keyword_set(silent) then	$
	           message,'Now creating indexed files',/INF
           dbindex,items
       endif

  endif

  dbclose

;  Mark successful completion, and return.

  status = 1
  return
  end
