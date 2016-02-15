Pro MID_UP_TABLE,table,ncol,nrow,data
;---------------------------------------------------------------------
;+
; NAME:
;	MID_UP_TABLE
;
; PURPOSE:
;	Open and update data from a MIDAS table.
;
; CALLING SEQUENCE:
;	MID_UP_TABLE,table,ncol,nrow,data
;
; INPUTS:
;	Table =  file name of MIDAS table or Logical Unit Number.  
;		* If a filename is given, the file will be opened and closed 
;		using a local LUN.  An extension -- not to be
;		supplied -- of .TBL is assumed.  No version number is
;		allowed: the most recent version is used.  
;		* If an LUN is given, the file associated with that LUN will be 
;		used.
;	Ncol =   number of columns in the array to be written into the 
;		MIDAS table.  Long integer (I*4).  Must be less than or
;		equal to the number of columns in the MIDAS table.
;		"Ncol" and "nrow" are checked for consistency only.
;	Nrow =   number of rows in the array to be written into the MIDAS 
;		table. Long integer (I*4).  Must be less than or equal to
;		the number of rows in the MIDAS table.
;		"Ncol" and "nrow" are checked for consistency only.
;	Data =   array of data to be written into the MIDAS table.  Floating 
;		(R*4).  Data is of dimensions nrow*ncol.  The select column 
;		in the MIDAS table is disregarded.
;
; RESTRICTIONS:
;	Dimensions of the array to be written into the MIDAS table must be
;	less than or equal to the corresponding dimensions of the latter.
;	Note that in the case of smaller dimensions, only part of the given
;	MIDAS table in altered (i.e. the 'upper left hand corner').
;	No changes in the MIDAS table descriptors are carried out.  No 
;	changes in the MIDAS table 'select' are carried out.
;	Real data handled only.
;
;	Midas extensions (.tbl) in lower case.
;
; OUTPUTS:
;	None.
;
; ALGORITHM:
;	We first consider the File Control Block of the MIDAS table file
;	to determine the start of descriptor information and the start of
;	the data. 
;	Next we consider the Descriptor Directory Entry for `tblcontr' (the 
;	number of columns and rows allocated; followed by the number of 
;	columns and rows in the actual table).
;	We check consistency between input array dimensions and the dimensions
;	of the MIDAS table into which the array is to be written.
;	Finally we carry out the writing of the data values.
;
; EXAMPLE:
;	Given two MIDAS tables, PSF.TBL and PSF2.TBL, of the same dimensions;
;	we wish to overwrite PSF2 with the contents of PSF.
;
;	MID_RD_TABLE,'PSF',NC,NR,DAT
;	DATTR = TRANSPOSE(DAT)
;	MID_UP_TABLE,'PSF2',NC,NR,DATTR
;
; AUTHORS:
;	FM - Fionn Murtagh, ST-ECF, Munich.
;
;MODIFICATION HISTORY:
;	OCT 1988  FM  Initial programming and decoding of MIDAS files.
;	MAY 1989  FM  Name of routine changed; allowance made for LUN instead
;		of table name; bug fixes.
;	FEB 1991  FM  Conversion to V. 2, IDL, Unix.
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;----------------------------------------------------------------

;-- Some error checking and precautions.

n = n_params(0)
if n ne 4 then begin
   print,'Calling sequence: MID_UP_TABLE,table,ncol,nrow,data'
   print,'where inputs are TABLE - a string or a number - NCOL,NROW,DATA'
   return
endif

info_table = size(table)
if (info_table[0] ne 0) then begin
   print,'Calling sequence: MID_UP_TABLE,table,ncol,nrow,data'
   print,'TABLE must be a scalar'
   return
endif

;-- The regular shutdown procedure for this routing can also handle the
;-- shutdown required after an error.

on_ioerror, clean_up

;------------------------------
;-- STEP 1: OPEN THE TABLE FILE
;------------------------------

;-- Open the ".tbl" file (update).  If TABLE was a string then allocate
;-- a LUN and open the file.  Otherwise assume the file is open and use
;-- TABLE as the LUN.

if (info_table[1] eq 7) then begin
   get_lun,lun1
   openu,lun1,table+'.tbl'
endif else begin
   lun1 = table
endelse

;----------------------------------------------
;-- STEP 2: DETERMINE THE BEGINNING OF THE DATA
;----------------------------------------------

;-- Map the first block into 512-length byte array called
;-- "fcb", i.e. "frame control block".

fcb = assoc(lun1,bytarr(512))


;-- Define the byte string "h" from "fcb", i.e. from the first
;-- element onwards (remember, IDL reckons from location 0). 

h = fcb[0]

;-- Pick up the value of "mainseg4" as a long integer (i.e. I*4) from
;-- location 185 of "h".  "Mainseg4" is a pointer to the first block of 
;-- data: we know this from studying file MID_DISK:[MIDAS.NEW.INCL]
;-- FCBCOM.INC.  This latter file defines the structure of the frame
;-- control block.

mainseg4 = long(h,184)

;-- If we knew the number of data (nrow*ncol) values required, i.e. the

;-- dimensions of the table array, we would now be in a position to 
;-- give the command: data=assoc(lun1,fltarr(nvals),mainseg4-1)
;-- to access the table values. 
;-- What we do is to map the data as a floating array of "nvals"
;-- values, starting at the beginning of block "mainseg4 -1".  The 
;-- reason why we subtract 1 from "mainseg4" is that this is a (number

;-- of blocks) offset, - i.e. the number of blocks to skip.
;-- Before giving such a command, we must determine the table dimensions.


;-------------------------------------
;-- STEP 2: DETERMINE TABLE DIMENSIONS 
;-------------------------------------

;-- We assume that the descriptor information begins in block 2; map this.
;-- The following statements (commented out) are examples of how we
;-- generally pick up information relative to the any given descriptor.
;-- We have descriptor name (ch*15), type (C,I,R or D; ch*5), no. of 
;-- bytes per descr. elt. (i*2), no. elts. in descr. (i*2), start block 
;-- of descr. (i*4) and offset within start block (i*2).  Note that this 
;-- comprises 30 bytes of information.  The specification of this 30-byte 
;-- "descriptor directory entry" comes from file MID_DISK:[MIDAS.NEW.INCL]
;-- DSCDIR.INC.

;-- Assume first that the "descriptor directory entry" is in block 2.

desc = assoc(lun1, bytarr(512),512)
d = desc[0]

;-- Look for descriptor "tblcontr".  The "name" field of the "descriptor
;-- directory entry" will necessarily start at one of the elements 46, 76,
;-- 106, ... of this block.  The block is 512 bytes long.  If we have to
;-- continue looking in block 3, the byte elements at which "name" can
;-- start are 14, 44, ...  We make the assumption that in practice the
;-- "tblcontr" descriptor will always be found in blocks 1 or 2.

for i = 46,496,30 do begin
    name = string(d[i:i+14])
    if (name eq 'TBLCONTR       ') then goto, jump1
endfor
;-- Still not found? Look in block 3.
desc = assoc(lun1, bytarr(512),1024)
d = desc[0]
for i = 14,494,30 do begin
    name = string(d[i:i+14])
    if (name eq 'TBLCONTR       ') then goto, jump1
endfor
print,'ERROR: cannot find descriptor TBLCONTR.'
return

jump1: 

;-- Assume usual values for "type" (I), "bytelem" (4) and "noelem" (8).
;-- Get values of "start" and "index".

start = long(d,i+24)
index = fix(d,i+28)

;-- This offset or index (to where the actual data associated with a
;-- descriptor is to be found) is in long integer units.  Multiplying by
;-- 4 gives the byte units.  Furthermore, each set of actual descriptor
;-- data is prefaced with three 4-byte numbers.  These are: no. of
;-- elts., block no. of next elt., index of next extension (see 
;-- file MID_DISK:[MIDAS.NEW.LIB.ST]MIDRDLDB.FOR for this).  Hence
;-- we allow for these three numbers by adding 12 bytes to the offset.

index = index*4 + 12

;-- Given this value, how many virtual (i.e. 512-byte long) blocks
;-- must we skip (in addition to the first, "fcb", block)?

skip = fix(index/512)
start = start + skip - 1

;-- Offset within this (virtual or 512-byte long) block.

index = index - skip*512

;-- Now map this block into memory and determine "ncol" and "nrow".

start=start*512
vals = assoc(lun1,bytarr(512),start)
v = vals[0]
;-- Now skip 2 long integers corresponding to the allocated nos. of 
;-- cols., rows.
nc = long(v,index+8)
nr = long(v,index+12)
print,'Number of columns, rows in MIDAS table:  ',nc,nr

;-- Check for consistency in dimensions between the MIDAS table and
;-- the array to be written into it.

if (nr lt nrow) or (nc lt ncol) then begin
   print,'ERROR: dimensions of data array too large.'
   return
endif

;--------------------------------
;-- STEP 4: FINALLY, MAP THE DATA
;--------------------------------

;-- Transpose data array?
;-- data = transpose(data)

offst = mainseg4-1
offst=offst*512
dat = assoc(lun1,fltarr(nr,nc+1),offst)  ; Allow for select col.
datarr = dat[0]           ; Table data as an array.
datarr[0,1] = data        ; Write data, starting col. beyond select col.
dat[0] = datarr           ; Map back.


;-- Step 5: clean up

clean_up:

;-- If the TABLE input was a string we need to clean up the file action.

if (info_table[1] eq 7) then begin
   free_lun,lun1
endif

return
end
