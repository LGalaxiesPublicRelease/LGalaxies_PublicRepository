pro dbcreate,name,newindex,newdb,maxitems,EXTERNAL=EXTERNAL
;+
; NAME: 
;       DBCREATE
; PURPOSE:      
;       Create a new data base (.dbf), index (.dbx) or description (.dbh) file
; EXPLANATION:
;       A database definition (.dbd) file must already exist.
;       The default directory must be a ZDBASE: directory.
;
; CALLING SEQUENCE:     
;       dbcreate, name,[ newindex, newdb, maxitems]  [,/EXTERNAL]  
;
; INPUTS:       
;       name- name of the data base (with no qualifier), scalar string. 
;               The description will be read from the file "NAME".dbd 
;
; OPTIONAL INPUTS:      
;       newindex - if non-zero then a new index file is created,
;               otherwise it is assumed that changes do not affect the
;               index file. (default=0)
;       newdb - if non-zero then a new data base file (.dbf) will
;               be created. Otherwise changes are assumed not to affect
;               the file's present format.
;       maxitems - maximum number of items in data base.
;               If not supplied then the number of items is
;               limited to 200.
;
; OUTPUTS:
;       NONE.
;
; OPTIONAL INPUT KEYWORD:       
;
;       external - If set, then the database is written with an external data
;               representation.  This allows the database files to be used on
;               any computer platform, e.g. through NFS mounts, but some
;               overhead is added to reading the files.  The default is to
;               write the data in the native format of the computer being used.
;
;               This keyword is only paid attention to if NEWDB or NEWINDEX
;               are nonzero.  Otherwise, the database is opened to find
;               out if it uses external representation or not.
;
;               Extreme caution should be used if this keyword is used with
;               only NEWINDEX set to a nonzero value.  This mode is allowed so
;               that databases written on machines which already use the
;               external data representation format, e.g. Sun workstations, to
;               be marked external so that other machines can read them.
;
; PROCEDURE CALLS:      
;       GETTOK(), FIND_WITH_DEF(), HOST_TO_IEEE, ZPARCHECK
;
; RESTRICTIONS: 
;       If newdb=0 is not specified, the changes to the .dbd file can
;       not alter the length of the records in the data base file.
;       and may not alter positions of current fields in the file.
;       permissible changes are:
;               1) utilization of spares to create a item or field
;               2) change in field name(s)
;               3) respecification of index items
;               4) changes in default print formats
;               5) change in data base title
;               6) changes in pointer specification to other data
;                       data bases
;
;       !priv must be 2 or greater to execute this routine.
;
; SIDE EFFECTS:  
;       data base description file ZDBASE:name.dbh is created
;       and optionally ZDBASE:name.dbf (data file) and
;       ZDBASE.dbx (index file) if it is a new data base.
;
; REVISION HISTORY:     
;       D. Lindler, GSFC/HRS, October 1987
;       Modified:  Version 1, William Thompson, GSFC, 29 March 1994
;                  Version 2, William Thompson, GSFC/CDS (ARC), 28 May 1994
;                  Added EXTERNAL keyword.
;       Version 4, William Thompson, GSFC, 3 November 1994
;                       Modified to allow ZDBASE to be a path string.
;       8/14/95  JKF/ACC - allow EXTERNAL data for newindex OR newdb modes.
;       Make sure all databases closed before starting W. Landsman June 1997
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Added new unsigned and 64 bit integer datatypes W. Landsman July 2001
;-
;----------------------------------------------------------
On_error,2                         ;Return to caller

if N_Params() LT 1 then begin
      print,'Syntax - dbcreate, name, [ newindex, newdb, maxitems ]'
      print,'         !PRIV must be 2 or greater to execute this routine'
      return
endif
;
; check privilege
;
if !priv LT 2 then  $
        message,'!PRIV must be 2 or greater to execute this routine'
;
; check parameters
;
zparcheck, 'DBCREATE', name, 1, 7, 0, 'Database Name'
if N_params() LT 2 then newindex = 0
if N_params() LT 3 then newdb = 0
if N_params() LT 4 then maxitems = 200
filename = strtrim(name,2)

 dbclose                         ;Close any databases already open
;
; Decide whether or not external data representation should be used.
;   8/14/95  JKF/ACC - allow EXTERNAL data for newindex OR newdb modes.
;
if ((newindex ne 0) or (newdb ne 0)) or $
                ((findfile(filename+'.dbh'))[0] eq '') then begin
        extern = keyword_set(external)
end else begin
        openr,tempunit,filename+'.dbh',/block,/get_lun
        point_lun,tempunit,119
        extern = 0b
        readu,tempunit,extern
        free_lun,tempunit
endelse
;
; set up data buffers
;
names = strarr(maxitems)                        ;names of items
numvals = intarr(maxitems)+1                    ;number of values
type = intarr(maxitems)                         ;data type
nbytes = intarr(maxitems)                       ;number of bytes in item
desc = strarr(maxitems)                         ;descriptions of items
sbyte = intarr(maxitems)                        ;starting byte position
format = strarr(maxitems)                       ;print formats
headers = strarr(3,maxitems)                    ;print headers
headers[*,*]='               '                  ;init headers
title = ''                                      ;data base title
index = intarr(maxitems)                        ;index type
pointers = strarr(maxitems)                     ;pointer array
npointers = 0
maxentries = 30000L
alloc = 100L
;
; first item is always entry number
;
names[0] = 'ENTRY'
type[0] = 3             ;longword integer
nbytes[0] = 4           ;four bytes
desc[0] = 'Entry or Record Number'
format[0] = 'I8'
headers[1,0] = 'ENTRY'
nitems = 1
nextbyte = 4            ;next byte position in record
;
; open .dbd file
;
get_lun, unit                   ;get free unit number
On_ioerror, BAD_IO              ;On I/O errors go to BAD_IO
openr, unit, find_with_def(filename+'.dbd', 'ZDBASE')
;
; read and process input data
;
block='TITLE'                           ;assume first block is title
inputst=''
while not eof(unit) do begin            ;loop on records in the file
;
; process next line of input
;
    readf,unit,inputst
    print,inputst
    st=gettok(inputst,';')
    if strtrim(st,2) eq '' then goto,next       ;skip blank lines
    if strmid(st,0,1) eq '#' then begin
        block=strupcase(strmid(st,1,strlen(st)-1));begin new block
        goto,next
    end
;
    case strtrim(block,2) of

        'TITLE' : title=st

        'MAXENTRIES' : maxentries=long(strtrim(st,2))

        'ITEMS' : begin
;
;               process statement in form
;                       <itemname> <datatype> <description>
;
                item_name=" "
                item_name=strupcase(gettok(st,' '))
                st = strtrim(st, 1)
                item_type = " "
                item_type=gettok(st,' ')
                st = strtrim(st, 1)
                desc[nitems]=st
                if item_name eq '' then $
                        message,'Invalid item name',/IOERROR
                names[nitems]=gettok(item_name,'(')
                if item_name ne '' then $               ;is it a vector
                        numvals[nitems]=fix(gettok(item_name,')')) 
                if item_type eq '' then $
                  message,'Item data type not supplied for item ' + $
                          strupcase(item_name),/IOERROR
                data_type=strmid(strupcase(gettok(item_type,'*')),0,1)
                num_bytes=item_type
                if num_bytes eq '' then num_bytes='4'
                if (data_type eq 'R') or (data_type eq 'I') or $
                   (data_type eq 'U') then $
                                data_type=data_type+num_bytes
                case data_type of
                        'B' : begin & idltype= 1 & nb=1 & ff='I6' & end
                        'L' : begin & idltype= 1 & nb=1 & ff='I6' & end
                        'I2': begin & idltype= 2 & nb=2 & ff='I7' & end
                        'I4': begin & idltype= 3 & nb=4 & ff='I11' & end
                        'I8': begin & idltype= 14 & nb=8 & ff='I22' & end
                        'R4': begin & idltype= 4 & nb=4 & ff='G12.6' & end
                        'R8': begin & idltype= 5 & nb=8 & ff='G20.12' & end
                        'U2': begin & idltype= 12 & nb=2 & ff='I7' & end
                        'U4': begin & idltype= 13 & nb=4 & ff='I11' & end
                        'U8': begin & idltype= 15 & nb=8 & ff='I22' & end
                        'C' : begin
                                idltype = 7
                                nb=fix(num_bytes)
                                ff='A'+num_bytes
                              end
                        else: message,'Invalid data type "'+ item_type+ $
                                       '" specified',/IOERROR
                endcase
                format[nitems]=ff                       ;default print format
                headers[1,nitems]=names[nitems] ;default print header
                type[nitems]=idltype            ;idl data type for item
                nbytes[nitems]=nb               ;number of bytes for item
                sbyte[nitems]=nextbyte          ;position in record for item
                nextbyte=nextbyte+nb*numvals[nitems] ;next byte position
                nitems=nitems+1
                end

        'FORMATS': begin
;
;                process strings in form:
;                       <item name> <format> <header1>,<header2>,<header3>
;
                item_name=" "
                item_name=strupcase(gettok(st,' '))
                item_no=0
                while item_no lt nitems do begin
                        if strtrim(names[item_no]) eq item_name then begin
                                st = strtrim(st, 1)
                                format[item_no]=gettok(st,' ')
                                if strtrim(st,2) ne '' then begin
                                        st = strtrim(st, 1)
                                        headers[0,item_no]=gettok(st,',')
                                        headers[1,item_no]=gettok(st,',')
                                        headers[2,item_no]=strtrim(st)
                                endif
                        endif
                        item_no=item_no+1
                endwhile
                end

        'POINTERS': begin
;
;               process record in form:
;                       <item name> <data base name>
;
                item_name=strupcase(gettok(st,' '))
                item_no=0
                while item_no lt nitems do begin
                        if strtrim(names[item_no]) eq item_name then $
                                pointers[item_no]=strupcase(strtrim(st, 1))
                        item_no=item_no+1
                endwhile
                endcase

        'INDEX': begin
;
;               process record of type:
;               <item name> <index type>
;
                item_name=strupcase(gettok(st,' '))
                st = strtrim(st, 1)
                indextype=gettok(st,' ')
                item_no=0
                while item_no lt nitems do begin
                        if strtrim(names[item_no]) eq item_name then begin
                            case strupcase(indextype) of
                                'INDEX' : index[item_no]=1
                                'SORTED': index[item_no]=2
                                'SORT'  : index[item_no]=3
                                'SORT/INDEX' : index[item_no]=4
                                else    : message,'Invalid index type',/IOERROR
                            endcase
                        endif
                        item_no=item_no+1
                endwhile
                end
        else : begin
                print,'DBCREATE-- invalid block specfication of ',block
                print,'   Valid values are #TITLE, #ITEMS, #FORMATS, #INDEX,'
                print,'   #MAXENTRIES or #POINTERS'
               end
        endcase
next:
endwhile; loop on records
;
; create data base descriptor record --------------------------------------
;
;       byte array of 120 values
;
;       bytes
;         0-18   data base name character*19
;         19-79  data base title character*61
;         80-81  number of items (integer*2)
;         82-83  record length of DBF file (integer*2)
;         84-118 values filled in by DBOPEN
;         119    equals 1 if keyword EXTERNAL is true.
;
totbytes=((nextbyte+3)/4*4)  ;make record length a multiple of 4
drec = bytarr(120)
drec[0:79]=32b                      ;blanks
drec[0] = byte(strupcase(filename))
drec[19] = byte(title)
drec[80] = byte(nitems,0,2)
drec[82] = byte(totbytes,0,2)
drec[119] = byte(extern)
;
; create item description records
;
;  irec(*,i) contains decription of item number i with following
;  byte assignments:
;       0-19    item name (character*20)
;       20-21   IDL data type (integet*2)
;       22-23   Number of values for item (1 for scalar) (integer*2)
;       24-25   Starting byte position i record (integer*2)
;       26-27   Number of bytes per data value (integer*2)
;       28      Index type
;       29-97   Item description
;       98-99   Field length of the print format
;       100     Pointer flag
;       101-119 Data base this item points to
;       120-125 Print format
;       126-170 Print headers
;       171-199 Added by DBOPEN
irec=bytarr(200,nitems)
rec=bytarr(200)
headers = strmid(headers,0,15)       ;Added 15-Sep-92
for i=0,nitems-1 do begin
        rec[0:19]=32b  &  rec[101:170]=32b    ;Default string values are blanks
        rec[29:87] = 32b
        rec[0]  = byte(names[i])
        rec[20] = byte(type[i],0,2)
        rec[22] = byte(numvals[i],0,2)
        rec[24] = byte(sbyte[i],0,2)
        rec[26] = byte(nbytes[i],0,2)
        rec[28] = index[i]
        rec[29] = byte(desc[i])
        if strtrim(pointers[i]) ne '' then rec[100]=1 else rec[100]=0
        rec[101]= byte(strupcase(pointers[i]))
        rec[120]= byte(format[i])
        ff=strtrim(format[i])
        flen=fix(gettok(strmid(ff,1,strlen(ff)-1),'.'))
        rec[98] = byte(flen,0,2)
        rec[126]= byte(headers[0,i]) > 32b    ;Modified Nov-91
        rec[141]= byte(headers[1,i]) > 32b
        rec[156]= byte(headers[2,i]) > 32b
        irec[0,i]=rec
end
;
; Make sure user is on ZDBASE and write description file
;
 close,unit
 if extern then $
        openw,unit,filename+'.dbh',/block else  $
        openw,unit,filename+'.dbh',/segmented
   
On_ioerror, NULL 
if extern then begin
        tmp = fix(drec,80,1) & host_to_ieee,tmp & drec[80] = byte(tmp,0,2)
        tmp = fix(drec,82,1) & host_to_ieee,tmp & drec[82] = byte(tmp,0,2)
;
        tmp = fix(irec[20:27,*],0,4,nitems)
        host_to_ieee,tmp
        irec[20,0] = byte(tmp,0,8,nitems)
;
        tmp = fix(irec[98:99,*],0,1,nitems)
        host_to_ieee,tmp
        irec[98,0] = byte(tmp,0,2,nitems)
;
        tmp = fix(irec[171:178,*],0,4,nitems)
        host_to_ieee,tmp
        irec[171,0] = byte(tmp,0,8,nitems)
endif
writeu, unit, drec
writeu, unit, irec
;
; if new data base create .dbf and .dbx files -----------------------------
;

if newdb then begin
    close,unit
    if !VERSION.OS EQ "vms" then $
         openw, unit, filename+'.dbf', totbytes, /NONE   $
    else openw, unit, filename+'.dbf'
    header = bytarr(totbytes)
    p = assoc(unit,header)
    p[0] = header
end

;
; determine if any indexed items
;
nindex = total(index GT 0)
;
; create empty index file if needed
;
if (nindex GT 0) and (newindex) then begin
        indexed = where(index GT 0)
;
; create header array
;       header=intarr(7,nindex)
;               header(i,*) contains values
;               i=0     item number
;               i=1     index type
;               i=2     idl data type for the item
;               i=3     starting block for header
;               i=4     starting block for data
;               i=5     starting block for indices (type 3)
;               i=6     starting block for unsorted data (type 4)
;
        nb = (maxentries+511)/512       ;number of 512 value groups
        nextblock = 1
        header = lonarr(7,nindex)
        for ii = 0, nindex-1 do begin
                item = indexed[ii]
                header[0,ii] = item
                header[1,ii] = index[item]
                header[2,ii] = type[item]
                data_blocks = nbytes[item]*nb
                if index[item] NE 1 $
                             then header_blocks = (nbytes[item]*nb+511)/512 $
                             else header_blocks = 0
                if (index[item] eq 3) or (index[item] EQ 4) then $
                                 index_blocks=(4*nb) else index_blocks=0
                if index[item] EQ 4 then unsort_blocks = data_blocks else $
                                                        unsort_blocks=0
                header[3,ii] = nextblock
                header[4,ii] = nextblock+header_blocks
                header[5,ii] = header[4,ii]+data_blocks
                header[6,ii] = header[5,ii]+index_blocks
                nextblock = header[6,ii]+unsort_blocks
        end
        totblocks = nextblock
        close, unit
        openw, unit, filename+'.dbx', 512, /BLOCK  
;
        p = assoc(unit,lonarr(2))
        tmp = [long(nindex),maxentries]
        if extern then host_to_ieee, tmp
        p[0] = tmp
;
        p = assoc(unit,lonarr(7,nindex),8)
        tmp = header
        if extern then host_to_ieee, tmp
        p[0] = tmp
endif
free_lun, unit
return
;
BAD_IO: free_lun,unit
print, !MSG_PREFIX+!ERR_STRING
print, !SYSERR_STRING
return
;
end
