pro dbcompress, dbname
;+
; NAME:
;       DBCOMPRESS
; PURPOSE:
;       Compress a .dbf database file after a call to DBDELETE
; EXPLANATION:
;       The procedure DBDELETE will remove specified entries from a database
;       but it will not free the unused space.     DBCOMPRESS will compress
;       the .dbf file so that it only contains valid entries.   
; CALLING SEQUENCE:
;       DBCOMPRESS, dbname
; INPUT PARAMETERS: 
;       dbname - Name of the database to be compressed, scalar string
; NOTES:
;       (1) Will not compress the index (.dbx) file.   The size of the .dbx file
;       is controlled by the MaxEntries value in the database definition 
;       (.dbd) file
;       (2) The updated .dbf file is written in the current directory.
;       This may need to be moved into the ZDBASE directory.
; PROCEDURE CALLS:
;       DBOPEN, DB_INFO(), FIND_WITH_DEF()
; REVISION HISTORY:
;       Written, W. Landsman      Raytheon STX        May 1998 
;       Converted to IDL V5.0 June 1998
;-
; Get the record length and number of entries
;
 if N_params() LT 1 then begin
      print,'Syntax - DBCOMPRESS, dbname'
      return
 endif
 dbopen,dbname
 len = db_info('length')
 len = len[0]
 N_entries = db_info('entries')
 N_entries = N_entries[0]

; Open the .dbf file directly (since DBRD won't let us read record 0)

 dbfname = find_with_def(dbname + '.dbf','ZDBASE')
 openr,lun1,dbfname,/GET_LUN
 a = assoc(lun1,bytarr(len))
 file_info = fstat(lun1)
 if file_info.size EQ (N_entries+1)*len then begin
	message,'No Compression needed for database ' + db_info('NAME',0),/INFO
	return
 end

;Open a temporary output file

 openw,lun2,'tmp.dbf',/GET_LUN
 b = assoc(lun2,bytarr(len))

;Copy all records up to N_entries (including record 0) into the output file.

 for i = 0L, N_entries do b[i] = a[i] 

 free_lun,lun1
 free_lun,lun2

; Make the temporary output file the new .dbf file

 case !VERSION.OS_FAMILY of
 'vms': spawn,'rename tmp.dbf ' + dbfname
 else: spawn,'mv tmp.dbf ' + dbfname
 endcase

 return
 end
