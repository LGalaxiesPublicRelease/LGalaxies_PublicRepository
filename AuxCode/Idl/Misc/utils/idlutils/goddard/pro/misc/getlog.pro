Function GETLOG,lname
;+
; NAME:
;     GETLOG
; PURPOSE:
;     Formats a logical directory for the given operating system.
;
; CALLING SEQUENCE:
;     result = GETLOG(lname)
;
; INPUTS:
;     lname   - the base name of the logical (without special characters).
;
; OUTPUTS:
;     Returns appropriate string.
;     Under VMS the logical is not translated since it may correspond to
;     multiple directories.    
;
; RESTRICTIONS:
;     Assumes that the directory logical will have meaning to the host
;     operating system.
; PROCEDURE:
;       The operating system in !version.os_family is checked. If it equals:
;
;               'vms'           then a ':' is appended.
;               'windows'       directory name is translated with GETENV()
;                               and a '\' is appended
;               'unix'          directory name is translated with GETENV()
;                               and a '/' is appended
;               'MacOS'         directory name is translated with GETENV()
;
; EXAMPLE:
;       Open the file 'stars.dbh' in the logical directory ZDBASE in an 
;       operating system independent way:
;               IDL> openr,1,getlog('ZDBASE') + 'stars.dbh'
;
; MODIFICATION HISTORY:
;       Written, JDNeill, May, 1990.
;       Modified, JDNeill,Sep, 1990 -- for unix return full path instead of
;               just environment variable name.
;       Modified, I. Freedman, HSTX April 1994 -- for MacOS return full path
;       Bug in CASE statement fixed.    JDO, HSTX, May 2 1994.
;       Added Windows compatibility  W. Landsman      September 1995
;       Converted to IDL V5.0   W. Landsman   September 1997
;-
;-----------------------------------------------------------------------------
        CASE !VERSION.OS_FAMILY OF
  'vms':     return, lname+':' 
  'MacOS':   return, getenv(strupcase(lname)) 
  'windows': return, strtrim( getenv(strupcase(lname)),2) + '\'
   ELSE:     return, getenv(strupcase(lname))+'/'   ; assume UNIX 
  ENDCASE

;
        end     ; getlog.pro
