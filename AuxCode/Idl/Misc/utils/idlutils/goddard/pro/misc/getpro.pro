pro getpro,proc_name            ;Obtain a copy of a procedure
;+
; NAME:
;     GETPRO
; PURPOSE:
;     Search !PATH for a procedure, and copy into user's working directory
; EXPLANATION:
;     Extract a procedure from an IDL Library or directory given in the 
;     !PATH  system variable and place it in the current default directory
;     (presumably to be edited by the user).  GETPRO can also be used to 
;     obtain a copy of the default startup file.
;
; CALLING SEQUENCE:
;     GETPRO, [ proc_name ]          ;Find PROC_NAME in !PATH and copy
;
; OPTIONAL INPUT:
;     proc_name - Character string giving the name of the IDL procedure or 
;               function.  Do not give an extension.   If omitted, 
;               the program will prompt for PROC_NAME.   
;
; OUTPUTS:
;     None.
;
; SIDE EFFECTS:
;      A file with the extension .pro and a name given by PROC_NAME will
;      be created on the user's directory.
;
; PROCEDURE:
;      The system variable !PATH is parsed into individual libraries or 
;      directories.   Each library or directory is then searched for the
;      procedure name.   When found, a SPAWN is used to extract or copy
;      the procedure into the user's directory.    If not found in !PATH,
;      then the ROUTINE_INFO() function is used to determine if it is an
;      intrinsic IDL procedure.  
;
; EXAMPLE:
;       Put a copy of the USER library procedure CURVEFIT on the current
;       directory
;
;       IDL> getpro, 'CURVEFIT'
;
; RESTRICTIONS:
;       User will be unable to obain source code for a native IDL function
;       or procedure, or for a FORTRAN or C routine added with CALL_EXTERNAL.
;       User must have write privilege to the current directory
;
;       This procedure is not used with Macintosh IDL. 
; PROCEDURE CALLS:
;       FDECOMP, ZPARCHECK
; REVISION HISTORY:
;      Written W. Landsman, STX Corp.   June 1990
;      Now use intrinsic EXPAND_PATH() command  W. Landsman November 1994
;      Use ROUTINE_NAMES() to check for intrinsic procs  W. Landsman July 95
;      Update for Windows/IDL     W. Landsman      September 95
;      Check if procedure is in current directory  W. Landsman  June 1997
;      Converted to IDL V5.0   W. Landsman   September 1997
;      Use ROUTINE_INFO instead of undocumented ROUTINE_NAMES W.L. October 1998
;-
  On_error,2                                     ;Return to caller on error
  os = !VERSION.OS                               ;Operating system

  IF os EQ 'MacOS' then message, $
   'This procedure not available for the Macintosh -- use the the FINDER menu'

  if N_params() EQ 0 then begin ;Prompt for procedure name?
        proc_name = ' ' 
        read,'Enter name of procedure you want a copy of: ',proc_name     
  
  endif else zparcheck, 'getpro', proc_name, 1, 7, 0, 'Procedure name'

  fdecomp, proc_name, disk, dir, name      ;Don't want file extensions
  name = strtrim( name, 2 )  
;                                    ;Make sure user has write privileges

  openw, lun, 'temp.tmp', /DELETE, /GET_LUN, ERROR=ERR
  if err NE 0 then begin               ;Problem writing a temporary file?
     cd,current=curdir   
     message,curdir + ' has insufficient privilege or file protection violation'
  endif
  free_lun,lun

;Set up separate copy commands for VMS and Unix

  case !VERSION.OS_FAMILY OF 

  'vms': begin
        dirsep = ''
        name = strupcase(name)
        end

  'Windows':dirsep = '\'
  else: dirsep = '/'

  ENDCASE   

  path_dir = expand_path(!PATH,/ARRAY, Count = N_dir)

;    Loop over each directory in !PATH until procedure name found

   for idir = 0, N_dir-1 do begin

     dir = path_dir[idir]

     if strmid(dir,0,1) EQ '@' then begin          ;VMS text Library?

        libname = strmid( dir,1,strlen(dir)-1 )         ;Remove the "@" symbol
        spawn,'library/extract='+name+'/out='+name+'.pro '+ $
                               libname,out,count=i
        if i EQ 0 then begin                           ;Success?
            message,name + '.PRO extracted from ' + libname,/INF
            return
        endif

   endif else begin                              ;Directory

        a = findfile(dir + dirsep + name +'.pro',COUNT=i)

        if I GE 1 then begin                     ;Found by FINDFILE?
        Case !VERSION.OS_FAMILY of
         'Windows': spawn, 'copy ' + dir + dirsep + name+ '.pro' + '*.*'
         'vms': spawn, 'copy ' + dir + dirsep + name+ '.pro' + ' *'
          else: spawn,'cp ' + a[0] + ' .'
        endcase
          message,'Procedure '+ NAME + ' copied from directory '+ dir,/INF
          return

        endif

    endelse
  endfor

; At this point !PATH has been searched and the procedure has not been found
; Check if the procedure is already in the current directory

  openr,lun, name + '.pro',/GET_LUN, ERROR = error
  if error EQ 0 then begin
        free_lun, lun
        message,name + '.pro already exists in the current directory',/INF
        return
  endif 


; Now check if it is an intrinsic IDL procedure or function.  

  funcnames = routine_info(/system,/func)
  test = where ( funcnames EQ strupcase(name), fcount)

  funcnames = routine_info(/system)
  test = where ( funcnames EQ strupcase(name), pcount)

  if (fcount EQ 0) and (pcount EQ 0) then begin   

     message,'Procedure '+NAME+' not found in the !PATH search string',/CONT
     message,'Check your spelling or search the individual directories',/INF

  endif else begin     

  if fcount GT 0 then $
       message,NAME + ' is an intrinsic IDL function',/CONT  $
  else message,NAME + ' is an intrinsic IDL procedure',/CONT
       message,'No source code is available',/INF

  endelse
  
  return
  
  end 
