function spec_dir,filename,extension
;+
; NAME:
;     SPEC_DIR
; PURPOSE:
;     Complete a file specification by appending the default disk or directory
;
; CALLING SEQUENCE:                      
;     File_spec = SPEC_DIR( filename, [ extension ] )
; INPUT:
;     filename - character string giving partial specification of a file
;               name.  Examples for different operating systems include the
;                       following:
;               VMS: '$1$DUA5:TEST.DAT','[.SUB]TEST'
;               Unix: 'pro/test.dat', '$IDL_HOME/test','~/subpro'
;               MacOS: ':Programs:test'
;               Windows: '\pro\test.dat','d:\pro\test'
;
; OPTIONAL INPUT:
;     exten - string giving a default file name extension to be used if
;             filename does not contain one.  Do not include the period.
;
; OUTPUT:
;     File_spec - Complete file specification using default disk or 
;               directory when necessary.  
;
; EXAMPLE:
;      IDL> a = spec_dir('test','dat')
;
;      is equivalent to the commands
;      IDL> cd, current=cdir
;      IDL> a = cdir + delim + 'test.dat'
;
;      where delim is the OS-dependent separator 
; METHOD:
;      SPEC_DIR() decomposes the file name using FDECOMP, and appends the 
;      default directory (obtained from the CD command) if necessary.   
;      Under VMS, SPEC_DIR() will also try to translate disk and directory 
;      logical names.
;
;      SPEC_DIR() does not check whether the constructed file name actually
;      exists.
; PROCEDURES CALLED:
;      EXPAND_TILDE(), FDECOMP
; REVISION HISTORY:
;      Written W. Landsman         STX         July, 1987
;      Added Unix compatibility, W.  Landsman, STX   August 1991
;      Added Windows and Macintosh compatibility   W. Landsman  September, 1995
;      Work for relative Unix directory            W. Landsman  May, 1997
;      Expand Unix tilde if necessary              W. Landsman  September 1997
;      Converted to IDL V5.0   W. Landsman   September 1997
;      Fix VMS call to TRNLOG()  W. Landsman       September 2000
;-
 On_error,2                                     ;Return to user

 unix = !VERSION.OS_FAMILY EQ 'unix'
 filname = filename
 if unix then if strpos(filname,'~') GE 0 then filname = expand_tilde(filname) 
 fdecomp,filname,disk,dir,name,ext             ;Decompose filename

 if (ext EQ '') and ( N_params() GT 1) then $   ;Use supplied default extension?
                    ext = extension

 environ = (unix) and (strmid(dir,0,1) EQ '$')

 if not environ then begin
 if (unix) and (strmid(dir,0,1) NE '/')  then begin
     cd,current=curdir
     dir = curdir + '/' + dir
 endif

 if (dir EQ '') and (!VERSION.OS NE "vms") and (not environ) then begin

    cd,current=dir
    if name NE '' then begin
        case !VERSION.OS_FAMILY of 
        'windows': dir = dir + '\'    ;Get current default directory
        'MacOS': 
         else: dir = dir + '/'
        endcase
    endif
 
 endif else begin

   if ( disk EQ '' ) or ( dir EQ '' ) then begin
     cd,current=defdir                          ;Get current default directory
     fdecomp,defdir,curdisk,curdir
     if disk EQ '' then disk = curdisk else begin
       if !VERSION.OS EQ "vms" then begin
         logname = strmid(disk,0,strpos(disk,':'))
         test = trnlog(logname,fname)
         if test then begin
            if strmid(fname,strlen(fname)-1,1) EQ ']' then begin
               if strmid(fname,strlen(fname)-2,1) NE '.' then $
                                  return,spec_dir(fname+name,ext)
            endif else return,spec_dir(fname+':'+dir+name,ext)
         endif
       endif
     endelse

    if dir eq '' then dir = curdir else if !VERSION.OS EQ 'vms' then begin
        if strpos(dir,'.') eq 1 then dir = $
           strmid(curdir,0,strlen(curdir)-1) + strmid(dir,1,strlen(dir)-1)
     endif

 endif
 
 endelse
 endif

 if ext ne '' then ext = '.'+ext

 if !VERSION.OS ne "vms" then return,dir+name+ext else  $             ;Unix
                               return,strupcase(disk+dir+name+ext)   ;VMS

end
