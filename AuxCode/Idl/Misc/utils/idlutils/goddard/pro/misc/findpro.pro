pro FindPro, Proc_Name, NoPrint=NoPrint, DirList=DirList, ProList=ProList
;+
; NAME:
;     FINDPRO
; PURPOSE:
;     Find all locations of a procedure in the IDL !PATH
; EXPLANATION:
;     FINDPRO searces for the procedure name (as a .pro or a .sav file) in all 
;     IDL libraries or directories given in the !PATH system variable.  
;               
; CALLING SEQUENCE:
;    FINDPRO, [ Proc_Name, /NoPrint, DirList = , ProList = ]
;
; OPTIONAL INPUT:
;     Proc_Name - Character string giving the name of the IDL procedure or 
;             function. Do not include the ".pro" extension. If Proc_Name is
;             omitted, the program will prompt for PROC_NAME.  "*" wildcards
;             are permitted.
;
; OPTIONAL KEYWORD INPUT:
;     /NoPrint - if set, then the file's path is not printed on the screen and
;             absolutely no error messages are printed on the screen.  If not
;             set, then - since the MESSAGE routine is used - error messages 
;             will be printed but the printing of informational messages
;             depends on the value of the !Quiet variable.
;
; OPTIONAL KEYWORD OUTPUTS:
;     DirList - The directories in which the file is located are returned in
;             the keyword as a string array.
;             If the procedure was found in a VMS text library, then the
;             full path and name of that library is returned and is prefixed
;             by an "@" sign as a flag that it is a library, not a directory.
;             If the procedure is an intrinsic IDL procedure, then the 
;             value of DirList = ['INTRINSIC'].
;             If the procedure is not found, the value of DirList = [''].
;     ProList - The list (full pathnames) of procedures found.  Useful if you
;             are looking for the name of a procedure using wildcards.
;
;     The order of the names in DirList and ProList is identical to the order
;     in which the procedure name appears in the !PATH
; PROCEDURE:
;     The system variable !PATH is parsed using EXPAND_PATH into individual 
;     libraries or directories.   Each library or directory is then 
;     searched for the procedure name.  If not found in !PATH, then the 
;     the name is compared with the list of intrinsic IDL procedures given
;     by the ROUTINE_INFO function. 
;
; EXAMPLE:
;     (1) Find the procedure CURVEFIT.  Assume for this example that the user
;     also has a copy of the CURVEFIT.PRO procedure in her home directory
;     on a Unix machine.
;
;       IDL> findpro, 'curvefit', DIRLIST=DirList
;       Procedure curvefit.pro found in directory  .
;       Procedure curvefit.pro found in directory  /home/idl/lib/userlib 
;       IDL> help, DirList
;       DIRLIST         STRING    = Array(2) 
;       IDL> help, DirList(0), DirList(1)
;       <Expression>    STRING    = '.'
;       <Expression>    STRING    = '/home/idl/lib/userlib' 
;
;     (2) Find all procedures in one's !path containing the characters "zoom" 
;
;       IDL> findpro,'*zoom*'
; RESTRICTIONS:
;       User will be unable to find a path for a native IDL function
;       or procedure, or for a FORTRAN or C routine added with CALL_EXTERNAL.
;       Remember that Unix is case sensitive, and most procedures will be in
;       lower case.
;
; PROCEDURES USED:
;       ZPARCHECK, FDECOMP, UNIQ()
; REVISION HISTORY:
;       Based on code extracted from the GETPRO procedure, J. Parker 1994
;       Use the intrinsic EXPAND_PATH function    W. Landsman Nov. 1994
;       Use ROUTINE_NAMES() to check for intrinsic procs   W. Landsman Jul 95
;       Added Macintosh, WINDOWS compatibility    W. Landsman   Sep. 95
;       Removed spurious first element in PROLIST  W. Landsman  March 1997
;       Don't include duplicate directories  in !PATH  WL   May 1997
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Use ROUTINE_INFO instead of undocumented ROUTINE_NAMES W.L. October 1998
;       Also check for save sets   W. Landsman  October 1999 
;       Force lower case check for VMS  W. Landsman January 2000 
;       Only return .pro or .sav files in PROLIST   W. Landsman  January 2002 
;       Force lower case check for .pro and .sav    D. Swain  September 2002 
;
;-
;/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

 On_error,2                           ;Return to caller on error
 OS = !VERSION.OS                     ;VMS or Unix operating system

 if (N_params() EQ 0) then begin      ;Prompt for procedure name?
   Proc_Name = ' ' 
   read,'Enter name of procedure for which you want the path: ',Proc_Name
 endif else zparcheck, 'getpro', Proc_Name, 1, 7, 0, 'Procedure name'

 NoPrint = keyword_set(NoPrint)
 DirList = strarr(1)
 ProList = strarr(1)

 fdecomp, Proc_Name, Disk, Dir, Name      ;Don't want file extensions
 Name = strtrim( Name, 2 )  

; Set up separate file and directory separators for VMS and Unix

 case !VERSION.OS_FAMILY of
 
 'vms':  begin   
           DirSep  = ''
           Name    = strupcase(Name)
           lname = name
           remchar,lname,'*'
        end
 'Windows': DirSep = '\'
 'MacOS': DirSep = ''                 ;Fixed 21-Sep-1995
 else: DirSep = '/'
 
 endcase

 pathdir = expand_path(!PATH,/ARRAY, Count = N_dir)
 cd,current = dir
; Remove duplicate directories  in !PATH but keep original order
 path_dir = [dir]
 for i = 0,N_dir -1 do begin
      test = where(path_dir EQ pathdir[i], Ndup)
      if Ndup EQ 0 then path_dir = [path_dir,pathdir[i]]
 endfor
 N_dir = N_elements(path_dir)

; Loop over each directory in !PATH until procedure name found

   for idir = 0, N_dir-1 do begin

   dir = path_dir[idir]
   if (strmid(Dir,0,1) eq '@') then begin          ;Text Library?

         LibName = strmid( Dir, 1, strlen(Dir)-1 )      ;Remove the "@" symbol
         spawn, 'library /list ' + LibName, List
         if strpos(name,'*') GE 0 then begin
                lfound = where( strpos(list, lname) GE 0, Nfound)
         endif else lfound = where( list EQ name, Nfound)

         if (Nfound GT 0) then begin
           for j=0,Nfound-1 do begin
            DirList = [DirList, Dir]
            ProList = [ProList, name]
            Mess = list[lfound[j]] + ' found in the library  ' + LibName
            message, Mess, /CONT, NOPRINT=NoPrint, /NOPREFIX, /NONAME
          endfor
         endif

   endif else begin                              ;Directory
     found = 0b
      ProsFound = findfile(Dir+DirSep+Name+'.*', COUNT=Nfile)

      if (Nfile ge 1) then begin                     ;Found by FINDFILE?
       for j = 0,nfile-1 do begin
         fdecomp, ProsFound[j], ddisk,ddir,fname,ext
         if ((strlowcase(ext) EQ 'pro') or (strlowcase(ext) EQ 'sav')) then begin
         	found = 1b
         	profound = prosfound[j]
         	if strlowcase(ext) EQ 'pro' then $  
         	  message,/Con,NoPrint = NoPrint,/NoPrefix, /Noname, $
         	 'Procedure ' + fname + ' found in directory  ' + disk + Dir $
         	else if strlowcase(ext) EQ 'sav' then $ 
     				message,/Con,NoPrint = NoPrint,/NoPrefix, /Noname, $
         		'Save set ' + fname + '.sav found in directory  ' + disk + Dir
         	endif
        endfor
         if found then begin
         DirList = [DirList, Dir]
         ProList = [ProList, ProFound]
         endif
       endif
   endelse

endfor

 if (N_elements(DirList) GT 1) then begin
        DirList = DirList[1:*]
        ProList = ProList[1:*]
 endif
      

; At this point !PATH has been searched.  If the procedure has not been found
; (nothing is in DirList) check if it is an intrinsic IDL procedure or function

 if (DirList[0] eq '') then begin

  funcnames = routine_info(/system,/func)
  test = where ( funcnames EQ strupcase(name), fcount)

  funcnames = routine_info(/system)
  test = where ( funcnames EQ strupcase(name), pcount)

   if (fcount EQ 0) and (pcount EQ 0) then begin
      if not(NoPrint) then begin
         message, 'Procedure '+Name+' not found in a !PATH directory.', /CONT
         message, 'Check your spelling or search individual directories.', /INF
      endif
   endif else begin 
      DirList = ['INTRINSIC']
      ProList = ['INTRINSIC']
      if not(NoPrint) then begin
         message, 'Procedure ' + Name + ' is an intrinsic IDL procedure.', /CONT
         message, 'No path information available.', /INF
      endif
   endelse

 endif
  
 return
 end   
