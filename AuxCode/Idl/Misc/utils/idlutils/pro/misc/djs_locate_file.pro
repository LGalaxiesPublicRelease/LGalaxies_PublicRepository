;-----------------------------------------------------------------------
;+
; NAME:
;   djs_locate_file()
;
; PURPOSE:
;   Locate full file name (including path) of a file, searching IDL paths
;
; CALLING SEQUENCE:
;   fullname = djs_locate_file( filename )
;
; INPUT:
;   filename:   File name to find somewhere in path, including any extensions
;
; OUTPUTS:
;   fullname:   File name of first located file (including full path),
;               or '' if no matches found
;
; PROCEDURES CALLED:
;   os_family()
;
; REVISION HISTORY:
;   Written by D. Schlegel, 27 May 1997, Durham
;   Modified version of GETPRO in Goddard library.
;-
;-----------------------------------------------------------------------
function djs_locate_file, filename

   os = !VERSION.OS                               ;Operating system

   ; Need 1 parameter
   if N_params() LT 1 then begin
      print, 'Syntax - fullname = djs_locate_file( filename )'
      return, ''
   endif

   case os of
   'vms': begin
      dirsep = ''
      filename = strupcase(name)
      end
   'windows': dirsep = '\'
   else: dirsep = '/'
   endcase   

   path_dir = expand_path(!PATH, /array, count=N_dir)

   ; Loop over each directory in !PATH until procedure name found
   for idir = 0, N_dir-1 do begin
      dir = path_dir[idir]
      fullname = dir + dirsep + filename
      matches = findfile(fullname, count=count)
      if ( count GT 0 ) then begin
         return, fullname
      endif
   endfor

   return, ''
end 
;-----------------------------------------------------------------------
