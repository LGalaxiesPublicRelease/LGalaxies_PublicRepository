;-----------------------------------------------------------------------
;+
; NAME:
;   djs_lockfile
;
; PURPOSE:
;   Test if a file is already "locked", and lock it if not.
;
; CALLING SEQUENCE:
;   res = djs_lockfile( filename, [ lun=, append= ] )
;
; INPUT:
;   filename:   File name
;
; OPTIONAL INPUTS:
;   lun:        If this argument exists, then open FILENAME for read/write
;               access and return the pointer (LUN number) for that file.
;               Do this only if we are able to lock the file.
;   append:     If set, then append to any file that already exists if
;               opening the file using LUN.  Ignored if the LUN argument
;               is not present.
;
; OUTPUTS:
;   res:        Return 0 if file already locked, or 1 if not in which case
;               we would have just locked it.
;
; COMMENTS:
;   For Unix systems running IDL 5.4 or later, we use the SPAWN command
;   to create a symbolic link from FILENAME.lock -> FILENAME.  This can
;   be done atomically, such that it is impossible for two processes
;   to build that same link at once.
;
;   For other operating systems or ealier versions of IDL (which do
;   not allow SPAWN to return an error), we use a lock file which
;   has a single byte written to it to indicate that FILENAME should
;   be locked (as determined by any subsequent calls to this routine).
;
;   For all OS-es, unlock files with DJS_UNLOCKFILE.
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   30-Apr-2000  Written by D. Schlegel, APO
;-
;-----------------------------------------------------------------------
function djs_lockfile, filename, lun=lun, append=append

   if (n_elements(filename) NE 1 OR size(filename,/tname) NE 'STRING') $
    then begin
      splog, 'FILENAME must be specified'
      return, 0
   endif

   lockfile = filename + '.lock'

   osfamily = !version.os_family
   if (osfamily EQ 'unix' AND !version.release LT '5.4') then $
    osfamily = 'other'

   case osfamily of
   'unix': begin
      spawn, '\ln -s ' + filename + ' ' + lockfile, $
       spawnres, spawnerr
      if (keyword_set(spawnerr)) then begin
         res = 0
      endif else begin
         res = 1
         if (arg_present(lun)) then $
          openw, lun, filename, /get_lun, append=append
      endelse
      end
   else: begin
      openw, olun, lockfile, /append, /get_lun
      if ((fstat(olun)).size EQ 0) then begin
         writeu, olun, 1B
         flush, olun ; Flush output immediately
         res = 1
         if (arg_present(lun)) then begin
            openw, lun, filename, /get_lun, append=append
         endif
      endif else begin
         res = 0
      endelse
      close, olun ; This will flush output to this file
      free_lun, olun
      end
   endcase

   return, res
end
;-----------------------------------------------------------------------
