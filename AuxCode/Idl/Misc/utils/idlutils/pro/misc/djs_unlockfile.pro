;-----------------------------------------------------------------------
;+
; NAME:
;   djs_unlockfile
;
; PURPOSE:
;   Unlock a file if locked with DJS_LOCKFILE().
;
; CALLING SEQUENCE:
;   djs_unlockfile, filename, [lun= ]
;
; INPUT:
;   filename:   File name
;
; OPTIONAL INPUTS:
;   lun:        If this argument exists, then close the file associated
;               with this LUN number.  This is useful if FILENAME has
;               been opened with DJS_LOCKFILE().
;
; OUTPUTS:
;
; COMMENTS:
;   We use a lock file, which is either a symbolic link or a file with
;   a single byte written to it, to indicate that FILENAME has been
;   locked by DJS_LOCKFILE().  This routine deletes that file.
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   30-Apr-2000  Written by D. Schlegel, APO
;-
;-----------------------------------------------------------------------
pro djs_unlockfile, filename, lun=lun

   ; Close the file before unlocking it
   if (arg_present(lun)) then begin
      close, lun
      free_lun, lun
   endif

   lockfile = filename + '.lock'

   osfamily = !version.os_family
   if (osfamily EQ 'unix' AND !version.release LT '5.4') then $
    osfamily = 'other'

   case osfamily of
   'unix': begin
      ; In this case, we've created a symlink, which must be
      ; deleted with the Unix command "rm".  Otherwise, using
      ; the IDL command OPENW,/DELETE would delete the lock
      ; file but also set the contents of the actual file to null!
      spawn, '\rm ' + lockfile, $
       spawnres, spawnerr
      end
   else: begin
      openw, olun, lockfile, /get_lun, /delete
      close, olun
      free_lun, olun
      end
   endcase

   return
end
;-----------------------------------------------------------------------
