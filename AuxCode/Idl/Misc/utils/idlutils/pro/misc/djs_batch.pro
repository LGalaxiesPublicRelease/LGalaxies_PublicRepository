;+
; NAME:
;   djs_batch
;
; PURPOSE:
;   Batch processing script for running jobs locally or across a network.
;
; CALLING SEQUENCE:
;   djs_batch, topdir, [localfile], [outfile], protocol, remotehost, $
;    remotedir, command, [ priority=, selecthost=, wtime= ]
;
; INPUTS:
;   topdir     - Local top-level directory for input and output files.
;                Also use this directory for remote hosts where REMOTEDIR
;                is not specified.
;   localfile  - Array of pointers to input files on local machine [NPROGRAM].
;                This input is optional.
;   outfile    - Array of pointers to output files created on remote machine
;                and copied to local machine upon completion [NPROGRAM]
;                This input is optional.
;   protocol   - List of protocols for remote hosts.  Valid values are:
;                'ssh', 'ssh1', 'ssh2', 'rsh', or ''.  One must set to
;                no protocol ('') if the remote host name is 'localhost'.
;                Otherwise, one must always set a protocol.
;   remotehost - List of remote hosts [NHOST]
;   remotedir  - List of remote directories; scalar or [NHOST]
;   command    - Command to execute to begin a job; scalar or [NPROGAM]
;
; OPTIONAL KEYWORDS:
;   priority   - Priority for each job, where the jobs with the largest
;                value are done first [NPROGRAM]
;   selecthost - If set, then assign each job to only a host that matches
;                the selected host per job [NPROGRAM]
;   wtime      - Sleep time between checking status of all jobs; default to
;                600 seconds.
;
; OUTPUTS:
;
; COMMENTS:
;   The file names will support wildcards.  For example, if you want to
;   return all files from the directory REMOTEDIR/abc on the remote machine,
;   then set OUTFILE = 'abc/*'.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;   batch_spawn
;   count_freelun()
;   create_program_list()
;   create_host_list()
;   batch_if_done()
;   batch_assign_job
;   batch_finish_job
;
; REVISION HISTORY:
;   17-Oct-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro batch_spawn, command, retval, pid=pid, unit=unit

   splog, command
   if (arg_present(retval)) then $
    spawn, command, retval, pid=pid $
   else if (arg_present(unit)) then $
    spawn, command, pid=pid, unit=unit $
   else $
    spawn, command, pid=pid

   return
end

;------------------------------------------------------------------------------
function create_program_list, localfile, outfile, command, $
 priority=priority, selecthost=selecthost

   nprog = n_elements(command)

   ftemp = create_struct( name='PROGLIST_STRUCT', $
    'PROGNAME', '', $
    'LOCALFILE', ptr_new(), $
    'OUTFILE', ptr_new(), $
    'COMMAND', '', $
    'PID', 0L, $
    'UNIT', 0L, $
    'PRIORITY', 1.d0, $
    'SELECTHOST', '', $
    'STATUS', 'UNASSIGNED' )

   proglist = replicate(ftemp, nprog)

   for iprog=0, nprog-1 do begin
      if (keyword_set(localfile)) then $
       proglist[iprog].localfile = localfile[iprog] ; (pointer)
      if (keyword_set(outfile)) then $
       proglist[iprog].outfile = outfile[iprog] ; (pointer)
      proglist[iprog].progname = 'JOB#' + strtrim(string(iprog),2)
   endfor

   proglist.command = command
   if (keyword_set(priority)) then proglist.priority = priority
   if (keyword_set(selecthost)) then proglist.selecthost = selecthost

   return, proglist
end

;------------------------------------------------------------------------------
function create_host_list, protocol, remotehost, remotedir, topdir

   ; Do not try to use more hosts than there are free file pointers (LUNs).
   ; This is because the SPAWN command uses GET_LUN to open pipes,
   ; and will crash if it tries to allocate too many pipes.
   ; Specifically, it will crash in the SPAWN command in the BATCH_SPAWN
   ; procedure.
   nhost = n_elements(remotehost)
   nlunmax = count_freelun()
   if (nhost GT nlunmax) then begin
      splog, 'WARNING: Trimming to use only the first ', nlunmax, ' hosts'
      nhost = nlunmax
   endif

   if (nhost EQ 0) then begin
      splog, 'REMOTEHOST not set.  Quitting.'
      return, ''
   endif
   if (total(strtrim(remotehost) EQ '') GT 0) then begin
      splog, 'Some entries for REMOTEHOST are empty strings.  Quitting'
      return, ''
   endif

   ftemp = create_struct( name='HOSTLIST_STRUCT', $
    'PROGNAME', '', $
    'REMOTEHOST', '', $
    'REMOTEDIR', '', $
    'PROTOCOL', '', $
    'CPSTRING', '', $
    'STATUS', 'IDLE' )

   hostlist = replicate(ftemp, nhost)

   hostlist[*].remotehost = remotehost[0:nhost-1]
   hostlist[*].remotedir = remotedir[0:nhost-1]
   hostlist[*].protocol = protocol[0:nhost-1]

   for ihost=0, nhost-1 do begin
      if (hostlist[ihost].remotehost EQ 'localhost') then begin
         if (keyword_set(hostlist[ihost].protocol)) then $
          message, 'It is invalid to specify a protocol for host=localhost'
      endif else begin
         if (NOT keyword_set(hostlist[ihost].protocol)) then $
          message, 'Protocol needed for remote host ' $
           + hostlist[ihost].remotehost
      endelse
   endfor

   for ihost=0, nhost-1 do begin
      if (keyword_set(hostlist[ihost].remotedir)) then begin
         case hostlist[ihost].protocol of
            'ssh' : hostlist[ihost].cpstring = 'scp'
            'ssh1': hostlist[ihost].cpstring = 'scp1'
            'ssh2': hostlist[ihost].cpstring = 'scp2'
            'rsh' : hostlist[ihost].cpstring = 'rcp'
            ''    : hostlist[ihost].cpstring = ''
            else  : message, 'Invalid protocol: ' + hostlist[ihost].protocol
         endcase
      endif else begin
         ; Retain CPSTRING=''
         hostlist[ihost].remotedir = topdir
      endelse
   endfor

   return, hostlist
end

;------------------------------------------------------------------------------
function batch_if_done, remotehost, remotedir, protocol, pid

   sq = "'"

   ; Test whether a remote job is done...
;   if (keyword_set(protocol)) then begin
;      prothost = protocol + ' ' + remotehost + ' '
;      batch_spawn, prothost + sq+' ps -p '+string(pid)+sq, retstring
;   endif else begin
;      batch_spawn, 'ps -p '+string(pid), retstring
;   endelse

   ; Test if the local parent job is done...
   batch_spawn, 'ps -p '+string(pid), retstring

   ; If there is only one line in the return string, then this process
   ; was not found with the 'ps' command, so the job must be done.
   ; Or if the parent process was a <zombie> or <defunct>, then it must be done.
   if (n_elements(retstring) EQ 1) then begin
      retval = 'DONE'
   endif else if ((strpos(retstring[1],'zombie'))[0] NE -1 $
               OR (strpos(retstring[1],'defunct'))[0] NE -1) then begin
      retval = 'DONE'
   endif else begin
      retval = 'NOTDONE'
   endelse

   splog, 'Status of job on ' + remotehost + ' = ' + retval

   return, retval
end

;------------------------------------------------------------------------------
pro batch_assign_job, ihost, iprog

   common com_batch, hostlist, proglist

   sq = "'"

   if (hostlist[ihost].status NE 'IDLE') then $
    message, 'Host is not idle'

   splog, ''
   splog, 'Assigning job ' + proglist[iprog].progname $
    + ' to host ' + hostlist[ihost].remotehost

   cpstring = hostlist[ihost].cpstring
   prothost = hostlist[ihost].protocol + ' ' + hostlist[ihost].remotehost + ' '

   ;----------
   ; Copy files to a remote machine if a remote directory is specified.
   ; Otherwise, the host is assumed to be able to see the local directory.

   if (keyword_set(cpstring)) then begin

      ; Create directories on remote machine for input files.
      ; Create all remote directories at once, then copy files into one
      ; directory at a time.
      if (keyword_set(proglist[iprog].localfile)) then begin
         allinput = djs_filepath(*proglist[iprog].localfile, $
          root_dir=hostlist[ihost].remotedir)
         junk = fileandpath(allinput, path=newdir)

         iuniq = uniq(newdir, sort(newdir))
         batch_spawn, prothost + 'mkdir -p ' $
          + string(newdir[iuniq]+' ',format='(99a)')

         for i=0, n_elements(iuniq)-1 do begin
            newdir1 = newdir[iuniq[i]]
            indx = where(newdir EQ newdir1)

            tmp1 = string((*proglist[iprog].localfile)[indx]+' ', $
             format='(99a )')
            batch_spawn, cpstring + ' ' + tmp1 + ' ' $
             + hostlist[ihost].remotehost + ':' + newdir1
         endfor
      endif

      ; Create directories on remote machine for output files
      ; (Not necessary if OUTFILE was not specified, and is a null pointer)
      if (keyword_set(proglist[iprog].outfile)) then begin
         alloutput = djs_filepath(*proglist[iprog].outfile, $
          root_dir=hostlist[ihost].remotedir)
         junk = fileandpath(alloutput, path=newdir)
         iuniq = uniq(newdir, sort(newdir))
         batch_spawn, prothost + 'mkdir -p ' $
          + string(newdir[iuniq]+' ',format='(99a)')
      endif

   endif else begin

      ; Only need to create local directories for output files
      ; (Not necessary if OUTFILE was not specified, and is a null pointer)
      if (keyword_set(proglist[iprog].outfile)) then begin
         junk = fileandpath(*proglist[iprog].outfile, path=newdir)
         iuniq = uniq(newdir, sort(newdir))
         batch_spawn, 'mkdir -p ' $
          + string(newdir[iuniq]+' ',format='(99a)')
      endif

   endelse

   ;----------
   ; Launch this job on a remote machine if a protocol is listed,
   ; otherwise launch locally (in which case, the host name should
   ; be "localhost").  Launch the job in the background.

   if (keyword_set(hostlist[ihost].protocol)) then begin
      batch_spawn, prothost + sq+'cd '+hostlist[ihost].remotedir+'; ' $
       +proglist[iprog].command+sq, pid=thispid, unit=thisunit
   endif else begin
      batch_spawn, proglist[iprog].command, pid=thispid, unit=thisunit
   endelse

   ;----------
   ; Save the process ID number and the unit (LUN) number of the pipe
   ; to this process.

splog, "THISPID " ,thispid
   proglist[iprog].pid = thispid
   proglist[iprog].unit = thisunit

   hostlist[ihost].status = 'BUSY'
   hostlist[ihost].progname = proglist[iprog].progname
   proglist[iprog].status = 'RUNNING'

   return
end

;------------------------------------------------------------------------------
pro batch_finish_job, ihost, iprog

   common com_batch, hostlist, proglist

   if (hostlist[ihost].status NE 'BUSY') then $
    message, 'Host is not busy'

   splog, 'Finishing job ' + proglist[iprog].progname $
    + ' on host ' + hostlist[ihost].remotehost

   ; The following will close the pipe to the parent process which is
   ; now a <zombie> or <defunct>.
   free_lun, proglist[iprog].unit

   cpstring = hostlist[ihost].cpstring
   prothost = hostlist[ihost].protocol + ' ' + hostlist[ihost].remotehost + ' '

   if (keyword_set(cpstring) AND keyword_set(proglist[iprog].outfile)) $
     then begin

      alloutput = djs_filepath(*proglist[iprog].outfile, $
       root_dir=hostlist[ihost].remotedir)

      ; Create directories on local machine for output files.
      ; Create all local directories at once, then copy files into one
      ; directory at a time.
      junk = fileandpath(*proglist[iprog].outfile, path=newdir)

      iuniq = uniq(newdir, sort(newdir))
      batch_spawn, 'mkdir -p ' $
       + string(newdir[iuniq]+' ',format='(99a)')
      for i=0, n_elements(iuniq)-1 do begin
         newdir1 = newdir[iuniq[i]]
         indx = where(newdir EQ newdir1)

         ; Copy output files from remote machine to local
         tmp1 = string(hostlist[ihost].remotehost+':' $
          +alloutput[indx]+' ',format='(99a )')
         batch_spawn, cpstring + ' ' + tmp1 + ' ' + newdir1

      endfor

      ; Remove remote output files
      batch_spawn, prothost + ' rm -f ' $
       + string(alloutput+' ',format='(99a )')

      ; Remove remote input files
      allinput = djs_filepath(*proglist[iprog].localfile, $
       root_dir=hostlist[ihost].remotedir)
      batch_spawn, prothost + ' rm -f ' $
       + string(allinput+' ',format='(99a )')

   endif

   hostlist[ihost].status = 'IDLE'
   hostlist[ihost].progname = ''
   proglist[iprog].status = 'DONE'
   proglist[iprog].pid = 0L

   return
end

;------------------------------------------------------------------------------
pro djs_batch, topdir, localfile, outfile, protocol, remotehost, remotedir, $
 command, priority=priority, selecthost=selecthost, wtime=wtime

   common com_batch, hostlist, proglist

   if (NOT keyword_set(wtime)) then wtime = 600 ; Default to wait 10 mins

   if (keyword_set(topdir)) then cd, topdir

   ;----------
   ; Create a list of programs to execute (and their status)

   proglist = create_program_list(localfile, outfile, command, $ 
    priority=priority, selecthost=selecthost)
   nprog = n_elements(proglist)
   splog, 'Number of batch programs = ', nprog
   for iprog=0, nprog-1 do $
    splog, proglist[iprog].progname + ' = ' + proglist[iprog].command

   ;----------
   ; Create a list of available remote hosts (and their status)

   hostlist = create_host_list(protocol, remotehost, remotedir, topdir)
   if (NOT keyword_set(hostlist)) then return
   nhost = n_elements(hostlist)
   splog, 'Number of hosts = ', nhost

   ;----------
   ; Find which programs are already done by looking at local files

;   for iprog=0, nprog-1 do begin
;      qdone = batch_if_done('', '', '', (*proglist[iprog].outfile)[0], $
;       proglist[iprog].endstring)
;      if (qdone EQ 'DONE') then proglist[iprog].status = 'DONE' $
;       else proglist[iprog].status = 'UNASSIGNED'
;   endfor

   ;---------------------------------------------------------------------------
   ; MAIN LOOP
   ;---------------------------------------------------------------------------

   ndone = -1
   while (ndone LT nprog) do begin

      ;----------
      ; Find any jobs that may have completed

      for ihost=0, nhost-1 do begin
         if (hostlist[ihost].status EQ 'BUSY') then begin
            j = (where(proglist.progname EQ hostlist[ihost].progname, ct))[0]
            if (ct NE 1) then $
             message, 'No or multiple program names for this host'
            qdone = batch_if_done(hostlist[ihost].remotehost, $
             hostlist[ihost].remotedir, hostlist[ihost].protocol, $
             proglist[j].pid)
            if (qdone EQ 'DONE') then $
             batch_finish_job, ihost, j
         endif
      endfor

      ;----------

      iunassign = where(proglist.status EQ 'UNASSIGNED', nunassign)
      irun = where(proglist.status EQ 'RUNNING', nrunning)
      idone = where(proglist.status EQ 'DONE', ndone)

      iidle = where(hostlist.status EQ 'IDLE', nidle)

      splog, 'Current time = ', systime()
      splog, 'Number of UNASSIGNED jobs = ', nunassign
      splog, 'Number of RUNNING jobs = ', nrunning
      splog, 'Number of DONE jobs = ', ndone

      ;----------
      ; Assign jobs, doing the highest priority jobs first.
      ; If SELECTHOST is specified for a given job, then that job
      ; can only be assigned to hosts by that name.

      if (nidle GT 0 AND nunassign GT 0) then begin
         for j=0, nidle-1 do begin ; Loop over available hosts
            ; k indexes the jobs that could be assigned to host index iidle[j]
            k = (where(proglist.status EQ 'UNASSIGNED' $
             AND (proglist.selecthost EQ '' $
                  OR proglist.selecthost EQ hostlist[iidle[j]].remotehost)))
            if (k[0] NE -1) then begin
               junk = max(proglist[k].priority, kmax)
               kbest = k
               batch_assign_job, iidle[j], k[kmax]
            endif
         endfor
      endif

      ;----------
      ; Sleep

      splog, 'Sleeping for ', wtime, ' seconds'
      if (ndone LT nprog) then wait, wtime

   endwhile

   splog, 'All jobs have completed at ', systime()

   return
end
;------------------------------------------------------------------------------
