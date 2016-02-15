;+
; NAME:  
;       unixfind
;
; PURPOSE: 
;       recursively search for filenames faster than IDL's FILE_SEARCH family.
;
; IDL's file_search, etc. actually do a stat() call on each examined file, 
; which can make them pretty painful  over NFS.
; This is a shortcut for a few simple cases. 
;
; Over-optimization? Well, here is the search which prodded me into
; writing this:
;
;   fl1 = unixfind('/u/dss/spectro_v5','spFluxcalib*')
;     7.2860990 seconds (after cache load)
;
;   fl2 = file_search('/u/dss/spectro_v5', 'spFluxcalib*')
;     280.20624 seconds (after cache load, right after the above)
;
; OK, OK, filling the cache can still be huge.
;
; CALLING SEQUENCE:
;       filelist = unixfind(rootdir, pattern, [/onlydirs, /onlyfiles,
;       maxdepth=N, /striproot, debug=debug])
;
; INPUTS:
;       rootdir   - the name of the root directory to search. No
;                   single quotes allowed. If this is a symbolic link
;                   it is followed.
;       pattern   - a "find"-compatible pattern. No single-quotes
;                   allowed.
;
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
;       /onlydirs   - limit search to directory names (i.e. "-type d")
;       /onlyfiles  - limit search to file names (i.e. "-type f")
;       maxdepth=N  - limit search to N directories below the root.
;       /striproot  - if set, remove the rootdir from the returned filenames.
;       debug       - if >0, print command which we spawn.
;                     if >1, do not execute it.
; OUTPUTS:
;       Return the list of matched files or directories.
;       Returns '' if no matches or on error.
;
; SIDE EFFECTS:
;       Spawn, of, the, unix, devil.
;
; RESTRICTIONS:
;       Input strings are passed to a unix command unchecked. Fools beware.
;       Few of the glitzy find options supported.
;       Lists of roots or patterns should be supported.
;
; EXAMPLE:
;   platedirs = unixfind('/u/dss/spectro', '[0-9]*', maxdepth=1, /onlydirs)
;   maps = unixfind('/u/dss/astrolog', 'plPlugMap*')
;
; MODIFICATION HISTORY:
;
;-
;
function unixfind,root,namematch,$
                  onlydirs=onlydirs,onlyfiles=onlyfiles,$
                  maxdepth=maxdepth,striproot=striproot,$
                  debug=debug

    ; Finesse top-level links
    if (strmid(root, 0, 1, /rev) ne '/') then $
      root = string(root, "/")

    ; Deal with this more intelligently should it come up.
    if (strpos(string(root,namematch), "'") GE 0) then begin
        print,"Cannot handle single-quotes. Fixme if you care."
        return,''
    end

    if (keyword_set(onlydirs) and keyword_set(onlyfiles)) then begin
        splog, "/onlydirs and /onlyfiles cannot both be set."
        return,''
    end

    opts = ['']

    ; Optionally limit the search depth
    if (keyword_set(maxdepth)) then opts = [opts, string(" -maxdepth ", maxdepth, format='(a,i2)')]

    ; Optionally limit to directories or files
    if (keyword_set(onlydirs)) then opts = [opts, "-type d"] $
    else if (keyword_set(onlyfiles)) then opts = [opts, "-type f"]

    optstring = strjoin(opts, ' ') ; "-a "?
    cmd = "find '" + root + "'" + optstring + " -name '" + namematch + "' -print"

    if (keyword_set(debug)) then begin
        if (debug GT 1) then begin 
            splog, "NOT spawning: ", cmd
            return,''
        end else splog, "spawning: ", cmd
    end

    spawn, cmd, files

    ; Optionally strip the root dir.
    if (keyword_set(striproot) and size(files, /dim) GT 0) then $
      files = strmid(files, strlen(root))
    
    return, files
end
