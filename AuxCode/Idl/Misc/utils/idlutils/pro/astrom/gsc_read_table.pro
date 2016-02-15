;+
; NAME:
;   gsc_read_table
;
; PURPOSE:
;   Read one table from the GSC (Guide star catalogue)
;
; CALLING SEQUENCE:
;   cat=gsc_read_table(fname, maglim=maglim, poserrlim=poserrlim)
;
; INPUTS:
;   fname      - file name(s) of catalogue FITS table (e.g. 0060.gsc )
;                      (can be array)
; KEYWORDS:
;   maglim     - array of [lomag, himag] magnitudes to keep
;   poserrlim  - discard stars with poserr higher than this limit
;	
; OUTPUTS:
;   cat        - structure returned by function call
;
; EXAMPLES:
;
; COMMENTS:
;   Note: these catalogues are written as F9.5 ASCII FITS tables.  So 
;   there is no reason to read them as double precision(?)
;
; REVISION HISTORY:
;   2002-Mar-27  Written by D. Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function gsc_read_table, fname, maglim=maglim, poserrlim=poserrlim

  cat1 = {ra: 0., dec:0., mag:0.}  ; no double precision...

  nfile = n_elements(fname)
  for i=0L, nfile-1 do begin 
     gsc = mrdfits(fname[i], 1, /silent)
     
     w = where((gsc.mag GT maglim[0]) AND (gsc.mag LT maglim[1]) AND $
               (gsc.pos_err LT poserrlim) AND $
               (gsc.class EQ 0) AND (gsc.dec_deg LT 90.001) AND $
               (gsc.dec_deg GT -30.), ct)
     
     if ct GT 0 then begin 
        cat  = replicate(cat1, ct)
        gscw = gsc[w]
        cat.ra  = gscw.ra_deg
        cat.dec = gscw.dec_deg
        cat.mag = gscw.mag
        
        catall = keyword_set(catall) ? [catall, cat] : cat
     endif 
  endfor 
  
  if NOT keyword_set(catall) then catall = 0 ;set default
  
  return, catall
end

