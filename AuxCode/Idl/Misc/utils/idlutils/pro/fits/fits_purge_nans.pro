;-----------------------------------------------------------------------
;+
; NAME:
;   fits_purge_nans
;
; PURPOSE:
;   Purge invalid (NaN) values from FITS headers
;
; CALLING SEQUENCE:
;   fits_purge_nans, hdr, [ /verbose ]
;
; INPUTS:
;   hdr       - FITS header
;
; OPTIONAL INPUTS:
;   verbose   - If set, then report which keywords are being disposed
;
; OUTPUTS:
;   hdr       - (Modified)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This procedure removes header keywords that are not strings,
;   but are returned by SXPAR() as NaN-valued (as determined by
;   the FINITE() function).
;
;   This procedure was written to deal with raw SDSS headers that
;   sometimes contain header lines like 'ALT     = NaN', where the
;   NaN is not contained in quotes.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   splog
;   sxpar
;
; REVISION HISTORY:
;   14-Jul-2004  Written by David Schlegel, Princeton.
;-
;-----------------------------------------------------------------------
pro fits_purge_nans, hdr, verbose=verbose

   ; Need at least 1 parameter
   if (N_params() LT 1) then begin
      print, 'Syntax - fits_purge_nans, hdr'
      return
   endif

   num = n_elements(hdr)
   qgood = bytarr(num) + 1B
   for i=0L, num-1 do begin
      thiskey = strmid(hdr[i],0,8)
      if (thiskey NE 'COMMENT ') then begin
         thisval = sxpar([hdr[i],'END'], thiskey)
         if (size(thisval,/tname) NE 'STRING') then begin
            if (finite(thisval) EQ 0) then begin
               if (keyword_set(verbose)) then $
                splog, 'Warning: Disposing of invalid FITS header keyword ' + thiskey
               qgood[i] = 0B
            endif
         endif
      endif
   endfor

   igood = where(qgood, ngood)
   if (ngood EQ 0) then hdr = 0 $
    else if (ngood LT num) then hdr = hdr[igood]

   return
end
;-----------------------------------------------------------------------
