;+
; NAME:
;   ucac_readindex()
;
; PURPOSE:
;   Return the indices for seeking within the raw UCAC data files.
;
; CALLING SEQUENCE:
;   uindex = ucac_readindex()
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;   uindex     - Structure with raw UCAC data indices
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   readfmt
;
; REVISION HISTORY:
;   27-May-2003  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function ucac_readindex

   common com_ucac, uindex

   if (keyword_set(uindex)) then return, uindex

   ;----------
   ; Check inputs

   ucac_dir = getenv('UCAC_DIR')
   if (NOT keyword_set(ucac_dir)) then begin
      print, 'Environment variable UCAC_DIR must be set!'
      return, 0
   endif

   ;----------
   ; Read the index file (if not already read and cached in memory)

; Should we read the binary version of this file instead of ASCII ???
   indexfile = filepath('u2index.txt', root_dir=ucac_dir, subdir='info')
   readfmt, indexfile, 'I6,I8,I9,I4,I4,F6.1,F5.1', $
    nsbin, naz, nat, zn, jj, dcmax, ramax, skipline=10
   uindex = replicate( create_struct( $
    'NSBIN', 0L, $
    'NAZ'  , 0L, $
    'NAT'  , 0L, $
    'ZN'   , 0L, $
    'JJ'   , 0L, $
    'DCMAX', 0d, $
    'RAMAX', 0d ), n_elements(nsbin))
   uindex.nsbin = nsbin
   uindex.naz = naz
   uindex.nat = nat
   uindex.zn = zn
   uindex.jj = jj
   uindex.dcmax = dcmax
   uindex.ramax = ramax

   return, uindex
end
;------------------------------------------------------------------------------
