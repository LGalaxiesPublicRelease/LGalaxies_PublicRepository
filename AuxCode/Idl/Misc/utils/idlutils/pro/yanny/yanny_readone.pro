;+
; NAME:
;   yanny_readone
;
; PURPOSE:
;   Read one structure from a Yanny parameter file into an IDL structure.
;
; CALLING SEQUENCE:
;   data = yanny_readone( filename, [ selectname, hdr=, enums=, structs=, $
;    /anonymous, stnames=, /quick, errcode= ] )
;
; INPUTS:
;   filename   - Input file name for Yanny parameter file
;
; OPTIONAL INPUTS:
;   selectname - Name of structure to select.  If not specified, then
;                the first structure is returned.
;
; OUTPUT:
;   data       - Selected structure, or 0 if not found.
;
; OPTIONAL OUTPUTS:
;   hdr        - Header lines in Yanny file, which are usually keyword pairs.
;   enums      - All "typedef enum" structures.
;   structs    - All "typedef struct" structures, which define the form
;                for all the PDATA structures.
;   anonymous  - If set, then all returned structures are anonymous; set this
;                keyword to avoid possible conflicts between named structures
;                that are actually different.
;   stnames    - Names of structures.  If /ANONYMOUS is not set, then this
;                will be equivalent to the IDL name of each structure in PDATA,
;                i.e. tag_names(PDATA[0],/structure_name) for the 1st one.
;                This keyword is useful for when /ANONYMOUS must be set to
;                deal with structures with the same name but different defns.
;                Note that this will contain all of the structure names, even
;                though this routine only returns that data from one of them.
;   quick      - This keyword is only for backwards compatability, and
;                has no effect.
;   errcode    - Returns as non-zero if there was an error reading the file.
;
; COMMENTS:
;   This is a simple wrapper for YANNY_READ.
;
; EXAMPLES:
;   Select the IDCOMMENT structure from the Yanny file idReport-52522.par.
;   Set the /ANONYMOUS flag to prevent conflicts with any other structures
;   of the same name.
;     IDL> data = yanny_readone('idReport-52522.par', 'IDCOMMENT', /anon)
;
; BUGS:
;
; PROCEDURES CALLED:
;   yanny_free
;   yanny_read
;
; REVISION HISTORY:
;   05-Oct-2002  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function yanny_readone, filename, selectname, $
 hdr=hdr, enums=enums, structs=structs, anonymous=anonymous, $
 stnames=stnames, quick=quick, errcode=errcode

   retval = 0

   yanny_read, filename, pp, $
    hdr=hdr, enums=enums, structs=structs, anonymous=anonymous, $
    stnames=stnames, quick=quick, errcode=errcode
   if (keyword_set(selectname)) then begin
      if (n_elements(selectname) GT 1) then begin
         print, 'SELECTNAME must be a scalar'
         yanny_free, pp
         return, retval
      endif
      if (keyword_set(stnames)) then begin
         i = where(stnames EQ selectname[0], ct)
         if (ct GT 0) then begin
            retval = *pp[i[0]]
         endif
      endif
   endif else begin
      if (keyword_set(pp)) then retval = *pp[0]
   endelse

   yanny_free, pp

   return, retval
end
;------------------------------------------------------------------------------
