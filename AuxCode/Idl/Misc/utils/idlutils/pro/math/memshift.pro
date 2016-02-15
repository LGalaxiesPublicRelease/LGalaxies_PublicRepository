;+
; NAME:
;   memshift
;
; PURPOSE:
;   Shift elements in an array, which can be a structure array.
;
; CALLING SEQUENCE:
;   memshift, array, isrc, idest, nmove
;
; INPUTS:
;   array      - Array of any type except string or pointer; structures
;                are allowed
;   isrc       - Starting position in memory to copy from
;   idest      - Starting position in memory to copy to
;   nmove      - Number of array elements to copy
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;   array      - (Modified.)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine is equivalent to the IDL command:
;     array[idest:idest+nmove-1] = array[isrc:isrc+nmove-1]
;   but is faster, at least on the Linux platform.  Note that the
;   memory in the source and destination can be overlapping.
;   The C code calls the Unix memmove() function.
;
;   If an attempt is made to copy out-of-bounds, then this procedure
;   intentionally crashes using the MESSAGE function.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   Dynamic link to memshift.o
;
; REVISION HISTORY:
;   18-Oct-2003  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro memshift, array, isrc1, idest1, nmove1

   num = n_elements(array)
   isrc = long(isrc1[0])
   idest = long(idest1[0])
   nmove = long(nmove1[0])
   if (isrc LT 0 OR idest LT 0 OR nmove LE 0 $
    OR isrc+nmove GT num OR idest+nmove GT num) then $
    message, 'Indices are out-of-bounds!'

   case size(array, /type) of
   1 : nper = 1 ; BYTE
   2 : nper = 2 ; INT
   3 : nper = 4 ; LONG
   4 : nper = 4 ; FLOAT
   5 : nper = 8 ; DOUBLE
   6 : nper = 8 ; complex FLOAT
   7 : return ; STRING
   8 : nper = n_tags(array, /length) ; STRUCTURE
   9 : nper = 16 ; complex DOUBLE
   12 : nper = 2 ; UINT
   13 : nper = 4 ; ULONG
   14 : nper = 8 ; LONG64
   15 : nper = 8 ; ULONG64
   else : return ; All other types are unsupported
   endcase

   soname = filepath('libmath.'+idlutils_so_ext(), $
    root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')
   retval = call_external(soname, 'memshift', array, $
    isrc*nper, idest*nper, nmove*nper)

   return
end
;------------------------------------------------------------------------------
