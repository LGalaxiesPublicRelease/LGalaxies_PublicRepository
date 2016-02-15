;------------------------------------------------------------------------------
;+
; NAME:
;   bh_rdfort
;
; PURPOSE:
;   Read 4*E(B-V) for Burstein & Heiles reddening maps.
;   Read values directly from Michael Strauss' Fortran data files.
;
;   This program returns identical results as from Michael Strauss'
;   program, but does not allow the interpolation option.
;
;   Regions with no data return values of -14, -99 or -396.  Valid
;   data is always > -0.22.
;
; CALLING SEQUENCE:
;   value = bh_rdfort( [ gall, galb, infile=infile, outfile=outfile, $
;    bhpath=bhpath ] )
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   gall:       Galactic longitude(s) in degrees
;   galb:       Galactic latitude(s) in degrees
;   infile:     If set, then read LCOORD and BCOORD from this file
;   outfile:    If set, then write results to this file
;   bhpath:     Path name for BH data files
;
; OUTPUTS:
;   value:      4*E(B-V) from Burstein & Heiles reddening maps
;
; PROCEDURES CALLED:
;   readcol
;
; REVISION HISTORY:
;   Written by D. Schlegel, 2 Oct 1997, Durham
;-
;------------------------------------------------------------------------------
; INTERNAL SUPPORT PROCEDURES:
;
; bh_rdmany()
;------------------------------------------------------------------------------
function bh_rdmany, b_index, l_index, filename, nrec
   l_int = long(l_index) - 1
   b_int = long(b_index) - 1
   num = N_elements(l_int)
   aext = lonarr(num)

   ; Read one element at a time from the ASSOC file.
   ; This is done because FILENAME cannot be called with more than one
   ; index at a time.
   for ii=0L, num-1 do begin
      aext[ii] = filename(b_int[ii]*nrec + l_int[ii])
   endfor

   return, aext
end
;------------------------------------------------------------------------------
function bh_rdfort, gall, galb, infile=infile, outfile=outfile, $
 bhpath=bhpath

   if ( (NOT keyword_set(gall) OR NOT keyword_set(galb) ) $
        AND NOT keyword_set(infile) ) then begin
      print, 'Must set either coordinates or infile'
      return, -1
   endif

   ; If INFILE is defined, then read galactic coordinates from that file
   if (keyword_set(infile)) then $
    readcol, infile, gall, galb, format='F,F'

   if NOT keyword_set(bhpath) then bhpath = ''
   
   nred = 1200
   nhi = 201
   DDTOR = !dpi / 180

   bhfile = ['redsouth.dat', 'rednorth.dat', $
             'hinorth.dat',  'hisouth.dat' ]

   ; Open the four data files, and associate them with variables
   ilun = intarr(4)
   for ifile=0, 3 do begin
      get_lun, itemp
      ilun[ifile] = itemp
      openr, ilun[ifile], bhpath+bhfile[ifile]
   endfor
   redsouth = assoc( ilun[0], lonarr(1))
   rednorth = assoc( ilun[1], lonarr(1))
   hinorth = assoc( ilun[2], lonarr(1))
   hisouth = assoc( ilun[3], lonarr(1))

   aext = fltarr( N_elements(gall) )

   indx = where(galb LT -62.)
   if (indx[0] NE -1) then begin
      b_index = 101.51 + sin(gall[indx]*DDTOR) * (90.+galb[indx]) / 0.3
      l_index  = 101.51 + cos(gall[indx]*DDTOR) * (90.+galb[indx]) / 0.3
      temp = bh_rdmany(b_index, l_index, hisouth, nhi)
      aext[indx] = -36.7 + 0.0357*temp
   endif

   indx = where(galb GE -62. AND galb LE -10.)
   if (indx[0] NE -1) then begin
      b_index = -(galb[indx] + 10.)/0.6 + 1.51
      l_index = gall[indx]/0.3 + 1.01 < 1200.
      aext[indx] = bh_rdmany(b_index, l_index, redsouth, nred)
   endif

   indx = where(galb GE 10. AND galb LT 62.)
   if (indx[0] NE -1) then begin
      b_index = (galb[indx]-10.)/0.6 + 1.51
      l_index = gall[indx]/0.3 + 1.01 < 1200.
      aext[indx] = bh_rdmany(b_index, l_index, rednorth, nred)
   endif

   indx = where(galb GE 62.)
   if (indx[0] NE -1) then begin
      b_index = 101.51 + sin(gall[indx]*DDTOR)*(90.-galb[indx])/0.3
      l_index = 101.51 + cos(gall[indx]*DDTOR)*(90.-galb[indx])/0.3
      temp = bh_rdmany(b_index, l_index, hinorth, nhi)
      aext[indx] = -36.7 + 0.0357*temp
   endif

   aext = long(aext)

   ; Changed as per e-mail from David Burstein to Saul Perlmutter 12jul96
   ; The constant was 0.005 in Strauss' original program, and the return
   ; value was also forced to be non-negative.
   value = 0.00005 + aext/250.

   ; Return -99 if within 10 degrees of the Galactic plane
   indx = where(galb GT -10. AND galb LT 10.)
   if (indx[0] NE -1) then begin
      value[indx] = -99.0
   endif

   ; Close input files
   for ifile=0, 3 do begin
      close, ilun[ifile]
      free_lun, ilun[ifile]
   endfor

   ; If OUTFILE is defined, then write to output file
   if (keyword_set(outfile)) then begin $
      get_lun, olun
      openw, olun, outfile
      printf, olun, format='(f12.5)', value
      close, olun
      free_lun, olun
   endif

   return, value
end
;------------------------------------------------------------------------------
