pro astrmfix, hdr, chip
;+
; NAME:
;       ASTRMFIX
; PURPOSE:
;       Calculate a rough HST WFPC or FOC astrometry solution
; EXPLANATION:
;       This program will calculate a rough HST WFPC or FOC astrometry solution
;       using the keyword PSANGLEV3 which gives the angle of the V3 axis of
;       HST.    Called by WFPCREAD.
;
; CALLING SEQUENCE:
;       AstrmFix, hdr, chip
;
; INPUT - OUTPUT:
;       hdr - FITS header (string array) from either WFPC or FOC.   Header will
;               be updated with rough astrometry 
;
; INPUT:        
;       chip - Scalar (typically 0-3) giving the WFPC chip to read.
;
; PROCEDURES CALLED:
;       EXTAST, SXPAR(), SXADDPAR
; HISTORY:
;       ??-???-???? Written by Eric W. Deutsch
;       22-OCT-1992 Changed all calculations to double precision. (E. Deutsch)
;       22-OCT-1992 Updated PC Pixel size of 0.04389 from WFPC IDT OV/SV manual(EWD)
;       22-OCT-1992 Updated WF Pixel size of 0.1016 from WFPC IDT OV/SV manual(EWD)
;       11-JAN-1993 Added warning message and changed CD001001... to CD1_1... (EWD)
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Remove calls to obsolete !ERR variable  W. Landsman   December 2000
;-

  if (n_params(0) lt 2) then begin
    print,'Call: ASTRMFIX,header_array,chip_number (1-8)'
    print,'e.g.> ASTRMFIX,h,3'
    return
    endif

  check=sxpar(hdr,'INSTRUME')
  if (strn(check) ne 'WFPC') and (strn(check) ne 'FOC') then begin
    print,'ASTRMFIX: Unrecognized INSTRUME in Header.  [END]'
    return
    endif
  V3word='PSANGLV3'
  angle=sxpar(hdr,'PSANGLV3',Count = N_angle)
  if (N_angle EQ 0) then begin
    V3word='PA_V3'
    angle = sxpar(hdr,'PA_V3', Count = N_pa_v3)
    if (N_pa_v3 EQ 0) then begin
      print,'ASTRMFIX: Keyword PSANGLV3 or PA_V3 must be in header.  [END]'
      return
      endif
    endif
  print,'HST V3 Position Angle['+V3word+']: ',angle

  if (strn(check) eq 'FOC') then begin
    cam=strn(check) & hand=-1
    check=sxpar(hdr,'OPTCRLY')
    pix_size=0 & base_angle=0
    if (strn(check) eq 'F96') then begin
      base_angle=55.26
      pix_size=.022
      endif
    if (strn(check) eq 'F48') then begin
      base_angle=-65.
      pix_size=.044
      endif
    if (pix_size eq 0) then begin
      print,'ASTRMFIX: Unrecognized Relay or missing keyword OPTCRLY.  [END]'
      return
      endif
    angle=base_angle+angle
    goto,CONVERT
    endif

  if (chip gt 4) then chip=chip-4
  if (chip lt 1) or (chip gt 4) then begin
    print,'ASTRMFIX: Chip number must be 1-8 or 1-4.  [END]'
    return
    endif

  cam=strn(sxpar(hdr,'CAMERA')) & hand=1
  if (cam ne 'WF') and (cam ne 'WFC') and (cam ne 'PC') then begin
    print,'ASTRMFIX: Unrecognized Camera type "',cam,'".  [STOP]'
    return
    endif

  if (cam eq 'PC') then begin
    base_angle=270
    pix_size=.04389d                    ; from WFPC IDL OV/SV report
  endif else begin
    base_angle=315
    pix_size=.1016                      ; from WFPC IDL OV/SV report
    endelse

  chip_rot=[3,0,1,2]
  base_angle=base_angle+chip_rot[chip-1]*90
  angle=base_angle-angle+180

CONVERT:
  flag=1
  while (flag eq 1) do begin
    flag=0
    if (angle lt 0) then begin
      angle=angle+360 & flag=1
      endif
    if (angle ge 360) then begin
      angle=angle-360 & flag=1
      endif
    endwhile

  extast,hdr,astr,noparams                      ; to do CD1_1 -> CD001001
                                                ;  if necessary
                                                ; out of date (at 1/11/93) but
                                                ; will the reverse become
                                                ; fashionable again?
  cd = dblarr(2,2)
  print,'North Points to angle: ',angle
  tmp='Left' & if (hand eq -1) then tmp='Right'
  print,'East Points 90 degrees to the ',tmp
  angle=angle/180.d*!dpi
  cd[0,0] = cos(angle+!dpi/2.*hand)
  cd[0,1] = sin(angle+!dpi/2.*hand)
  cd[1,0] = cos(angle)
  cd[1,1] = sin(angle)
  cd = cd*pix_size/3600.d

  message,'Astrometry in header modified!!',/INF

  sxaddpar,hdr,'CD1_1',cd[0,0],' dRA/dX  Calculated from '+V3word
  sxaddpar,hdr,'CD1_2',cd[0,1],' dRA/dY  Calculated from '+V3word
  sxaddpar,hdr,'CD2_1',cd[1,0],' dDEC/dX Calculated from '+V3word
  sxaddpar,hdr,'CD2_2',cd[1,1],' dDEC/dY Calculated from '+V3word

  return
end
