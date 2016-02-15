pro extast,hdr,astr,noparams, alt=alt
;+
; NAME:
;     EXTAST
; PURPOSE:
;     Extract ASTrometry parameters from a FITS image header.
; EXPLANATION:
;     The astrometry in the header can be in either CD (Coordinate
;     description) format, or CROTA and CDELT (AIPS-type) format.
;     However, the output astrometry will always be in CD format.
;
; CALLING SEQUENCE:
;     EXTAST, hdr, [ astr, noparams, ALT= ]   
;
; INPUT:
;     HDR - variable containing the FITS header (string array)
;
; OUTPUTS:
;     ASTR - Anonymous structure containing astrometry info from the FITS 
;             header ASTR always contains the following tags (even though 
;             some projections do not require all the parameters)
;      .CD   -  2 x 2 array containing the astrometry parameters CD1_1 CD1_2
;               in DEGREES/PIXEL                                 CD2_1 CD2_2
;      .CDELT - 2 element vector giving physical increment at reference pixel
;      .CRPIX - 2 element vector giving X and Y coordinates of reference pixel
;               (def = NAXIS/2) in FITS convention (first pixel is 1,1)
;      .CRVAL - 2 element double precision vector giving R.A. and DEC of 
;             reference pixel in DEGREES
;      .CTYPE - 2 element string vector giving projection types, default
;             ['RA---TAN','DEC--TAN']
;      .LONGPOLE - scalar giving native longitude of the celestial pole 
;             (default = 180 for zenithal projections) 
;      .LATPOLE - scalar giving native latitude of the celestial pole default=0)
;      .PV2 - Vector of projection parameter associated with latitude axis
;             PV2 will have up to 21 elements for the ZPN projection, up to 3 
;             for the SIN projection and no more than 2 for any other 
;             projection  
;      .DISTORT - optional substructure specifying any distortion parameters
;                 currently implemented only for "SIP" (Spitzer Imaging 
;                 Polynomial) distortion parameters
;
;       NOPARAMS -  Scalar indicating the results of EXTAST
;             -1 = Failure - Header missing astrometry parameters
;             1 = Success - Header contains CROTA + CDELT (AIPS-type) astrometry
;             2 = Success - Header contains CDn_m astrometry, rec.    
;             3 = Success - Header contains PCn_m + CDELT astrometry. 
;             4 = Success - Header contains ST  Guide Star Survey astrometry
;                           (see gsssextast.pro )
; OPTIONAL INPUT KEYWORDS:
;       ALT -  single character 'A' through 'Z' or ' ' specifying an alternate 
;              astrometry system present in the FITS header.    The default is
;              to use the primary astrometry or ALT = ' '.   If /ALT is set, 
;              then this is equivalent to ALT = 'A'.   See Section 3.3 of 
;              Greisen & Calabretta (2002, A&A, 395, 1061) for information about
;              alternate astrometry keywords.
; PROCEDURE:
;       EXTAST checks for astrometry parameters in the following order:
;
;       (1) the CD matrix PC1_1,PC1_2...plus CDELT*, CRPIX and CRVAL
;       (3) the CD matrix CD1_1,CD1_2... plus CRPIX and CRVAL.   
;       (3) CROTA2 (or CROTA1) and CDELT plus CRPIX and CRVAL.
;
;       All three forms are valid FITS according to the paper "Representations 
;       of World Coordinates in FITS by Greisen and Calabretta (2002, A&A, 395,
;       1061 http://www.aoc.nrao.edu/~egreisen) although form (1) is preferred/
;
; NOTES:
;       An anonymous structure is created to avoid structure definition
;       conflicts.    This is needed because some projection systems
;       require additional dimensions (i.e. spherical cube
;       projections require a specification of the cube face).
;
; PROCEDURES CALLED:
;      GSSSEXTAST, ZPARCHECK
; REVISION HISTORY
;      Written by B. Boothman 4/15/86
;      Accept CD001001 keywords               1-3-88
;      Accept CD1_1, CD2_1... keywords    W. Landsman    Nov. 92
;      Recognize GSSS FITS header         W. Landsman    June 94
;      Converted to IDL V5.0   W. Landsman   September 1997
;      Get correct sign, when converting CDELT* to CD matrix for right-handed
;      coordinate system                  W. Landsman   November 1998
;      Consistent conversion between CROTA and CD matrix  October 2000
;      CTYPE = 'PIXEL' means no astrometry params  W. Landsman January 2001
;      Don't choke if only 1 CTYPE value given W. Landsman  August 2001
;      Recognize PC00n00m keywords again (sigh...)  W. Landsman December 2001
;      Recognize GSSS in ctype also       D. Finkbeiner Jan 2002
;      Introduce ALT keyword              W. Landsman June 2003
;      Fix error introduced June 2003 where free-format values would be
;      truncated if more than 20 characters.  W. Landsman Aug 2003
;      Further fix to free-format values -- slash need not be present Sep 2003
;      Default value of LATPOLE is 90.0  W. Landsman February 2004
;      Allow for distortion substructure, currently implemented only for
;          SIP (Spitzer Imaging Polynomial)   W. Landsman February 2004 
;      Correct LONGPOLE computation if CTYPE = ['*DEC','*RA'] W. L. Feb. 2004
;      Assume since V5.3 (vector STRMID)  W. Landsman Feb 2004
;      Yet another fix to free-format values   W. Landsman April 2004
;      Introduce PV2 tag to replace PROJP1, PROJP2.. etc.  W. Landsman May 2004
;      Convert NCP projection to generalized SIN   W. Landsman Aug 2004
;-
; On_error,2

 if ( N_params() LT 2 ) then begin
     print,'Syntax - EXTAST, hdr, astr, [ noparams, ALT = ]'
     return
 endif

 proj0 = ['CYP','CEA','CAR','MER','SFL','PAR','MOL','AIT','BON','PCO', $
          'TSC','CSC','QSC']
 radeg = 180.0D0/!DPI
 keyword = strtrim(strmid( hdr, 0, 8), 2)

; Extract values from the FITS header.   This is either up to the first slash
; (free format) or first space

 space = strpos( hdr, ' ', 10) + 1
 slash = strpos( hdr, '/', 10)  > space
 
 N = N_elements(hdr)
 len = (slash -10) > 20
 len = reform(len,1,N)
 lvalue = strtrim(strmid(hdr, 10, len),2)
 zparcheck,'EXTAST',hdr,1,7,1,'FITS image header'   ;Make sure valid header
 noparams = -1                                    ;Assume no astrometry to start

 if N_elements(alt) EQ 0 then alt = '' else if (alt EQ '1') then alt = 'A' $
    else alt = strupcase(alt)
 ctype = ['','']
 l = where(keyword EQ 'CTYPE1'+alt,  N_ctype1)
 if N_ctype1 GT 0 then ctype[0] = lvalue[l[0]]
 l = where(keyword EQ 'CTYPE2'+alt,  N_ctype2)
 if N_ctype2 GT 0 then ctype[1] = lvalue[l[0]]
 remchar,ctype,"'"
 ctype = strtrim(ctype,2)

; If the standard CTYPE* astrometry keywords not found, then check if the
; ST guidestar astrometry is present

 check_gsss = (N_ctype1 EQ 0)
 if N_ctype1 GE 1  then check_gsss = (strmid(ctype[0], 5, 3) EQ 'GSS')

 if check_gsss then begin

        l = where(keyword EQ 'PPO1'+alt,  N_ppo1)
        if N_ppo1 EQ 1 then begin 
                gsssextast, hdr, astr, gsssparams
                if gsssparams EQ 0 then noparams = 4
                return
        endif
        ctype = ['RA---TAN','DEC--TAN']
  endif

  if (ctype[0] EQ 'PIXEL') then return
  if N_ctype2 EQ 1 then if (ctype[1] EQ 'PIXEL') then return

 crval = dblarr(2)

 l = where(keyword EQ 'CRVAL1'+alt,  N_crval1)
 if N_crval1 GT 0 then crval[0] = lvalue[l[0]]
 l = where(keyword EQ 'CRVAL2'+alt,  N_crval2)
 if N_crval2 GT 0 then crval[1] = lvalue[l[0]]
 if (N_crval1 EQ 0) or (N_crval2 EQ 0) then return  


 crpix = fltarr(2)
 l = where(keyword EQ 'CRPIX1'+alt,  N_crpix1)
 if N_crpix1 GT 0 then crpix[0] = lvalue[l[0]]
 l = where(keyword EQ 'CRPIX2'+alt,  N_crpix2)
 if N_crpix2 GT 0 then crpix[1] = lvalue[l[0]]
 if (N_crpix1 EQ 0) or (N_crpix2 EQ 0) then return  


 cd = dblarr(2,2)
cdelt = [1.0d,1.0d]

 l = where(keyword EQ 'PC1_1' + alt,  N_pc11) 
 if N_PC11 GT 0 then begin 
        cd[0,0]  = lvalue[l]
        l = where(keyword EQ 'PC1_2' + alt,  N_pc12) 
        if N_pc12 GT 0 then cd[0,1]  = lvalue[l[0]]
        l = where(keyword EQ 'PC2_1' + alt,  N_pc21) 
        if N_pc21 GT 0 then cd[1,0]  = lvalue[l[0]]
        l = where(keyword EQ 'PC2_2' + alt,  N_pc22) 
        if N_pc22 GT 0 then cd[1,1]  = lvalue[l[0]]
         l = where(keyword EQ 'CDELT1' + alt,  N_cdelt1) 
        if N_cdelt1 GT 0 then cdelt[0]  = lvalue[l[0]]
        l = where(keyword EQ 'CDELT2' + alt,  N_cdelt2) 
        if N_cdelt2 GT 0 then cdelt[1]  = lvalue[l[0]]
        noparams = 3
 endif else begin 

    l = where(keyword EQ 'CD1_1' + alt,  N_cd11) 
     if N_CD11 GT 0 then begin        ;If CD parameters don't exist, try CROTA
        cd[0,0]  = strtrim(lvalue[l],2)
        l = where(keyword EQ 'CD1_2' + alt,  N_cd12) 
        if N_cd12 GT 0 then cd[0,1]  = lvalue[l[0]]
        l = where(keyword EQ 'CD2_1' + alt,  N_cd21) 
        if N_cd21 GT 0 then cd[1,0]  = lvalue[l[0]]
        l = where(keyword EQ 'CD2_2' + alt,  N_cd22) 
        if N_cd22 GT 0 then cd[1,1]  = lvalue[l[0]]
        noparams = 2
    endif else begin

; Now get rotation, first try CROTA2, if not found try CROTA1, if that
; not found assume North-up.   Then convert to CD matrix - see Section 5 in
; Greisen and Calabretta

        l = where(keyword EQ 'CDELT1' + alt,  N_cdelt1) 
        if N_cdelt1 GT 0 then cdelt[0]  =lvalue[l[0]]
        l = where(keyword EQ 'CDELT2' + alt,  N_cdelt2) 
        if N_cdelt2 GT 0 then cdelt[1]  = lvalue[l[0]]
        if (N_cdelt1 EQ 0) or (N_Cdelt2 EQ 0) then return   ;Must have CDELT1 and CDELT2

        l = where(keyword EQ 'CROTA2' + alt,  N_crota) 
        if N_Crota EQ 0 then $
            l = where(keyword EQ 'CROTA1' + alt,  N_crota) 
        if N_crota EQ 0 then crota = 0.0d else $
                             crota = double(lvalue[l[0]])/RADEG
        cd = [ [cos(crota), -sin(crota)],[sin(crota), cos(crota)] ] 
 
       noparams = 1           ;Signal AIPS-type astrometry found
     
  endelse
  endelse

  proj = strmid( ctype[0], 5, 3)
  case proj of 
 'ZPN': npv = 21
 'SZP': npv = 3
 else:  npv = 2
  endcase

  index = proj EQ 'ZPN' ? strtrim(indgen(npv),2) : strtrim(indgen(npv)+1,2)
      pv2 = dblarr(npv)
      for i=0,npv-1 do begin 
      l = where(keyword EQ 'PV2_' + index[i] + alt,  N_pv2)
      if N_pv2 GT 0 then pv2[i] = lvalue[l[0]] 
      endfor
 
          
  l = where(keyword EQ 'PV1_3' + alt,  N_pv1_3)
  if N_pv1_3 GT 0 then  longpole = double(lvalue[l[0]]) else begin
      l = where(keyword EQ 'LONPOLE' + alt,  N_lonpole)
      if N_lonpole GT 0 then  longpole = double(lvalue[l[0]]) 
  endelse

; If LONPOLE (or PV1_3) is not defined in the header, then we must determine 
; its default value.    This depends on the value of theta0 (the native
; longitude of the fiducial point) of the particular projection)

  conic = (proj EQ 'COP') or (proj EQ 'COE') or (proj EQ 'COD') or $
          (proj EQ 'COO')


  if N_elements(longpole) EQ 0 then  begin 
    if conic then begin 
      if N_pv2 EQ 0 then message, $
     'ERROR -- Conic projections require a PV2_1 keyword in FITS header' else $
      theta0 = PV2[0]
    endif else if (proj EQ 'ZAP') or (proj EQ 'SZP') or (proj EQ 'TAN') or $
          (proj EQ 'STG') or (proj EQ 'SIN') or (proj EQ 'ARC') or $
          (proj EQ 'ZPN') or (proj EQ 'ZEA') or (proj EQ 'AIR') then begin
       theta0 = 90.0
    endif else theta0 = 0. 
    celcoord = strmid(ctype[1],0,4)
;Check to see if RA and DEC are reversed in CRVAL
    if (celcoord EQ 'RA--') or (celcoord EQ 'GLON') or (celcoord EQ 'ELON') $
           then cellat = crval[0] else cellat = crval[1]
    if cellat GE theta0 then longpole = 0.0 else longpole = 180.0
 endif

  l = where(keyword EQ 'LATPOLE' + alt,  N_latpole)
  if N_latpole GT 0 then  latpole = double(lvalue[l[0]]) else latpole = 90.0d


; Convert NCP projection to generalized SIN projection (see Section 6.1.2 of 
; Calabretta and Greisen (2002)

  if proj EQ 'NCP' then begin
       ctype = repstr(ctype,'NCP','SIN')
       PV2 = [0., 1/tan(crval[1]/radeg) ]
       longpole = 180.0
  endif 

; Note that the dimensions and datatype of each tag must be explicit, so that
; there is no conflict with structure definitions from different FITS headers

  ASTR = {CD: cd, CDELT: cdelt, $
                CRPIX: crpix, CRVAL:crval, $
                CTYPE: string(ctype), LONGPOLE: double( longpole[0]),  $
                LATPOLE: double(latpole[0]), PV2: pv2 }

; Check for any distortion keywords

  if strlen(ctype[0]) GE 12 then begin
       distort_flag = strmid(ctype[0],9,3)
       case distort_flag of 
       'SIP': begin
              l = where(keyword EQ 'A_ORDER',  N) 
              if N GT 0 then a_order  = lvalue[l[0]] else a_order = 0
              l = where(keyword EQ 'B_ORDER',  N) 
              if N GT 0 then b_order  = lvalue[l[0]] else b_order = 0
              l = where(keyword EQ 'AP_ORDER',  N) 
              if N GT 0 then ap_order  = lvalue[l[0]] else ap_order = 0
              l = where(keyword EQ 'BP_ORDER',  N) 
              if N GT 0 then bp_order  = lvalue[l[0]] else bp_order = 0
  a = fltarr(a_order+1,a_order+1) & b = fltarr(b_order+1,b_order+1) 
  ap = fltarr(ap_order+1,ap_order+1) &  bp = fltarr(bp_order+1,bp_order+1)

  for i=0, a_order do begin
    for j=0, a_order do begin
     l = where(keyword EQ 'A_' + strtrim(i,2) + '_' + strtrim(j,2), N)
     if N GT 0 then a[i,j] = lvalue[l[0]]
  endfor & endfor

   for i=0, b_order  do begin
    for j=0, b_order do begin
     l = where(keyword EQ 'B_' + strtrim(i,2) + '_' + strtrim(j,2), N)
     if N GT 0 then b[i,j] = lvalue[l[0]]
  endfor & endfor

   for i=0, bp_order do begin
    for j=0, bp_order do begin
     l = where(keyword EQ 'BP_' + strtrim(i,2) + '_' + strtrim(j,2), N)
     if N GT 0 then bp[i,j] = lvalue[l[0]]
  endfor & endfor

    for i=0, ap_order do begin
    for j=0, ap_order do begin
     l = where(keyword EQ 'AP_' + strtrim(i,2) + '_' + strtrim(j,2), N)
     if N GT 0 then ap[i,j] = lvalue[l[0]]
  endfor & endfor
   
  distort = {name:distort_flag, a:a, b:b, ap:ap, bp:bp}
  astr = create_struct(temporary(astr), 'distort', distort)
  end
  else: message,/con,'Unrecognized distortion acronym: ' + distort_flag 
  endcase
  endif
  return
  end
