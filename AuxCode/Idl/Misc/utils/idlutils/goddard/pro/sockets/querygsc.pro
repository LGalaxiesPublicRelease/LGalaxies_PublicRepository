function Querygsc, target, dis,magrange = magrange, HOURS = hours
;+
; NAME: 
;   QUERYGSC
;
; PURPOSE: 
;   Query the Guide Star Catalog (GSC V2.2) at STScI by position
; 
; EXPLANATION:
;   Uses the IDL SOCKET command to query the GSC 2.2 database over the Web.    
;   Requires IDL V5.4 or later.
; 
;   The GSC  all-sky catalog will contain an estimated 2 billion objects
;   and will be complete to a magnitude of at least J=18 and as faint as J=21 at
;   high galactic latitudes. Using the observations in different bandpasses at 
;   different epochs allows the computation of both colors and proper motions. 
;   These data are in an object-oriented database at 
;   http://www-gsss.stsci.edu/support/data_access.htm.  The final version 
;   (GSC 2.3),  expected to be released in 2004, will also contain proper 
;   motions.   
;
; CALLING SEQUENCE: 
;     info = QueryGSC(targetname_or_coords, [ dis, Magrange =, /HOURS] )
;
; INPUTS: 
;      TARGETNAME_OR_COORDS - Either a scalar string giving a target name, 
;          (with J2000 coordinates determined by SIMBAD), or a 2-element
;          numeric vector giving the J2000 right ascension in *degrees* (or 
;          hours if /HOURS is set) and the target declination in degrees.
;
; OPTIONAL INPUT:
;    dis - Search radius in arcminutes to search around specified target
;          Default is 5 arcminutes
;
; OPTIONAL INPUT KEYWORDS:
;
;    /HOURS - If set, then the right ascension is both input and output (in
;             the info .ra tag) in hours instead of degrees
;
;    MAGRANGE - two element vector giving the magnitude range (on either the
;           F plate or the J plate) to search for  GSC stars.   
;           Default is [0,30]
;
; OUTPUTS: 
;   info - IDL structure containing information on the GSC stars within the 
;          specified distance of the specified center.   There are (currently)
;          23 tags in this structure  -- for further information see
;           http://www-gsss.stsci.edu/gsc/gsc2/gsc22_release_notes.htm  
;
;          .GSCID2 - GSC 2.2 identification number
;          .RA,.DEC - Position in degrees (double precision).   RA is given in
;                   hours if the /HOURS keyword is set.
;          .RAERR, .DECERR - uncertainty (in arcseconds) in the RA and Dec
;          .EPOCH - mean epoch of the observation
;          .RAPM,DECPM - RA and Dec proper motion (mas/year) 
;          .RAPMERR,DECPMERR - Uncertainty RA and Dec proper motion (mas/year) 
;          .FMAG, .FMAGERR - magnitude and error in photographic F
;          .JMAG, .JMAGERR - magnitude and error in photographic J
;          .VMAG, .VMAGERR - V magnitude and error 
;          .NMAG, .NMAGERR - magnitude and error
;          .A - semi-major axis in pixels
;          .E - eccentricity of extended objects
;          .PA - Position angle of extended objects in degrees
;          .C - classification (0-5): 0-star, 1-galaxy, 2-blend, 3-nonstar,
;                                     4-unclassified, 5-defect
;          .STATUS -10 digit field  used to encode more detailed information 
;              about the properties of the catalog object.   For more info, see 
;http://www-gsss.stsci.edu/gsc/gsc2/gsc22_release_notes.htm#SourceStatusFlagCodes
;
; EXAMPLE: 
;          Plot a histogram of the photographic J magnitudes of all GSC 2.2 
;          stars within 10 arcminutes of the center of the globular cluster M13 
;
;          IDL> info = querygsc('M13',10)
;          IDL> plothist,info.jmag,xran=[10,20]
;
; PROCEDURES USED:
;          QUERYSIMBAD, RADEC, WEBGET()
;
; MINIMUM IDL VERSION
;         V5.4  (uses SOCKET)
; MODIFICATION HISTORY: 
;         Written by W. Landsman  SSAI  August 2002
;         Fixed parsing of RA and Dec  W. Landsman September 2002
;
;-
  if N_params() LT 2 then begin
       print,'Syntax - info = QueryGSC(targetname_or_coord, dis,'
       print,'                        [/Hours, MagRange=[m1,m2]} )'
       print,'   RA (degrees), Dec (degrees) -- search coordinates of center)'
       print,'  dis -- search radius in arcminutes'
       if N_elements(info) GT 0 then return,info else return, -1
  endif
  if N_elements(dis) EQ 0 then dis = 5
    if N_elements(target) EQ 2 then begin
      ra = float(target[0])
      dec = float(target[1])
  endif else begin
       QuerySimbad, target, ra,dec, Found = Found
       if found EQ 0 then message,'Target name ' + target + $
                 ' could not be translated by SIMBAD'
  endelse  
   radec,ra,dec,hr,mn,sc,deg,dmn,dsc,hours=keyword_set(hours)
  ;;
  QueryURL = "http://www-gsss.stsci.edu/cgi-bin/gsc22query.exe?ra=" + $
  string(hr,'(i2.2)') +'%3A' + string(mn,'(i2.2)') + '%3A'  + strtrim(sc,2) + $          
  '&dec=' + strtrim(string(deg,'(i3.2)'),2) +'%3A' + $
   string(dmn,'(i2.2)') + '%3A' + strtrim(dsc,2) + '&r2=' + strtrim(dis,2) 
  if N_elements(magrange) EQ 2 then begin
  QueryURL = QueryURL + '&m1=' + strtrim(magrange[0],2) + $
                        '&m2=' + strtrim(magrange[1],2)
  endif
  QueryURL = QueryURL + '&submit2=Submit+Request'         
  ;;  
  Result = webget(QueryURL)
;
  t = result.text

  nstar = N_elements(t) -3
  if strmid(t[0],0,5) NE 'Usage' and nstar GT 0 THEN BEGIN
  headers = strsplit(t[0],/extract)
  info = create_struct(Name='gsc',headers, '',0.0d,0.0d,0.0,0.0,0.0, $
   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0ULL)
  info = replicate(info,nstar)
  for i=0,nstar-1 do begin
      temp = strsplit(t[i+2],/extract)
      for j=0,22 do info[i].(j) = temp[j]
  endfor
   ENDIF ELSE BEGIN 
      message, 'No objects returned by server. The server answered:', /info
      print, Result.Text
      if N_elements(info) GT 0 then return, info else return, -1
  ENDELSE 
  if keyword_set(hours) then info.ra = info.ra/15.0d 
 return,info
END 
  
