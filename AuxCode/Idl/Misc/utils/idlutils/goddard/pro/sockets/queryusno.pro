function Queryusno, target, dis,magrange = magrange, HOURS = hours
;+
; NAME: 
;   QUERYUSNO
;
; PURPOSE: 
;   Query the USNO-A2.0 Catalog at the ESO/ST-ECF Archive by position
; 
; EXPLANATION:
;   Uses the IDL SOCKET command to query the USNO-A2.0 database over the Web.    
;   Requires IDL V5.4 or later.
; 
;   With the introduction of QUERYVIZIER this routine became mostly obsolete
;   as the newer USNO-B1 catalog can be accessed from QUERYVIZIER.
;
;   USNO-A2.0 contains entries for over a half billion stars (526,230,881, to 
;   be exact!) which were detected in the digitized images of three photographic
;   sky surveys. For the entire northern sky and the southern sky down to 
;   declinations of -30°, all the photographic plates were part of the original
;   Palomar Optical Sky Survey (POSS-I).  Photographs were taken on blue- and 
;   red-sensitive emulsions. Only those stars which were detected in both colors
;   were included in the USNO-A2.0 catalog. The rest of the southern sky was 
;   covered by the Science Research Council (SRC)-J survey and the European 
;   Southern Observatory (ESO)-R survey.  Only stars appearing in both 
;   colors were accepted for the final catalogue.   Coordinates are J2000 
;   at the epoch of the mean of the blue and red exposure. 
;
; CALLING SEQUENCE: 
;     info = QueryUSNO(targetname_or_coords, [ dis, Magrange =, /HOURS] )
;
; INPUTS: 
;      TARGETNAME_OR_COORDS - Either a scalar string giving a target name, 
;          (with J2000 coordinates determined by SIMBAD), or a 2-element
;          numeric vector giving the J2000 right ascension in *degrees* and
;          the target declination in degrees.
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
;    MAGRANGE - two element vector giving the magnitude range (on the
;           r plate) to search for  stars.   Default is to return all stars
;
; OUTPUTS: 
;   info - IDL structure containing information on the USNO-A2 stars within the 
;          specified distance of the specified center.   There are (currently)
;          5 tags in this structure  -- for further information see
;           http://ftp.nofs.navy.mil/projects/pmm/readme.html
;
;          .ID - USNO-A2.0 identification number
;          .RA,.DEC - Position in degrees (double precision).   RA is given in
;                   hours if the /HOURS keyword is set.
;          .r_mag, .b_mag - magnitudes on the red and blue plates
;
; EXAMPLE: 
;          Plot a histogram of the photographic r magnitudes of all USNO-A2 
;          stars within 10 arcminutes of the center of the globular cluster M13 
;
;          IDL> info = queryusno('M13',10)
;          IDL> plothist,info.r_mag,xran=[10,20]
;
; PROCEDURES USED:
;          QuerySIMBAD, RADEC, WEBGET()
;
; MODIFICATION HISTORY: 
;         Written by W. Landsman  SSAI  September 2002
;
;-
  if N_params() LT 2 then begin
       print,'Syntax - info = QueryUSNO(targetname_or_coord, dis,'
       print,'                        [/Hours, MagRange=[m1,m2]} )'
       print,'   RA (degrees), Dec (degrees) -- search coordinates of center)'
       print,'  dis -- search radius in arcminutes'
       return,info
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
  QueryURL = "http://archive.eso.org/" + $
      'skycat/servers/usnoa_res?catalogue=usnoa&epoch=2000.0&chart=0&ra=' +$
  string(hr,'(i2.2)')  + '+' + string(mn,'(i2.2)')  + '+' + strtrim(sc,2) + $          
  '&dec=' + strtrim(string(deg,'(i3.2)'),2) + '+' + $
   string(dmn,'(i2.2)') + '+' + strtrim(dsc,2) + '&radmax=' + strtrim(dis,2) 
  if N_elements(magrange) EQ 2 then begin
  QueryURL = QueryURL + '&magbright=' + strtrim(magrange[0],2) + $
                        '&magfaint=' + strtrim(magrange[1],2)
  endif

 QueryURL = QueryURL + '&format=2&sort=d%27'         
  ;;  print
  Result = webget(QueryURL)
;

  t = result.text
  N = N_elements(t)

  if N GT 51 then begin 

     t = t[28:N-25]
     nstar = N_elements(t)
     t[0] = strmid(t[0],4,80)
   info = create_struct(Name='usno','ID',' ','ra',0.0d,'dec',0.0d,'r_mag',0.0, $
   'b_mag',0.0)
  info = replicate(info,nstar)
  for i=0,nstar-1 do begin
      temp = strsplit(t[i],/extract)
      info[i].id = temp[1]
      info[i].ra = temp[2]
      info[i].dec = temp[3]
      info[i].r_mag = temp[4]
      info[i].b_mag = temp[5]
  endfor
  if keyword_set(hours) then info.ra = info.ra/15.0d 
   ENDIF ELSE BEGIN 
      message, 'No objects returned by server', /info
      if N_elements(info) GT 0 then return, info else return,-1
  ENDELSE 
 return,info
END 
  
