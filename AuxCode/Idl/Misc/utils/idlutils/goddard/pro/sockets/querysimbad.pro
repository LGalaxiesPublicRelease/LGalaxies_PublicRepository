PRO QuerySimbad, name, ra, de, id, Found = found, NED = ned, ERRMSG = errmsg
;+
; NAME: 
;   QUERYSIMBAD
;
; PURPOSE: 
;   Query the SIMBAD or NED astronomical name resolver to obtain coordinates
;
; EXPLANATION: 
;   Uses the IDL SOCKET command to query either the SIMBAD or NED nameserver 
;   over the Web to return J2000 coordinates.     Requires IDL V5.4 or later.
;
;   For details on the SIMBAD service, see http://simbad.u-strasbg.fr/Simbad 
;   and for the NED service, see http://ned.ipac.caltech.edu/
;
; CALLING SEQUENCE: 
;    QuerySimbad, name, ra, dec, [ id, Found=, /NED, ERRMSG=]
;
; INPUTS: 
;    name - a scalar string containing the target name in SIMBAD (or NED)
;           nomenclature. For SIMBAD details see
;           http://vizier.u-strasbg.fr/cgi-bin/Dic-Simbad .
;
; OUTPUTS: 
;     ra - the right ascension of the target in J2000.0 in *degrees* 
;     dec - declination of the target in degrees
;
; OPTIONAL INPUT KEYWORD:
;     /NED - if set, then nameserver of the NASA Extragalactic database is
;           used to resolve the name and return coordinates.   Note that
;           /NED cannot be used with Galactic objects
; OPTIONAL OUTPUT: 
;     id - the primary SIMBAD (or NED) ID of the target, scalar string
;
; OPTIONAL KEYWORD OUTPUT:
;     found - set to 1 if the translation was successful, or to 0 if the
;           the object name could not be translated by SIMBAD or NED
;     Errmsg - if supplied, then any error messages are returned in this
;            keyword, rather than being printed at the terminal.   May be either
;            a scalar or array.
;            
; EXAMPLES:
;     (1) Find the J2000 coordinates for the ultracompact HII region
;         G45.45+0.06 
;
;      IDL> QuerySimbad,'GAL045.45+00.06', ra, dec
;      IDL> print, adstring(ra,dec,1)
;           ===>19 14 20.77  +11 09  3.6
; PROCEDURES USED:
;       REPSTR(), WEBGET()
;
; MINIMUM IDL VERSION:
;    V5.4 (uses SOCKET)
; MODIFICATION HISTORY: 
;     Written by M. Feldt, Heidelberg, Oct 2001   <mfeldt@mpia.de>
;     Minor updates, W. Landsman   August 2002
;     Added option to use NED server, better parsing of SIMBAD names such as 
;          IRAS F10190+5349    W. Landsman  March 2003
;     Turn off extended name search for NED server, fix negative declination
;     with /NED    W. Landsman  April 2003
;
;-
  if N_params() LT 3 then begin
       print,'Syntax - QuerySimbad, name, ra, dec, [ id, Found=, /NED, ERRMSG=]'
       print,'   Input - object name, scalar string'
       print,'   Output -  Ra, dec of object'
       return
  endif
  ;;
  printerr = not arg_present(errmsg)
  object = repstr(name,'+','%2B')
    if keyword_set(NED) then begin
 object = repstr(strcompress(object),' ','+')
 QueryURL = "http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?objname=" + $
              strtrim(object,2) + '&img_stamp=NO&list_limit=0&extend=no' 
  endif else begin
 object = repstr(strcompress(object),' ','%20')

  QueryURL = "http://archive.eso.org/skycat/servers/sim-server?&o=" + $
              strcompress(object,/remove)
  endelse
  ;;

  Result = webget(QueryURL)
  found = 0
  ;;
  if keyword_set(NED) THEN BEGIN
      if (strmid(result.text[5],1,3) NE 'PRE') and $
         (N_Elements(result.text) GE 16) THEN BEGIN
            found = 1
             t = result.text[15]
            headpos = strpos(t,'A>')
            t = strmid(t,headpos+2,80)
             id = strtrim( strmid(t,0,32),2)
            hr = fix(strmid(t,33,2))
            mn = fix(strmid(t,36,2))
            sc = float(strmid(t,39,4))
            ra = ten(hr,mn,sc)*15.0d
            dsgn = strmid(t,45,1)
            deg = fix(strmid(t,46,2))
            dmn = fix(strmid(t,49,2))            
            dsc = fix(strmid(t,52,2))
            de = ten(deg,dmn,dsc)
            if dsgn EQ '-' then de = -de

      endif else begin

      errmsg = 'No objects returned by NED server'
      if printerr then  message, errmsg, /info
     endelse
 endif else begin

  IF strmid(Result.Text[0], 0, 2) EQ 'Id' THEN BEGIN 
      found = 1
      ;;
      ;; prepare the result fields
      ;;
      ra = dblarr(N_Elements(Result.Text)-3)
      de =  ra
      id = strarr(N_Elements(Result.Text)-3) 
      ;;
      ;; decode the result
      ;;
      FOR ii=2, N_Elements(Result.Text)-2 DO BEGIN
          TheseFields = strsplit(Result.Text[ii], string(9B), /extract)
          id[ii-2] = TheseFields[0]
          ra[ii-2] = float(TheseFields[1])
          de[ii-2] = float(TheseFields[2])
      ENDFOR 
      ;;
      ;; ready for return
      ;;
      IF N_Elements(Result.Text) EQ 4 THEN BEGIN 
          ra = ra[0] ; do not return single-element arrays
          de = de[0]
          id = id[0]
      ENDIF 
      return
  ENDIF ELSE BEGIN 
      errmsg = ['No objects returned by SIMBAD.   The server answered:', $
                 result.text]
      if printerr then begin
         message, errmsg[0], /info
         print, Result.Text
     endif
  ENDELSE
  ENDELSE 
END 
  
