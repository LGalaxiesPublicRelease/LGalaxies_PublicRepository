function Queryvizier, catalog, target, dis, VERBOSE=verbose, CANADA = canada, $
               CONSTRAINT = constraint
;+
; NAME: 
;   QUERYVIZIER
;
; PURPOSE: 
;   Query any catalog in the Vizier database by position
; 
; EXPLANATION:
;   Uses the IDL SOCKET command to provide a positional query of any catalog 
;   in the the Vizier (http://vizier.u-strasbg.fr/) database over the Web.    
;   Requires IDL V5.4 or later.
; 
;    
; CALLING SEQUENCE: 
;     info = QueryVizier(catalog, targetname_or_coords, [ dis )
;
; INPUTS: 
;      CATALOG - Scalar string giving the name of the VIZIER catalog to be
;            searched.    The complete list of catalog names is available at
;            http://vizier.u-strasbg.fr/vizier/cats/U.htx . 
;
;            Popular VIZIER catalogs include  
;            '2MASS-PSC' - 2MASS point source catalog (2003)
;            'GSC2.2' - Version 2.2 of the HST Guide Star Catalog (2003)
;            'USNO-B1' - Verson B1 of the US Naval Observatory catalog (2003)
;            'NVSS'  - NRAO VLA Sky Survey (1998)
;            'B/DENIS/DENIS' - 2nd Deep Near Infrared Survey of southern Sky
;            'I/259/TYC2' - Tycho-2 main catalog (2000)
;            'I/239/HIP_MAIN' - Hipparcos main catalog (1997)
;
;          Note that some names will prompt a search of multiple catalogs
;          and QUERYVIZIER will only return the result of the first search.
;          Thus, setting catalog to "HIPPARCOS" will search all catalogs 
;          associated with the Hipparcos mission, and return results for the
;          first catalog found.    To specifically search the Hipparcos or
;          Tycho main catalogs use the VIZIER catalog names listed above
;                             
;      TARGETNAME_OR_COORDS - Either a scalar string giving a target name, 
;          (with J2000 coordinates determined by SIMBAD), or a 2-element
;          numeric vector giving the J2000 right ascension in *degrees* and 
;          the target declination in degrees.
;
; OPTIONAL INPUT:
;    dis - scalar or 2-element vector.   If one value is supplied then this
;          is the search radius in arcminutes.     If two values are supplied
;          then this is the width (i.e., in longitude direction) and height
;          of the search box.   Default is a radius search with radius of
;          5 arcminutes
;
; OUTPUTS: 
;   info - Named IDL structure containing information on the catalog stars 
;          within the specified distance of the specified center.   The
;          structure name is taken from the searched VIZIER catalog (which may
;          differ from the input 'CATALOG' parameter if multiple catalogs are 
;          searched).  The structure tag names are identical with the VIZIER 
;          catalog column names, with the exception of an occasional underscore
;          addition, if necessary to convert the column name to a valid 
;          structure tag.    The VIZIER Web  page should consulted for the 
;          column names and their meaning for each particular catalog..
;           
;          If the tagname is numeric and the catalog field is blank then either
;          NaN  (if floating) or -1 (if integer) is placed in the tag.
;
;          If no sources are found within the specified radius, or an
;          error occurs in the query then -1 is returned. 
; OPTIONAL KEYWORDS:
;          /CANADA - By default, the query is sent to the main VIZIER site in
;            Strasbourg, France.   If /CANADA is set then the VIZIER site
;            at the Canadian Astronomical Data Center (CADC) is used instead.
;            Note that not all Vizier sites have the option to return
;            tab-separated values (TSV) which is required by this program.
;
;          CONSTRAINT - string giving additional nonpositional numeric 
;            constraints on the entries to be selected.     For example, when 
;            in the GSC2.2  catalog, to only select sources with Rmag < 16 set 
;            Constraint = 'Rmag < 16'.    Multiple constraints can be 
;            separated by commas.    Use '!=' for "not equal", '<=' for smaller
;            or equal, ">=" for greater than or equal.  See the complete list
;            of operators at  
;                 http://vizier.u-strasbg.fr/doc/asu.html#AnnexQual
;            For this keyword only, **THE COLUMN NAME IS CASE SENSITIVE** and 
;            must be written exactly as displayed on the VIZIER Web page.  
;            Thus for the GSC2.2 catalog one must use 'Rmag' and not 'rmag' or
;            'RMAG'.
;         
;           /VERBOSE - If set then the query sent to the VIZIER site is
;               displayed, along with the returned title(s) of found catalog(s)
; EXAMPLES: 
;          (1) Plot a histogram of the J magnitudes of all 2MASS point sources 
;          stars within 10 arcminutes of the center of the globular cluster M13 
;
;          IDL>  info = queryvizier('2MASS-PSC','m13',10)
;          IDL> plothist,info.jmag,xran=[10,20]
;
;          (2)  Find the brightest R magnitude GSC2.2 source within 3' of the 
;               J2000 position ra = 10:12:34, dec = -23:34:35
;          
;          IDL> str = queryvizier('GSC2.2',[ten(10,12,34)*15,ten(-23,34,35)],3)
;          IDL> print,min(str.rmag,/NAN)
;
;          (3) Find sources with V < 19 in the Magellanic Clouds Photometric 
;              Survey (Zaritsky+, 2002) within 5 arc minutes of  the position 
;              00:47:34 -73:06:27
;
;              Checking the VIZIER Web page we find that this catalog is
;          IDL>  catname =  'J/AJ/123/855/table1'
;          IDL>  ra = ten(0,47,34)*15 & dec = ten(-73,6,27)
;          IDL> str = queryvizier(catname, [ra,dec], 5, constra='Vmag<19')
;
; PROCEDURES USED:
;          GETTOK(),IDL_VALIDNAME()(if prior to V6.0), RADEC, REMCHAR, REPSTR(),
;          WEBGET()
; TO DO:
;       (1) Allow specification of output sorting
; MODIFICATION HISTORY: 
;         Written by W. Landsman  SSAI  October 2003
;         Give structure name returned by VIZIER not that given by user
;                    W. Landsman   Feburary 2004 
;         Don't assume same format for all found sources W. L. March 2004
;         Added CONSTRAINT keyword for non-positional constraints WL July 2004
;         Remove use of EXECUTE() statement WL June 2005
;-
  if N_params() LT 3 then begin
       print,'Syntax - info = QueryVizier(catalog, targetname_or_coord, dis,'
       print,'                            [ /VERBOSE, /CANADA, CONSTRAINT= ]'
       print,'                       '
       print,'  Coordinates (if supplied) should be J2000 RA (degrees) and Dec'
       print,'  dis -- search radius or box in arcminutes'
       if N_elements(info) GT 0 then return,info else return, -1
  endif

if keyword_set(canada) then root = "http://vizier.hia.nrc.ca" $
                       else  root = "http://vizier.u-strasbg.fr" 
 
  if N_elements(catalog) EQ 0 then $
            message,'ERROR - A catalog name must be supplied as a keyword'
  targname = 0b
 if N_elements(dis) EQ 0 then dis = 5
 if N_elements(dis) EQ 2 then $
    search = "&-c.bm=" + strtrim(dis[0],2) + '/' + strtrim(dis[1],2) else $
    search = "&-c.rm=" + strtrim(dis,2) 
    if N_elements(target) EQ 2 then begin
      ra = float(target[0])
      dec = float(target[1])
      RADEC,ra,dec,hr,mn,sc,deg,dmn,dsc,hours=keyword_set(hours)
   endif else begin 
       object = repstr(target,'+','%2B')
        object = repstr(strcompress(object),' ','+')
       targname = 1b 
  endelse
 source = catalog

; Add any addition constraints to the search.    Convert an URL special 
; special characters in the constraint string.

 if N_elements(constraint) EQ 0 then constraint = '' 
 if strlen(constraint) GT 0 then begin
     urlconstrain = strtrim(constraint,2)
     urlconstrain = repstr(urlconstrain, ',','&')
     urlconstrain = repstr(urlconstrain, '<','=%3C')
     urlconstrain = repstr(urlconstrain, '>','=%3E')
     urlconstrain = repstr(urlconstrain, '!','=!')
     search = search + '&' + urlconstrain
 endif
 ;
 if targname then $
  QueryURL = $
   root + "/cgi-bin/asu-tsv/?-source=" + source + $
     "&-c=" + object + search + '&-out.max=unlimited' else $
  queryURL = $
   root + "/cgi-bin/asu-tsv/?-source=" + source + $
       "&-c.ra=" + strtrim(ra,2) + '&-c.dec=' + strtrim(dec,2) + $
       search + '&-out.max=unlimited'
   ;;  
 if keyword_set(verbose) then message,queryurl,/inf

  Result = webget(QueryURL)

;
  t = strtrim(result.text,2)
  keyword = strtrim(strmid(t,0,7),2)

  lcol = where(keyword EQ "#Column", Nfound)
  if Nfound EQ 0 then begin
       if max(strpos(strlowcase(t),'errors')) GE 0 then begin 
            message,'ERROR - Unsuccessful VIZIER query',/CON 
            print,t
       endif else message,'No sources found within specified radius',/INF
       return,-1
  endif
; Check to see if more than one catalog has been searched
  if N_elements(lcol) GT 1 then begin
      diff  = lcol[1:*] - lcol
      g = where(diff GT 1, Ng)
      if Ng GT 0 then lcol = lcol[0:g[0]]
  endif

  if keyword_set(verbose) then begin
    titcol = where(keyword EQ '#Title:', Ntit)    
    if Ntit GT 0 then message,/inform, $
        strtrim(strmid(t[titcol[0]],8),2)
  endif

  trow = t[lcol]
  dum = gettok(trow,' ')
  colname = gettok(trow,' ')
  fmt = gettok(trow,' ')
  remchar,fmt,'('
  remchar,fmt,')'
  remchar,colname,')'
  val = fix(strmid(fmt,1,4))
  for i=0,N_elements(colname)-1 do $
         colname[i] = IDL_VALIDNAME(colname[i],/convert_all)

;Check if any Warnings or fatal errors in the VIZIER output
   badflag = strmid(keyword,0,5)
   warn = where(badflag EQ '#++++', Nwarn)
   if Nwarn GT 0 then for i=0,Nwarn-1 do $
        message,'Warning: ' + strtrim(t[warn[i]],2),/info
   
   fatal = where(badflag EQ '#****', Nfatal)
   if Nfatal GT 0 then for i=0,Nfatal-1 do $
        message,'Error: ' + strtrim(t[fatal[i]],2),/info


   table = where(keyword EQ '#Table', Ntable)
   if Ntable GT 0 then begin
          source = strmid(t(table[0]),7)
          remchar,source,':'
   endif 

  source = IDL_VALIDNAME(source, /convert_all)

 ntag = N_elements(colname)
 for i=0,Ntag-1 do begin
 case strmid(fmt[i],0,1) of 
  'A': cval = ' '
  'I': if val[i] LE 4 then cval = 0 else $
                           cval = 0L
  'F': if val[i] LE 7 then  cval = 0. else $ 
                           cval = 0.0d
   else: message,'ERROR - unrecognized format'
 
  endcase
   if i EQ 0 then   info = create_struct(colname[0], cval) else $
   info =  create_struct(temporary(info), colname[i],cval)
 endfor
   info = create_struct(temporary(info), name=source)


  i0 = max(lcol) + 5  
  iend = where( t[i0:*] EQ '', Nend)
  if Nend EQ  0  then iend = N_elements(t) else iend = iend[0] + i0
  nstar = iend - i0 
    info = replicate(info, nstar)

; Find positions of tab characters 
  t = t[i0:iend-1]

  for j=0,Ntag-1 do begin
      x = strtrim( gettok(t,string(9b),/exact ),2)
       dtype = size(info[0].(j),/type)
      if dtype NE 7 then begin
      bad = where(strlen(x) EQ 0, Nbad)
     if (Nbad GT 0) then $
           if (dtype EQ 4) or (dtype EQ 5) then x[bad] = 'NaN' $
                                           else x[bad] = -1
      endif
      info.(j) = x 
   endfor

 return,info
END 
  
