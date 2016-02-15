PRO QueryDSS, target, Image,  Header, IMSIZE=ImSIze, ESO=eso, STSCI=stsci, $
              NED=ned, SURVEY = survey, OUTFILE = outfile
;+
; NAME: 
;   QueryDSS
;
; PURPOSE: 
;    Query the digital sky survey (DSS) on-line at  the ESO or STSCI servers
;
; EXPLANATION: 
;     The script can query the DSS survey and retrieve an image and FITS 
;     header either from the European Southern Observatory (ESO) or the 
;     Space Telescope Science Institute (STScI) servers.
;     See http://archive.eso.org/dss/dss and/or 
;     http://archive.stsci.edu/dss/index.html for details.
;
; CALLING SEQUENCE: 
;      QueryDSS, targetname_or_coords, Im, Hdr, [IMSIZE= , /ESO, /STSCI ]
;
; INPUTS:
;      TARGETNAME_OR_COORDS - Either a scalar string giving a target name, 
;          (with J2000 coordinates determined by SIMBAD (default) or NED), or 
;          a 2-element numeric vector giving the J2000 right ascension in 
;          *degrees* and the target declination in degrees.
;
; OPTIONAL INPUTS: 
;          None
;
; OPTIONAL KEYWORD PARAMETERS: 
;     ImSize - Numeric scalar giving size of the image to be retrieved in 
;                 arcminutes.    Default is 10 arcminute.
;
;     /ESO - Use the ESO server for image retrieval.    Default is to use
;            the STScI server if user is in the Western hemisphere, and 
;            otherwise to use the ESO server.
;
;     /NED - Query the Nasa Extragalactic Database (NED) for the
;            target's coordinates.  The default is to use Simbad for
;            the target search.
;
;     /STSCI - Use the STSCI server for image retrieval.  Default is to use
;            the STScI server if user is in the Western hemisphere, and 
;            otherwise to use the ESO server.    
;
;     SURVEY - Scalar string specifying which survey to retrieve.  
;          Possible values are 
;          '1' - First generation (red), this is the default
;          '2b' - Second generation blue
;          '2r' - Second generation red
;          '2i' - Second generation near-infrared
; 
;      Note that 2nd generation images may not be available for all regions
;      of the sky.   Also note that the first two letters of the 'REGION'
;      keyword in the FITS header gives the bandpass 'XP' - Red IIIaF, 
;      'XJ' - Blue IIIaJ, 'XF' - Near-IR IVN
;
; OUTPUTS: 
;       Im - The image returned by the server. If there is an error, this 
;             contains a single 0.
;
;       Hdr - The FITS header of the image. Empty string in case of errors.
;
;       If the OutFile keyword is set then no outputs are returned (only the
;       file is written).
; SIDE EFFECTS: 
;     If Im and Hdr exist in advance,  they are overwritten.
;
; RESTRICTIONS: 
;      Relies on a working network connection. 
;
; PROCEDURE: 
;      Construct a query-url,  call WEBGET() and sort out the server's 
;      answer.
;
; EXAMPLE:           
;      Retrieve an 10'  image surrounding the ultracompact HII region
;       G45.45+0.06.   Obtain the 2nd generation blue image.
;
;       IDL> QueryDSS, 'GAL045.45+00.06', image, header, survey = '2b'
;       IDL> tvscl, image
;       IDL> hprint, header
;       IDL> writefits,'dss_image.fits', image, header
; Note that the coordinates could have been specified directly, rather than
; giving the target name.
;       IDL> QueryDSS, [288.587, 11.1510], image, header,survey='2b'
;
; To write a file directly to disk, use the OutFile keyword
;
;       IDL> QueryDSS, [288.587, 11.1510], survey='2b', out='gal045_2b.fits'
;   
; PROCEDURES CALLED:
;       QUERYSIMBAD, WEBGET()
; MINIMUM IDL VERSION:
;       V5.4  (uses SOCKET)
; MODIFICATION HISTORY: 
;       Written by M. Feldt, Heidelberg, Oct 2001 <mfeldt@mpia.de>
;       Option to supply target name instead of coords  W. Landsman Aug. 2002
;       Added OUTFILE, /NED keywords W. Landsman   April 2003
;       Don't abort on Simbad failure W. landsman/J. Brauher  June 2003
;
;-
  if N_params() LT 1 then begin
      print,'Syntax - QueryDSS, TargetName_or_coords, image, header'
      print,"           [Imsize= ,/ESO, /STScI, Survey = ['1','2b','2r','2i'] "
      print,'            /NED, OutFile = ]'
      return
   endif
;; Is the user in the Western Hemisphere?
   if not keyword_set(ESO) and not keyword_set(STScI) then begin
       timezone =  (systime(/JULIAN,/UTC)- systime(/julian))*24
       stsci = (timezone GE 4) and (timezone LE 10)
   endif
  ;;
  if N_elements(target) EQ 2 then begin
      ra = float(target[0])
      dec = float(target[1])
  endif else begin
       QuerySimbad, target, ra,dec, NED= ned, Found = Found
       if found EQ 0 then begin 
             message,/inf,'Target name ' + target + $
                 ' could not be translated by SIMBAD'
             return
       endif
  endelse  
  if not keyword_set(ESO) then eso =  1b-keyword_Set(stsci) 
  IF NOT Keyword_Set(ImSize) THEN ImSize = 10
  Equinox = 'J2000'
  ;;
  ;;
 if N_elements(survey) EQ 0 then survey = '1'
 dss = strlowcase(strtrim(strmid(survey,0,2),2))
 if ESO then begin
  case dss of 
  '1': dss = 'DSS1'
  '2b': dss = 'DSS2-blue'
  '2r': dss = 'DSS2-red'
  '2i': dss = 'DSS2-infrared'
  else: message,'Unrecognized Survey - should be 1, 2b, 2r or 2i'
 endcase
 endif

  IF eso THEN $ 
    QueryURL=strcompress("http://archive.eso.org/dss/dss/image?ra="+$
                       string(RA)+$
                       "&dec="+$
                       string(DEC)+$
                       "&x="+$
                       string(ImSize)+$
                       "&y="+$
                       string(ImSize)+$
                       "&Sky-Survey="+dss +"&mime-type=download-fits", /remove) $
  ELSE $
    QueryURL=strcompress("http://archive.stsci.edu/cgi-bin/dss_search?ra="+$
                         string(RA)+$
                         "& dec="+$
                         string(DEC)+$
                         "& equinox="+$
                         Equinox +$
                         "& height="+$
                         string(ImSize) +$
                         "&generation=" + dss +$                       
                         "& width="+$
                         string(ImSize)+$
                         "& format=FITS", /remove)
  ;;

  if keyword_set(OutFile) then begin
      if not ESO then dss = 'DSS' + dss
      message,'Writing ' + dss + ' FITS file ' + outfile,/inf
      Result = webget(QueryURL, copyfile= outfile)
      return
  endif
  Result = webget(QueryURL)
  Image = Result.Image
  Header = Result.ImageHeader
  ;;
  ;; error ?
  ;;
  IF N_Elements(Image) NE 1 THEN return
  message, 'Problem retrieving your image! The server answered:', /info
  print, Result.Text
END 
