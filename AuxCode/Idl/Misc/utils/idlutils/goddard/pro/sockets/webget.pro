;+
; NAME: 
;    WEBGET()
;
; PURPOSE: 
;    Use the IDL SOCKET procedure to get data from http servers
;
; EXPLANATION: 
;     WEBGET() can access http servers - even from behind a firewall - 
;     and perform simple downloads. Currently, text and FITS files can be 
;     accessed.    Requires IDL V5.4 or later on Unix or Windows, V5.6 on
;     Macintosh
;
; CALLING SEQUENCE: 
;      a=webget(URL)
;
; INPUTS: 
;      URL - scalar string giving a fully qualified url of the form
;          'http://server.eso.org/path/file.html'.    WEBGET() can
;          also use other valid URLs that contain 'GET'-codes.
;
; OPTIONAL INPUT KEYWORD PARAMETERS: 
;       COPYFILE - if set to a valid filename (file must have write permission),
;            the data contents of the web server's answer is copied to that 
;            file. 
;       /SILENT - If set, the information error messages are suppressed
; OUTPUTS: A structure with the following fields:
;
;            .Header - the HTTP header sent by the server
;
;            .Text   - The text part of the downloaded file. If the
;                     content type of the file was not of class
;                     'text',  this will be an empty string.
;
;            .ImageHeader - Header file of a FITS-image. FITS images
;                          are read when the content type is
;                          'image/fits' or 'application/octet-stream'
;                          (for dss-access). If the file is not a FITS
;                          image,  this will be an empty string.
;
;            .Image - The FITS image read from the server. If the file
;                    did not contain a FITS image,  this will be zero.
;
;
; RESTRICTIONS: 
;     The mime-type recognition is extremely limited. Only the content-type is 
;     determined. Any text-file  will be stored in out.Text. The only other 
;     category which can be fetched is FITS files,  which will be stored in 
;     out.Image and out.ImageHeader.
;
;     PROXY: If you are behind a firewall and have to access the net through a 
;         Web proxy,  set the environment variable 'http_proxy' to point to 
;         your proxy server and port, e.g. 
;         'setenv http_proxy=http://web-proxy.mpia-hd.mpg.de:3128'
;
;               The URL *MUST* begin with "http://".
;
; PROCEDURE: 
;     Open a socket to the webserver and download the header. After deciding 
;     whether it is text or binary, either store the text or try to read a 
;     FITS file.
;
; EXAMPLE: 
;      IDL> a=webget('http://www.mpia.de/index.html')
;      IDL> print,a.Text
;      or
;
;          > PointingRA=0.0
;          > PointingDE=30.0
;          > QueryURL = strcompress("http://archive.eso.org/dss/dss/image?ra="+$
;          >                          string(PointingRA)+$
;          >                          "&dec="+$
;          >                          string(PointingDE)+$
;          >                          "&x=10&y=10&Sky-Survey=DSS1&mime-type=download-fits", $
;          >                          /remove)
;          > a=webget(QueryURL)
;          > tvscl,a.Image
;          > print,a.ImageHead
;
;
; MINIMUM IDL VERSION:
;     V5.4  (uses SOCKET)
; MODIFICATION HISTORY: 
;     Written by M. Feldt, Heidelberg, Oct 2001 <mfeldt@mpia.de>
;     Use /swap_if_little_endian keyword to SOCKET  W. Landsman August 2002
;     Less restrictive search on Content-Type   W. Landsman   April 2003
;
;-

PRO MimeType,  Header, Class, Type, Length
;;
;; MIME type recognition
;
  Class = 'text'
  Type = 'simple'               ; in case no information found...    
  def = strupcase(strmid(header,0,13))
  g = where(def EQ 'CONTENT-TYPE:', Ng)
  if Ng GT 0 then begin
       ClassAndType = strmid(Header[g[0]], 14, strlen(Header[g[0]])-1)
       Class = (strsplit(ClassAndType, '/', /extract))[0]
       Type = (strsplit(ClassAndType, '/', /extract))[1]
  ENDIF 
  def = strupcase(strmid(header,0,15))
  g = where(def EQ 'CONTENT-LENGTH:', Ng)
  if Ng GT 0 then $
         Length = long(strmid(Header[g[0]], 15, strlen(Header[g[0]])-1))
  return
END 

FUNCTION webget,  url,  SILENT=silent, COPYFILE=copyfile
  ;;
  ;;
  ;; sockets supported in unix & windows since V5.4, Macintosh since V5.6

   proc = routine_info(/system)
   g  = where(proc EQ 'SOCKET', Ng)
   If Ng EQ 0 THEN BEGIN 
      IF NOT Keyword_set(silent) THEN $
        dummy=Dialog_Message('Sorry,  web-access not supported in '+!version, /error)
      return, ''
  ENDIF
  ;;
  ;; define the result fields
  ;;
  Header = strarr(256)
  Data = strarr(256)
  Image = 0
  ImageHeader = ''
  ;;
  ;; open the connection and request the file
  ;;
  Proxy = getenv('http_proxy')
  IF Proxy NE '' THEN BEGIN 
      ;;
      ;; sort out proxy name
      ;;
      LastColon = StrPos(Proxy, ':', /Reverse_Search)
      ProxyPort = fix(StrMid(Proxy, LastColon+1, StrLen(Proxy)))
      ProxyServer = StrMid(Proxy, 7, LastColon-7)
      ;; open the connection and send the 'GET' command
      ProtocolString = " HTTP/1.0 User-Agent: IDL/"+!version.release
      socket, unit, ProxyServer,  ProxyPort, /get_lun, /swap_if_little_endian
      printf, unit, 'GET '+url+ProtocolString
      printf, unit, ''          ; a second carriage return is needed by proxies
  ENDIF ELSE BEGIN 
      ;;
      ;; same thing easier without proxy
      ;;
      slash1 = StrPos(strmid(url, 7, StrLen(url)), '/')
      Server = StrMid(url, 7, slash1 )
      purl = strmid(url,slash1+7, StrLen(url))
      Port = 80
      socket, unit, Server,  Port, /get_lun,/swap_if_little_endian
      printf, unit, 'GET '+purl +  ' HTTP/1.0'  
      printf, unit, 'HTTP/1.0 User-Agent:  IDL ' + !VERSION.RELEASE + ' on ' + $
                 !VERSION.OS + '/' + !VERSION.ARCH
      printf, unit, ''

  ENDELSE 

  LinesRead = 0
  text = 'xxx'
  ;;
  ;; now read the header
  ;;
On_IOERROR, done
  WHILE  text NE '' do begin
      readf, unit, text
      Header[LinesRead] = text
      LinesRead = LinesRead+1
      IF LinesRead MOD 256 EQ 0 THEN $
        Header=[Header, StrArr(256)]
  ENDWHILE 
DONE: On_IOERROR, NULL
  ;;
  Header = Header[0:LinesRead-1]
  MimeType, Header, Class,  Type, Length; analyze the header
  ;;
  IF Keyword_Set(CopyFile) THEN BEGIN
      openw, wunit, CopyFile, /get_lun
      aaa = bytarr(Length,/nozero)
      readu, unit, aaa
      writeu, wunit, aaa
      free_lun, wunit
      free_lun, unit
      return, 1
  ENDIF 
  ;;
  text = '' ;initialize text fields
  LinesRead = 0l
  ;;

  CASE Class OF 
      'text': BEGIN 
          ;;
          ;; read anything of class 'text'
          WHILE  eof(unit) EQ 0 do begin
              readf, unit, text
              Data[LinesRead] = text
              LinesRead = LinesRead+1
              IF LinesRead MOD 256 EQ 0 THEN $
                Data=[Data, StrArr(256)]
          ENDWHILE 
          Data = Data[0:LinesRead-1]
      END 
      'image':BEGIN
          CASE Type OF
              'x-fits':BEGIN ; answer from stsci server
                  ;;
                  Image = readfits(unit, ImageHeader)
              END 
          ENDCASE 
      END 
      'application':BEGIN 
          CASE Type OF
              'octet-stream':BEGIN ; try reading a FITS file because ESO 
                                   ; answers this way
                  Image = readfits(unit, ImageHeader)
               END 
          ENDCASE 
      END 
  ENDCASE 

  IF LinesRead EQ 0 THEN Data = ''
  free_lun, unit
  return, {Header:Header, Text:Data, ImageHeader:ImageHeader,  Image: Image}
END





    


