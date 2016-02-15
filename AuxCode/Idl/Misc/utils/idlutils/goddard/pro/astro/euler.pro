PRO EULER,AI,BI,AO,BO,SELECT, FK4 = FK4, SELECT = select1
;+
; NAME:
;     EULER
; PURPOSE:
;     Transform between Galactic, celestial, and ecliptic coordinates.
; EXPLANATION:
;     Use the procedure ASTRO to use this routine interactively
;
; CALLING SEQUENCE:
;      EULER, AI, BI, AO, BO, [ SELECT, /FK4, SELECT = ] 
;
; INPUTS:
;       AI - Input Longitude in DEGREES, scalar or vector.  If only two 
;               parameters are supplied, then  AI and BI will be modified to 
;               contain the output longitude and latitude.
;       BI - Input Latitude in DEGREES
;
; OPTIONAL INPUT:
;       SELECT - Integer (1-6) specifying type of coordinate transformation.  
;
;      SELECT   From          To        |   SELECT      From            To
;       1     RA-Dec (2000)  Galactic   |     4       Ecliptic      RA-Dec    
;       2     Galactic       RA-DEC     |     5       Ecliptic      Galactic  
;       3     RA-Dec         Ecliptic   |     6       Galactic      Ecliptic  
;
;      If not supplied as a parameter or keyword, then EULER will prompt for 
;      the value of SELECT
;      Celestial coordinates (RA, Dec) should be given in equinox J2000 
;      unless the /FK4 keyword is set.
; OUTPUTS:
;       AO - Output Longitude in DEGREES
;       BO - Output Latitude in DEGREES
;
; INPUT KEYWORD:
;       /FK4 - If this keyword is set and non-zero, then input and output 
;             celestial and ecliptic coordinates should be given in equinox 
;             B1950.
;       /SELECT  - The coordinate conversion integer (1-6) may alternatively be 
;              specified as a keyword
; NOTES:
;       EULER was changed in December 1998 to use J2000 coordinates as the 
;       default, ** and may be incompatible with earlier versions***.
; REVISION HISTORY:
;       Written W. Landsman,  February 1987
;       Adapted from Fortran by Daryl Yentis NRL
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Made J2000 the default, added /FK4 keyword  W. Landsman December 1998
;       Add option to specify SELECT as a keyword W. Landsman March 2003
;-
 On_error,2

 npar = N_params()
 if npar LT 2 then begin
    print,'Syntax - EULER, AI, BI, A0, B0, [ SELECT, /FK4, SELECT= ]'
    print,'    AI,BI - Input longitude,latitude in degrees'
    print,'    AO,BO - Output longitude, latitude in degrees'
    print,'    SELECT - Scalar (1-6) specifying transformation type'
    return
 endif

  twopi   =   2.0d*!DPI
  fourpi  =   4.0d*!DPI
  deg_to_rad = 180.0d/!DPI

;   J2000 coordinate conversions are based on the following constants
;   (see the Hipparcos explanatory supplement).
;  eps = 23.4392911111d              Obliquity of the ecliptic
;  alphaG = 192.85948d               Right Ascension of Galactic North Pole
;  deltaG = 27.12825d                Declination of Galactic North Pole
;  lomega = 32.93192d                Galactic longitude of celestial equator  
;  alphaE = 180.02322d              Ecliptic longitude of Galactic North Pole
;  deltaE = 29.811438523d            Ecliptic latitude of Galactic North Pole
;  Eomega  = 6.3839743d              Galactic longitude of ecliptic equator              

  if keyword_set(FK4) then begin 

  equinox = '(B1950)' 
  psi   = [ 0.57595865315D, 4.9261918136D,  $
            0.00000000000D, 0.0000000000D,  $  
            0.11129056012D, 4.7005372834D]     
  stheta =[ 0.88781538514D,-0.88781538514D, $
            0.39788119938D,-0.39788119938D, $
            0.86766174755D,-0.86766174755D]    
  ctheta =[ 0.46019978478D, 0.46019978478D, $
            0.91743694670D, 0.91743694670D, $
            0.49715499774D, 0.49715499774D]    
   phi  = [ 4.9261918136D,  0.57595865315D, $
            0.0000000000D, 0.00000000000D, $
	    4.7005372834d, 0.11129056012d]


 endif else begin 

  equinox = '(J2000)'
  psi   = [ 0.57477043300D, 4.9368292465D,  $
            0.00000000000D, 0.0000000000D,  $  
            0.11142137093D, 4.71279419371D]     
  stheta =[ 0.88998808748D,-0.88998808748D, $
            0.39777715593D,-0.39777715593D, $
            0.86766622025D,-0.86766622025D]    
  ctheta =[ 0.45598377618D, 0.45598377618D, $
            0.91748206207D, 0.91748206207D, $
            0.49714719172D, 0.49714719172D]    
   phi  = [ 4.9368292465D,  0.57477043300D, $
            0.0000000000D, 0.00000000000D, $
            4.71279419371d, 0.11142137093d]

 endelse
;
 if N_elements(select) EQ 0 then $
          if N_elements(select1) EQ 1 then select=select1
 if N_elements(select) EQ 0 then begin
        print,' '
        print,' 1 RA-DEC ' + equinox + ' to Galactic'
        print,' 2 Galactic       to RA-DEC' + equinox
        print,' 3 RA-DEC ' + equinox + ' to Ecliptic'
        print,' 4 Ecliptic       to RA-DEC' + equinox
        print,' 5 Ecliptic       to Galactic'
        print,' 6 Galactic       to Ecliptic'
;
        select = 0
        read,'Enter selection: ',select
 endif

 I  = select - 1                         ; IDL offset
 a  = ai/deg_to_rad - phi[i]
 b = bi/deg_to_rad
 sb = sin(b) &	cb = cos(b)
 cbsa = cb * sin(a)
 b  = -stheta[i] * cbsa + ctheta[i] * sb
 bo    = asin(b<1.0d)*deg_to_rad
;
 a =  atan( ctheta[i] * cbsa + stheta[i] * sb, cb * cos(a) )
 ao = ( (a+psi[i]+fourpi) mod twopi) * deg_to_rad


 if ( npar EQ 2 ) then begin
	ai = ao & bi=bo
 endif

 return
 end
