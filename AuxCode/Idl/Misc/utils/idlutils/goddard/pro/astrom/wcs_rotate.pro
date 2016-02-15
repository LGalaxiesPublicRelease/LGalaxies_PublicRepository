;+
; NAME:
;       WCS_ROTATE 
;
; PURPOSE:
;       Rotate between standard (e.g. celestial) and native coordinates
; EXPLANATION:
;       Computes a spherical coordinate rotation between native coordinates 
;       and  standard celestial coordinate system (celestial, Galactic, or
;       ecliptic).   Applies the equations in Appendix A of the paper 
;       "Representation of Celestial Coordinates in FITS" by Calabretta 
;       Greisen (2002, A&A, 395, 1077).    Also see 
;       http://www.aoc.nrao.edu/~egreisen
;
; CATEGORY:
;       Mapping and Auxiliary FITS Routine
;
; CALLING SEQUENCE:
;       WCS_ROTATE, longitude, latitude, phi, theta, crval, 
;               [LONGPOLE = , LATPOLE = , /REVERSE, /ORIGIN ]
;
; INPUT PARAMETERS:
;       crval - 2 element vector containing standard system coordinates (the 
;               longitude and latitude) of the reference point
;
; INPUT OR OUTPUT PARAMETERS
;       longitude - longitude of data, scalar or vector, in degrees, in the
;               standard celestial coordinate system
;       latitude - latitude of data, same number of elements as longitude, 
;               in degrees
;       phi - longitude of data in the native system, in degrees, scalar or
;               vector
;       theta - latitude of data in the native system, in degrees, scalar or
;               vector
;
;       If the keyword(REVERSE) is set then phi and theta are input parameters
;       and longitude and latitude are computed.    Otherwise, longitude and
;       latitude are input parameters and phi and theta are computed.
;
; OPTIONAL KEYWORD INPUT PARAMETERS:
;
;      ORIGIN - If this keyword is set and non-zero, then the reference point
;               given by CRVAL in the native system is assumed to be at the
;               origin of the coordinates, rather than at the North Pole.
;               ORIGIN should be set for cylindrical projections (Cylindrical
;               perspective-CYP, Cartesian - CAR, Mercator - MER, Cylindrical
;               Equal area - CEA) and conventional projections (Bonne's equal
;               area - BON, Polyconic - PCO, Sinusoidal - GLS, Parabolic - PAR,
;               Aitoff - AIT, Mollweide - MOL, COBE quadrilateralized sphere -
;               CSC, Quadrilateralized Spherical Cube - QSC, and Tangential
;               Spherical Cube - TSC)
;
;       LONGPOLE - native longitude of standard system's North Pole, default
;               for a Zenithal system is 180 degrees
;       LATPOLE -  native latitude of the standard system's North Pole
;       /REVERSE - if set then phi and theta are input parameters and longitude
;                  and latitude are computed.    By default, longitude and
;                  latitude are input parameters and phi and theta are computed.
; REVISION HISTORY:
;       Written    W. Landsman               December, 1994
;       Fixed error in finding North Pole if /ORIGIN and LONGPOLE NE 180
;       Xiaoyi Wu and W. Landsman,   March, 1996
;       Fixed implementation of March 96 error, J. Thieler,  April 1996
;       Updated to IDL V5.0   W. Landsman    December 1997
;       Fixed determination of alpha_p if /ORIGIN and LONGPOLE EQ 180
;               W. Landsman    May 1998
;       Ensure argument of ASIN() is -1<x<-1 after roundoff 
;               W. Landsman/R. Arendt  June 2002
;       Call WCS_GETPOLE, accept LATPOLE keyword, update cylindrical coords
;               W. Landsman  June 2003 
;       Don't attempt to rotate NaN values   W. Landsman  May 2004
;       
;-

pro wcs_rotate, longitude, latitude, phi, theta, crval, LONGPOLE = longpole, $
          LATPOLE = latpole,REVERSE=reverse, ORIGIN = origin, THETA0 = theta0


; check to see that enough parameters (at least 4) were sent
 if (N_params() lt 5) then begin
    print,'Syntax - WCS_ROTATE, longitude, latitude, phi, theta, crval'
    print,'               LATPOLE =,  LONGPOLE = , /REVERSE, /ORIGIN' 
    return
 endif 

 ; DEFINE ANGLE CONSTANTS 
 pi = !DPI
 pi2 = pi/2.d0
 radeg = 1.8d2/pi

 if keyword_set( REVERSE) then begin
        if min([ N_elements(phi), N_elements(theta) ]) EQ 0 then          $
        message,'ERROR - Native Coordinates (phi,theta) not defined'    
 endif else begin
        if min([ N_elements(longitude), N_elements(latitude) ]) EQ 0 then $ 
        message, 'ERROR - Celestial Coordinates (long,lat) not defined' 
 endelse

; Longpole is the longitude in the native system of the North Pole in the
; standard system (default = 180 degrees).

 if N_elements(longpole) eq 0 then longpole = 1.8d2
 phi_p = double(longpole)/radeg
 sp = sin(phi_p)
 cp = cos(phi_p)

; If Theta0 = 90 then CRVAL gives the coordinates of the origin in the
; native system.   This must be converted (using Eq. 7 in Greisen & Calabretta
; with theta0 = 0) to give the coordinates of the North pole (alpha_p, delta_p)

 if theta0 EQ 90 then begin
        alpha_p = double(crval[0])/radeg
        delta_p = double(crval[1])/radeg
 endif else WCS_GETPOLE, crval, longpole, theta0, alpha_p, delta_p, $
            LATPOLE = latpole
    
; compute useful quantities relating to reference angles
  sa = sin(alpha_p)
  ca = cos(alpha_p)
  sd = sin(delta_p)
  cd = cos(delta_p)

; calculate rotation matrix 

  r = [ [-sa*sp - ca*cp*sd,  ca*sp - sa*cp*sd, cp*cd ] , $
        [ sa*cp - ca*sp*sd, -ca*cp - sa*sp*sd, sp*cd ] , $
        [ ca*cd           ,  sa*cd           , sd    ] ]

; solve the set of equations for each datum point

 if keyword_set(REVERSE) then begin
        latitude = phi
        longitude = theta
        g = where( finite(phi) and finite(theta), Ng )
        if Ng EQ 0 then return
        phi1 = double(phi[g])/radeg
        theta1 = double(theta[g])/radeg
        r = transpose(r)
 endif else begin
        phi = longitude
        phi1 = double(longitude)/radeg
        theta1 = double(latitude)/radeg
 endelse

; define the right-hand side of the equations

 l = cos(theta1)*cos(phi1)
 m = cos(theta1)*sin(phi1)
 n = sin(theta1)

; find solution to the system of equations and put it in b
; Can't use matrix notation in case l,m,n are vectors

 b0 = r[0,0]*l + r[1,0]*m + r[2,0]*n
 b1 = r[0,1]*l + r[1,1]*m + r[2,1]*n
 b2 = (r[0,2]*l + r[1,2]*m + r[2,2]*n) > (-1) < 1 ;Account for possible roundoff

; use b0,b1,b2 to compute "native" latitude and longitude

 if keyword_set(REVERSE) then begin
        latitude[g] = asin(b2)*radeg
        longitude[g] = atan( b1, b0)*radeg
 endif else begin
        theta = asin(b2)*radeg
        phi = atan( b1, b0)*radeg
 endelse

 return
 end
