;+
; NAME:
;       WCS_GETPOLE 
;
; PURPOSE:
;       Compute the coordinates of the native pole for a non-polar projection
; EXPLANATION:
;       For non-polar (cylindrical or conic) projections, the native pole is
;       not at the reference point, and WCS_GETPOLE is used to determine the
;       position of the native pole.    See section 2.4 of the paper 
;       "Representation of Celestial Coordinates in FITS" by Calabretta 
;       Greisen (2002, A&A, 395, 1077, also available at  
;       http://www.aoc.nrao.edu/~egreisen    Called by WCS_ROTATE
;
; CALLING SEQUENCE:
;       WCS_GETPOLE,  crval, lonpole, theta0, alpha_p, delta_p, LATPOLE= ]
;
; INPUT PARAMETERS:
;       crval - 2 element vector containing standard system coordinates (the 
;               longitude and latitude) of the reference point in degrees
;       lonpole - native longitude of the celestial North Pole (degrees)
;       theta0 - native latitude of the fiducial point
; OUTPUT PARAMETERS:
;       alpha_p, delta_p - celestial longitude and latitude of the native pole
;               (degrees)
; OPTIONAL KEYWORD INPUT PARAMETERS:
;       LATPOLE - native latitude of the celestial North Pole (degrees)
; REVISION HISTORY:
;       Written    W. Landsman               June, 2003
;       Fix calculation when theta0 is not 0 or 90     February 2004
;-

pro WCS_GETPOLE, crval, lonpole, theta0, alpha_p, delta_p, LATPOLE = latpole
          
; check to see that enough parameters (at least 4) were sent
 if (N_params() lt 5) then begin
    print,'Syntax - WCS_GETPOLE,  crval, lonpole, theta0 = ,alpha_p, delta_p, '
    print,'                [LATPOLE=, Phi_0 =]' 
    return
 endif 

 ; DEFINE ANGLE CONSTANTS 
 pi = !DPI
 pi2 = pi/2.d0
 radeg = 1.8d2/pi
 alpha_0 = double(crval[0])/radeg
 delta_0 = double(crval[1])/radeg

  if theta0 EQ 90 then begin
     alpha_p = alpha_0
     delta_p = delta_0
     return
 endif

; Longpole is the longitude in the native system of the North Pole in the
; standard system (default = 180 degrees).

 phi_p = double(lonpole)/radeg
 sp = sin(phi_p)
 cp = cos(phi_p)
  sd = sin(delta_0)
 cd = cos(delta_0)
 tand = tan(delta_0)

; If /ORIGIN is set then CRVAL gives the coordinates of the origin in the
; native system.   This must be converted (using Eq. 7 in Greisen & Calabretta
; with theta0 = 0) to give the coordinates of the North pole (alpha_p, delta_p)

 if (theta0 EQ 0.0) then begin
        if (delta_0 EQ 0) and (lonpole EQ 90.0d) then delta_p = latpole else $
        delta_p = acos( sd/cp)               ;Updated May 98
        if (latpole NE 90) then  $
            if abs(latpole + delta_p) LT abs(latpole - delta_p) then  $
            delta_p = -delta_p  
        if (lonpole EQ 1.8d2) or (cd EQ 0.0) then alpha_p = alpha_0 else $
                alpha_p = alpha_0  - atan(sp/cd, -tan(delta_p)*tand )
 endif else begin                ;General case for arbitary theta0
        ctheta = cos(theta0/RADEG)
        stheta = sin(theta0/RADEG)
        term1 = atan(stheta, ctheta*cp ) 
        term2 = acos( sd/( sqrt(1.0d - ctheta^2*sp^2)  ))
        if term2 EQ 0.0 then delta_p = term1 else begin
           delta_p1 = abs( (term1 + term2)*radeg)
           delta_p2 = abs((term1 - term2)*radeg)
           case 1 of 
           (delta_p1 GT 90) and (delta_p2 GT 90):message,'No valid solution'
           (delta_p1 LE 90) and (delta_p2 GT 90): delta_p = term1 + term2
           (delta_p1 GT 90) and (delta_p2 LE 90): delta_p = term1 - term2
           else: begin             ;Two valid solutions
                 delta_p1 = (term1 + term2)*radeg
                 delta_p2 = (term1 - term2)*radeg
                 if abs(latpole-delta_p1) LT abs(latpole - delta_p2) then $
                       delta_p = term1+term2 else delta_p = term1 - term2
                 end
           endcase
           if (cd EQ 0.0) then alpha_p = alpha_0 else begin
              sdelt = sin(delta_p)
              if (sdelt EQ 1) then alpha_p = alpha_0 - phi_p - !DPI else $
              if (sdelt EQ -1) then alpha_p = alpha_0 -phi_p else $
              alpha_p = alpha_0 - $
               atan( (stheta-sin(delta_p)*sd)/(cos(delta_p)*cd), sp*ctheta/cd )
           endelse
         endelse
 endelse 

 return
 end
