;-------------------------------------------------------------
;+
; NAME:
;       REDSHIFT
; PURPOSE:
;       Interactively converts between redshift, Recession velocity, & Distance
; EXPLANATION:
;       This simple program assumes a linear Hubble law and no cosmological
;       constant.    For more general and precise conversions use the program
;       lumdist.pro 
; CALLING SEQUENCE:
;       redshift, [h, /HELP]
; INPUTS:
;       h = optional Hubble constant (def = 50 km/s/Mpc).      in
; OUTPUTS:
;       Results are displayed at the terminal screen
; NOTES:
;       Note: H may be changed at any time by typing h=new_value.
;       Also displays angular size equivalence and photometric information.
;
; MODIFICATION HISTORY:
;       R. Sterner. 17 July, 1987.
;       Johns Hopkins University Applied Physics Laboratory.
;       RES 7 Jan, 1988 --- added H0.
;
; Copyright (C) 1987, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;       Converted to IDL V5.0   W. Landsman   September 1997
;-
;-------------------------------------------------------------
 
        PRO REDSHIFT, H0, help=hlp
 
        if keyword_set(hlp) then begin
          print,' Converts between redshift, Recession velocity, and Distance'
          print,' redshift, [h]'
          print,'   h = optional Hubble constant (def = 50 km/s/Mpc).      in'
          print,' Note: H may be changed at any time by typing h=new_value.'
          print,'   Also displays angular size equivalence and photometric '+$
            'information.'
          return
        endif
 
        C = 2.9979E5    ; km/s.
        H = 50          ; km/s/Mpc.
        IF N_PARAMS(0) GT 0 THEN H = H0
 
        PRINT,' '
        PRINT,' ---==< Redshift >==---'
        PRINT,' Converts between Redshift, Recession velocity, Distance.'
LOOP:   PRINT,' '
        PRINT,' Enter Redshift as:'
        PRINT,'   Z = xxx'
        PRINT,' Enter Recession Velocity as:'
        PRINT,'   V = xxx (km/s)'
        PRINT,' Enter Distance as:'
        PRINT,'   D = xxx (Mpc)'
        PRINT,' To change Hubble constant (current value = '+STRTRIM(H,2)+'):'
        PRINT,'   H = xxx (km/s/Mpc)'
        PRINT,' '
        TXT = ''
        READ, ' Entry: ', TXT
        IF TXT EQ '' THEN RETURN
        TXT = STRUPCASE(TXT)
        TXT = REPCHR(TXT,'=')
        W = GETWRD(TXT, 0)
        X = GETWRD(TXT, 1)
 
        CASE W OF
'H':    BEGIN
          H = X + 0.
          GOTO, LOOP
        END
'Z':    BEGIN
          Z = X + 0.
          V = C*((Z+1)^2-1)/((Z+1)^2+1)
          D = V/H
        END
'V':    BEGIN
          V = X + 0.
          Z = SQRT((1+V/C)/(1-V/C)) - 1
          D = V/H
        END
'D':    BEGIN
          D = X + 0.
          V = D*H
          Z = SQRT((1+V/C)/(1-V/C)) - 1
        END
else:   begin
          print,' Example:  z=.329'
          print,' '
          goto, loop
        end
        ENDCASE
 
        PRINT,' '
        PRINT,' For H = '+STRTRIM(H,2)+' km/s/Mpc:'
        PRINT,' Redshift, Z = ',strtrim(Z,2)
        PRINT,' Recession Velocity, V = '+strtrim(V,2)+' km/s'
        PRINT,' Distance, D = '+strtrim(D,2)+' Mpc ('+strtrim(D*3.258,2)+$
          ' million lt-yrs)'
        PRINT,' '
        print,' REDSHIFT EFFECTS'
        print,' Angular size/Linear size:'
        print,'   Angular size is increased by a factor of (1+z) = '+$
          strtrim((1.+z),2)
        print,'     due to recession velocity.'
        t = d*tan((1./60.)/(1.+z)/!radeg)*1000.         ; kpc.
        print,"   So 1' in the sky at this distance is a linear distance of "+$
          strtrim(t,2)+" kpc,"
        print,'   So 1" in the sky at this distance is a linear distance of '+$
          strtrim(t/60.,2)+' kpc.'
        print,' Photometry:'
        print,'   Distance Modulus = m-M = '+strtrim(5.*alog10(d*1e6/10.), 2)
        print,'   Rest wavelengths increased by (1+z) so for the UBVRI system'
        print,'     U shifts from 360 nm to '+strtrim(360.*(1+z),2)+' nm'
        print,'     B shifts from 420 nm to '+strtrim(420.*(1+z),2)+' nm'
        print,'     V shifts from 520 nm to '+strtrim(520.*(1+z),2)+' nm'
        print,'     R shifts from 680 nm to '+strtrim(680.*(1+z),2)+' nm'
        print,'     I shifts from 825 nm to '+strtrim(825.*(1+z),2)+' nm'
        print,' Press any key to continue'
        k = get_kbrd(1)
        print,' '
        GOTO, LOOP
 
        END
