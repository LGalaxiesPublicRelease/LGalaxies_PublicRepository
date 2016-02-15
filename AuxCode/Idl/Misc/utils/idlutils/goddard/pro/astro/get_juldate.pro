pro get_juldate,jd
;+
; NAME:
;    GET_JULDATE
; PURPOSE:
;     Return the current Julian Date
;
; EXPLANATION:
;     This procedure became partially obsolete with the introduction of the
;     /JULIAN keyword to the intrinsic SYSTIME function in IDL V5.2.   Note
;     however, that SYSTIME(/JULIAN) always returns the *local* time, whereas
;     for most machines, GET_JULDATE returns  Universal Time (i.e. Greenwich
;     mean time.)
;
;     In V5.4, GET_JULDATE became completely obsolete with the introduction
;     of the /UTC keyword to SYSTIME().   So GET_JULDATE,jd is equivalent to
;     jd = SYSTIME(/JULIAN,/UTC).
;
; CALLING SEQUENCE:
;       GET_JULDATE,jd
;
; INPUTS:
;       None
;
; OUTPUTS:
;       jd = Current Julian Date, double precision scalar
;
; EXAMPLE:
;       Return the current hour, day, month and year as integers
;
;       IDL> GET_JULDATE, JD                  ;Get current Julian date
;       IDL> DAYCNV, JD, YR, MON, DAY, HOURS  ;Convert to hour,day month & year
;
; METHOD:
;       The systime(1) function is used to obtain the number of days after 
;       1-JAN-1970.     The offset to Julian days is then computed.
;
;       WARNING!   This procedure assumes that systime(1) returns the value
;       of Universal Time (UT).    This appears to be true for most Unix
;       workstations and DOS machines, but not for VMS or Macintoshes, 
;       for which systime(1) returns the local time.     Users
;       may need to add the difference between UT and local time to the value
;       of JD, depending on the particular installation.
;
; REVISION HISTORY:
;       Written Wayne Landsman                March, 1991
;       Converted to IDL V5.0   W. Landsman   September 1997
;-
 if N_Params() LT 1 then begin
     Print,'Syntax - GET_JULDATE, JD'
     return
 endif

 x = systime(1)             ;Number of seconds elapsed since 1-jan-1970
 days = x/24./3600.         ;Number of days elapsed since 1-jan-1970
 jd0 = 2440587.5d           ;Julian Date of 1-jan-1970
 jd = jd0 + days            ;Return current Julian date

 end
