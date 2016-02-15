PRO srcor,x1in,y1in,x2in,y2in,dcr,ind1,ind2,option=option,magnitude=magnitude,$
   spherical=spherical
;+
; NAME:
;       SRCOR
; PURPOSE:
;       Correlate the source positions found on two lists.
; CALLING SEQUENCE:
;       srcor,x1in,ylin,x2in,y2in,dcr,ind1,ind2
; INPUTS:
;       x1in,y1in - First set of x and y coordinates.  The program
;                   marches through this list element by element,
;                   looking in list 2 for the closest match.  So, the program
;                   will run faster if this is the shorter of the two lists.
;                   Unless you use the option or magnitude keyword, there is
;                   nothing to guarantee unique matches.  
;       x2in,y2in - Second set of x and y coordinates.  This list is
;                   searched in its entirety every time one element of list 1
;                   is processed.
;       dcr - Critical radius outside which correlations are rejected;
;             but see 'option' below.
; OPTIONAL KEYWORD INPUT:
;       option - Changes behavior of program and description of output
;                lists slightly, as follows: 
;       OPTION=0 or left out
;             Same as older versions of SRCOR.  The closest match from list2
;             is found for each element of list 1, but if the distance is
;             greater than DCR, the match is thrown out.  Thus the index
;             of that element within list 1 will not appear in the IND1 output
;             array.
;       OPTION=1
;             Forces the output mapping to be one-to-one.  OPTION=0 results,
;             in general, in a many-to-one mapping from list 1 to list 2.
;             Under OPTION=1, a further processing step is performed to
;             keep only the minimum-distance match, whenever an entry from
;             list 1 appears more than once in the initial mapping.
;       OPTION=2
;             Same as OPTION=1, except the critical distance parameter DCR
;             is ignored.  I.e., the closest object is retrieved from list 2
;             for each object in list 1 WITHOUT a critical-radius criterion,
;             then the clean-up of duplicates is done as under OPTION=1.
;       magnitude
;             An array of stellar magnitudes corresponding to x1in and y1in.  
;             If this is supplied, then the brightest star from list 1
;             within the selected distance of the star in list 2 is taken.
;             The option keyword is ignored in this case.
;       spherical
;             If SPHERICAL=1, it is assumed that the input arrays are in
;             celestial coordinates (RA and Dec), with x1in and x2in in
;             decimal hours and y1in and y2in in decimal degrees.  If
;             SPHERICAL=2 then it is assumed that the input arrays are in
;             longitude and latitude with x1in,x2in,y1in,y2in in decimal
;             degrees.  In both cases, the critial radius dcr is in
;             *arcseconds*.  Calculations of spherical distances are made
;             with the gcirc program.
; OUTPUTS:
;       ind1 - index of matched stars in first list
;       ind2 - index of matched stars in second list
; COMMON BLOCKS:
;       none
; SIDE EFFECTS:
;       none
; METHOD:
;       See under keyword LEVEL above.
; REVISON HISTORY:
;       Adapted from UIT procedure  J.Wm.Parker, SwRI 29 July 1997
;       Converted to IDL V5.0   W. Landsman   September 1997
;       
;-
;
 ON_Error,2   ; Return if error (incl. non-info message)

;;;
;   If not enough parameters, then print out the syntax.
;
IF N_params() lt 7 THEN BEGIN
  print,'SRCOR calling sequence: '
  print,'srcor,x1in,y1in,x2in,y2in,dcr,ind1,ind2 [,option={0, 1, or 2}] $'
  print,'      [,magnitude=mag_list_1, spherical={1 or 2}]'
  RETURN
ENDIF

;;;
;   Keywords.
;
IF not keyword_set(option) THEN option=0
message,/info,'Option code = '+strtrim(option,2)
IF (option lt 0) or (option gt 2) THEN MESSAGE,'Invalid option code.'

SphereFlag = keyword_set(Spherical)

;;;
;   Store the input variables into internal arrays that we can manipulate and
; modify.
;
x1=float(x1in)
y1=float(y1in)
x2=float(x2in)
y2=float(y2in)

;;;
;   If the Spherical keyword is set, then convert the input values (degrees
; and maybe hours) into radians, so GCIRC doesn't have to make this calculation
; each time it is called in the FOR loop.  Also convert the critical radius
; (which is in arcsec, so convert by 3600.) to radians
;
if SphereFlag then begin
   dcr2 = dcr
   if (Spherical eq 1) then XScale = 15.0 else XScale = 1.0 
   d2r  = !DPI/180.0d0
   x1 = x1 * (XScale * d2r)
   y1 = y1 * d2r
   x2 = x2 * (XScale * d2r)
   y2 = y2 * d2r
   dcr2 = dcr2 * (d2r / 3600.)
endif else dcr2=dcr^2


;;;
;   Set up some other variables.
;
n1=n_elements(x1) & message,/info,strtrim(n1,2)+' sources in list 1'
n2=n_elements(x2) & message,/info,strtrim(n2,2)+' sources in list 2'
nmch=0
ind1=-1L & ind2=-1L

;;;
;   The main loop.  Step through each index of list 1, look for matches in 2.
;
FOR i=0L,n1-1 DO BEGIN
   xx=x1[i] & yy=y1[i] 
   if SphereFlag then gcirc,0,xx,yy,x2,y2,d2 else d2=(xx-x2)^2+(yy-y2)^2
   dmch=min(d2,m)
   IF (option eq 2) or (dmch le dcr2) THEN BEGIN
      nmch=nmch+1
      IF nmch eq 1 THEN BEGIN 
         ind1=long(i) 
         ind2=long(m)
      ENDIF ELSE BEGIN
         ind1=[ind1,i]
         ind2=[ind2,m]
      ENDELSE
   ENDIF
ENDFOR
message,/info,strtrim(nmch,2)+' matches found.'

;;;
;   Modify the matches depending on input options.
;
use_mag = (n_elements(magnitude) ge 1)
IF (option eq 0) and (not use_mag) THEN RETURN
IF use_mag THEN BEGIN
   message,/info,'Cleaning up output list using magnitudes.'
ENDIF ELSE BEGIN
   IF option eq 1 then message,/info,'Cleaning up output list (option = 1).'
   IF option eq 2 then message,/info,'Cleaning up output list (option = 2).'
ENDELSE

FOR i=0L,max(ind2) DO BEGIN
   csave = n_elements(ind2)
   ww = where(ind2 eq i,count) ; All but one of the list in WW must
                               ; eventually be removed.
   IF count gt 1 THEN BEGIN
      IF use_mag THEN BEGIN
         dummy = min(magnitude[ind1[ww]],m)
      ENDIF ELSE BEGIN
         xx=x2[i] & yy=y2[i]
         if SphereFlag then gcirc,0,xx,yy,x1[ind1[ww]],y1[ind1[ww]],d2 else $
                            d2=(xx-x1[ind1[ww]])^2+(yy-y1[ind1[ww]])^2
         IF n_elements(d2) ne count THEN MESSAGE,'Logic error 1'
         dummy = min(d2,m)
      ENDELSE
      remove,m,ww              ; Delete the minimum element
                               ; from the deletion list itself.

      remove,ww,ind1,ind2      ; Now delete the deletion list from
                               ; the original index arrays.
      IF n_elements(ind2) ne (csave-count+1) THEN MESSAGE,'Logic error 2'
      IF n_elements(ind1) ne (csave-count+1) THEN MESSAGE,'Logic error 3'
      IF n_elements(ind2) ne n_elements(ind1) THEN MESSAGE,'Logic error 4'
   ENDIF
ENDFOR

message,/info,strtrim(n_elements(ind1),2)+' left in list 1'
message,/info,strtrim(n_elements(ind2),2)+' left in list 2'

;
RETURN
end
