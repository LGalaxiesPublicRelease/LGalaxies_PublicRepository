pro gal_uvw, u, v, w, distance = distance, LSR = lsr, ra=ra,dec=dec, $
     pmra = pmra, pmdec=pmdec, vrad = vrad, plx  = plx
;+
; NAME:
;     GAL_UVW
; PURPOSE:
;     Calculate the Galactic space velocity (U,V,W) of star  
; EXPLANATION:
;     Calculates the Galactic space velocity U, V, W of star given its 
;     (1) coordinates, (2) proper motion, (3) distance (or parallax), and 
;     (4) radial velocity.
; CALLING SEQUENCE:
;     GAL_UVW, U, V, W, [/LSR, RA=, DEC=, PMRA= ,PMDEC=, VRAD= , DISTANCE= 
;              PLX= ]
; OUTPUT PARAMETERS:
;      U - Velocity (km/s) positive toward the Galactic *anti*center
;      V - Velocity (km/s) positive in the direction of Galactic rotation
;      W - Velocity (km/s) positive toward the North Galactic Pole 
; REQUIRED INPUT KEYWORDS:
;      User must supply a position, proper motion,radial velocity and distance
;      (or parallax).    Either scalars or vectors can be supplied.
;     (1) Position:
;      RA - Right Ascension in *Degrees*
;      Dec - Declination in *Degrees*
;     (2) Proper Motion
;      PMRA = Proper motion in RA in arc units (typically milli-arcseconds/yr)
;      PMDEC = Proper motion in Declination (typically mas/yr)
;     (3) Radial Velocity
;      VRAD = radial velocity in km/s
;     (4) Distance or Parallax
;      DISTANCE - distance in parsecs 
;                 or
;      PLX - parallax with same distance units as proper motion measurements
;            typically milliarcseconds (mas)
;
; OPTIONAL INPUT KEYWORD:
;      /LSR - If this keyword is set, then the output velocities will be 
;             corrected for the solar motion (U,V,W)_Sun = (-9,+12,+7)  
;             (Mihalas & Binney, 1981) to the local standard of rest
;  EXAMPLE:
;      (1) Compute the U,V,W coordinates for the halo star HD 6755.  
;          Use values from Hipparcos catalog, and correct to the LSR
;      ra = ten(1,9,42.3)*15.    & dec = ten(61,32,49.5)
;      pmra = 627.89  &  pmdec = 77.84         ;mas/yr
;      dis = 144    &  vrad = -321.4
;      gal_uvw,u,v,w,ra=ra,dec=dec,pmra=pmra,pmdec=pmdec,vrad=vrad,dis=dis,/lsr
;          ===>  u=154  v = -493  w = 97        ;km/s
;
;      (2) Use the Hipparcos Input and Output Catalog IDL databases (see 
;      http://idlastro.gsfc.nasa.gov/ftp/zdbase/) to obtain space velocities
;      for all stars within 10 pc with radial velocities > 10 km/s
;
;      dbopen,'hipparcos,hic'      ;Need Hipparcos output and input catalogs
;      list = dbfind('plx>100,vrad>10')      ;Plx > 100 mas, Vrad > 10 km/s
;      dbext,list,'pmra,pmdec,vrad,ra,dec,plx',pmra,pmdec,vrad,ra,dec,plx
;      ra = ra*15.                 ;Need right ascension in degrees
;      GAL_UVW,u,v,w,ra=ra,dec=dec,pmra=pmra,pmdec=pmdec,vrad=vrad,plx = plx 
;      forprint,u,v,w              ;Display results
; METHOD:
;      Follows the general outline of Johnson & Soderblom (1987, AJ, 93,864)
;      except that U is positive outward toward the Galactic *anti*center, and 
;      the J2000 transformation matrix to Galactic coordinates is taken from  
;      the introduction to the Hipparcos catalog.   
; REVISION HISTORY:
;      Written, W. Landsman                       December   2000
;-
 if N_Params() EQ 0 then begin
       print,'Syntax - GAL_UVW, U, V, W, [/LSR, RA=, DEC=, PMRA= ,PMDEC=, VRAD='
       print,'                  Distance=, PLX='
       print,'         U, V, W - output Galactic space velocities (km/s)'
       return
 endif
 
 Nra = N_elements(ra)
 if (nra EQ 0) or (N_elements(dec) EQ 0) then message, $
     'ERROR - The RA, Dec (J2000) position keywords must be supplied (degrees)'
 if (N_elements(vrad) LT Nra) then message, $
      'ERROR - A Radial Velocity (km/s) must be supplied for each star'
 if (N_elements(pmra) LT Nra) or (N_elements(pmdec) LT Nra) then message, $ 
      'ERROR - A proper motion must be supplied for each star'
 if N_elements(distance) GT 0 then begin 
       bad = where(distance LE 0, Nbad)
       if Nbad GT 0 then message,'ERROR - All distances must be > 0'
       plx = 1e3/distance          ;Parallax in milli-arcseconds
 endif else begin
       if N_elements(plx) EQ 0 then message, $
             'ERROR - Either a parallax or distance must be specified'
       bad = where(plx LE 0.0, Nbad)
       if Nbad GT 0 then message,'ERROR - Parallaxes must be > 0'
 endelse

 cosd = cos(dec/!radeg)
 sind = sin(dec/!radeg)
 cosa = cos(ra/!RADEG)
 sina = sin(ra/!RADEG)

 u = fltarr(Nra) & v = u & w = u

 k = 4.74047     ;Equivalent of 1 A.U/yr in km/s   
 t = [ [ 0.0548755604, +0.4941094279, -0.8676661490], $ 
      [ 0.8734370902, -0.4448296300,-0.1980763734], $
      [ 0.4838350155, 0.7469822445, +0.4559837762] ]

 for i = 0,Nra -1 do begin
 a = [ [cosa[i]*cosd[i],sina[i]*cosd[i] ,sind[i] ], [-sina[i], cosa[i], 0],  $
          [-cosa[i]*sind[i],-sina[i]*sind[i],cosd[i]] ]
 b = t#a
 vec = [vrad[i], k*pmra[i]/plx[i], k*pmdec[i]/plx[i] ]
 starvel = b#vec
 if keyword_set(lsr) then starvel = starvel + [-9,12,7]
 u[i] = starvel[0]
 v[i] = starvel[1]
 w[i] = starvel[2]
 endfor

 sz = size(ra)
 if sz(0) EQ 0 then begin
     u = u[0]  & v = v[0]   & w = w[0]
 endif

 return
 end
