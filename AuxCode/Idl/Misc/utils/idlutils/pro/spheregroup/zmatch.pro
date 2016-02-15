;+
; NAME:
;   zmatch
; PURPOSE:
;   Match objects in ra, dec, and z
; CALLING SEQUENCE:
;   zmatch, z1, ra1, dec1, z2, ra2, dec2, mproj, mz, $
;      match1, match2, dproj=, dz=, nmatch=, maxmatch=
; INPUTS:
;   z1          - redshifts 
;   ra1         - ra coordinates in degrees (N-dimensional array)
;   dec1        - dec coordinates in degrees (N-dimensional array)
;   z2          - redshifts 
;   ra2         - ra coordinates in degrees (N-dimensional array)
;   dec2        - dec coordinates in degrees (N-dimensional array)
;   mproj       - projected distance which defines a match (redshift units)
;   mz          - redshift distance which defines a match (redshift units)
; OPTIONAL INPUTS:
;   maxmatch    - maximum number of matches per object (default to 1)
;                 0 means unlimited
; OUTPUTS:
;   dproj      - projected distance for each match
;   dz         - redshift distance for each match
;   match1     - List of indices of matches in list 1; -1 if no matches
;   match2     - List of indices of matches in list 2; -1 if no matches
;   nmatch     - number of matches
; COMMENTS:
;   Calls matchnd. 
;   Bases "closest" match on projected distance.
; REVISION HISTORY:
;   Begun 2007-07-12 MRB (NYU)
;-
;------------------------------------------------------------------------------
pro zmatch, z1, ra1, dec1, z2, ra2, dec2, mproj, mz, match1, $
            match2, dproj=dproj, dz=dz, nmatch=nmatch, $
            maxmatch=maxmatch

; Need at least 3 parameters
if (N_params() LT 10) then begin
    print, 'Syntax - zmatch, z1, ra1, dec1, z2, ra2, dec2, match_angle, match_z, match1, match2, $'
    print, ' dangle, dz'
    return
endif

if(n_elements(maxmatch) eq 0) then maxmatch=1L

theta1=(90.D)-dec1
theta2=(90.D)-dec2

x1=angles_to_x(ra1, theta1)
x2=angles_to_x(ra2, theta2)
for i=0L,2L do $
  x1[i,*]=x1[i,*]*z1
for i=0L,2L do $
  x2[i,*]=x2[i,*]*z2

;; first find all POSSIBLE matches
matchlength=sqrt(mproj^2+mz^2)*1.1
splog, 'starting matchnd ...'
matchnd, x1, x2, matchlength, m1=m1, m2=m2, d12=d12, $
  nmatch=nmatch, maxmatch=0L, nd=3L

if(nmatch eq 0) then begin
    match1=-1L
    match2=-1L
    dproj=-1L
    dz=-1L
    return
endif

;; then check if they pass both criteria
splog, 'sorting out projections ...'
isort=sort(m1)
iuniq=uniq(m1[isort])
istart=0L
diffz=dblarr(n_elements(m1))
diffproj=dblarr(n_elements(m1))
closerank=lonarr(n_elements(m1))
for i=0L, n_elements(iuniq)-1L do begin
    iend=iuniq[i]
    icurr=isort[istart:iend]

    ;; find parallel difference
    mid=(x1[*,m1[icurr]]+x2[*,m2[icurr]])*0.5
    midnorm=sqrt(total(mid^2,1))
    middir=mid
    for j=0L,2L do $
      middir[j,*]=middir[j,*]/midnorm
    diff=(x1[*,m1[icurr]]-x2[*,m2[icurr]])
    diffz[icurr]=total(middir*diff,1)
    
    ;; find projected difference
    diffpar=middir
    for j=0L,2L do $
      diffpar[j,*]=diffz[icurr]*diffpar[j,*]
    diffperp=diff-diffpar
    diffproj[icurr]=sqrt(total(diffperp^2,1))

    closerank[icurr[sort(diffproj[icurr])]]=lindgen(n_elements(icurr))
    
    istart=iend+1L
endfor

imatch=where(abs(diffz) lt mz AND abs(diffproj) lt mproj AND $
             (closerank lt maxmatch OR maxmatch eq 0), nmatch)
if(nmatch eq 0) then begin
    match1=-1L
    match2=-1L
    dproj=-1L
    dz=-1L
    return
endif 

match1=m1[imatch]
match2=m2[imatch]
dproj=diffproj[imatch]
dz=abs(diffz[imatch])

return
end
;------------------------------------------------------------------------------
