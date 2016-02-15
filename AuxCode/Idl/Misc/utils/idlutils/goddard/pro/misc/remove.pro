pro remove,index, v1, v2, v3, v4, v5, v6, v7
;+
; NAME:
;       REMOVE
; PURPOSE:
;       Contract a vector or up to 7 vectors by removing specified elements   
; CALLING SEQUENCE:
;       REMOVE, index, v1,[ v2, v3, v4, v5, v6, v7]     
; INPUTS:
;       INDEX - scalar or vector giving the index number of elements to
;               be removed from vectors.  Duplicate entries in index are
;               ignored.    An error will occur if one attempts to remove
;               all the elements of a vector.
;
; INPUT-OUTPUT:
;       v1 - Vector or array.  Elements specifed by INDEX will be 
;               removed from v1.  Upon return v1 will contain
;               N fewer elements, where N is the number of values in
;               INDEX.
;
; OPTIONAL INPUT-OUTPUTS:
;       v2,v3,...v7 - additional vectors containing
;               the same number of elements as v1.  These will be
;               contracted in the same manner as v1.
;
; EXAMPLES:
;       (1) If INDEX = [2,4,6,4] and V = [1,3,4,3,2,5,7,3] then after the call
;
;               IDL> remove,index,v      
;
;       V will contain the values [1,3,3,5,3]
;
;       (2) Suppose one has a wavelength vector W, and three associated flux
;       vectors F1, F2, and F3.    Remove all points where a quality vector,
;       EPS is negative
;
;               IDL> bad = where( EPS LT 0, Nbad)
;               IDL> if Nbad GT 0 then remove, bad, w, f1, f2, f3
;
; METHOD:
;       If more than one element is to be removed, then HISTOGRAM is used
;       to generate a 'keep' subscripting vector.    To minimize the length of 
;       the subscripting vector, it is only computed between the minimum and 
;       maximum values of the index.   Therefore, the slowest case of REMOVE
;       is when both the first and last element are removed.
;
; REVISION HISTORY:
;       Written W. Landsman        ST Systems Co.       April 28, 1988
;       Cleaned up code          W. Landsman            September, 1992
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Major rewrite for improved speed   W. Landsman    April 2000
;-
 On_error,2

 npar = N_params()
 if npar LT 2 then begin
      print,'Syntax - remove, index, v1, [v2, v3, v4, v5, v6, v7]'
      return
 endif

 npts = N_elements(v1)

 max_index = max(index, MIN = min_index)

 if ( min_index LT 0 ) or (max_index GT npts-1) then message, $
             'ERROR - Index vector is out of range'

 if ( max_index Eq min_index ) then begin 
     Ngood = 0  
    if npts EQ 1 then message, $ 
         'ERROR - Cannot delete all elements from a vector'
  endif else begin  ;Remove only 1 element?
         

;  Begin case where more than 1 element is to be removed.   Use HISTOGRAM
;  to determine then indices to keep

 nhist = max_index - min_index +1 

 hist = histogram( index)      ;Find unique index values to remove
 keep = where( hist EQ 0, Ngood ) + min_index

 if ngood EQ 0 then begin 
    if ( npts LE nhist ) then message, $
          'ERROR - Cannot delete all elements from a vector'
  endif 
 endelse

 imin = min_index - 1
 imax = max_index + 1
 i0 = (min_index EQ 0) + 2*(max_index EQ npts-1)       

 case i0 of 
 3: begin
    v1 = v1[keep]
    if Npar GE 3 then v2 = v2[keep]
    if Npar GE 4 then v3 = v3[keep]
    if Npar GE 5 then v4 = v4[keep]
    if Npar GE 6 then v5 = v5[keep]
    if Npar GE 7 then v6 = v6[keep]
    if Npar GE 8 then v7 = v7[keep]
    end

 1: begin
    if Ngood GT 0 then $
        v1 = [v1[keep],v1[imax:*] ] else v1 = v1[imax:*]
    if Npar GE 3 then if Ngood GT 0 then  $
        v2 = [v2[keep],v2[imax:*] ] else v2 = v2[imax:*]
    if Npar GE 4 then if NGood GT 0 then $
        v3 = [v3[keep],v3[imax:*] ] else v3 = v3[imax:*]
    if Npar GE 5 then if NGood GT 0 then $
        v4 = [v4[keep],v4[imax:*] ] else v4 = v4[imax:*]
    if Npar GE 6 then if NGood GT 0 then $
        v5 = [v5[keep],v5[imax:*] ] else v5 = v5[imax:*]
    if Npar GE 7 then if NGood GT 0 then $
        v6 = [v6[keep],v6[imax:*] ] else v6 = v6[imax:*]
    if Npar GE 8 then if NGood GT 0 then $
        v7 = [v7[keep],v7[imax:*] ] else v7 = v7[imax:*]
    end

  2: begin 
     if NGood GT 0 then $
        v1 = [v1[0:imin], v1[keep] ] else v1 = v1[0:imin]
     if Npar GE 3 then if NGood GT 0 then $
        v2 = [v2[0:imin], v2[keep] ] else v2 = v2[0:imin]
     if Npar GE 4 then if Ngood GT 0 then $
        v3 = [v3[0:imin], v3[keep] ] else v3 = v3[0:imin]
     if Npar GE 5 then if Ngood GT 0 then $
        v4 = [v4[0:imin], v4[keep] ] else v4 = v4[0:imin]
     if Npar GE 6 then if Ngood GT 0 then $
        v5 = [v5[0:imin], v5[keep] ] else v5 = v5[0:imin]
     if Npar GE 7 then if Ngood GT 0 then $
        v6 = [v6[0:imin], v6[keep] ] else v6 = v6[0:imin]
     if Npar GE 8 then if Ngood GT 0 then $
        v7 = [v7[0:imin], v7[keep] ] else v7 = v7[0:imin]
   end

 0: begin
    if NGood GT 0 then $
         v1 = [v1[0:imin], v1[keep], v1[imax:*] ] else $
         v1 = [v1[0:imin], v1[imax:*] ] 
    if Npar GE 3 then if NGood GT 0 then $
         v2 = [v2[0:imin], v2[keep], v2[imax:*] ] else $
         v2 = [v2[0:imin], v2[imax:*] ] 
    if Npar GE 4 then if NGood GT 0 then $
         v3 = [v3[0:imin], v3[keep], v3[imax:*] ] else $
         v3 = [v3[0:imin], v3[imax:*] ] 
    if Npar GE 5 then if NGood GT 0 then $
         v4 = [v4[0:imin], v4[keep], v4[imax:*] ] else $
         v4 = [v4[0:imin], v4[imax:*] ] 
    if Npar GE 6 then if NGood GT 0 then $
         v5 = [v5[0:imin], v5[keep], v5[imax:*] ] else $
         v5 = [v5[0:imin], v5[imax:*] ] 
    if Npar GE 7 then if NGood GT 0 then $
         v6 = [v6[0:imin], v6[keep], v6[imax:*] ] else $
         v6 = [v6[0:imin], v6[imax:*] ] 
    if Npar GE 8 then if NGood GT 0 then $
         v7 = [v7[0:imin], v7[keep], v7[imax:*] ] else $
         v7 = [v7[0:imin], v7[imax:*] ] 

    end
 endcase

    
 
 return
 end
