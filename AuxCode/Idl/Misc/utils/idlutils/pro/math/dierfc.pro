;+
; NAME:
;   dierfc
;
; PURPOSE:
;       Inverse of the Complementary Error Function "erfc^{-1}(x)" 
;
; CALLING SEQUENCE:
;   result = dierfc( input )
;
; INPUTS:
;   input      - Arbitrary array of values from 0 to 2.
;               (positive values returned for inputs between 0 and 1)
;               exact 0 return NaN,    
;
; OUTPUTS:
;   result     - The output array of type double, with range from
;                  -infinity to +infinity.  
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The results outside of -20 < results < +20 may lack desired accuracy 
;
; EXAMPLES:
;    inverse = dierfc([0.0,0.0027,0.0456,1.0d,1.6827,1.9])
;    sigma = -sqrt(2.0) * inverse*sqrt(2.0)
;    print, sigma, format='(6f10.4)'
;      -Infinity   -3.0000   -1.9991    0.0000    1.0000    1.6449
;
; COPYRIGHT:
;    Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
;    You may use, copy, modify this code for any purpose and
;    without fee. You may distribute this ORIGINAL package.
;
; REVISION HISTORY:
;   11-Jun-2002   Adapted by S. Burles, MIT
;-
;------------------------------------------------------------------------------
function dierfc, y

    IF N_params() LT 1 then begin
       print, 'Syntax: res = dierfc(input)'
       print, 'Input should be GT 0 and  LT 2'
       return, 0
    endif

    above = y GT 1
    z = y + above * 2.0 * (1.0 - y)
    inf = z EQ 0
    z = z + inf
 
;    if (y GT 1) then z = 2 - y

    w = 0.916461398268964d - alog(z)
    u = sqrt(w)
    s = (alog(u) + 0.488826640273108d) / w
    t = 1 / (u + 0.231729200323405d)
    x = u * (1 - s * (s * 0.124610454613712d + 0.5d)) - $
        ((((-0.0728846765585675d * t + 0.269999308670029d) * t + $
        0.150689047360223d) * t + 0.116065025341614d) * t + $
        0.499999303439796d) * t
    t = 3.97886080735226d / (x + 3.97886080735226d)
    u = t - 0.5
    s = (((((((((0.00112648096188977922d * u + $
        1.05739299623423047d-4) * u - 0.00351287146129100025d) * u - $
        7.71708358954120939d-4) * u + 0.00685649426074558612d) * u + $
        0.00339721910367775861d) * u - 0.011274916933250487d) * u - $
        0.0118598117047771104d) * u + 0.0142961988697898018d) * u + $
        0.0346494207789099922d) * u + 0.00220995927012179067d
    s = ((((((((((((s * u - 0.0743424357241784861d) * u - $
        0.105872177941595488d) * u + 0.0147297938331485121d) * u + $
        0.316847638520135944d) * u + 0.713657635868730364d) * u + $
        1.05375024970847138d) * u + 1.21448730779995237d) * u + $
        1.16374581931560831d) * u + 0.956464974744799006d) * u + $
        0.686265948274097816d) * u + 0.434397492331430115d) * u + $
        0.244044510593190935d) * t - $
        z * exp(x * x - 0.120782237635245222d)
    x = x + s * (x * s + 1)

    ans = x - 2.0 * above * x

    placeinf = where(inf)
    if placeinf[0] NE -1 then begin
      ans[placeinf] = (1.0 - 2.0*above[placeinf]) * !VALUES.F_INFINITY
    endif

    return, ans

end
