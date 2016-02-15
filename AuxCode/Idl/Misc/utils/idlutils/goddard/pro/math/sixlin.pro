pro sixlin,xx,yy,a,siga,b,sigb
;+
; NAME:
;       SIXLIN
; PURPOSE:
;       Compute linear regression coefficients by six different methods.
; EXPLANATION:
;       Adapted from the FORTRAN program (Rev. 1.1)  supplied by Isobe, 
;       Feigelson, Akritas, and Babu Ap. J. Vol. 364, p. 104 (1990).   
;       Suggested when there is no understanding about the nature of the 
;       scatter about a linear relation, and NOT when the errors in the 
;       variable are calculable.
;
; CALLING SEQUENCE:
;       SIXLIN, xx, yy, a, siga, b, sigb   
;
; INPUTS:
;       XX - vector of X values
;       YY - vector of Y values, same number of elements as XX
;
; OUTPUTS:
;       A - Vector of 6 Y intercept coefficients
;       SIGA - Vector of standard deviations of 6 Y intercepts
;       B - Vector of 6 slope coefficients
;       SIGB - Vector of standard deviations of slope coefficients
;
;       The output variables are computed using linear regression for each of 
;       the following 6 cases:
;               (0) Ordinary Least Squares (OLS) Y vs. X
;               (1) Ordinary Least Squares  X vs. Y
;               (2) Ordinary Least Squares Bisector
;               (3) Orthogonal Reduced Major Axis
;               (4) Reduced Major-Axis 
;               (5) Mean ordinary Least Squares
;
; NOTES:
;       Isobe et al. make the following recommendations
;
;       (1) If the different linear regression methods yield similar results
;               then quoting OLS(Y|X) is probably the most familiar.
;
;       (2) If the linear relation is to be used to predict Y vs. X then
;               OLS(Y|X) should be used.   
;
;       (3) If the goal is to determine the functional relationship between
;               X and Y then the OLS bisector is recommended.
;
; REVISION HISTORY:
;       Written   Wayne Landsman          February, 1991         
;       Corrected sigma calculations      February, 1992
;       Converted to IDL V5.0   W. Landsman   September 1997
;-
 On_error, 2                                   ;Return to Caller

 if N_params() LT 5 then begin   
    print,'Syntax - SIXLIN, xx, yy, a, siga, b, sigb'   
    return
  endif

 a = dblarr(6) & b=a & siga = a & sigb =a
 x = double(xx)      ;Keep input X and Y vectors unmodified
 y = double(yy)
 rn = N_elements(x)

 if rn LT 2 then $
    message,'Input X and Y vectors must contain at least 2 data points'

 if rn NE N_elements(y) then $
    message,'Input X and Y vectors must contain equal number of data points'

; Compute averages and sums

 xavg = total(x)/rn
 yavg = total(y)/rn
 x = x - xavg
 y = y - yavg
 sxx = total(x^2)
 syy = total(y^2)
 sxy = total(x*y)
 if sxy EQ 0. then $
      message,'SXY is zero, SIXLIN is terminated'
 if sxy LT 0. then sign = -1.0 else sign = 1.0

; Compute the slope coefficients

 b[0] = sxy / sxx
 b[1] = syy / sxy
 b[2] = (b[0]*b[1] - 1.D + sqrt((1.D + b[0]^2)*(1.D +b[1]^2)))/(b[0] + b[1] )
 b[3] = 0.5 * ( b[1] - 1.D/b[0] + sign*sqrt(4.0D + (b[1]-1.0/b[0])^2))
 b[4] = sign*sqrt( b[0]*b[1] )
 b[5] = 0.5 * ( b[0] + b[1] )

; Compute Intercept Coefficients

 a = yavg - b*xavg

;  Prepare for computation of variances

 gam1 = b[2] / ( (b[0] + b[1]) *   $
         sqrt( (1.D + b[0]^2)*(1.D + b[1]^2)) )
 gam2 = b[3] / (sqrt( 4.D*b[0]^2 + ( b[0]*b[1] - 1.D)^2))
 sum1 = total( ( x*( y - b[0]*x ) )^2)
 sum2 = total( ( y*( y - b[1]*x ) )^2)
 sum3 = total( x * y * ( y - b[0]*x) * (y - b[1]*x ) )
 cov = sum3 / ( b[0]*sxx^2 )

; Compute variances of the slope coefficients

 sigb[0] = sum1 / sxx^2
 sigb[1] = sum2 / sxy^2
 sigb[2] = (gam1^2) * ( ( (1.D + b[1]^2) ^2 )*sigb[0] +  $
                  2.D*(1.D + b[0]^2) * (1.D + b[1]^2)*cov + $
                  (  (1.D + b[0]^2)^2)*sigb[1] )
 sigb[3] = (gam2^2)*( sigb[0]/b[0]^2 + 2.D*cov + b[0]^2*sigb[1] )
 sigb[4] = 0.25*(b[1]*sigb[1]/b[1] + $
                     2.D*cov + b[0]*sigb[1]/b[1] )
 sigb[5] = 0.25*(sigb[0] + 2.D*cov + sigb[1] )

; Compute variances of the intercept coefficients

 siga[0] = total( ( ( y - b[0]*x) * (1.D - rn*xavg*x/sxx) )^2 )
 siga[1] = total( ( ( y - b[1]*x) * (1.D - rn*xavg*y/sxy) )^2 ) 
 siga[2] = total( ( (x * (y - b[0]*x) * (1.D + b[1]^2) / sxx + $
                  y * (y - b[1]*x) * (1.D + b[0]^2) / sxy)*  $
                  gam1 * xavg * rn - y + b[2] * x) ^ 2)
 siga[3] = total( ( ( x * ( y - b[0]*x) / sxx + $
                   y * ( y - b[1]*x) * b[0]^2/ sxy) * gam2 * $
                   xavg * rn / sqrt( b[0]^2) - y + b[3]*x) ^ 2 )
 siga[4] = total( ( ( x * ( y - b[0] * x) * sqrt( b[1] / b[0] ) / sxx + $
                   y * ( y - b[1] * x) * sqrt( b[0] / b[1] ) / sxy) * $
                  0.5 * rn * xavg - y + b[4] * x)^2 )
 siga[5] = total( ( (x * ( y - b[0] * x) / sxx +  $
                  y * ( y - b[1] * x) / sxy)*    $
                  0.5 * rn * xavg - y + b[5]*x )^2 )

; Convert variances to standard deviation

 sigb = sqrt(sigb)
 siga = sqrt(siga)/rn
 
 return
 end
