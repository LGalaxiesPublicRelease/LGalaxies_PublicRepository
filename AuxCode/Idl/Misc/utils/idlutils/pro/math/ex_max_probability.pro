;+
; NAME:
;   ex_max_probability
; PURPOSE:
;   Return probabilities given ex_max results
; INPUTS:
;   point       [d,N] array of data points - N vectors of dimension d
;   amp         [M] array of gaussian amplitudes
;   mean        [d,M] array of gaussian mean vectors
;   var         [d,d,M] array of gaussian variance matrices
; OPTIONAL INPUTS:
; OUTPUTS:
; OPTIONAL OUTPUTS:
;   probability [N,M] array of probabilities of point i in gaussian j
; BUGS:
;   should be called by ex_max for consistency (would involve
;   including entropy calc here)
; DEPENDENCIES:
;   idlutils
; REVISION HISTORY:
;   2002-Nov-20  Blaton
;-
pro ex_max_probability, point,amp,mean,var,probability=probability

; check dimensions
  ngauss= n_elements(amp)                      ; M
  dimen= n_elements(var)/n_elements(mean)
  ndata= n_elements(point)/dimen 
  splog, ndata,' data points,',dimen,' dimensions,',ngauss,' gaussians'

; cram inputs into correct format
  point= reform(double(point),dimen,ndata)
  amp= reform(double(amp),ngauss)
  mean= reform(double(mean),dimen,ngauss)
  var= reform(double(var),dimen,dimen,ngauss)

; allocate space for probabilities and updated parameters
  probability= reform(dblarr(ndata*ngauss),ndata,ngauss)
  exponent= reform(dblarr(ndata*ngauss),ndata,ngauss)

; compute (un-normalized) probabilities
    for j=0L,ngauss-1 do begin
      invvar= invert(var[*,*,j],/double)
      if dimen GT 1 then $
        normamp= amp[j]*sqrt(determ(invvar,/double))/(2.0*!PI)^(dimen/2) $
      else $
        normamp= amp[j]*sqrt(invvar)/(2.0*!PI)^(dimen/2)
      for i=0L,ndata-1 do begin
        delta= point[*,i]-mean[*,j]
        exponent[i,j]= delta#invvar#delta
      endfor
      probability[*,j]= normamp*exp(-0.5*exponent[*,j])
    endfor

; normalize probabilities
    probability= probability/(total(probability,2,/double)#(dblarr(ngauss)+1))

end
