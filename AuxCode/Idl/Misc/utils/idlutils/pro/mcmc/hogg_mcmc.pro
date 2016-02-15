;+
; NAME:
;   hogg_mcmc
; PURPOSE:
;   Make a Markov Chain Monte Carlo chain.
; INPUTS:
;   seed       - random number seed for randomu()
;   step_func  - function that takes a step in parameter space
;   like_func  - function that computes the likelihood
;   nstep      - number of links
;   pars       - initial parameters (can be an array or structure)
; KEYWORDS:
;   log        - like_func returns log_e(likelihood), not straight
;                likelihood
; OUTPUTS:
;   pars       - array of parameters, sorted from most to least likely
;   like       - array of likelihoods for the pars
; BUGS:
; REVISION HISTORY:
;   2005-03-31  started - Hogg
;-
pro hogg_mcmc, seed,step_func,like_func,nstep,pars,like,log=log
npars= n_elements(pars)
oldpars= pars
pars= replicate(pars[0],npars,nstep)
like= dblarr(nstep)
for ii=0L,nstep-1L do begin
    hogg_mcmc_step, seed,oldpars,oldlike,step_func,like_func,newpars,newlike, $
      log=log
    pars[*,ii]= newpars
    like[ii]= newlike
    oldpars= newpars
    oldlike= newlike
    print, $
      strjoin(strarr(21)+string(byte(8)),''), $
      'HOGG_MCMC: ',100L*ii/nstep,' percent', $
      format= '($,A,A,I2.2,A)'
endfor
print, strjoin(strarr(21)+string(byte(8)),'')+'HOGG_MCMC: done      '
sindx= reverse(sort(like))
pars= pars[*,sindx]
like= like[sindx]
return
end
