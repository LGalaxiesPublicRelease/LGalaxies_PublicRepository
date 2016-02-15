;+
; NAME:
;   hogg_mcmc_step
; PURPOSE:
;   Make one step in a Markov Chain Monte Carlo.
; INPUTS:
;   seed       - random number seed for randomu()
;   pars       - initial parameters (can be an array or structure)
;   like       - initial likelihood
;   step_func  - function that takes a step in parameter space
;   like_func  - function that computes the likelihood
; KEYWORDS:
;   log        - like_func returns log_e(likelihood), not straight
;                likelihood
; OUTPUTS:
;   newpars    - new parameters
;   newlike    - new likelihood
; BUGS:
; REVISION HISTORY:
;   2005-03-31  started - Hogg
;   2006-02-28  HUGE bug fixed - Masjedi
;-
pro hogg_mcmc_step, seed,pars,like,step_func,like_func,newpars,newlike,log=log
if (NOT keyword_set(like)) then like= call_function(like_func,pars)
newpars= call_function(step_func,seed,pars)
newlike= call_function(like_func,newpars)
if keyword_set(log) then likeratio= newlike-like $
else likeratio= newlike/like
randomnumber= randomu(seed)
if keyword_set(log) then randomnumber= alog(randomnumber)
if NOT ((newlike GT like) OR $
        (randomnumber LT likeratio)) then begin
    newpars= pars
    newlike= like
endif
return
end
