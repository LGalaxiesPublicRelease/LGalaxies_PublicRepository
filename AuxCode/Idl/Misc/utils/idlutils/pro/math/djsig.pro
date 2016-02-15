; wrapper for djs_iterstat
; D. Finkbeiner 14 Oct 1999
; DPF 27 Jun 2000 - added keywords

function djsig, x, sigrej=sigrej, maxiter=maxiter

  djs_iterstat, x, sigma=sigma, sigrej=sigrej, maxiter=maxiter

  return,sigma
end
