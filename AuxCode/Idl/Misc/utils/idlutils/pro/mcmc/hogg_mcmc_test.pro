;+
; NAME:
;   hogg_mcmc_test
; PURPOSE:
;   Test the hogg_mcmc code.
; REVISION HISTORY:
;   2005-03-31  started - Hogg
;-
function hogg_mcmc_test_like, pars
common hogg_mcmc_test_block, xx,yy,yy_ivar
mean= pars[0]
chisq= total((yy-mean)*yy_ivar*(yy-mean),/double)
like= -0.5*(chisq-n_elements(yy))
return, like
end

function hogg_mcmc_test_step, seed,pars
common hogg_mcmc_test_block, xx,yy,yy_ivar
npars= n_elements(pars)
newpars= pars+randomn(seed,npars)
return, newpars
end

pro hogg_mcmc_test_setup
common hogg_mcmc_test_block, xx,yy,yy_ivar
ndata= 1
mean= 10.0
sigma= fltarr(ndata)+1.0
seed= -10
yy= mean+sigma*randomn(seed,ndata)
yy_ivar= 1.0/sigma^2
return
end

function hogg_mcmc_test_like_2d, pars
common hogg_mcmc_test_block, xx,yy,yy_ivar
model= pars[0]+pars[1]*xx
chisq= total((yy-model)*yy_ivar*(yy-model),/double)
like= -0.5*(chisq-n_elements(yy))
return, like
end

pro hogg_mcmc_test_setup_2d
common hogg_mcmc_test_block, xx,yy,yy_ivar
ndata= 40
aa= -1.0
bb= 10.0
xx= 5.0*randomu(seed,ndata)
sigma= fltarr(ndata)+1.0
seed= -10
yy= aa*xx+bb+sigma*randomn(seed,ndata)
yy_ivar= 1.0/sigma^2
return
end

pro hogg_mcmc_test
common hogg_mcmc_test_block, xx,yy,yy_ivar
set_plot, 'PS'
hogg_plot_defaults

nstep= 100000
seed= -1L
hogg_mcmc_test_setup
pars= yy[0]
hogg_mcmc, seed,'hogg_mcmc_test_step','hogg_mcmc_test_like',nstep,pars,like, $
  /log

!P.TITLE= '1d test'
yyrange= 10.0+3.0*[-1,1]
nyy= n_elements(yy)
plot, dindgen(nyy),yy,psym=6, $
  xrange=[0.0,nyy]-0.5, $
  yrange=yyrange
djs_oploterr, dindgen(nyy),yy,yerr=1.0/sqrt(yy_ivar),psym=6
xxfit= [-nyy,2*nyy]
oplot, xxfit,pars[0]+0.0*xxfit,psym=0
plot, pars,exp(like),psym=1,symsize=0.01, $
  xrange=yyrange, $
  yrange=[-0.1,1.1]*max(exp(like))
oplot, xxfit,0.0*xxfit,psym=0
hogg_plothist, pars,xrange=yyrange,npix=100

nstep= 3000
seed= -1L
hogg_mcmc_test_setup_2d
pars= dblarr(2)
pars[1]= (yy[0]-yy[1])/(xx[0]-xx[1])
pars[0]= yy[0]-xx[0]*pars[1]
hogg_mcmc, seed,'hogg_mcmc_test_step','hogg_mcmc_test_like_2d',nstep,pars, $
  like,/log

!P.TITLE= '2d test'
xxrange= [0,5]
yyrange= [0,13]
p0range= [7,13]
p1range= [-2,0]
plot, xx,yy,psym=6, $
  xrange=xxrange, $
  yrange=yyrange
djs_oploterr, xx,yy,yerr=1.0/sqrt(yy_ivar),psym=6
xxfit= xxrange
oplot, xxfit,pars[0,0]+pars[1,0]*xxfit,psym=0
plot, pars[0,*],pars[1,*],psym=1,symsize=0.01, $
  xrange= p0range, $
  yrange= p1range
oplot, xxfit,0.0*xxfit,psym=0
hogg_plothist, pars[0,*],npix=100,xrange=p0range
hogg_plothist, pars[1,*],npix=100,xrange=p1range
device,/close
return
end
