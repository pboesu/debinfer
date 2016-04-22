### revised version of the inference function

##' Bayesian inference for a deterministic differential equation model
##'  (with models solved via an ode or dde solver in the function specified with argument
##' "de.model") with an observation model, and with observation error
##' variances specified in sds.
##'
##'
##' @title de_mcmc
##' @param N
##' @param data
##' @param all.params
##' @param inits
##' @param sim
##' @param samp.p
##' @param Tmax
##' @param cnt
##' @param burnin
##' @param plot
##' @param sizestep
##' @param w.t
##' @param which
##' @param test
##' @param my.par
##' @param myswitch
##' @param mymap
##' @return
##' @author Philipp Boersch-Supan
##' de_mcmc(N = 5000, data = N_obs, sim = logistic_model, obs.model = logistic_obs_model, r, K, sigma_obs, ...
de_mcmc <- function(N, data, de.model, obs.model, all.params,
                   Tmax, cnt=10, burnin=0.1,
                   plot=TRUE, sizestep=0.01, w.t=1, which=2,
                   test=TRUE, my.par=c(1,1),  myswitch=NULL,
                   mymap=NULL, ...)

{ #right now this is just a wrapper for the old function, reassigning inputs from the debinfer_parlist object to dde_mcmc
  p.start = list()
  w.p =
  params
  inits
  sds
  hyper = list
  pdfs =
  prop.sd
  data.times

  samps <- deb_mcmc(N=N, p.start=list(r=0.5, K=5, loglogsd.N=-2), data=data, w.p=c("r", "K", "loglogsd.N"), params=parms,
           inits=c(N=0.1), sim=de.model, sds=list(N=0.01), hyper=list(r=list(mean=0, sd=1),K=list(meanlog=1, sdlog=1), loglogsd.N=list(mean=-2, sd=1)),
           pdfs = list(r='norm', K='lnorm', loglogsd.N='norm'),
           prop.sd=c(r=0.005, K=0.1, loglogsd.N=0.5),
           Tmax=Tmax, cnt=cnt, burnin=burnin, plot=plot, sizestep=sizestep, which = which,
           data.times = N_obs$time, obs.model=obs.model)
  return(samps)

}



