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
##'
##' @export
##'
de_mcmc <- function(N, data, de.model, obs.model, all.params,
                   Tmax, data.times, cnt=10, burnin=0.1,
                   plot=TRUE, sizestep=0.01, which=1,
                   myswitch=NULL,
                   mymap=NULL, verbose =FALSE, ...)

{ #right now this is just a wrapper for the old function, reassigning inputs from the debinfer_parlist object to dde_mcmc
  ##remapping here
  # get all parameter names
  p.names <- sapply(all.params, function(x) x$name)
  is.free <- !sapply(all.params, function(x) x$fixed)
  is.init <- sapply(all.params, function(x) x$var.type)=="init"

  #inits are matched by order in deSolve --> how to ensure they are in the correct order?? add init.order variable to parameter declaration?



  #get all start values
  p.start <- lapply(all.params, function(x) x$value)[is.free]
  names(p.start) = p.names[is.free]
  # get the parameters that are to be estimated
  w.p <- p.names[is.free]
  #check what this is needed for, except for the "true" likelihood calculation
  params <- unlist(lapply(all.params, function(x) x$value))
  names(params) <-  p.names
  #initial values for DE (no ordering!)
  inits <- sapply(all.params, function(x) x$value)[is.init]
  names(inits) <- p.names[is.init]

 # sds <- NULL# is this obsolete if it is not used in the obs model?
  hyper = lapply(all.params, function(x) x$hyper)[is.free]
  names(hyper) <- p.names[is.free]
  pdfs = lapply(all.params, function(x) x$prior)[is.free]
  names(pdfs) = p.names[is.free]

  prop.sd <- sapply(all.params[is.free], function(x) x$prop.var)
  names(prop.sd) <- p.names[is.free]



  mcmc_samples <- deb_mcmc(N=N, p.start=p.start, data=data, w.p=w.p, params=parms,
                           inits=inits, sim=de.model, sds=NULL,
                           hyper=hyper,
                           pdfs = pdfs,
                           prop.sd=prop.sd,
                           Tmax=Tmax, cnt=cnt, burnin=burnin, plot=plot, sizestep=sizestep, which = which,
                           data.times = data.times, obs.model=obs.model, verbose = verbose, ...)
  return(mcmc_samples)

}



