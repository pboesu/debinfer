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



  mcmc_samples <- deb_mcmc(N=N, p.start=p.start, data=data, w.p=w.p, params=params,
                           inits=inits, sim=de.model, sds=NULL,
                           hyper=hyper,
                           pdfs = pdfs,
                           prop.sd=prop.sd,
                           Tmax=Tmax, cnt=cnt, burnin=burnin, plot=plot, sizestep=sizestep, which = which,
                           data.times = data.times, obs.model=obs.model, verbose = verbose, ...)
  return(mcmc_samples)

}


###########################################
###### code from chytrid branch here ######
###########################################


##' Bayesian inference for a deterministic DEB model (with models
##' solved via an ode solver in the function specified with argument
##' "sim") with an observation model, and with observation error
##' variances specified in sds.
##'
##'
##' @title de_mcmc_rev
##' @param N
##' @param data
##' @param all.params
##' @param ref.params an optional named vector containing a set of reference parameters, e.g. the true paramaters underlying a simulated data set
##' @param ref.inits  an optional named vector containing a set of reference initial values, e.g. the true initial values underlying a simulated data set
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
de_mcmc_rev <- function(N, data, de.model, obs.model, all.params, ref.params=NULL, ref.inits=NULL,
                              Tmax, data.times, cnt=10, burnin=0.1,
                              plot=TRUE, sizestep=0.01, which=1,
                              myswitch=NULL,
                              mymap=NULL, verbose =FALSE, ...)
{
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


  ## first we calculate a few lengths, etc, that we use for the for
  ## loops later.

  l<-length(samp.p)
  np<-length(all.params)

  ltot<-0
  for(i in 1:l) ltot<-ltot+length(samp.p[[i]]$params)

  ## for testing, here is code that calculates (and prints out) the
  ## posterior prob of the reference parameters, which can be passed in
  ## through ref.params

  if(!is.null(ref.params) && !is.null(ref.inits)){
    sim.ref<-make.states(sim, ref.params, ref.inits, Tmax, which=which, sizestep, w.t,
                         myswitch=myswitch, mymap=mymap, ...)
    prob.ref<-log.post.params(ref.params, data, samp.p, sim.ref)
    print(paste("(unnormalized) posterior probability of the reference parameters= ",
                prob.old, sep=""))
  }

  ## build an mcmc object to hold the posterior samples. These include
  ## "samples" even for the parameters that are being held fixed, to
  ## make the code more straightforward and generic, as well as being
  ## useful for debugging. The samps structure will also keep track of
  ## the log posterior probability of that particular sample.

  samps <- mcmc(matrix(numeric((np+1)*N), ncol=np+1, nrow=N, dimnames = list(NULL, c(p.names, "lpost"))))

  return(samps)
}
