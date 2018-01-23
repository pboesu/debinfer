## This file holds the functions needed to simulate the DE
## model forward in time.


#' solve_de
#'
#' @import deSolve
#' @import PBSddesolve
#' @param sim function; solver compatible specification of the DE
#' @param params numeric; named vector of parameter values
#' @param inits numeric; initial values. Must be in the same order as specified within sim!
#' @param Tmax numeric; maximum timestep
#' @param numsteps numeric
#' @param solver Choice of solver to use 1 or "ode" = deSolve::ode, 2 or "dde" = PBSddesolve::dde, 3 or "dede" = deSolve::dede
#' @param sizestep for solver
#' @param method solver method
#' @param verbose passed to deSolve::ode
#' @param data.times numeric a vector of times at which the ODE is to be evaluated. Defaults to NULL. If value is supplied it takes precedence over any value supplied to \code{numsteps} or \code{sizesteps}.
#' @param ... additional arguments to solver
#'
#'
#' @return integrated ode object. Data structure depends on the employed solver.
#'
#' @export
solve_de<-function(sim, params, inits, Tmax, numsteps=10000,
                    solver="ode", sizestep=NULL, verbose=FALSE, data.times=NULL, method = "lsoda", ...){
    #TODO: recalculate inits from parameters if necessary
    if( !(solver %in% c(1,2,3) || solver %in% c("ode","dde","dede") )){
      stop(paste('invalid value [', solver, '] for argument "solver". Must be one of "ode","dde","dede". To set a specific solver algorithm use the "method" argument (see ?ode or ?dede).'))
      }
    if(!is.null(data.times)){
      #this is fragile. really the data should be in a class that ensures proper times, no missing data etc. pp. Also this now assumes observations at identical times for all observed variables.
      times <- data.times #this shouldn't be done every time the solver is called. solver times should be set up at the start of the mcmc procedure and then passed in through times argument
    } else { # currently there's no proper handling of multiple non-NULL arguments/all NULL arguments
      if(is.null(sizestep)) times<- seq(0, Tmax, length=numsteps)
      if(is.null(numsteps))  times<- seq(0, Tmax, by=sizestep)
    }
      if(solver == 1 || solver == "ode"){
        #require(deSolve)
        out <- try(ode(inits, times, sim, parms=params, verbose=verbose , method=method, ...))
      }
      if(solver == 2 || solver == "dde"){
        #require(PBSddesolve)
        #on.exit(freeglobaldata())
        out <- try(PBSddesolve::dde(inits, times, sim, parms=params, ...))
      }
      if(solver == 3 || solver == "dede"){
        #require(PBSddesolve)
        #on.exit(freeglobaldata())
        out <- try(dede(y=inits, times=times, func=sim, parms=params))
      }
  return(out)
}




#' post_sim
#'
#' @import deSolve
#' @import PBSddesolve
#' @param x debinfer_result object
#' @param n number of simulations
#' @param output character, "sims", "all", "HDI"
#' @param burnin integer, number of samples to discard from the start of the mcmc chain
#' @param times numeric a vector of times at which the ODE is to be evaluated. Defaults to NULL.
#' @param prob A numeric scalar in the interval (0,1) giving the target probability content of the intervals. The nominal probability content of the intervals is the multiple of 1/nrow(obj) nearest to prob.
#' @param ... additional arguments to solver
#' @return a post_sim object containing a list of de solutions or summaries thereof
#' @import coda
#' @import plyr
#' @export
post_sim<-function(x, n=100, times, output = "all" , burnin = NULL, prob = 0.95, ...)
{
#sample parameter values
  if(!is.null(burnin)) {
    #get number of samples of mcmc object
    mcmc_end <- attr(x$samples, "mcpar")[2]
  x$samples <- window(x$samples, burnin, mcmc_end)
  }
  #reshape and relabel in case there is only one free paramater (mcmc object is vector, not matrix)
  if(!is.matrix(x$samples)) {x$samples <- as.matrix(x$samples)
  colnames(x$samples) <- names(freeparams(x))
  }
  samps <- as.matrix(x$samples[sample(nrow(x$samples), size=n, replace=FALSE),])
  colnames(samps) <- colnames(x$samples)

  #apply restore_and_solve over the samples
  if (n==1){
    sims <- restore_and_solve(x = x, samp = samps, times = times, ...)
  } else {
    sims <- plyr::alply(samps, 1, function(samp) restore_and_solve(x = x, samp = samp, times = times, ...) )
  }
  class(sims) <- "post_sim"

  if (output == "HDI" | output=="all"){
    newlist <- reshape_post_sim(sims)
    HDI <- llply(newlist[2:length(newlist)], function(x) HPDinterval(as.mcmc(x), prob = prob))
    time <- newlist$time
    medianlist <- llply(newlist[2:length(newlist)], function(x) apply(x,2,median))
  }
  if (output == "sims") return(sims)
  if (output == "HDI") return(HDI)
  if (output == "all"){
    out <- list(sims=sims, HDI = HDI, median=medianlist, time = time)
    class(out) <- "post_sim_list"
    return(out)
  }
  }


#' an internal function to recombine inits and pars
#'
#' @param x debinfer_result object
#' @param samp sample
#' @param times times
#' @param ... additional arguments to solver
#'
#' @return solve_de output
#' @export
restore_and_solve <- function(x, samp, times, ...){
  params <- depars(x)
  inits <- deinits(x)
  #paste in sample values
  for (i in names(samp)){
    if (i %in% names(params)) params[i]<-samp[i]
    if (i %in% names(inits)) inits[i]<-samp[i]
  }
  #TODO recalc inits from pars if any init is designated as such
  #check if any initial values need to be calculated
  var_types <- vapply(x$all.params, function(y) y$var.type, character(1))
  if (any(var_types == "initfunc")){
    initfunc_idx <- which(var_types == "initfunc")
    initfunc_name <- names(x$all.params[initfunc_idx])
    initfunc <- x$all.params[[initfunc_idx]]$initfunc
    #recalculate inits
    inits <- initfunc(inits, params)
  }
  #solve DE model
  soln <- solve_de(sim = x$de.model, params = params, inits = inits, data.times = times, solver = x$solver, ...)
  return(soln)
}





