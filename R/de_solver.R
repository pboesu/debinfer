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
#' @param data.times numeric a vector of times at which the ODE is to be evaluated. Defaults to NULL.
#' @param ... additional arguments to solver
#' If value is supplied it takes precendence over any value supplied to \code{numsteps} or \code{sizesteps}.
#'
#' @return integrated ode - describe structure
#'
#'
#' @examples example
#' @export
solve_de<-function(sim, params, inits, Tmax, numsteps=10000,
                    solver="ode", sizestep=NULL, verbose=FALSE, data.times=NULL, method = "lsoda", ...){

    if(!is.null(data.times)){
      #this is fragile. really the data should be in a class that ensures proper times, no missing data etc. pp. Also this now assumes observations at identical times for all observed variables.
      times <- data.times #this shouldn't be done every time the solver is called. solver times should be set up at the start of the mcmc procedure and then passed in through times argument
    } else { # currently there's no proper handling of multiple non-NULL arguments/all NULL arguments
      if(is.null(sizestep)) times<- seq(0, Tmax, length=numsteps)
      if(is.null(numsteps))  times<- seq(0, Tmax, by=sizestep)
    }

    if(solver == 1 || solver == "ode"){
    #require(deSolve)
    out<-ode(inits, times, sim, parms=params, verbose=verbose , method=method, ...)
  }
  if(solver == 2 || solver == "dde"){
    #require(PBSddesolve)
    #on.exit(freeglobaldata())
    out<-PBSddesolve::dde(inits, times, sim, parms=params, ...)
  }
  if(solver == 3 || solver == "dede"){
    #require(PBSddesolve)
    #on.exit(freeglobaldata())
    out <- dede(y=inits, times=times, func=sim, parms=params)
    out <- as.data.frame(out)
  }
  return(out)
}




#' post_sim
#'
#' @import deSolve
#' @import PBSddesolve
#' @param x debinfer_result object
#' @param n number of simulations
#' @param output character, "median", "mult", hdi
#' @param times numeric a vector of times at which the ODE is to be evaluated. Defaults to NULL.
#' @param ... additional arguments to solver
#'
#'
#' @return integrated ode - describe structure
#'
#'
#' @examples example
#' @export
post_sim<-function(result, n,
                   sizestep=NULL, verbose=FALSE, data.times=NULL, method = "lsoda", ...)
{

}





