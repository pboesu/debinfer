## This file holds many of the functions needed to simulate the DEB
## model forward in time. It needs the package PBSddesolve to function
## properly.


##install PBS solve with
##sudo R CMD INSTALL PBSddesolve_1.08.9_bobby.tar.gz
#library(PBSddesolve)


#' solve_de
#'
#' @import deSolve
#' @import PBSddesolve
#' @param sim function; solver compatible specification of the DE
#' @param params numeric; named vector of parameter values
#' @param inits numeric; initial values. Must be in the same order as specified within sim!
#' @param Tmax numeric; maximum timestep
#' @param numsteps numeric
#' @param which interger. whihch solver to use 1= deSolve, 2 = PBSddesolve
#' @param sizestep for solver
#' @param verbose passed to deSolve::ode
#' @param data.times numeric a vector of times at which the ODE is to be evaluated. Defaults to NULL.
#' If value is supplied it takes precendence over any value supplied to \code{numsteps} or \code{sizesteps}.
#'
#' @return integrated ode - describe structure
#'
#'
#' @examples example
#' @export
solve_de<-function(sim, params, inits, Tmax=400, numsteps=10000,
                    which=1, sizestep=NULL, verbose=FALSE, data.times=NULL){

    if(!is.null(data.times)){
      #this is fragile. really the data should be in a class that ensures proper times, no missing data etc. pp. Also this now assumes observations at identical times for all observed variables.
      times <- data.times
    } else { # currently there's no proper handling of multiple non-NULL arguments/all NULL arguments
      if(is.null(sizestep)) times<- seq(0, Tmax, length=numsteps)
      if(is.null(numsteps))  times<- seq(0, Tmax, by=sizestep)
    }

    if(which==1){
    #require(deSolve)
    out<-ode(inits, times, sim, parms=params, verbose=verbose , method='lsoda')
  }
  if(which==2){
    #require(PBSddesolve)
    #on.exit(freeglobaldata())
    out<-dde(inits, times, sim, parms=params)
  }
  return(out)
}



#' test_sim
#'
#' @param sim
#' @param params
#' @param inits
#' @param ylim
#' @param Tmax
#' @param numsteps
#' @param which
#' @param scale
#'
#' @return simulated ODE and a plot
#'
#'
#' @examples example
test_sim <- function(sim, params, inits, ylim=c(0,600), Tmax=200,
                  numsteps=10000, which=2, scale=10){

  out<-solve_de(sim, params, inits, Tmax, numsteps, which)

  plot.DEB(out, scale)
  return(out)

}


#' Title
#'
#' @param out
#' @param scale
#' @param scaled.length
#'
#' @return plot
#'
#'
#' @examples plot
plot.DEB<-function(out, scale=100, scaled.length=TRUE){

  par(mfrow=c(2,2))
   plot(out[,1],out[,2], type="l", lty=5, col="blue",
        lwd=2, xlab="time", ylab="Food")
   plot(out[,1],out[,3], type="l", col="red",  lty=4, lwd=2, xlab="time",
        ylab="scaled reserves")
  if(scaled.length){
    plot(out[,1],out[,4], type="l", col="green",lwd=2, xlab="time",
         ylab="scaled length")
  }
  else{
    plot(out[,1],out[,4], type="l", col="green",lwd=2, xlab="time",
         ylab="length")
  }
   plot(out[,1],out[,6]/scale, type="l", col=1, xlab="time",
        ylab=paste("maturity & reproduction/", scale, sep=""))
   lines(out[,1],out[,5], lty=2, col=2)
}


#' Title
#'
#' @param out
#'
#' @return plots
#'
#'
#' @examples plot
plot.DEB.red<-function(out){

  par(mfrow=c(2,1))
   plot(out$t,out$L, lty=5, col="blue",
        lwd=2, xlab="time", ylab="Length", pch=19)
   plot(out$t,out$Negg, col="red",  lty=4, lwd=2, xlab="time",
        ylab="Number of Eggs", pch=19)
}




#' Title
#'
#' The following functions extract time series observations from the
#' forward simulation. In particular, it extracts only those data
#' points which are zero mod w.t
#'
#' @param data
#' @param w.t
#' @param Tmax
#'
#' @return data
#'
#'
#' @examples examples
extract.data <- function(data, w.t=1, Tmax){
  #need to assert here, that the time vector is actually evaluated at integer values, which may not be the case when it is created using the numsteps argument in solve.DEB
  ww<-which(data[,'time']%%w.t==0) #function breaks here when handling deSolve output
  ##print(ww)
  if(length(ww)!=(Tmax+1)){
    print("something is wrong")
    break
  }
  else data<-data[ww,]
  return(data)

}



#' Title
#'
#' This function takes data from a forward simulation,
#' extracts a subset of the data points using extract.data,
#' and adds observational noise using add.noise
#'
#' @param dt
#' @param sds
#' @param params
#' @param w.t
#' @param Tmax
#'
#' @return dt
#'
#'
#' @examples examples
make.obs<-function(dt, sds, params, w.t, Tmax, ode.pars){

  ##print(w.t)
  dt<-extract.data(dt, w.t, Tmax)
  #print(head(data))
  #return(data)
  dt<-add.noise(dt, sds, params, ode.pars)
  return(dt)

}



#' Title
#'
#'This function takes the simulator, params, inits, etc, and runs a
#'forward simulation, and then extracts a subset of the points from
#'the simulator (determined by w.t), without noise
#'
#' @param sim
#' @param params
#' @param inits
#' @param Tmax
#' @param which
#' @param sizestep
#' @param w.t
#' @param data.times numeric passed on to \code{solve.DEB}
#'
#' @return a subset of the points from the simulator (determined by w.t), without noise
#'
#'
# #' @examples EXAMPLES
make.states<-function(sim, params, inits, Tmax, which=1, sizestep=0.01, w.t=1, data.times=NULL){

  ##dt<-0
  ##print(dt)
  dt<-solve_de(sim, params, inits, Tmax, numsteps=NULL, which, sizestep, data.times=data.times)

  #this is a dirty solution to the situation that the ode solver terminates prematurely. all remaining timepoints are set to 0, which should badly impact the likelihood (although it might still be better than numerically valid but "bad" solutions that don't end prematurely). Ideally the MCMC sampler should reject a sample if the ode solver stops prematurely

  if (max(dt[,"time"]) < Tmax) {
    print(paste("integration failed prematurely at t = ", max(dt[,"time"])))
    if (!is.null(data.times)){
      new.time <- data.times
    } else {new.time <- seq(0,Tmax,by=sizestep)} #breaks when sizestep is NULL and numstep is specified instead
    #0 as "penalized value" for premature termination yields weird results
    fake.solution <- matrix(1e99,nrow = length(new.time),ncol = dim(dt)[2])
    dimnames(fake.solution) <- dimnames(dt)
    fake.solution[,'time'] <- new.time
    matched <- which(fake.solution[,'time'] %in% dt[,'time'])
    fake.solution[matched,2:ncol(fake.solution)]<-dt[matched,2:ncol(dt)]
    dt <- fake.solution
  }

  if (is.null(data.times)) dt<-extract.data(dt, w.t, Tmax) #only evaluate this when no specific time points are given. again this is a dirty fix. all of this should be dealt with inside the solver function, so that solve.DEB always returns a well formed data structure containing all expected timesteps and then NAs/penalized values/NaNs as necessary.

  return(as.list(as.data.frame(dt))) # this returns a list of the named vectors

}


