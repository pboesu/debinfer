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
make.states<-function(sim, params, inits, Tmax, which=1, sizestep=0.01, w.t=1, data.times=NULL, method = "lsoda", ...){

  ##dt<-0
  ##print(dt)
  dt<-solve_de(sim, params, inits, Tmax, numsteps=NULL, which, sizestep, data.times=data.times, method = method, ...)

  #this is a dirty solution to the situation that the ode solver terminates prematurely. all remaining timepoints are set to 0, which should badly impact the likelihood (although it might still be better than numerically valid but "bad" solutions that don't end prematurely). Ideally the MCMC sampler should reject a sample if the ode solver stops prematurely

  if (max(dt[,"time"]) < Tmax) {
    print(paste("integration failed prematurely at t = ", max(dt[,"time"])))
    if (!is.null(data.times)){
      new.time <- data.times
    } else {new.time <- seq(0,Tmax,by=sizestep)} #breaks when sizestep is NULL and numstep is specified instead
    #0 as "penalized value" for premature termination yields weird results
    fake.solution <- matrix(NA,nrow = length(new.time),ncol = dim(dt)[2])
    dimnames(fake.solution) <- dimnames(dt)
    fake.solution[,'time'] <- new.time
    matched <- which(fake.solution[,'time'] %in% dt[,'time'])
    fake.solution[matched,2:ncol(fake.solution)]<-dt[matched,2:ncol(dt)]
    dt <- fake.solution
  }

  if (is.null(data.times)) dt<-extract.data(dt, w.t, Tmax) #only evaluate this when no specific time points are given. again this is a dirty fix. all of this should be dealt with inside the solver function, so that solve.DEB always returns a well formed data structure containing all expected timesteps and then NAs/penalized values/NaNs as necessary.

  return(as.list(as.data.frame(dt))) # this returns a list of the named vectors

}



