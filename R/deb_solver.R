## This file holds many of the functions needed to simulate the DEB
## model forward in time. It needs the package PBSddesolve to function
## properly.


##install PBS solve with
##sudo R CMD INSTALL PBSddesolve_1.08.9_bobby.tar.gz
#library(PBSddesolve)


#' Title
#'
#' @import deSolve
#' @param sim
#' @param params
#' @param inits
#' @param Tmax
#' @param numsteps
#' @param which
#' @param sizestep
#' @param verbose passed to deSolve::ode
#'
#' @return integrated ode - describe structure
#'
#'
#' @examples example
solve.DEB<-function(sim, params, inits, Tmax=400, numsteps=10000,
                    which=1, sizestep=NULL, verbose=TRUE){

  if(is.null(sizestep)) times<- seq(0, Tmax, length=numsteps)
  if(is.null(numsteps))  times<- seq(0, Tmax, by=sizestep)
  if(which==1){
    require(deSolve)
    out<-ode(inits, times, sim, parms=params, verbose=verbose ,method='lsoda')
  }
  if(which==2){
    require(PBSddesolve)
    #on.exit(freeglobaldata())
    out<-dde(inits, times, sim, parms=params)
  }
  return(out)
}



#' Title
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
test.sim <- function(sim, params, inits, ylim=c(0,600), Tmax=200,
                  numsteps=10000, which=2, scale=10){

  out<-solve.DEB(sim, params, inits, Tmax, numsteps, which)

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
extract.data<-function(data, w.t=1, Tmax){

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
#'
#' @return a subset of the points from the simulator (determined by w.t), without noise
#'
#'
# #' @examples EXAMPLES
make.states<-function(sim=DEB1, params, inits, Tmax, which=2, sizestep=0.01, w.t=1){

  ##dt<-0
  ##print(dt)
  dt<-solve.DEB(sim, params, inits, Tmax, numsteps=NULL, which, sizestep)

  dt<-extract.data(dt, w.t, Tmax)

  return(list(t=dt$t, l=dt$y3, n=dt$y5))

}


