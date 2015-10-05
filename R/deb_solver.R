## This file holds many of the functions needed to simulate the DEB
## model forward in time. It needs the package PBSddesolve to function
## properly.


##install PBS solve with
##sudo R CMD INSTALL PBSddesolve_1.08.9_bobby.tar.gz
library(PBSddesolve)


solve.DEB<-function(sim, params, inits, Tmax=400, numsteps=10000,
                    which=1, sizestep=NULL){

  if(is.null(sizestep)) times<- seq(0, Tmax, length=numsteps)
  if(is.null(numsteps))  times<- seq(0, Tmax, by=sizestep)
  if(which==1){
    require(odesolve)
    out<-lsoda(inits, times, sim, parms=params, verbose=TRUE)
  }
  if(which==2){
    require(PBSddesolve)
    #on.exit(freeglobaldata())
    out<-dde(inits, times, sim, parms=params)
  }
  return(out)
}



test.sim <- function(sim, params, inits, ylim=c(0,600), Tmax=200,
                  numsteps=10000, which=2, scale=10){

  out<-solve.DEB(sim, params, inits, Tmax, numsteps, which)

  plot.DEB(out, scale)
  return(out)
 
}


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


plot.DEB.red<-function(out){
  par(mfrow=c(2,1))
   plot(out$t,out$L, lty=5, col="blue",
        lwd=2, xlab="time", ylab="Length", pch=19)
   plot(out$t,out$Negg, col="red",  lty=4, lwd=2, xlab="time",
        ylab="Number of Eggs", pch=19)
}




## The following functions extract time series observations from the
## forward simulation. In particular, it extracts only those data
## points which are zero mod w.t

extract.data<-function(data, w.t=1, Tmax){

  ww<-which(data$time%%w.t==0)
  ##print(ww)
  if(length(ww)!=(Tmax+1)){
    print("something is wrong")
    break
  }
  else data<-data[ww,]
  return(data)
  
}

## This function takes data from a forward simulation, extracts a
## subset of the data points using extract.data, and adds
## observational noise using add.noise

make.obs<-function(dt, sds, params, w.t, Tmax){

  ##print(w.t)
  dt<-extract.data(dt, w.t, Tmax)
  #print(head(data))
  #return(data)
  dt<-add.noise(dt, sds, params)
  return(dt)
  
}


## This function takes the simulator, params, inits, etc, and runs a
## forward simulation, and then extracts a subset of the points from
## the simulator (determined by w.t), without noise

make.states<-function(sim=DEB1, params, inits, Tmax, which=2, sizestep=0.01, w.t=1){

  ##dt<-0
  ##print(dt)
  dt<-solve.DEB(sim, params, inits, Tmax, numsteps=NULL, which, sizestep)
  
  dt<-extract.data(dt, w.t, Tmax)

  return(list(t=dt$t, l=dt$y3, n=dt$y5))
  
}


