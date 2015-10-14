## functions for determining the posterior and prior probs are
## contained here:
## source("deb_post_prior.R")

## Bayesian inference for a deterministic DEB model (with models
## solved via an ode solver in the function specified with argument
## "sim") with an observation model, and with observation error
## variances specified in sds.



#' deb.mcmc
#'
#' @param N
#' @param p.start
#' @param data
#' @param w.p
#' @param params
#' @param inits
#' @param sim
#' @param sds
#' @param hyper
#' @param prop.sd tuning parameters, that is the standard deviation of the proposal distribution for each parameter
#' @param Tmax
#' @param cnt
#' @param burnin
#' @param plot
#' @param sizestep
#' @param w.t
#' @param which which ode solver to use. 1=desolve 2=PBSdesolve
#' @param data.times numeric vector of timepoints at which to solve the DEB model and evaluate the likelihoods. must match the timepoints for which observations are available
#' @param free.inits string, name of function to evaluate
#'
#' @return returns sth
#'
#'
#' @examples example
deb.mcmc<-function(N, p.start, data, w.p, params, inits, sim=DEB1,
                   sds, hyper, prop.sd, Tmax, cnt, burnin=0.1,
                   plot=TRUE, sizestep=0.01, w.t=1, which=1, data.times=NULL, free.inits=NULL)
{

  ## first determine the number of parameters to be infered
  np<-length(p.start)
  if(np!=length(w.p)) {
    print("not properly initializing the mcmc")
    break
  }

  ## build a data frame to hold the posterior samples
  if(np==1) samps<-data.frame(matrix(0, ncol=np+1, nrow=N+1))
  else samps<-data.frame(matrix(0, ncol=np, nrow=N+1))
  names(samps)<-w.p

  ## initialize the data frame with the starting values passed into
  ## the function as p.start
  if(np==1) samps[1,]<-c(p.start, NULL)
  else samps[1,]<-p.start

  ## sample parameters (current and proposed) are stored in p.old and
  ## p.new, respectively. Here, p.old is initialized with the values
  ## passed as params, and then the values for the parameters that are
  ## going to be infered are changed to those in samps.

  p.old<-params

  ## for testing, I'll run things and see what the posterior prob of
  ## the real params are
  sim.old<-make.states(sim, p.old, inits, Tmax, which=which, sizestep, w.t, data.times=data.times)
  prob.old<-log.post.params(params, w.p, data, p.old, hyper, sim.old, sds)
  print(paste(Sys.time()," posterior likelihood of the real parameters= ", prob.old, sep=""))

  for(k in 1:np) p.old[w.p[k]]<-samps[1,k]

  ## run the data simulation to make the underlying states (for
  ## determining the likelihoods) using the parameters stored in
  ## p.old.
  if (!is.null(free.inits)){
    inits <- do.call(free.inits, args=list(from.pars = p.old), quote=T)
  }
  sim.old<-make.states(sim, p.old, inits, Tmax, which=which, sizestep, w.t, data.times=data.times)

  ## check the posterior probability to make sure you have reasonable
  ## starting values, and to initialize prob.old
  prob.old<-log.post.params(samps[1,], w.p, data, p.old, hyper, sim.old, sds)
  print(paste(Sys.time()," initial posterior likelihood = ", prob.old, sep=""))

  if(!is.finite(prob.old)){
    print("bad starting values")
    break
  }


  for(i in 1:N){
    if(i%%cnt == 0){
      print(paste(Sys.time(), " sample number", i, sep=" "))
      ##for(j in 1:np) print(paste(w.p[j], "=", samps[i,j], sep=" "))
      if(plot){
        if(np>1 ) par(mfrow=c(ceiling(np/3),3), bty="n")
        for(j in 1:np) plot(samps[0:i,j], type="l", main=w.p[j])
      }
    }

    samps[i+1,]<-samps[i,]

    for(k in 1:np){
      p.new<-p.old

      u<-runif(1)
      q<-propose(samps[i,k], prop.sd[w.p[k]])##, i)

      ## add something here to automatically reject if the proposed
      ## parameter value is outside some limit?

      #diagnostic
      #print(paste("proposing", w.p[k], " = ", q$b))

      p.new[w.p[k]]<-q$b
      samps[i+1,k]<-q$b
      #recalculate initial values from parameters
      if (!is.null(free.inits)){
        inits.new <- do.call(free.inits, args=list(from.pars = p.new), quote=T)
        #diagnostic
        #print(inits.new)
      }

      sim.new<-make.states(sim, p.new, inits.new, Tmax, which=which, sizestep, w.t, data.times=data.times)

      ## currently only calculating prob.old outside the loop, and
      ## setting prob.old<-prob.new if we accpt the draw. This cuts
      ## down on total number of calculations
      ## prob.old<-log.post.params(samps[i,k], w.p, data, p.old,
      ## hyper, sim.old, sds)
      prob.new<-log.post.params(samps[i+1,], w.p, data, p.new, hyper, sim.new, sds)
      if(i%%cnt==0) print(paste("prob.old = ", prob.old, "; prob.new = ", prob.new, sep=""))
      if(is.finite(prob.new) && is.finite(prob.old)){
        A<-exp( prob.new + q$lbak - prob.old - q$lfwd )
      }
      else{
        A<-0
        print("whoops! must have proposed outside the correct region")
      }
      if( u > A ) samps[i+1,k]<-samps[i,k]
      else{
        sim.old<-sim.new #is this assignment necessary?
        p.old<-p.new
        prob.old<-prob.new
      }
    }

  }

  ##plot(b[100:N],type="l")
  lim<-min(1, burnin*N)
  samps <- samps[lim:(N+1),]

  return(list(samps=samps))

}



#' propose
#'
#' @param b
#' @param sd
#'
#' @return list of b.new and lfwd and lbak
#'
#'
#' @examples example
propose<-function(b, sd)##, i, freq=50, size=50 )##l=5, h=6)
{
  ##b.new<-runif(1,	l/h*b, h/l*b)
  ##fwd<-dunif(b.new, l/h*b, h/l*b)
  ##bak<-dunif(b, l/h*b.new, h/l*b.new)

  ##if(i%%freq==0) sd<-sd*size
  b.new<-rnorm(1, b, sd=sd)
  fwd<-dnorm(b.new, b, sd=sd)
  bak<-dnorm(b, b.new, sd=sd)

  return(list(b=b.new, lfwd=log(fwd), lbak=log(bak)))
}










