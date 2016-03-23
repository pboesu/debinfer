## functions for determining the posterior and prior probs are
## contained here:
## source("deb_post_prior.R")

## Bayesian inference for a deterministic DEB model (with models
## solved via an ode solver in the function specified with argument
## "sim") with an observation model, and with observation error
## variances specified in sds.



#' deb_mcmc
#'
#' @param N integer; Number of MCMC iterations
#' @param p.start numeric vector; start values for free parameters
#' @param data the data; need to formalize expected format
#' @param w.p character vector; names of free parameters
#' @param params named numeric; vector of all parameters, i.e. free and fixed
#' @param inits named numeric; vector of initial values of the state variable(s)
#' @param sim function; deSolve compatible specification of the differential equation system for which parameters are to be estimated
#' @param sds numeric; observation error variances
#' @param hyper; list of lists(?) containing the hyperparameters for each free parameter
#' @param pdfs; character vector containing names of prior pdfs
#' @param prop.sd numeric vector; tuning parameters, that is the standard deviation of the proposal distribution for each parameter
#' @param Tmax numeric; maximum timestep for the solver
#' @param cnt numeric; interval at which to refresh MCMC chain plots
#' @param burnin ??
#' @param plot logical; plot snapshots of MCMC chains? Useful for tuning.
#' @param sizestep; sizestep of solver [maybe change the way arguments are passed onto the solver using lists of options, and/or ellipses]
#' @param w.t ???
#' @param which which ode solver to use. 1=desolve 2=PBSdesolve
#' @param data.times numeric vector of timepoints at which to solve the DEB model and evaluate the likelihoods. must match the timepoints for which observations are available
#' @param free.inits string, name of function to evaluate
#'
#' @return returns sth
#'
#'
#' @export
deb_mcmc<-function(N, p.start, data, w.p, params, inits, sim=DEB1,
                   sds, hyper, pdfs, prop.sd, Tmax, cnt, burnin=0.1,
                   plot=TRUE, sizestep=0.01, w.t=1, which=1, data.times=NULL,
                   free.inits=NULL, obs.model)
{
  #check observation model function
  if (!is.function(obs.model)) stop("obs.model must be a function")
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
  prob.old<-log_post_params(samp=params, w.p=w.p, data=data, p=p.old, pdfs=pdfs, hyper=hyper, sim.data=sim.old, sds=sds, obs.model=obs.model)
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
  prob.old <- log_post_params(samp=samps[1,], w.p=w.p, data=data, p=p.old, pdfs=pdfs, hyper=hyper, sim.data=sim.old, sds=sds, obs.model=obs.model)
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
      } else {inits.new <- inits}

      sim.new<-make.states(sim, p.new, inits.new, Tmax, which=which, sizestep, w.t, data.times=data.times)

      ## currently only calculating prob.old outside the loop, and
      ## setting prob.old<-prob.new if we accpt the draw. This cuts
      ## down on total number of calculations
      ## prob.old<-log_post_params(samps[i,k], w.p, data, p.old,
      ## hyper, sim.old, sds)
      prob.new <- log_post_params(samp=samps[i+1,], w.p=w.p, data=data, p=p.new, pdfs=pdfs, hyper=hyper, sim.data=sim.new, sds=sds, obs.model=obs.model)
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
  #this piece of the function doesn't make sense  to me
  #essentially the burnin argument is ignored
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




#' log_prior_params
#'
#'
#'
#'
#' @export
log_prior_params<-function(samp, pdfs, w.p, hyper){
  lp<-0
  len<-length(w.p)
  if(len==1){
    ##print(paste("w.p =", w.p, "samp =", samp[w.p[1]], sep=" "))
    lp<-list(NULL)
    names(lp)<-w.p
  }
  else{
    lp<-data.frame(matrix(0, nrow=1, ncol=len))
    names(lp)<-w.p
  }

  ##print(c(as.numeric(samp), w.p, hyper[1]))
  for(i in 1:len){
    p<-w.p[i]
    s<-as.numeric(samp[p])

    lp[[p]] <- logd_prior(s, pdfs[[p]], hypers=hyper[[p]])

  }

  return(lp)

}

#' log_post_params
#'
#'
#'
#'
#' @export
log_post_params <- function(samp, w.p, data, p, pdfs, hyper, sim.data, sds, verbose.lik=FALSE, obs.model){

  log_data <- obs.model

  llik <- log_data(data=data, sim.data=sim.data, samp=samp, sds=sds)

  #if(length(w.p)==1) lprior<-as.numeric(log_prior_params(samp, w.p, hyper))
  #else {
    lprior<-sum(log_prior_params(samp, pdfs, w.p, hyper))
    ##if(!is.finite(lprior)) break
  #}
  ##print(c(b, lik, prior))

  if(is.na(llik)) break
  if(is.na(lprior)) break

  return( llik + lprior )
}

