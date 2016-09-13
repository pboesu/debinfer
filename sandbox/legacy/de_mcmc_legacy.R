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
#' @param method ode method, see ?deSolve::ode
#'
#' @return returns sth
#'
#'
#' @export
deb_mcmc<-function(N, p.start, data, w.p, params, inits, sim=DEB1,
                   sds, hyper, pdfs, prop.sd, Tmax, cnt, burnin=0.1,
                   plot=TRUE, sizestep=0.01, w.t=1, which=1, data.times=NULL,
                   free.inits=NULL, obs.model, verbose=FALSE, method = "lsoda")
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
  sim.old<-make.states(sim, p.old, inits, Tmax, which=which, sizestep, w.t, data.times=data.times, method=method)
  prob.old<-log_post_params(samp=params, w.p=w.p, data=data, p=p.old, pdfs=pdfs, hyper=hyper, sim.data=sim.old, sds=sds, obs.model=obs.model)
  print(paste(Sys.time()," posterior likelihood of the real parameters= ", prob.old, sep=""))

  for(k in seq_len(np)) p.old[w.p[k]]<-samps[1,k]

  ## run the data simulation to make the underlying states (for
  ## determining the likelihoods) using the parameters stored in
  ## p.old.
  if (!is.null(free.inits)){
    inits <- do.call(free.inits, args=list(from.pars = p.old), quote=T)
  }
  sim.old<-make.states(sim, p.old, inits, Tmax, which=which, sizestep, w.t, data.times=data.times, method=method)

  ## check the posterior probability to make sure you have reasonable
  ## starting values, and to initialize prob.old
  prob.old <- log_post_params(samp=samps[1,], w.p=w.p, data=data, p=p.old, pdfs=pdfs, hyper=hyper, sim.data=sim.old, sds=sds, obs.model=obs.model)
  print(paste(Sys.time()," initial posterior likelihood = ", prob.old, sep=""))

  if(!is.finite(prob.old)){
    print("bad starting values")
    break
  }


  for(i in seq_len(N)){
    if(i%%cnt == 0){
      print(paste(Sys.time(), " sample number", i, sep=" "))
      ##for(j in 1:np) print(paste(w.p[j], "=", samps[i,j], sep=" "))
      if(verbose) print(samps[i,])
      if(plot){
        if(np>1 ) par(mfrow=c(ceiling(np/3),3), bty="n")
        for(j in seq_len(np)) plot(samps[0:i,j], type="l", main=w.p[j])
      }
    }

    samps[i+1,]<-samps[i,]

    for(k in seq_len(np)){
      p.new<-p.old

      u<-runif(1)
      q<-propose(b=samps[i,k], sd=prop.sd[w.p[k]])##, i)

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

      sim.new<-make.states(sim, p.new, inits.new, Tmax, which=which, sizestep, w.t, data.times=data.times, method=method)

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








###########################################
###### code from chytrid branch here ######
###########################################


##' Bayesian inference for a deterministic DEB model (with models
##' solved via an ode solver in the function specified with argument
##' "sim") with an observation model, and with observation error
##' variances specified in sds.
##'
##'
##' @title dde_mcmc
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
dde_mcmc<-function(N, data, all.params, inits, sim=DEB1, samp.p,
                   Tmax, cnt=10, burnin=0.1,
                   plot=TRUE, sizestep=0.01, w.t=1, which=2,
                   test=TRUE, my.par=c(1,1),  myswitch=NULL,
                   mymap=NULL, ...)
{

  ## first we calculate a few lengths, etc, that we use for the for
  ## loops later.

  l<-length(samp.p)
  np<-length(all.params)

  ltot<-0
  for(i in seq_len(l)) ltot<-ltot+length(samp.p[[i]]$params)

  ## for testing, here is code that calculates (and prints out) the
  ## posterior prob of the "real" parameters, which can be passed in
  ## through all.params

  if(test){
    sim.old<-make.states(sim, all.params, inits, Tmax, which=which, sizestep, w.t,
                         myswitch=myswitch, mymap=mymap, ...)
    prob.old<-log.post.params(all.params, data, samp.p, sim.old)
    print(paste("(unnormalized) posterior probability of the real parameters= ",
                prob.old, sep=""))
  }

  ## build a data frame to hold the posterior samples. These include
  ## "samples" even for the parameters that are being held fixed, to
  ## make the code more straightforward and generic, as well as being
  ## useful for debugging. The samps structure will also keep track of
  ## the log posterior probability of that particular sample.

  samps<-data.frame(matrix(0, ncol=np+1, nrow=N))
  names(samps)<-c(names(all.params), "lpost")

  ## initialize the samps structure and p
  samps[1,seq_len(np)]<-p<-all.params

  ## initialize the data frame with the starting values passed into
  ## the function in the structure samp.p
  for(i in seq_len(l)){
    samps[1, samp.p[[i]]$params]<-samp.p[[i]]$start
    p[samp.p[[i]]$params]<-samp.p[[i]]$start
  }


  ## run the data simulation to make the underlying states (for
  ## determining the likelihoods) using the parameters stored in p.

  sim.old<-make.states(sim, p, inits, Tmax, which=which, sizestep, w.t,
                       myswitch=myswitch, mymap=mymap, ...)

  ## check the posterior probability to make sure you have reasonable
  ## starting values, and to initialize prob.old
  samps$lpost[1]<-log.post.params(samps[1,], data, samp.p, sim.old)

  print(paste("initial posterior probability = ", samps$lpost[1], sep=""))

  if(!is.finite(samps$lpost[1])) stop("bad starting values")

  ## now we begin the MCMC

  out<-list(s=samps[1,], p=p, sim.old=sim.old)

  for(i in 2:N){
    ## printing and plotting output so we can watch the progress
    if(i%%cnt == 0){
      print(paste("sample number", i, sep=" "))
      if(plot) plot.output(i, samps, samp.p, l, ltot, my.par)
    }

    ## the meat of the MCMC is found in the function update.samps (see below)

    out<-update.sample(samps[i-1,], samp.p, data, sim, inits, out,
                       Tmax, sizestep, w.t, l, which, i, cnt,
                       myswitch=myswitch, mymap=mymap)
    samps[i,]<-out$s
    if(test){
      if(-samps$lpost[i-1]+samps$lpost[i-1]<=-10){
        stop("we've had a really large swing in the posterior prob")
      }
    }

  }

  ##plot(b[100:N],type="l")
  lim<-min(1, round(burnin*N))
  samps <- samps[lim:N,]

  return(list(samps=samps))

}



##' @title update_sample
##' @param samps
##' @param samp.p
##' @param data
##' @param sim
##' @param inits
##' @param out
##' @param Tmax
##' @param sizestep
##' @param w.t
##' @param l
##' @param which
##' @param i
##' @param cnt
##' @param myswitch
##' @param mymap
##' @param test
##' @return
update_sample<-function(samps, samp.p, data, sim, inits, out, Tmax, sizestep,
                        w.t, l, which, i, cnt,  myswitch=NULL, mymap=NULL, test=TRUE, ...)
{
  ## read in some bits
  s<-samps
  sim.old<-out$sim.old
  p<-out$p

  #randomize updating order
  x<-seq_len(l)
  s.x<-sample(x)

  for(k in s.x){

    s.new<-s
    p.new<-p

    q<-propose_params(s, samp.p[[k]])

    ## automatically reject if the proposed parameter value is
    ## outside of the reasonable limits (i.e. < 0 )
    zeros<-0
    zeros<-check.zeros(samp.p[[k]], q$b)
    if(zeros){
      if(length(samp.p[[k]]$params)>=2){
        print("proposed one of the joint params outside of the range. Moving on.")
      }
      else print(paste("proposed ", samp.p[[k]]$params, " outside of the range. moving on", sep=""))
      next
    }
    ## write the proposed params into the p.new and s.new.

    for(j in seq_along(samp.p[[k]]$params)){
      ww<-samp.p[[k]]$params[j]
      p.new[ww]<-s.new[ww]<-q$b[j]
    }

    ## simulate the dynamics forward with the new parameters
    sim.new<-make.states(sim, p.new, inits, Tmax, which=which, sizestep, w.t,
                         myswitch=myswitch, mymap=mymap, ...)

    ## The posteriorprob of the previous sample is saved as
    ## s$lpost. If we accept a draw, we will set s$lpost<-s.new$lpost

    s.new$lpost<-log.post.params(s.new, data, samp.p, sim.new)

    if(is.finite(s.new$lpost) && is.finite(s$lpost)){
      A<-exp( s.new$lpost + q$lbak - s$lpost - q$lfwd )
    }
    else{
      A<-0
      print("whoops! must have proposed outside the correct region")
    }

    ## print some output so we can follow the progress
    if(i%%cnt==0){
      print(paste("proposing " , samp.p[[k]]$params, ": prob.old = ",
                  signif(s$lpost, digits=5),
                  "; prob.new = ", signif(s.new$lpost, digits=5), "; A = ",
                  signif(A, digits=5),
                  sep=""))
    }

    ## take a draw from a unif distrib, and use it to accept/reject
    u<-runif(1)
    if( u < A ){ ## accept
      sim.old<-sim.new
      p<-p.new
      s<-s.new
    }
  }

  return(list(s=s, p=p, sim.old=sim.old))

}

##' utility function to check for negative proposals when parameter is not allowed to be negative (e.g. because it's a variance)
##'
##'
##' @title check_zeros
##' @param s.p
##' @param q.b
##' @return
check_zeros<-function(s.p, q.b){
  z<-0
  for(j in seq_along(s.p$params)){
    if(s.p$params[j]=="l.M.HP") next
    if(q.b[j]<0) z<-1
    if(s.p$params[j]=="fs"){
      if(q.b[j]>1) z<-1
    }
  }
  return(z)
}



##' single proposal function
##'
##' individual paramater proposal
##' @title propose_params
##' @param samps
##' @param s.p
##' @return
propose_params<-function(samps, s.p)
{
  if(length(s.p$params)==1){
    ##print(paste(s.p$params, " proposing single ", sep=" "))
    q<-propose_single(samps, s.p)
  }
  else{
    ##print(paste(s.p$params, " proposing jointly ", sep=" "))
    q<-propose_joint(samps, s.p)
  }
  return(q)

}


##' propose a parameter individually
##'
##' @title propose_single
##' @param samps
##' @param s.p
##' @return
propose_single<-function(samps, s.p)##, i, freq=50, size=50 )##l=5, h=6)
{ ## I'm feeding in the variance, so I need to take the square root....

  b<-as.numeric(samps[s.p$params])
  var<-s.p$var
  type<-s.p$type
  hyps<-s.p$hyp

  if(type=="rw"){
    if(length(var)>1){
      l<-var[1]
      h<-var[2]
      b.new<-runif(1, l/h*b, h/l*b)
      lfwd<-dunif(b.new, l/h*b, h/l*b, log=TRUE)
      lbak<-dunif(b, l/h*b.new, h/l*b.new, log=TRUE)
    }
    else{
      sd<-sqrt(var)
      b.new<-rnorm(1, b, sd=sd)
      lfwd<-dnorm(b.new, b, sd=sd, log=TRUE)
      lbak<-dnorm(b, b.new, sd=sd, log=TRUE)
    }
    return(list(b=b.new, lfwd=lfwd, lbak=lbak))
  }
  else if(type=="ind"){
    out<-prior_draw(b, hyps, s.p$params)
    return(out)
  }

}


##' joint proposal function
##'
##' Function to jointly propose parameters using a multivariate normal proposal distribution
##' @title propose_joint
##' @param samp
##' @param samp.p
##' @return
propose_joint<-function(samp, samp.p){


  b<-NULL
  if(samp.p$type=="rw"){
    b<-as.numeric(samp[samp.p$params])
    sigma<-samp.p$var

    b.new<-rmvnorm(1, mean=b, sigma=sigma)
    lfwd<-dmvnorm(b.new, b, sigma, log=TRUE)
    lbak<-dmvnorm(b, b.new, sigma, log=TRUE)
  }
  else if(samp.p$type=="ind"){
    if(is.null(samp.p$mean)) stop("not enough info for the independence sampler")
    mean<-as.numeric(samp.p$mean)
    b<-as.numeric(samp[samp.p$params])
    sigma<-samp.p$var

    b.new<-rmvnorm(1, mean=mean, sigma=sigma)
    lfwd<-dmvnorm(b.new, mean, sigma, log=TRUE)
    lbak<-dmvnorm(b, mean, sigma, log=TRUE)
  }

  ##print(c(b, b.new, lfwd, lbak))


  ##samp[s]<-b.new

  ##stop()
  return(list(b=b.new, lfwd=lfwd, lbak=lbak))

}

##' a plotting function
##'
##' @title plot_output
##' @param i
##' @param samps
##' @param samp.p
##' @param l
##' @param ltot
##' @param my.par
##' @param plot.post
##' @return
plot_output<-function(i, samps, samp.p, l, ltot, my.par=c(2,4), plot.post=TRUE){

  if( ltot > 1 ) par(mfrow=my.par, bty="n")
  for( j in seq_len(l)){
    ww<-samp.p[[j]]$params
    for(k in seq_along(ww)){
      plot(samps[seq_len(i-1),ww[k]], type="l", xlab="sample", ylab=ww[k])
    }
  }
  if(plot.post){
    my.min<-20
    if(i>my.min){
      x<-seq(my.min, (i-1), by=1)
      plot(x, samps$lpost[x], type="l", xlab="sample", ylab="log posterior prob")
    }
  }
}


##' draw from prior
##'
##' details of this function
##' @title prior_draw
##' @param b
##' @param hyp
##' @param p
##' @return
prior_draw<-function(b, hyp, p){

  param1<-hyp[1]
  param2<-hyp[2]

  gams<-c("sd.Z", "alpha", "sr",  "ds", "muz", "eta")
  norms<-c("Tmin", "Tmax")
  betas<-c("fs")

  if( p %in% gams ){
    b.new<-rgamma(1, shape=param1, rate=param2)
    lfwd<-dgamma(b.new, shape=param1, rate=param2, log=TRUE)
    lbak<-dgamma(b, shape=param1, rate=param2, log=TRUE)
  }
  else if( p %in% norms ){
    b.new<-rnorm(1, mean=param1, sd=param2)
    lfwd<-dnorm(b.new, mean=param1, sd=param2, log=TRUE)
    lbak<-dnorm(b, mean=param1, sd=param2, log=TRUE)
  }
  else if( p %in% betas ){
    b.new<-rbeta(1,  shape1=param1, shape2=param2)
    lfwd<-dbeta(b.new,  shape1=param1, shape2=param2, log=TRUE)
    lbak<-dbeta(b, shape1=param1, shape2=param2, log=TRUE)
  }

  return(list(b=b.new, lfwd=lfwd, lbak=lbak))
}



##' Bayesian inference for a deterministic differential equation model
##'  (with models solved via an ode or dde solver in the function specified with argument
##' "de.model") with an observation model, and with observation error
##' variances specified in sds.
##'
##'
##' @title de_mcmc old
##'
##'  this is just a wrapper for the old function, reassigning inputs from the debinfer_parlist object to dde_mcmc
##'
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
##'
##' @export
##'
de_mcmc_old <- function(N, data, de.model, obs.model, all.params,
                        Tmax, data.times, cnt=10, burnin=0.1,
                        plot=TRUE, sizestep=0.01, which=1,
                        myswitch=NULL,
                        mymap=NULL, verbose =FALSE, method = "lsoda", ...)

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
                           data.times = data.times, obs.model=obs.model, verbose = verbose, method = method, ...)
  return(mcmc_samples)

}


