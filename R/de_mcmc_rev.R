### revised version of the inference function


##' Bayesian inference for a deterministic DE model (with models
##' solved via an DE solver) with an observation model.
##'
##'
##' @title de_mcmc
##' @import coda
##' @import PBSddesolve
##' @param N integer, number of MCMC iterations
##' @param data data.frame of time course observations to fit the model to. The observations must be ordered ascending by time.
##' @param all.params debinfer_parlist containing all model, MCMC, and observation
##' @param ref.params an optional named vector containing a set of reference parameters, e.g. the true paramaters underlying a simulated data set
##' @param ref.inits  an optional named vector containing a set of reference initial values, e.g. the true initial values underlying a simulated data set
##' @param de.model a function defining a DE model, compliant with the solvers in deSolve or PBSddesolve
##' @param obs.model a function defining an observation model. Must be a function with arguments 'data', 'sim.data', 'samp'.
##' @param Tmax maximum timestep for solver
##' @param cnt integer interval at which to print and possibly plot information on the current state of the MCMC chain
##' @param plot logical, plot traces for all parameters at the interval defined by \code{cnt}
##' @param sizestep timestep for solver to return values at, only used if data.times is missing
##' @param solver the solver to use. 1 or "ode" = deSolve::ode; 2 or "dde" = PBSddesolve::dde; 3 or "dede" = deSolve::dde
##' @param data.times time points for which observations are available
##' @param verbose.mcmc logical display MCMC progress messages
##' @param verbose logical display verbose solver output
##' @param ... further arguments to the solver
##' @return a debinfer_result object containing input parameters, data and MCMC samples
##' @importFrom methods formalArgs
##' @export
de_mcmc <- function(N, data, de.model, obs.model, all.params, ref.params=NULL, ref.inits=NULL,
                              Tmax, data.times, cnt=10,
                              plot=TRUE, sizestep=0.01, solver="ode",
                              verbose.mcmc =TRUE, verbose = FALSE, ...)
{
  #check models
  if(!is.function(obs.model)) stop("obs.model must be a function")
  if(!identical(formalArgs(obs.model), c("data", "sim.data", "samp" ))) stop("obs.model must be a function with arguments 'data', 'sim.data', 'samp'")
  if(!(is.function(de.model) || is.character(de.model))) stop("de.model must be a function or function name (character)")
  if(!inherits(all.params, "debinfer_parlist")) stop("all.params must be of class debinfer_parlist")

  #identify cov matrices
  is.par <- vapply(all.params, class, character(1)) == "debinfer_par"
  is.cov <- vapply(all.params, class, character(1)) == "debinfer_cov"
  #subset into parameters and covariance matrices
  if (any(is.cov)){
    cov.matrices <- all.params[is.cov]
  } else {
    cov.matrices <- NULL
  }
  all.params <- all.params[is.par]
  #identify free parameters and inits
  p.names <- vapply(all.params, function(x) x$name, character(1))
  is.free <- !vapply(all.params, function(x) x$fixed, logical(1))
  is.init <- vapply(all.params, function(x) x$var.type, character(1))=="init"
  is.de <- vapply(all.params, function(x) x$var.type, character(1))=="de"
  #identify free parameters with joint proposal and split into blocks
  is.single <- vapply(all.params, function(x) is.null(x$joint) , logical(1)) & is.free
  joint.blocks <- unique(vapply(all.params[!is.single & is.free], function(x) x$joint, character(1)))

  #get all start values
  p.start <- lapply(all.params, function(x) x$value)[is.free]
  names(p.start) = p.names[is.free]
  # get the parameters that are to be estimated
  w.p <- p.names[is.free]
  #this is needed for the reference likelihood calculation
  params <- unlist(lapply(all.params, function(x) x$value))
  names(params) <-  p.names
  #initial values for DE (no re-ordering!)
  inits <- vapply(all.params, function(x) x$value, numeric(1))[is.init]
  names(inits) <- p.names[is.init]
  #inits are matched by order in deSolve. inform user of input order
  message(paste("Order of initial conditions is ", paste(names(inits), collapse = ", ")))

  hyper = lapply(all.params, function(x) x$hyper)[is.free]
  names(hyper) <- p.names[is.free]
  pdfs = lapply(all.params, function(x) x$prior)[is.free]
  names(pdfs) = p.names[is.free]

  prop.sd <- sapply(all.params[is.free], function(x) x$prop.var)
  names(prop.sd) <- p.names[is.free]


  ## first we calculate a few lengths, etc, that we use for the for
  ## loops later.

  n.free <- sum(is.free)
  n.joints <- length(joint.blocks)
  np<-length(all.params)


  #This only seems to be used in the trace plotting subroutine
  #ltot<-0
  #for(i in 1:l) ltot<-ltot+length(samp.p[[i]]$params)

  ## for testing, here is code that calculates (and prints out) the
  ## posterior prob of the reference parameters, which can be passed in
  ## through ref.params

  if(!is.null(ref.params) && !is.null(ref.inits)){
    sim.ref<-solve_de(sim = de.model, params = ref.params, inits = ref.inits, Tmax = Tmax, solver=solver, sizestep = sizestep, data.times = data.times, ...)
    prob.ref<-log_post_params(samp = ref.params, data = data, sim.data = sim.ref, w.p = w.p, obs.model = obs.model, pdfs = pdfs, hyper = hyper)
    message(paste("(unnormalized) posterior probability of the reference parameters= ",
                prob.ref, sep=""))
  }

  ## build an mcmc object to hold the posterior samples. These include
  ## "samples" even for the parameters that are being held fixed, to
  ## make the code more straightforward and generic, as well as being
  ## useful for debugging. The samps structure will also keep track of
  ## the log posterior probability of that particular sample.

  samps <- coda::mcmc(matrix(numeric((np+1)*N), ncol=np+1, nrow=N, dimnames = list(NULL, c(p.names, "lpost"))))

  #initialize the mcmc structure with the starting values
  samps[1,seq_len(np)] <- params

  ## run the data simulation to make the underlying states (for
  ## determining the likelihoods) using the parameters stored in p.

  sim.start<-solve_de(sim = de.model, params = params[is.de], inits = inits, Tmax = Tmax, solver=solver, sizestep = sizestep, data.times = data.times, ...)

  ## check that solver provides simulation values for all observations
  if(inherits(sim.start, "try-error")) {
    stop("solver failed on start values")
    } else {
       if(!all(data.times %in% sim.start[,"time"])) stop("solver times do not cover all data times")
    }


  ## check the posterior probability to make sure you have reasonable
  ## starting values, and to initialize prob.old
  try.lpost <- try(log_post_params(samp = params, data = data, sim.data = sim.start, obs.model = obs.model, pdfs = pdfs, hyper = hyper, w.p = w.p))
  if(inherits(try.lpost, "try-error")) {
    stop("log likelihood calculation failed with current starting values. Observation model may try to use undeclared parameters or model outputs.")
  } else {
    samps[1,"lpost"] <- try.lpost
  }
  if(!is.finite(samps[1,"lpost"])) stop("Infinite log likelihood with current starting values")

  message(paste("initial posterior probability = ", samps[1,"lpost"], sep=""))



  ## now we begin the MCMC

  out <- list(s=samps[1,], p=params, inits=inits, sim.old=sim.start)

  for(i in 2:N){

    ## the meat of the MCMC is found in the function update.samps (see below)
    #if(verbose.mcmc) message(i)
    out <- update_sample_rev(samps = samps[i-1,], samp.p = all.params[is.free], cov.mats = cov.matrices, data = data, sim = de.model, out = out,
                       Tmax = Tmax, sizestep = sizestep, data.times = data.times, l=n.free + n.joints, solver=solver, i=i, cnt=cnt, w.p = w.p,
                       obs.model = obs.model, pdfs = pdfs, hyper = hyper, verbose.mcmc = verbose.mcmc, verbose = verbose, is.de=is.de,
                       is.single = is.single, joint.blocks = joint.blocks, ...)
    samps[i,] <- out$s #make sure order is matched
#     if(test){
#       if(-samps$lpost[i-1]+samps$lpost[i-1]<=-10){
#         stop("we've had a really large swing in the posterior prob")
#       }
#     }

    ## printing and plotting output so we can watch the progress
    if(i%%cnt == 0){
      if (verbose.mcmc) message(paste("sample number", i, sep=" "))
      if (plot)         plot(window(samps[, is.free], 1, i), density=FALSE, ask=FALSE)#use window.mcmc
    }

  }

  #make a data structure that returns all parts of the analysis
  result <- list(
    de.model = de.model,
    obs.model = obs.model,
    all.params = all.params,
    ref.params = ref.params,
    ref.inits = ref.inits,
    iter = N,
    data = data,
    samples = coda::as.mcmc(samps[,w.p]),
    solver = solver,
    lpost = samps[,"lpost"]
  )
  result <- structure(result, class="debinfer_result")
  return(result)
}

##' update_sample_rev
##'
##' This is the workhorse of the MCMC algorithm
##'
##' @param samps row vector of samples from the previous mcmc iteration
##' @param samp.p the parlist created by setup_debinfer
##' @param cov.mats the covariance matrices
##' @param data the observation
##' @param sim the de.model
##' @param out list containing the initial or previous update i.e. list(s=samps[i-1,], inits=inits, p=params, sim.old=sim.start)
##' @param Tmax maximum timestep for solver
##' @param sizestep sizestep for solver when not using data.times
##' @param l number of parameters to be proposed
##' @param solver solver choice
##' @param i current MCMC iteration
##' @param cnt interval for printing/plotting information on chains
##' @param data.times times with observations
##' @param obs.model function containing obs model
##' @param pdfs names of prior pdfs
##' @param hyper list of hyperparameters
##' @param w.p names of free parameters
##' @param verbose.mcmc logical print MCMC progress messages
##' @param verbose logical, print additional information from solver
##' @param is.de logical, parameter is an input for the solver
##' @param is.single parameter is to be proposed individually
##' @param joint.blocks names of joint blocks
##' @param ... further arguments to solver
##' @export
update_sample_rev<-function(samps, samp.p, cov.mats, data, sim, out, Tmax, sizestep,
                        data.times, l, solver, i, cnt, obs.model, pdfs, hyper, w.p, verbose.mcmc, verbose, is.de,
                        is.single, joint.blocks, ...)
{
  ## read in some bits
  s<-samps
  sim.old<-out$sim.old
  p<-out$p
  ints <- out$inits

  singles <- names(is.single)[is.single]
  singles.blocks <- c(singles,joint.blocks)

  #randomize updating order
  #x<-seq_along(singles) # this is the sum of not joint parameters and joint blocks
  #then get the number of joint blocks and add the joint blocks
  s.x<-sample(singles.blocks) #resample order

  for(k in s.x){

    s.new<-s #sample vector
    p.new<-p #parameter vector
    i.new<-ints #initial value vector

    #pick

    if(k %in% singles){#if k is a single par
      ##print(paste(s.p$params, " proposing single ", sep=" "))
      q<-propose_single_rev(samps = s, s.p = samp.p[[k]])
    }
    else {
     if (k %in% joint.blocks){
       ##print(paste(s.p$params, " proposing jointly ", sep=" "))
       q<-propose_joint_rev(samps = s, s.ps = samp.p, cov.mat = cov.mats[[k]])
     } else {stop("Parameter is neither in a joint block nor set up for individual proposals. This should not happen.")}

    }

    ## automatically reject if the proposed parameter value is
    ## outside of the prior support
    if (k %in% singles){
      qprior <- logd_prior(q$b, pdfs[[k]], hypers=hyper[[k]])
    } else {
      if (k %in% joint.blocks){
        qprior <- sapply(dimnames(cov.mats[[k]]$sigma)[[1]], function(x) logd_prior(q$b[x], pdfs[[x]], hypers=hyper[[x]]), USE.NAMES = FALSE)##TIDY UP
        #if(verbose.mcmc) message(paste("assessing logd_prior for join proposal", dimnames(cov.mats[[k]]$sigma)[[1]], q$b, "; qprior = ", qprior))
      }
    }

    if (all(is.finite(qprior))){
      ## write the proposed params into the p.new and s.new.

      #for(j in 1:length(samp.p[[k]]$name)){#this will need to be able to handle joint proposals
        #ww<-samp.p[[k]]$name
      if (k %in% singles) jj <- k else jj <- dimnames(cov.mats[[k]]$sigma)[[1]] ##TIDY UP
      for (j in jj){
        if (samp.p[[j]]$var.type== "de" || samp.p[[j]]$var.type == "obs")  p.new[j]<-s.new[j]<-q$b[j]
        if (samp.p[[j]]$var.type== "init") i.new[j]<-s.new[j]<-q$b[j]
      }
      #}

      ## simulate the dynamics forward with the new parameters, but only if parameter in question is not an observation parameter
      if (all(vapply(jj, function(x) samp.p[[x]]$var.type, character(1)) == "obs")){ ##TIDY UP
        sim.new <- sim.old
        #if(verbose.mcmc)message(paste("keeping simulation",jj))#keep using last available de solution
      } else { #compute new solution
       sim.new<-solve_de(sim = sim , params = p.new[is.de], inits = i.new, Tmax = Tmax, solver=solver, sizestep = sizestep, data.times = data.times, ...)
       #if(verbose.mcmc)message(paste("renewing simulation",jj))#keep using last available de solution
      }



      ## The posteriorprob of the previous sample is saved as
      ## s$lpost. If we accept a draw, we will set s$lpost<-s.new$lpost

      ##calculate posterior likelihood, but only if solver did not fail
        if (!inherits(sim.new, "try-error")){
          s.new["lpost"] <- log_post_params(samp = s.new, data = data, sim.data = sim.new, obs.model = obs.model, pdfs = pdfs, hyper = hyper, w.p = w.p)
        } else {
          if (verbose) message(paste("Solver failed with current state =", s, "and proposal", paste(samp.p[[k]]$name,"=",q$b , collapse = ", ")))
          s.new["lpost"] <- -Inf
          }

      if(is.finite(s.new["lpost"]) && is.finite(s["lpost"])){
        A<-exp( s.new["lpost"] + q$lbak - s["lpost"] - q$lfwd )
      } else {
        A<-0
        if (verbose.mcmc) message(paste("posterior not finite for proposal", paste(samp.p[[k]]$name,"=",q$b , collapse = ", ")))
      }
    } else {
      A<-0
      if (verbose.mcmc) message(paste("proposal outside prior support for", paste(samp.p[[k]]$name,"=",q$b , collapse = ", ")))
    }

      ## print some output so we can follow the progress
      if(verbose.mcmc && i%%cnt==0  ){
        message(paste("proposing " , samp.p[[k]]$name, ": prob.old = ",
                    signif(s["lpost"], digits=5),
                    "; prob.new = ", signif(s.new["lpost"], digits=5), "; A = ",
                    signif(A, digits=5),
                    sep=""))
      }

      ## take a draw from a unif distrib, and use it to accept/reject
      u<-runif(1)
      if( u < A ){ ## accept
        sim.old<-sim.new
        p<-p.new
        s<-s.new
        ints <-i.new
      }
  }

  return(list(s=s, p=p, inits=ints, sim.old=sim.old))

}





##' propose a parameter individually
##'
##' @title propose_single_rev
##' @param samps current sample of the MCMC chain
##' @param s.p debinfer_par object representing the parameter that is to be proposed
##' @import stats
propose_single_rev<-function(samps, s.p)
{ ## I'm feeding in the variance, so I need to take the square root....

  b<-samps[s.p$name]
  var<-s.p$prop.var
  type<-s.p$samp.type
  hyps<-s.p$hypers

  if(type=="rw-unif"){
      l<-var[1]
      h<-var[2]
      b.new<-runif(1, l/h*b, h/l*b)
      names(b.new) <- names(b)
      lfwd<-dunif(b.new, l/h*b, h/l*b, log=TRUE)
      lbak<-dunif(b, l/h*b.new, h/l*b.new, log=TRUE)
      return(list(b=b.new, lfwd=lfwd, lbak=lbak))
    }
  if(type=="rw"){
      sd<-sqrt(var)
      b.new<-rnorm(1, b, sd=sd)
      names(b.new) <- names(b)
      lfwd<-dnorm(b.new, b, sd=sd, log=TRUE)
      lbak<-dnorm(b, b.new, sd=sd, log=TRUE)
      return(list(b=b.new, lfwd=lfwd, lbak=lbak))
    }
 if(type=="ind"){
    out<-prior_draw_rev(b, hyps, s.p$prior)
    return(out)
    }
 if(type=="rw-ref"){
   l.bound <- s.p$bounds[1]
   u.bound <- s.p$bounds[2]
   sd<-sqrt(var)
   b.new <- rnorm(1, b, sd=sd)
   names(b.new) <- names(b)
   while(b.new > u.bound || b.new < l.bound){
     if(b.new > u.bound) b.new <- 2*u.bound - b.new; #print(b.new)
     if(b.new < l.bound) b.new <- 2*l.bound - b.new; #print(b.new)
   }
   lfwd<-dnorm(b.new, b, sd=sd, log=TRUE)
   lbak<-dnorm(b, b.new, sd=sd, log=TRUE)
   return(list(b=b.new, lfwd=lfwd, lbak=lbak))
 }
  #for generalised reflecting sampler: calculate support of prior using qprior
  #ideally do this in setup_mcmc, to avoid doing the calculation for every proposal
  #use support as bounds and then use reflction as described in Yang and Rodriguez 2013
  #i.e. x < a --> x = 2a-x; x > b --> x = 2b-x

}

##' joint proposal function
##'
##' Function to jointly propose parameters using a multivariate normal proposal distribution
##' @title propose_joint
##' @param samps current sample of the MCMC chain
##' @param s.ps debinfer_parlist object representing the parameters that are to be proposed
##' @param cov.mat debinfer_cov object; covariance matrix for the proposal
##' @import stats
##' @import mvtnorm
propose_joint_rev<-function(samps, s.ps, cov.mat){


  b<-NULL
  joint.pars <- dimnames(cov.mat$sigma)[[1]]
  if(cov.mat$samp.type == "rw"){
    b<-samps[joint.pars]

    b.new<-rmvnorm(1, mean=b, sigma=cov.mat$sigma, method="svd")
    lfwd<-dmvnorm(b.new, b, cov.mat$sigma, log=TRUE)
    lbak<-dmvnorm(b, b.new, cov.mat$sigma, log=TRUE)
  }
  else if(cov.mat$samp.type =="ind"){
    stop("multivariate independence sampler not yet implemented")
    # if(is.null(s.p$mean)) stop("not enough info for the independence sampler")
    # mean<-as.numeric(s.p$mean)
    # b<-as.numeric(samps[s.p$params])
    # sigma<-s.p$var
    #
    # b.new<-rmvnorm(1, mean=mean, sigma=sigma)
    # lfwd<-dmvnorm(b.new, mean, sigma, log=TRUE)
    # lbak<-dmvnorm(b, mean, sigma, log=TRUE)
  } else stop(paste("unknown sampler type for multivariate block", k))

  ##print(c(b, b.new, lfwd, lbak))


  ##samp[s]<-b.new

  ##stop()
  return(list(b=b.new[1,], lfwd=lfwd, lbak=lbak))

}

##' draw from prior
##'
##' @title prior_draw_rev
##' @param b current value of a parameter
##' @param hypers list of hyper parameters, named appropriately for the corresponding prior.pdf
##' @param prior.pdf string name of probability distribution following base R conventions, or those of additionally loaded packages
##' @import stats
prior_draw_rev<-function(b, hypers, prior.pdf){
    #assemble random number generator function name
    rand <- paste("r", prior.pdf, sep="")
    #assemble arguments
    rand.args <- c(n=1, hypers)
    #draw from prior
    b.new<- do.call(rand, rand.args)
    #assemble density function name
    dens <- paste("d", prior.pdf, sep="")
    lfwd<-do.call(dens, c(x=b.new,  hypers, log=TRUE))
    lbak<-do.call(dens, c(x=unname(b),  hypers, log=TRUE))

  return(list(b=b.new, lfwd=lfwd, lbak=lbak))
}

#' log_post_params
#'
#' evaluate the likelihood given the data, the current deterministic model solution and the observation model
#' @param samp named numeric; current sample
#' @param w.p character; parameter names
#' @param data data
#' @param pdfs character, prior pdf names
#' @param hyper list, hyper parameters for the priors
#' @param sim.data solver output
#' @param obs.model function containing the observation model
#'
#' @export
log_post_params <- function(samp, w.p, data, pdfs, hyper, sim.data, obs.model){

  log_data <- obs.model

  llik <- log_data(data=data, sim.data=sim.data, samp=samp)

  #if(length(w.p)==1) lprior<-as.numeric(log_prior_params(samp, w.p, hyper))
  #else {
  lprior<-sum(log_prior_params(samp, pdfs, w.p, hyper))
  ##if(!is.finite(lprior)) break
  #}
  ##print(c(b, lik, prior))

  if(is.na(llik)) break #break doesn't work when it's inside a function
  if(is.na(lprior)) break

  return( llik + lprior )
}


#' log_prior_params
#'
#' evaluate prior density at current parameter values
#' @param samp named numeric; current sample
#' @param w.p character; parameter names
#' @param pdfs character, prior pdf names
#' @param hyper list of named hyper parameters for the priors
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
  for(i in seq_len(len)){
    p<-w.p[i]
    s<-as.numeric(samp[p])

    lp[[p]] <- logd_prior(s, pdfs[[p]], hypers=hyper[[p]])

  }

  return(lp)

}
