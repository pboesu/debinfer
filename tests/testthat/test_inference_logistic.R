library(deBInfer)
context("Testing an entire inference procedure using a simple logistic model")

test_that("Inference on simulated data w/ unknown obs. error returns simulation parameters. DEB version.", {
  #testthat::skip_on_cran()
  # define logistic model
  logistic_model <- function(time, y, parms) {
    with(as.list(c(y, parms)), {
      dN <- r * N * (1 - N / K)
      list(dN)
    })
  }
  expect_true(exists("logistic_model"))
  # set initial value for simulation
  y <- c(N = 0.1)
  # set parameter values
  parms <- c(r = 0.1, K = 10)
  # set simulation time points
  times <- seq(0, 120, 1)
  # solve ODE
  out <- ode(y, times, logistic_model, parms, method='lsoda')
  # sample from simulated data
  set.seed(143)
  N_obs <- as.data.frame(out[c(1,runif(35, 0, nrow(out))),]) #force include the first time-point (t=0)
  # add lognormal noise
  parms['loglogsd.N'] <- -4.6
  N_obs$N_noisy <- rlnorm(nrow(N_obs), log(N_obs$N),exp(parms['loglogsd.N']))
  # observations must be ordered for solver to work
  N_obs <- N_obs[order(N_obs$time),]

  # define an observation model
  # NB: lognormal errors are not great really for the obs model - should be changed to sth that actually allows N to be zero instead of using epsilon correction
  logistic_obs_model<-function(data, sim.data, sds, samp){

    llik.N<-sum(dlnorm(data$N_noisy, meanlog=log(sim.data$N+1e-6), sdlog=exp(samp[['loglogsd.N']]), log=TRUE))

    llik<-llik.N

    return(llik)
  }


  # do inference with deBInfer
  # MCMC iterations
  iter = 5000
  # define burnin
  burnin = 2000
  # inference call
  mcmc_samples <- deb_mcmc(N=iter, p.start=list(r=0.5, K=5, loglogsd.N=-2), data=N_obs, w.p=c("r", "K", "loglogsd.N"), params=parms,
                           inits=c(N=0.1), sim=logistic_model, sds=list(N=0.01), hyper=list(r=list(mean=0, sd=1),
                           K=list(meanlog=1, sdlog=1), loglogsd.N=list(mean=-2, sd=1)), pdfs = list(r='norm', K='lnorm', loglogsd.N='norm'),
                           prop.sd=c(r=0.005, K=0.1, loglogsd.N=0.5),
                           Tmax=max(N_obs$time), cnt=iter, burnin=200, plot=FALSE, sizestep=0.1, which = 1,
                           data.times = N_obs$time, obs.model=logistic_obs_model)
  #add more tests here checking the integrity & contents of the returned data structure
  #check accuracy of estimation (threshold is 0.5% of true parameter value)
  expect_equal(unname(colMeans(mcmc_samples$samps[burnin:iter,])/parms),c(1,1,1),tolerance = 1e-2)
})

test_that("Inference on simulated data w/ known obs. error returns simulation parameters. de_mcmc wrapper version.", {
  #testthat::skip_on_cran()
  # define logistic model
  logistic_model <- function(time, y, parms) {
    with(as.list(c(y, parms)), {
      dN <- r * N * (1 - N / K)
      list(dN)
    })
  }
  expect_true(exists("logistic_model"))
  # set initial value for simulation
  y <- c(N = 0.1)
  # set parameter values
  parms <- c(r = 0.1, K = 10)
  # set simulation time points
  times <- seq(0, 120, 1)
  # solve ODE
  out <- ode(y, times, logistic_model, parms, method='lsoda')
  # sample from simulated data
  set.seed(143)
  N_obs <- as.data.frame(out[c(1,runif(35, 0, nrow(out))),]) #force include the first time-point (t=0)
  # add lognormal noise
  parms['loglogsd.N'] <- -4.6
  N_obs$N_noisy <- rlnorm(nrow(N_obs), log(N_obs$N),exp(parms['loglogsd.N']))
  # observations must be ordered for solver to work
  N_obs <- N_obs[order(N_obs$time),]

  # define an observation model
  # NB: lognormal errors are not great really for the obs model - should be changed to sth that actually allows N to be zero instead of using epsilon correction
  logistic_obs_model<-function(data, sim.data, sds, samp){

    llik.N<-sum(dlnorm(data$N_noisy, meanlog=log(sim.data$N+1e-6), sdlog=exp(samp[['loglogsd.N']]), log=TRUE))

    llik<-llik.N

    return(llik)
  }

  library(deBInfer)
  r <- debinfer_par(name = "r", var.type = "de", fixed = FALSE,
                    value = 0.5, prior="norm", hypers=list(mean = 0, sd = 1),
                    prop.var=0.005, samp.type="rw")

  K <- debinfer_par(name = "K", var.type = "de", fixed = FALSE,
                    value = 5, prior="lnorm", hypers=list(meanlog = 1, sdlog = 1),
                    prop.var=0.1, samp.type="rw")

  loglogsd.N <- debinfer_par(name = "loglogsd.N", var.type = "obs", fixed = FALSE,
                             value = -10, prior="norm", hypers=list(mean = 0, sd = 1),
                             prop.var=0.5, samp.type="rw")

  #we also need to provide an initial condition for the DE
  N <- debinfer_par(name = "N", var.type = "init", fixed = TRUE,
                    value = 0.1)


  mcmc.pars <- setup_debinfer(r, K, loglogsd.N, N)

  # do inference with deBInfer
  # MCMC iterations
  iter = 5000
  # define burnin
  burnin = 2000
  # inference call

  mcmc_samples <- de_mcmc(N = iter, data=N_obs, de.model=logistic_model, obs.model=logistic_obs_model, all.params=mcmc.pars,
                   Tmax = max(N_obs$time), data.times=N_obs$time, cnt=iter+1, burnin=0.1,
                   plot=FALSE, sizestep=0.1, which=1)
  #add more tests here checking the integrity & contents of the returned data structure
  #check accuracy of estimation (threshold is 0.5% of true parameter value)
  expect_equal(unname(colMeans(mcmc_samples$samps[burnin:iter,])/parms),c(1,1,1),tolerance = 2e-2)
})



test_that("Inference on simulated data w/ known obs. error returns simulation parameters. DEB version.", {
  #testthat::skip_on_cran()
  # define logistic model
  logistic_model <- function(time, y, parms) {
    with(as.list(c(y, parms)), {
      dN <- r * N * (1 - N / K)
      list(dN)
    })
  }
  expect_true(exists("logistic_model"))
  # set initial value for simulation
  y <- c(N = 0.1)
  # set parameter values
  parms <- c(r = 0.1, K = 10)
  # set simulation time points
  times <- seq(0, 120, 1)
  # solve ODE
  out <- ode(y, times, logistic_model, parms, method='lsoda')
  # sample from simulated data
  set.seed(143)
  N_obs <- as.data.frame(out[c(1,runif(35, 0, nrow(out))),]) #force include the first time-point (t=0)
  # add lognormal noise
  N_obs$N_noisy <- rlnorm(nrow(N_obs), log(N_obs$N),0.01)
  # observations must be ordered for solver to work
  N_obs <- N_obs[order(N_obs$time),]

  # define an observation model
  # NB: lognormal errors are not great really for the obs model - should be changed to sth that actually allows N to be zero instead of using epsilon correction
  logistic_obs_model<-function(data, sim.data, sds, samp){

    llik.N<-sum(dlnorm(data$N_noisy, meanlog=log(sim.data$N + 1e-6), sdlog=sds$N, log=TRUE))

    llik<-llik.N

    return(llik)
  }


  # do inference with deBInfer
  # MCMC iterations
  iter = 5000
  # define burnin
  burnin = 2000
  # inference call
  mcmc_samples <- deb_mcmc(N=iter, p.start=list(r=0.5, K=5), data=N_obs, w.p=c("r", "K"), params=parms,
                           inits=c(N=0.1), sim=logistic_model, sds=list(N=0.01), hyper=list(r=list(mean=0, sd=1),
                                                                                            K=list(meanlog=1, sdlog=1)), pdfs = list(r='norm', K='lnorm'), prop.sd=c(r=0.001, K=0.1),
                           Tmax=max(N_obs$time), cnt=iter, burnin=200, plot=FALSE, sizestep=0.1, which = 1,
                           data.times = N_obs$time, obs.model=logistic_obs_model)
  #add more tests here checking the integrity & contents of the returned data structure
  #check accuracy of estimation (threshold is 0.5% of true parameter value)
  expect_equal(unname(colMeans(mcmc_samples$samps[burnin:iter,])/parms),c(1,1),tolerance = 5e-3)
})

