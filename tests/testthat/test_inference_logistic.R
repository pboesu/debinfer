library(deBInfer)
context("Testing an entire inference procedure using a simple logistic model")


# define logistic model
model <- function (time, y, parms) {
  with(as.list(c(y, parms)), {
    dN <- r * N * (1 - N / K)
    list(dN)
  })
}

test_that("Inference on simulated data returns simulation parameters.", {
  #testthat::skip_on_cran()

  expect_true(existsFunction("model"))
  # set initial value for simulation
  y <- c(N = 0.1)
  # set parameter values
  parms <- c(r = 0.1, K = 10)
  # set simulation time points
  times <- seq(0, 120, 1)
  # solve ODE
  out <- ode(y, times, model, parms, method='lsoda')
  # sample from simulated data
  set.seed(143)
  N_obs <- as.data.frame(out[c(1,runif(35, 0, nrow(out))),]) #force include the first time-point (t=0)
  # add lognormal noise
  N_obs$N_noisy <- rlnorm(nrow(N_obs), log(N_obs$N),0.01)
  # observations must be ordered for solver to work
  N_obs <- N_obs[order(N_obs$time),]

  # define an observation model
  # NB: lognormal errors are not great really for the obs model - should be changed to sth that actually allows N to be zero instead of using epsilon correction
  log.post.params<-function(samp, w.p, data, p, pdfs, hyper, sim.data, sds, verbose.lik=FALSE){

    llik.N<-sum(dlnorm(data$N_noisy, meanlog=log(sim.data$N+0.000000001), sdlog=sds$N, log=TRUE))

    llik<-llik.N

    if(length(w.p)==1) lprior<-as.numeric(log_prior_params(samp, w.p, hyper))
    else {
      lprior<-sum(log_prior_params(samp, pdfs, w.p, hyper))
      ##if(!is.finite(lprior)) break
    }
    ##print(c(b, lik, prior))

    if(is.na(llik)) break
    if(is.na(lprior)) break

    return( llik + lprior )
  }
  expect_true(existsFunction("log.post.params"))


  # do inference with deBInfer
  # MCMC iterations
  iter = 10000
  # define burnin
  burnin = 2000
  # inference call
  mcmc_samples <- deb_mcmc(N=iter, p.start=list(r=0.5, K=5), data=N_obs, w.p=c("r", "K"), params=parms,
                           inits=c(N=0.1), sim=model, sds=list(N=0.01), hyper=list(r=list(mean=0, sd=1),
                           K=list(meanlog=1, sdlog=1)), pdfs = list(r='norm', K='lnorm'), prop.sd=c(r=0.001, K=0.1),
                           Tmax=max(N_obs$time), cnt=iter+1, burnin=200, plot=FALSE, sizestep=0.1, which = 1,
                           data.times = N_obs$time)
  expect_equal(colMeans(mcmc_samples$samps[burnin:iter,]),parms,tolerance = 1e-2)
})
