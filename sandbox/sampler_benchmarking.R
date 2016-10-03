#benchmark bounded sampler

#testthat::skip_on_cran()
# define logistic model
logistic_model <- function(time, y, parms) {
  with(as.list(c(y, parms)), {
    dN <- r * N * (1 - N / K)
    list(dN)
  })
}

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
parms['logsd.N'] <- 0.01
N_obs$N_noisy <- rlnorm(nrow(N_obs), log(N_obs$N),(parms['logsd.N']))
# observations must be ordered for solver to work
N_obs <- N_obs[order(N_obs$time),]

# define an observation model
# NB: lognormal errors are not great really for the obs model - should be changed to sth that actually allows N to be zero instead of using epsilon correction
logistic_obs_model<-function(data, sim.data, samp){

  llik.N<-sum(dlnorm(data$N_noisy, meanlog=log(sim.data[,"N"]+1e-6), sdlog=samp[['logsd.N']], log=TRUE))

  llik<-llik.N

  return(llik)
}

r <- debinfer_par(name = "r", var.type = "de", fixed = FALSE,
                  value = 0.5, prior="norm", hypers=list(mean = 0, sd = 1),
                  prop.var=1e-5, samp.type="rw")

K <- debinfer_par(name = "K", var.type = "de", fixed = FALSE,
                  value = 5, prior="lnorm", hypers=list(meanlog = 1, sdlog = 1),
                  prop.var=0.1, samp.type="rw")

logsd.N <- debinfer_par(name = "logsd.N", var.type = "obs", fixed = FALSE,
                        value = 1, prior="lnorm", hypers=list(meanlog = 0, sdlog = 1),
                        prop.var=c(1,2), samp.type="rw-unif")

logsd.N.ref <- debinfer_par(name = "logsd.N", var.type = "obs", fixed = FALSE,
                        value = 1, prior="lnorm", hypers=list(meanlog = 0, sdlog = 1),
                        prop.var=0.001, samp.type="rw-ref")

#we also need to provide an initial condition for the DE
N <- debinfer_par(name = "N", var.type = "init", fixed = TRUE,
                  value = 0.1)


mcmc.pars <- setup_debinfer(r, K, logsd.N, N)
mcmc.pars.ref <- setup_debinfer(r, K, logsd.N.ref, N)
# do inference with deBInfer
# MCMC iterations
iter = 5000
# define burnin
burnin = 2000
# inference call

bm <- microbenchmark::microbenchmark(
rw_unif = de_mcmc(N = iter, data=N_obs, de.model=logistic_model, obs.model=logistic_obs_model, all.params=mcmc.pars,
                        Tmax = max(N_obs$time), data.times=N_obs$time, cnt=iter+1,
                        plot=FALSE, sizestep=0.1, solver=1),
rw_ref = de_mcmc(N = iter, data=N_obs, de.model=logistic_model, obs.model=logistic_obs_model, all.params=mcmc.pars.ref,
                        Tmax = max(N_obs$time), data.times=N_obs$time, cnt=iter+1,
                        plot=FALSE, sizestep=0.1, solver=1),
times=10)

rw_unif <- de_mcmc(N = iter, data=N_obs, de.model=logistic_model, obs.model=logistic_obs_model, all.params=mcmc.pars,
                  Tmax = max(N_obs$time), data.times=N_obs$time, cnt=iter+1,
                  plot=FALSE, sizestep=0.1, solver=1)
rw_ref <- de_mcmc(N = iter, data=N_obs, de.model=logistic_model, obs.model=logistic_obs_model, all.params=mcmc.pars.ref,
                 Tmax = max(N_obs$time), data.times=N_obs$time, cnt=iter+1,
                 plot=FALSE, sizestep=0.1, solver=1)

plot(window(rw_unif$samples, 3300, 5000))
coda::effectiveSize(window(rw_unif$samples, 3300, 5000))
plot(window(rw_ref$samples, 3300, 5000))
coda::effectiveSize(window(rw_ref$samples, 3300, 5000))
