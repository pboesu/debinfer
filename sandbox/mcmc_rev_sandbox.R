#sandbox for building up the new mcmc function
## ----install1, eval=FALSE------------------------------------------------
## install.packages("devtools")

## ----install1b-----------------------------------------------------------
#Load the devtools package.
library(devtools)

## ----install2------------------------------------------------------------
#install_github("pboesu/debinfer")
#reload deBinfer to ensure most current version is loaded
if ("package:deBInfer" %in% search()) detach("package:deBInfer", unload=TRUE)
library(deBInfer)

## ----ode-def, message=FALSE,warning=FALSE--------------------------------
library(deSolve)
logistic_model <- function (time, y, parms) {
  with(as.list(c(y, parms)), {#maybe not do the list coercion, and also the with() statement to save time?
    dN <- r * N * (1 - N / K)
    list(dN)
  })
}
y <- c(N = 0.1)
parms <- c(r = 0.1, K = 10)
times <- seq(0, 120, 1)
out <- ode(y, times, logistic_model, parms, method='lsoda')

## ---- echo=FALSE---------------------------------------------------------
plot(out)

## ------------------------------------------------------------------------
set.seed(143)
N_obs <- as.data.frame(out[c(1,runif(35, 0, nrow(out))),]) #force include the first time-point (t=0)


## ------------------------------------------------------------------------
# add lognormal noise
parms['loglogsd.N'] <- -4.6
N_obs$N_noisy <- rlnorm(nrow(N_obs), log(N_obs$N),exp(parms['loglogsd.N']))
#observations must be ordered for solver to work
N_obs <- N_obs[order(N_obs$time),]

plot(N_obs$time, N_obs$N, ylim=c(0, max(N_obs$N,N_obs$N_noisy)))
points(N_obs$time, N_obs$N_noisy, col="red")

out_obs <- ode(y, c(0,N_obs$time), logistic_model, parms, method='lsoda')
plot(out_obs)
lines(out)

## ----obs-model-----------------------------------------------------------
# the observation model
logistic_obs_model<-function(data, sim.data, sds, samp){
  llik.N<-sum(dlnorm(data$N_noisy, meanlog=log(sim.data$N + 1e-6), #here sim.data is the output of make.states, i.e. a list
                     sdlog=exp(samp[['loglogsd.N']]), log=TRUE)
  )

  llik<-llik.N

  return(llik)
}
# the observation model
logistic_obs_model_rev<-function(data, sim.data, sds, samp){

  llik.N<-sum(dlnorm(data$N_noisy, meanlog=log(sim.data[,"N"] + 1e-6), #here sim.data is the output from ode/dde i.e. a dataframe or matrix
                     sdlog=exp(samp[['loglogsd.N']]), log=TRUE)
  )

  llik<-llik.N

  return(llik)
}


## ----pars, results="hide", message=FALSE---------------------------------
library(deBInfer)
r <- debinfer_par(name = "r", var.type = "de", fixed = FALSE,
                  value = 0.5, prior="norm", hypers=list(mean = 0, sd = 1),
                  prop.var=0.00003, samp.type="rw")

K <- debinfer_par(name = "K", var.type = "de", fixed = FALSE,
                  value = 5, prior="lnorm", hypers=list(meanlog = 1, sdlog = 1),
                  prop.var=0.01, samp.type="rw")

loglogsd.N <- debinfer_par(name = "loglogsd.N", var.type = "obs", fixed = FALSE,
                           value = -2, prior="norm", hypers=list(mean = 0, sd = 1),
                           prop.var=0.5, samp.type="rw")


## ----inits---------------------------------------------------------------
N <- debinfer_par(name = "N", var.type = "init", fixed = TRUE, value = 0.1)

## ----setup---------------------------------------------------------------
mcmc.pars <- setup_debinfer(r, K, loglogsd.N, N)

## ----deBinfer, results="hide"--------------------------------------------
# do inference with deBInfer
# MCMC iterations
iter = 50

#original mcmc function

#library(profvis)
#library(microbenchmark)
ode_noplot <- microbenchmark::microbenchmark(
 mcmc_samples <- de_mcmc(N = iter, data=N_obs, de.model=logistic_model,
                         obs.model=logistic_obs_model, all.params=mcmc.pars,
                         Tmax = max(N_obs$time), data.times=N_obs$time, cnt=iter %/% 10,
                         burnin=0.1, plot=FALSE, sizestep=0.1, which=1),

#revised mcmc function
mcmc_rev <- de_mcmc_rev(N = iter, data=N_obs, de.model=logistic_model,
                        obs.model=logistic_obs_model_rev, all.params=mcmc.pars,
                        Tmax = max(N_obs$time), data.times=N_obs$time, cnt=iter %/% 10,
                        burnin=0.1, plot=FALSE, sizestep=0.1, which=1,
                        ref.params = parms, ref.inits = y),
times = 10)

ode_plot <- microbenchmark::microbenchmark(
  mcmc_samples <- de_mcmc(N = iter, data=N_obs, de.model=logistic_model,
                          obs.model=logistic_obs_model, all.params=mcmc.pars,
                          Tmax = max(N_obs$time), data.times=N_obs$time, cnt=iter %/% 10,
                          burnin=0.1, plot=TRUE, sizestep=0.1, which=1),

  #revised mcmc function
  mcmc_rev <- de_mcmc_rev(N = iter, data=N_obs, de.model=logistic_model,
                          obs.model=logistic_obs_model_rev, all.params=mcmc.pars,
                          Tmax = max(N_obs$time), data.times=N_obs$time, cnt=iter %/% 10,
                          burnin=0.1, plot=TRUE, sizestep=0.1, which=1,
                          ref.params = parms, ref.inits = y),
  times = 10)

