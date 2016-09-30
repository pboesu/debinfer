## ----opts, echo=FALSE, results='hide'------------------------------------
knitr::opts_chunk$set(dev='png', dpi=150)

## ----install1, eval=FALSE------------------------------------------------
#  install.packages("devtools")

## ----install1b-----------------------------------------------------------
#Load the devtools package.
library(devtools)

## ----install2, eval=FALSE, results='hide', message=FALSE, warning=FALSE----
#  install_github("pboesu/debinfer")

## ----loadlib, message=FALSE----------------------------------------------
library(deBInfer)

## ----ode-def, message=FALSE,warning=FALSE--------------------------------
library(deSolve)
logistic_model <- function (time, y, parms) {
    with(as.list(c(y, parms)), {
      dN <- r * N * (1 - N / K)
      list(dN)
      })
}

## ----pars-def------------------------------------------------------------
y <- c(N = 0.1)
parms <- c(r = 0.1, K = 10)
times <- seq(0, 120, 1)
out <- ode(y, times, logistic_model, parms, method="lsoda")

## ---- echo=FALSE---------------------------------------------------------
plot(out)

## ------------------------------------------------------------------------
set.seed(143)
#force include the first time-point (t=0)
N_obs <- as.data.frame(out[c(1,runif(35, 0, nrow(out))),]) 


## ------------------------------------------------------------------------
# add lognormal noise
  parms["sdlog.N"] <- 0.05
  N_obs$N_noisy <- rlnorm(nrow(N_obs), log(N_obs$N),parms["sdlog.N"])
#observations must be ordered for solver to work
N_obs <- N_obs[order(N_obs$time),]
#Plot true and noisy observations
plot(N_obs$time, N_obs$N, ylim=c(0, max(N_obs$N,N_obs$N_noisy)))
points(N_obs$time, N_obs$N_noisy, col="red")

## ----obs-model-----------------------------------------------------------
  # the observation model
  logistic_obs_model <- function(data, sim.data, samp){

    llik.N <- sum(dlnorm(data$N_noisy, meanlog = log(sim.data[,"N"] + 1e-6), 
                       sdlog = samp[["sdlog.N"]], log = TRUE))
    return(llik.N)
  }


## ----pars, results="hide", message=FALSE---------------------------------
library(deBInfer)
r <- debinfer_par(name = "r", var.type = "de", fixed = FALSE,
                value = 0.5, prior = "norm", hypers = list(mean = 0, sd = 1),
                prop.var = 0.0001, samp.type="rw")

K <- debinfer_par(name = "K", var.type = "de", fixed = FALSE,
                value = 5, prior = "lnorm", hypers = list(meanlog = 1, sdlog = 1),
                prop.var = 0.1, samp.type = "rw")

sdlog.N <- debinfer_par(name = "sdlog.N", var.type = "obs", fixed = FALSE,
                value = 0.05, prior = "lnorm", hypers = list(meanlog = 0, sdlog = 1),
                prop.var = c(3,4), samp.type = "rw-unif")

## ----inits---------------------------------------------------------------
N <- debinfer_par(name = "N", var.type = "init", fixed = TRUE, value = 0.1)

## ----setup---------------------------------------------------------------
mcmc.pars <- setup_debinfer(r, K, sdlog.N, N)

## ----deBinfer, results="hide", message=FALSE-----------------------------
# do inference with deBInfer
  # MCMC iterations
  iter <- 5000
  # inference call
  mcmc_samples <- de_mcmc(N = iter, data = N_obs, de.model = logistic_model, 
                          obs.model = logistic_obs_model, all.params = mcmc.pars,
                          Tmax = max(N_obs$time), data.times = N_obs$time, cnt = 500, 
                          plot = FALSE, verbose.mcmc = FALSE, solver = "ode")


## ----message=FALSE, warning=FALSE,fig.width = 8, fig.height = 8----------
plot(mcmc_samples)

## ------------------------------------------------------------------------
burnin <- 250
pairs(mcmc_samples, burnin = burnin, scatter=TRUE, trend=TRUE)

## ----post-dens, fig.height=5---------------------------------------------
par(mfrow = c(1,3))
post_prior_densplot(mcmc_samples, burnin = burnin, param = "r")
abline(v = 0.1, col = "red", lty = 2)
post_prior_densplot(mcmc_samples, burnin = burnin, param = "K")
abline(v = 10, col = "red", lty = 2)
post_prior_densplot(mcmc_samples, burnin = burnin, param = "sdlog.N")
abline(v = 0.05, col = "red", lty = 2)

## ----post-sims-----------------------------------------------------------
post_traj <- post_sim(mcmc_samples, n=500, times=0:100, burnin=burnin, output = "all", prob = 0.95)

## ----post-sims-plot------------------------------------------------------
#median and HDI
plot(post_traj, plot.type = "medianHDI", lty = c(2,1), lwd = 3, col = c("red","grey20"),
  panel.first = lines(out, col = "darkblue", lwd = 2)
  )
  legend("topleft", legend = c("posterior median", "95% HDI", "true model"),
  lty = c(2,1,1), lwd = c(3,2,2), col = c("red","grey20","darkblue"),
  bty = "n")

## ----post-sims-ensemble--------------------------------------------------
plot(post_traj, plot.type = "ensemble", col = "#FF000040")

## ----epsilon-sims, eval=FALSE--------------------------------------------
#  #reformulate the observation model such that epsilon is an observation parameter
#  logistic_obs_model_eps <- function(data, sim.data, samp){
#  
#    llik.N <- sum(dlnorm(data$N_noisy, meanlog = log(sim.data[,"N"] + samp[["epsilon"]]),
#                         sdlog = samp[["sdlog.N"]], log = TRUE))
#    return(llik.N)
#  }
#  
#  #declare epsilon
#  epsilon <- debinfer_par(name = "epsilon", var.type = "obs", fixed = TRUE,
#                            value = 1e-6)
#  
#  #set up MCMC parameters
#  mcmc.pars <- setup_debinfer(r, K, sdlog.N, epsilon, N)

## ----epsilon-sims-ctd, eval=FALSE----------------------------------------
#  # define a range of epsilon values
#  eps_range <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)
#  
#  # conduct inference for each value of eps_range
#  eps_sims <- lapply(eps_range, function(x){
#    #reassign value of epsilon
#    mcmc.pars$epsilon$value <- x
#    eps_samples <- de_mcmc(N = 10000, data = N_obs, de.model = logistic_model,
#                        obs.model = logistic_obs_model_eps, all.params = mcmc.pars,
#                        Tmax = max(N_obs$time), data.times = N_obs$time, cnt = 500,
#                        plot = FALSE, solver = "ode")
#    return(eps_samples)
#  })
#  
#  # add names to the inference runs
#  names(eps_sims) <- as.character(eps_range)
#  
#  # collate the samples in a data.frame for easy plotting while omitting
#  # a burnin of 1000 samples
#  eps_samps <- plyr::ldply(eps_sims, function(x)(x$samples[1000:10000,]), .id = "epsilon")
#  
#  # plot sensitivity analysis
#  par(mfrow = c(1,3), mar = c(3,3,2,2), mgp = c(1.7,0.5,0))
#  beanplot::beanplot(eps_samps$r ~ eps_samps$epsilon, log = "", what = c(0,1,1,0),
#                     col = "darkgrey", bw = 2e-4, maxwidth = 0.6, xlab = "epsilon",
#                     ylab = "r")
#  legend("bottomleft", legend = c("posterior mean", "true value"), lty = c(1,3),
#         col = c("black", "red"), lwd = 2, bty = "n")
#  abline(h = 0.1, col = "red", lwd = 2, lty = 3)
#  beanplot::beanplot(eps_samps$K ~ eps_samps$epsilon, log = "", what = c(0,1,1,0),
#                     col = "darkgrey", bw = 2e-2, maxwidth = 0.6, xlab = "epsilon"
#                     , ylab = "K")
#  abline(h = 10, col = "red", lwd = 2, lty = 3)
#  beanplot::beanplot(eps_samps$sdlog.N ~ eps_samps$epsilon, log="", what=c(0,1,1,0),
#                     col = "darkgrey", bw = 2e-3, maxwidth = 0.6, xlab = "epsilon",
#                     ylab = "sdlog.N")
#  abline(h = 0.05, col = "red", lwd = 2, lty = 3)

