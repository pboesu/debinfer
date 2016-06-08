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
  parms['sdlog.N'] <- 0.05
  N_obs$N_noisy <- rlnorm(nrow(N_obs), log(N_obs$N),parms['sdlog.N'])
#observations must be ordered for solver to work
N_obs <- N_obs[order(N_obs$time),]

plot(N_obs$time, N_obs$N, ylim=c(0, max(N_obs$N,N_obs$N_noisy)))
points(N_obs$time, N_obs$N_noisy, col="red")

out_obs <- ode(y, c(0,N_obs$time), logistic_model, parms, method='lsoda')
plot(out_obs, type="p")
lines(out)

## ----obs-model-----------------------------------------------------------
  # the observation model
  logistic_obs_model<-function(data, sim.data, samp){

    llik.N<-sum(dlnorm(data$N_noisy, meanlog=log(sim.data[,"N"] + 1e-6), 
                       sdlog=samp[['sdlog.N']], log=TRUE)
                )

    llik<-llik.N

    return(llik)
  }


## ----pars, results="hide", message=FALSE---------------------------------
library(deBInfer)
r <- debinfer_par(name = "r", var.type = "de", fixed = FALSE,
                value = 0.5, prior="norm", hypers=list(mean = 0, sd = 1),
                prop.var=0.005, samp.type="rw")

K <- debinfer_par(name = "K", var.type = "de", fixed = FALSE,
                value = 5, prior="lnorm", hypers=list(meanlog = 1, sdlog = 1),
                prop.var=0.1, samp.type="rw")

sdlog.N <- debinfer_par(name = "sdlog.N", var.type = "obs", fixed = FALSE,
                value = 0.1, prior="lnorm", hypers=list(meanlog = 0, sdlog = 1),
                prop.var=c(3,4), samp.type="rw-unif")

## ----inits---------------------------------------------------------------
N <- debinfer_par(name = "N", var.type = "init", fixed = TRUE, value = 0.1)

## ----setup---------------------------------------------------------------
mcmc.pars <- setup_debinfer(r, K, sdlog.N, N)

## ----deBinfer, results="hide"--------------------------------------------
# do inference with deBInfer
  # MCMC iterations
  iter = 5000
  # inference call
  mcmc_samples <- de_mcmc(N = iter, data=N_obs, de.model=logistic_model, 
                          obs.model=logistic_obs_model, all.params=mcmc.pars,
                          Tmax = max(N_obs$time), data.times=N_obs$time, cnt=iter+1, 
                          plot=FALSE, solver="ode")


## ----message=FALSE, warning=FALSE,fig.width = 8, fig.height = 8----------
plot(mcmc_samples)

## ------------------------------------------------------------------------
burnin = 200
pairs(mcmc_samples, burnin = burnin, scatter=TRUE, trend=TRUE)
post_prior_densplot(mcmc_samples, burnin = burnin)
summary(mcmc_samples)

## ----post-sims-----------------------------------------------------------
post_traj <- post_sim(mcmc_samples, n=200, times=0:100, burnin=burnin, output = 'all')

## ----post-sims-plot------------------------------------------------------
#median and HDI
plot(post_traj, plot.type = "medianHDI")
legend("topleft", legend=c("posterior median", "95% HDI"), lty=1, col=c("red","grey"))


## ----post-sims-ensemble--------------------------------------------------
plot(post_traj, plot.type = "ensemble", col = "#FF000040")

