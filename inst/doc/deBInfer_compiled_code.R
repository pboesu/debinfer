## ----opts, echo=FALSE, results='hide'------------------------------------
knitr::opts_chunk$set(dev='png')

## ----r-model-------------------------------------------------------------
modelR <- function(t, Y, parameters) {
  with (as.list(parameters),{

    dy1 = -k1*Y[1] + k2*Y[2]*Y[3]
    dy3 = k3*Y[2]*Y[2]
    dy2 = -dy1 - dy3

    list(c(dy1, dy2, dy3))
  })
}

## ----r-jacobian----------------------------------------------------------
jacR <- function (t, Y, parameters) {
  with (as.list(parameters),{

    PD[1,1] <- -k1
    PD[1,2] <- k2*Y[3]
    PD[1,3] <- k2*Y[2]
    PD[2,1] <- k1
    PD[2,3] <- -PD[1,3]
    PD[3,2] <- k3*Y[2]
    PD[2,2] <- -PD[1,2] - PD[3,2]

    return(PD)
  })
}

## ----run-r-model---------------------------------------------------------
library(deSolve)
parms <- c(k1 = 0.04, k2 = 1e4, k3=3e7)
Y     <- c(1.0, 0.0, 0.0)
times <- seq(0, 0.75*10^4, length.out = 1000)
PD    <- matrix(nrow = 3, ncol = 3, data = 0)
outR   <- ode(Y, times, modelR, parms = parms, jacfunc = jacR)
plot(outR)

## ----compile-C-----------------------------------------------------------
system("R CMD SHLIB mymod.c")

## ----load-so-------------------------------------------------------------
dyn.load(paste("mymod", .Platform$dynlib.ext, sep = ""))

## ----run-c-model---------------------------------------------------------
parms <- c(k1 = 0.04, k2 = 1e4, k3=3e7)
Y     <- c(y1 = 1.0, y2 = 0.0, y3 = 0.0)
times <- c(0, 0.4*10^(0:11) )

out <- ode(Y, times, func = "derivs", parms = parms,
           jacfunc = "jac", dllname = "mymod",
           initfunc = "initmod", nout = 1, outnames = "Sum")
plot(out)

## ----r-c-speed-----------------------------------------------------------
comp <- microbenchmark::microbenchmark(
  C = ode(Y, times, func = "derivs", parms = parms,
           jacfunc = "jac", dllname = "mymod",
           initfunc = "initmod", nout = 1, outnames = "Sum"),
  R = ode(Y, times, modelR, parms = parms, jacfunc = jacR)
)
print(comp, unit = "ms")

## ------------------------------------------------------------------------
set.seed(143)
#force include the first time-point (t=0)
obs <- as.data.frame(outR[c(1,runif(100, 0, nrow(outR))),]) 
obs <- obs[order(obs$time),]
plot(obs$time, obs[,2], type='o')
plot(obs$time, obs[,3], type='o')
plot(obs$time, obs[,4], type='o')

## ----obs-model-----------------------------------------------------------
  # the observation model
  obs_model <- function(data, sim.data, samp){

    llik.y1 <- sum(dnorm(obs[,2], mean = sim.data[,2], sd = samp[['sd.y1']], log = TRUE))
    llik.y2 <- sum(dnorm(obs[,3], mean = sim.data[,3], sd = samp[['sd.y2']], log = TRUE))
    llik.y3 <- sum(dnorm(obs[,4], mean = sim.data[,4], sd = samp[['sd.y3']], log = TRUE))
    return(llik.y1 + llik.y2 + llik.y3)
  }


## ----pars, results="hide", message=FALSE---------------------------------
library(deBInfer)
k1 <- debinfer_par(name = "k1", var.type = "de", fixed = FALSE,
                value = 0.1, prior = "norm", hypers = list(mean = 0, sd = 1),
                prop.var = 0.00001, samp.type="rw")

k2 <- debinfer_par(name = "k2", var.type = "de", fixed = FALSE,
                value = 5000, prior = "norm", hypers = list(mean = 5000, sd = 500),
                prop.var = 10, samp.type="rw")
k3 <- debinfer_par(name = "k3", var.type = "de", fixed = FALSE,
                value = 1e7, prior = "norm", hypers = list(mean = 1e7, sd = 1e6),
                prop.var = 5e6, samp.type="rw")

sd.y1 <- debinfer_par(name = "sd.y1", var.type = "obs", fixed = FALSE,
                value = 0.05, prior = "lnorm", hypers = list(meanlog = 0, sdlog = 1),
                prop.var = c(3,4), samp.type = "rw-unif")
sd.y2 <- debinfer_par(name = "sd.y2", var.type = "obs", fixed = FALSE,
                value = 0.05, prior = "lnorm", hypers = list(meanlog = 0, sdlog = 1),
                prop.var = c(3,4), samp.type = "rw-unif")
sd.y3 <- debinfer_par(name = "sd.y3", var.type = "obs", fixed = FALSE,
                value = 0.05, prior = "lnorm", hypers = list(meanlog = 0, sdlog = 1),
                prop.var = c(3,4), samp.type = "rw-unif")

## ----inits---------------------------------------------------------------
y1 <- debinfer_par(name = "y1", var.type = "init", fixed = TRUE, value = 1)
y2 <- debinfer_par(name = "y2", var.type = "init", fixed = TRUE, value = 0)
y3 <- debinfer_par(name = "y3", var.type = "init", fixed = TRUE, value = 0)

## ----setup---------------------------------------------------------------
mcmc.pars <- setup_debinfer(k1, k2, k3, sd.y1, sd.y2, sd.y3, y1, y2, y3)

## ----deBinfer------------------------------------------------------------
# do inference with deBInfer
  # MCMC iterations
iter <- 1000
  # inference call
Rt <- system.time(mcmc_samples <- de_mcmc(N = iter, data = obs, de.model = modelR, 
                          obs.model = obs_model, all.params = mcmc.pars,
                          Tmax = max(obs$time), data.times = obs$time, cnt = 500, 
                          plot = FALSE, solver = "ode", verbose.mcmc = FALSE))

Ct <- system.time(mcmc_samplesC <- de_mcmc(N = iter, data = obs, de.model = "derivs",
           jacfunc = "jac", dllname = "mymod",
           initfunc = "initmod", nout = 1, outnames = "Sum",
                          obs.model = obs_model, all.params = mcmc.pars,
                          Tmax = max(obs$time), data.times = obs$time, cnt = 50, 
                          plot = FALSE, solver = "ode", verbose.mcmc = FALSE))
print(Rt)
print(Ct)

## ----post-sim------------------------------------------------------------
system.time(post_traj <- post_sim(mcmc_samples, n=100, times=0:8000, burnin=100, output = 'all', prob = 0.95))
system.time(post_trajC <- post_sim(mcmc_samplesC, n=100, times=0:8000, burnin=100, output = 'all', prob = 0.95, jacfunc = "jac", dllname = "mymod",initfunc = "initmod", nout = 1, outnames = "Sum"))
par(mfrow=c(1,3))
plot(post_traj$time, post_traj$median$y1, type='l',lwd=2)
lines(post_trajC$time, post_trajC$median$y1, col="red", lty=3, lwd=4)
points(obs$time, obs[,2])
plot(post_traj$time, post_traj$median$y2, type='l')
lines(post_trajC$time, post_trajC$median$y2, col="red", lty=3, lwd=4)
points(obs$time, obs[,3])
plot(post_traj$time, post_traj$median$y3, type='l')
lines(post_trajC$time, post_trajC$median$y3, col="red", lty=3, lwd=4)
points(obs$time, obs[,4])


