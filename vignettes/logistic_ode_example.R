## ----install1, eval=FALSE------------------------------------------------
#  install.packages("devtools")

## ----install1b-----------------------------------------------------------
#Load the devtools package.
library(devtools)

## ----install2, eval=FALSE------------------------------------------------
#  install_github("pboesu/debinfer")

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
  logistic_obs_model<-function(data, sim.data, samp){

    llik.N<-sum(dlnorm(data$N_noisy, meanlog=log(sim.data[,"N"] + 1e-6),
                       sdlog=(samp[['loglogsd.N']]), log=TRUE)
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
                value = 1, prior="lnorm", hypers=list(meanlog = 1, sdlog = 1),
                prop.var=0.1, samp.type="rw")

loglogsd.N <- debinfer_par(name = "loglogsd.N", var.type = "obs", fixed = FALSE,
                value = 0.005, prior="lnorm", hypers=list(mean = 1, sd = 1),
                prop.var=0.1, samp.type="rw")

## ----inits---------------------------------------------------------------
N <- debinfer_par(name = "N", var.type = "init", fixed = TRUE, value = 0.1)

## ----setup---------------------------------------------------------------
mcmc.pars <- setup_debinfer(r, K, loglogsd.N, N)

## ----deBinfer, results="hide"--------------------------------------------
# do inference with deBInfer
  # MCMC iterations
  iter = 5000
  # inference call
  mcmc_samples <- de_mcmc(N = iter, data=N_obs, de.model=logistic_model,
                          obs.model=logistic_obs_model, all.params=mcmc.pars,
                          Tmax = max(N_obs$time), data.times=N_obs$time, cnt=100,
                          plot=TRUE, sizestep=0.1, which=1, verbose=TRUE)




## ----message=FALSE, warning=FALSE,fig.width = 8, fig.height = 8----------
burnin = 1500
pairs(mcmc_samples, burnin = burnin, scatter=TRUE, trend=TRUE)

plot(mcmc_samples$samples)
summary(mcmc_samples$samples)
#plot(N_obs$time,N_obs$N_noisy)

posterior.mean.sim <- solve_de(logistic_model, params=c(r=mean(mcmc_samples$samples[,"r"][burnin:iter]),K=mean(mcmc_samples$samples[,"K"][burnin:iter])), Tmax=max(N_obs$time), inits=c(N=0.1))
plot(posterior.mean.sim)
points(N_obs$time,N_obs$N_noisy)

## ------------------------------------------------------------------------
library(plyr)
simlist <- llply(sample(nrow(mcmc_samples$samples[burnin:iter,]), 1000),
                 function(x)solve_de(logistic_model, mcmc_samples$samples[burnin:iter,][x,],
                                     inits=c(N = 0.1), Tmax=120, numsteps=100))

sima <- simplify2array(simlist)
dim(sima)

simCI <- aaply(sima, .margins = c(1,2), quantile, probs=c(0.025,0.5,0.975))
for (p in 2:ncol(simlist[[1]])){
  plot(simCI[,,3][,c(1,p)], main="", type='n', ylim=c(0,11))

  #for (i in 1:length(simlist)){
  #  lines(simlist[[i]][,c(1,p)],lty=2,col='grey')
  #}
  for (i in 1:3){
    lines(simCI[,,i][,c(1,p)], lty=c(2,1,2)[i], col=c("red"))
  }
  if (dimnames(simlist[[1]])[[2]][p] == "N") points(N_obs$time,N_obs$N_noisy)
}

## ------------------------------------------------------------------------
#str(out)
#lines(posterior.mean.sim, ylim=c(0,15))

#l_ply(simlist, function(x)  lines(x[,1], x[,2]))
#lines(out[,1],out[,2], col='blue')

## ----pars2---------------------------------------------------------------
#we could also use  asymmetric uniform proposals to ensure only positive values of logsd.N are sampled

logsd.Nunif <-  debinfer_par(name = "logsd.N", var.type = "obs", fixed = FALSE,
                value = 0.01, prior="norm", hypers=list(mean = 0, sd = 1),
                prop.var=c(1,2), samp.type="rw-unif")

#in that case we need to adjust the observation model slightly

# the observation model
  logistic_obs_model_unif<-function(data, sim.data, samp){

    llik.N<-sum(dlnorm(data$N_noisy, meanlog=log(sim.data[,"N"] + 1e-6),
                       sdlog=samp[['logsd.N']], log=TRUE)
                )

    llik<-llik.N

    return(llik)
  }

  mcmc.pars_unif <- setup_debinfer(r, K, logsd.Nunif, N)

  mcmc_samples_unif <- de_mcmc(N = iter, data=N_obs, de.model=logistic_model,
                          obs.model=logistic_obs_model_unif, all.params=mcmc.pars_unif,
                          Tmax = max(N_obs$time), data.times=N_obs$time, cnt=iter %/% 10,
                          plot=TRUE, sizestep=0.1, which=1)


