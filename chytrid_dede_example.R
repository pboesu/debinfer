## ----install1, eval=FALSE------------------------------------------------
## install.packages("devtools")

## ----install1b-----------------------------------------------------------
#Load the devtools package.
library(devtools)

## ----install2, eval=FALSE, results='hide', message=FALSE, warning=FALSE----
## install_github("pboesu/debinfer")

## ----loadlib, message=FALSE----------------------------------------------
library(deBInfer)

## ----dde-def-------------------------------------------------------------
#dede version
CSZ.dede<-function(t,y,p){

  sr    <-p["sr"]
  fs    <-p["fs"]
  ds    <-p["ds"]
  eta   <-p["eta"]
  Tmin  <-p["Tmin"]
  ##Tmax  <-p["Tmax"]
  muz <-p["muz"]

  Rs<-Ms<-0
  lag1<-lag2<-0

  if (t>Tmin){
    lag1<-lagvalue(t-Tmin)
    Rs <- sr*fs*lag1[1]
  }

  phiZ <- eta*y[2]
  dy1 <- -(muz+sr)*y[1]
  dy2 <- Rs - Ms - ds*y[2]
  dy3 <- phiZ - (muz+sr)*y[3]

  if(y[1]<0) dy1<-0
  if(y[2]<0){
    dy2 <- Rs - Ms
    dy3 <- -(muz+sr)*y[3]
  }
  if(y[3]<0){
    dy3 <- dy3+(muz+sr)*y[3]
  }

  list(c(dy1,dy2,dy3))
}

## ----data----------------------------------------------------------------
#load chytrid data
data(chytrid)
#have a look at the variables
head(chytrid)
#plot the data
plot(chytrid, xlab='Time (days)', ylab='Zoospores x 10e4', xlim=c(0,10))

## ----obs-model-----------------------------------------------------------
# observation model
chytrid_obs_model<-function(data, sim.data, samp){

  ec<-0.01
  llik.Z<-0
  for(i in unique(data$time)){
    try(llik.Z<-llik.Z + sum(dpois(data$count[data$time==i], lambda=(sim.data[,'Z'][sim.data[,'time']==i]+ec), log=TRUE)))
  }
  llik<-llik.Z
  return(llik)
}

## ----vars----------------------------------------------------------------
sr <- debinfer_par(name = "sr", var.type = "de", fixed = FALSE,
                   value = 2, prior="gamma", hypers=list(shape = 5, rate = 1),
                   prop.var=c(3,4), samp.type="rw-unif")

fs <- debinfer_par(name = "fs", var.type = "de", fixed = FALSE,
                   value = 0.5, prior="beta", hypers=list(shape1 = 1, shape2 = 1),
                   prop.var=0.01, samp.type="ind")

ds <- debinfer_par(name = "ds", var.type = "de", fixed = FALSE,
                   value = 2, prior="gamma", hypers=list(shape = 1, rate = 1),
                   prop.var=0.1, samp.type="rw")

muz <- debinfer_par(name = "muz", var.type = "de", fixed = FALSE,
                    value = 1, prior="gamma", hypers=list(shape = 5, rate = 1),
                    prop.var=c(4,5), samp.type="rw-unif")

eta <- debinfer_par(name = "eta", var.type = "de", fixed = FALSE,
                    value = 10, prior="gamma", hypers=list(shape = 1, rate = 0.25),
                    prop.var=5, samp.type="rw")

Tmin <- debinfer_par(name = "Tmin", var.type = "de", fixed = FALSE,
                     value = 3, prior="unif", hypers=list(min = 2, max = 6),
                     prop.var=0.05, samp.type="rw")


# ----inits---------------------------------------------------------------
C <- debinfer_par(name = "C", var.type = "init", fixed = TRUE, value = 120)
S <- debinfer_par(name = "S", var.type = "init", fixed = TRUE, value = 0)
Z <- debinfer_par(name = "Z", var.type = "init", fixed = TRUE, value = 0)

## ----mcmc-setup----------------------------------------------------------
# ----setup---------------------------------------------------------------
mcmc.pars <- setup_debinfer(sr, fs, ds, muz, eta, Tmin, C, S, Z)

## ----de_mcmc-------------------------------------------------------------
## ----deBinfer, results="hide"--------------------------------------------
# do inference with deBInfer
# MCMC iterations
iter = 5000
# inference call
dede_rev <- de_mcmc(N = iter, data=chytrid, de.model=CSZ.dede,
                               obs.model=chytrid_obs_model, all.params=mcmc.pars,
                               Tmax = max(chytrid$time), data.times=c(0,chytrid$time), cnt=iter %/% 10,
                               plot=FALSE, sizestep=0.1, solver="dede", verbose = TRUE)

## ----message=FALSE, warning=FALSE,fig.width = 8, fig.height = 8----------
plot(dede_rev, ask=FALSE)

## ------------------------------------------------------------------------
burnin = 1
pairs(dede_rev, burnin = burnin, scatter=TRUE, trend=TRUE)
post_prior_densplot(dede_rev, burnin = burnin)
summary(dede_rev)

## ----post-sims-----------------------------------------------------------
post_traj <- post_sim(dede_rev, n=20, times=0:100, burnin=burnin, output = 'all')

## ----post-sims-plot------------------------------------------------------
#median and HDI
par(mfrow=c(1,3))
plot(post_traj, plot.type = "medianHDI", auto.layout = FALSE)
legend("topright", legend=c("posterior median", "95% HDI"), lty=1, col=c("red","grey"))


## ----post-sims-ensemble--------------------------------------------------
plot(post_traj, plot.type = "ensemble", col = "#FF000040")

