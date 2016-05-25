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

## ----de-def, message=FALSE,warning=FALSE---------------------------------
#dde version
CSZ.dde<-function(t,y,p){

  sr    <-p["sr"]
  fs    <-p["fs"]
  ds    <-p["ds"]
  eta   <-p["eta"]
  Tmin  <-p["Tmin"]
  muz <-p["muz"]

  Rs<-Ms<-0
  lag1<-lag2<-0

  if (t>Tmin){
    lag1<-pastvalue(t-Tmin)
    Rs <- sr*fs*lag1[1]
  }

  phiZ <- eta*y[2]
  dC <- -(muz+sr)*y[1]
  dS <- Rs - Ms - ds*y[2]
  dZ <- phiZ - (muz+sr)*y[3]

  if(y[1]<0) dC<-0
  if(y[2]<0){
    dS <- Rs - Ms
    dZ <- -(muz+sr)*y[3]
  }
  if(y[3]<0){
    dZ <- dZ+(muz+sr)*y[3]
  }

  list(c(dC,dS,dZ))
}

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
  ##if (t>(Tmin+Tmax)){
  ##  lag2<-pastvalue(t-(Tmin+Tmax))
  ##  Ms <- sr*fs*exp(-ds*Tmax)*lag2[1]
  ##}

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


## ----data, echo=FALSE----------------------------------------------------
library(deBInfer)
#load chytrid data
data(chytrid)
#knitr::kable(chytrid, caption='Replicated counts of chytrid fungus zoospores')
plot(chytrid, xlab='Time (days)', ylab='Zoospores x 10e4', xlim=c(0,10))

## ----obs-model-----------------------------------------------------------
## observation model
chytrid_obs_model<-function(data, sim.data, sds, samp){

  #z.temp<-sim.data$Z

  ##ec<-0.00001
  ec<-0.01
  llik.Z<-0
  for(i in unique(data$time)){
    llik.Z<-llik.Z + sum(dpois(data$count[data$time==i], lambda=(sim.data$Z[sim.data$time==i]+ec), log=TRUE))
  }
  llik<-llik.Z
  return(llik)
}

## ----pars, results="hide", message=FALSE---------------------------------
sr <- debinfer_par(name = "sr", var.type = "de", fixed = FALSE,
                   value = 2, prior="gamma", hypers=list(shape = 5, rate = 1),
                   prop.var=0.5, samp.type="rw")

fs <- debinfer_par(name = "fs", var.type = "de", fixed = FALSE,
                   value = 0.5, prior="beta", hypers=list(shape1 = 1, shape2 = 1),
                   prop.var=0.05, samp.type="ind")

ds <- debinfer_par(name = "ds", var.type = "de", fixed = FALSE,
                   value = 2, prior="gamma", hypers=list(shape = 1, rate = 1),
                   prop.var=0.1, samp.type="rw")

muz <- debinfer_par(name = "muz", var.type = "de", fixed = FALSE,
                    value = 1, prior="gamma", hypers=list(shape = 5, rate = 1),
                    prop.var=0.1, samp.type="rw")

eta <- debinfer_par(name = "eta", var.type = "de", fixed = FALSE,
                    value = 10, prior="gamma", hypers=list(shape = 1, rate = 0.25),
                    prop.var=5, samp.type="rw")

Tmin <- debinfer_par(name = "Tmin", var.type = "de", fixed = FALSE,
                     value = 3, prior="unif", hypers=list(min = 2, max = 6),
                     prop.var=0.2, samp.type="rw")


## ----inits---------------------------------------------------------------
C <- debinfer_par(name = "C", var.type = "init", fixed = TRUE, value = 120)
S <- debinfer_par(name = "S", var.type = "init", fixed = TRUE, value = 0)
Z <- debinfer_par(name = "Z", var.type = "init", fixed = TRUE, value = 0)

## ----setup---------------------------------------------------------------
mcmc.pars <- setup_debinfer(sr, fs, ds, muz, eta, Tmin, C, S, Z)

## ----deBinfer, results="hide"--------------------------------------------
# do inference with deBInfer
# MCMC iterations
iter = 500
# inference call
dde_plot <- microbenchmark::microbenchmark(
dde_old = de_mcmc(N = iter, data=chytrid, de.model=CSZ.dde,
                        obs.model=chytrid_obs_model, all.params=mcmc.pars,
                        Tmax = max(chytrid$time), data.times=c(0,chytrid$time), cnt=50,
                        burnin=0.1, plot=TRUE, sizestep=0.1, which=2, verbose = TRUE),


dde_rev = de_mcmc_rev(N = iter, data=chytrid, de.model=CSZ.dde,
                      obs.model=chytrid_obs_model, all.params=mcmc.pars,
                      Tmax = max(chytrid$time), data.times=c(0,chytrid$time), cnt=50,
                      burnin=0.1, plot=TRUE, sizestep=0.1, which=2, verbose = TRUE),

dede_rev <- de_mcmc_rev(N = iter, data=chytrid, de.model=CSZ.dede,
                               obs.model=chytrid_obs_model, all.params=mcmc.pars,
                               Tmax = max(chytrid$time), data.times=c(0,chytrid$time), cnt=50,
                               burnin=0.1, plot=TRUE, sizestep=0.1, which="dede", verbose = TRUE),
times = 10)

## ----solve-dde-----------------------------------------------------------
pars <- colMeans(mcmc_samples$samps)
names(pars) <- names(mcmc_samples$samps)

voylesfit <- dde(y=c(C=120,S=0,Z=0), times = seq(0,10,by=0.02), func=CSZ.dde, parms=pars )


#Sim CIs
library(coda)
library(plyr)
parsamps <- mcmc_samples$samps[sample(nrow(mcmc_samples$samps), size=500),]



siml <- alply(parsamps, 1, function(x){spars <- as.numeric(x)
names(spars) <- names(parsamps)
dde(y=c(120,0,0), times = seq(0,10,by=0.02), func=CSZ.dde, parms=spars)}, .progress='text')


plot(chytrid, xlim=c(0,9))
l_ply(siml, .fun=function(x){lines(x$time, x$y1, col='red3')})

l_ply(siml, .fun=function(x){lines(x$time, x$y2, col='blue3')})

l_ply(siml, .fun=function(x){lines(x$time, x$y3, col='green3')})

#figure
library(viridis)

#set colors
colours <- viridis(9)[c(1,4,7)]


with(voylesfit,{
  lwd=3
  plot(time,C, type='l', ylim=c(0,300),col=colours[1],lwd=lwd,ylab="Counts x 10^4", xlab="Time (days)")
  lines(time,Z,col=colours[2],lwd=lwd)
  lines(time,S,col=colours[3],lwd=lwd, lty=2)
  legend("topright", legend=c("C","Z","S", "Z_obs"), lwd=3, col=c(colours[c(1,2,3)],"black"), lty=c(1,1,2,NA), pch=c(NA,NA,NA,1))
  points(data_obs)
}
)


