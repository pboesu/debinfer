## here's code to run a sample of our nifty mcmc
## run with: nohup R CMD batch run_mcmc.R &



source("R/deb_mcmc.R")
source("R/deb_solver.R")
source("R/DEB_daphnia.r")

## need to source this one after deb_mcmc.R to tell it to use this
## versio of the prior/posterior calcs instead of the versions in
## deb_post_prior.R
source("R/deb_post_prior_daph.R")

##library("coda")

params<-setparams.DEB()
inits<-setinits.DEB()
inits[2]<-PHI(params["X.h"], inits[1])
Tmax<-50
sds<-list(L=0.5, Negg=50)


ss<-0.01
data<-solve.DEB(DEB.daphnia, params, inits, Tmax=Tmax, numsteps=NULL,
                which=1, sizestep=ss)

plot.DEB(data, scaled.length=FALSE, scale=100)
data<-make.obs(data, sds, params, w.t=1, Tmax=Tmax)

plot.DEB.red(data)

## p.start is a list giving the initial guesses for the params we want
## to infer since we're not currently proposing X.h, the initial
## condition for e0 is set, so we don't need to reset the "inits"
## inits<-setinits.DEB()
hyper<-make.hypers()
w.p<-c("kap", "g")
p.start<-c(0.5, 0.8)
prop.sd<-c(kap=0.004, g=0.004)

N<-200

samps<-deb.mcmc(N, p.start, data, w.p, params, inits, sim=DEB.daphnia, sds, hyper, prop.sd, Tmax, cnt=10, burnin=0, plot=TRUE, sizestep=ss)
##out<-mcmc(N=N, p.start=p.start, data, params, inits, sim=DEB1, sds, Tmax, burnin=0, cnt=50)

samps<-samps$samps


save(samps, file="test_samples.Rsave")

