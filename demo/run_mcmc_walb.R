## here's code to run a sample of our nifty mcmc
## run with: nohup R CMD batch run_mcmc.R &



source("R/deb_mcmc.R")
source("R/deb_solver.R")
source("R/DEB_walb.r")

## need to source this one after deb_mcmc.R to tell it to use this
## versio of the prior/posterior calcs instead of the versions in
## deb_post_prior.R
source("R/deb_post_prior_walb.R")

##library("coda")

params<-setparams.DEB.walb()
inits<-setinits.DEB.walb()

Tmax<-500
sds<-list(L=0.5, Ww=250)


ss<-0.5
system.time(data<-solve.DEB(DEB.walb, params, inits, Tmax=Tmax, numsteps=NULL,
                which=1, sizestep=ss, verbose = FALSE))

#data is now a deSolve object
plot(data)
#plot.DEB(data, scaled.length=FALSE, scale=100) # breaks with walb model
data<-make.obs(data, sds, params = c( shape = 0.721, gamma = 630), w.t=1, Tmax=Tmax, ode.pars=setparams.DEB.walb())#params are simply copied from daphnia model. doesn't make sense for walb
#data is now a list of 3, t, L, Ww
#plot.DEB.red(data)# works but not sensible for walb DEB
plot(data$t,data$L)
plot(data$t,data$Ww)
## p.start is a list giving the initial guesses for the params we want
## to infer since we're not currently proposing X.h, the initial
## condition for e0 is set, so we don't need to reset the "inits"
## inits<-setinits.DEB()
hyper<-make.hypers()
w.p<-c("f_slope", "f_intercept") #name the parameters that are to be estimated??

p.start<-c(0.5, 0.8) #what is this?
prop.sd<-c(f_slope=0.004, f_intercept=0.004)#what is this?

N<-200

samps<-deb.mcmc(N, p.start, data, w.p, params, inits, sim=DEB.walb, sds, hyper, prop.sd, Tmax, cnt=10, burnin=0, plot=TRUE, sizestep=ss)
##out<-mcmc(N=N, p.start=p.start, data, params, inits, sim=DEB1, sds, Tmax, burnin=0, cnt=50)

samps<-samps$samps


save(samps, file="test_samples.Rsave")

