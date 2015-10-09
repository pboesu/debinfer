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

Tmax<-250
sds<-list(L=0.5, Ww=250)


ss<-0.5
#solve ODE using initial paramters
old.data<-solve.DEB(DEB.walb, params, inits, Tmax=Tmax, numsteps=NULL,
                which=1, sizestep=ss, verbose = FALSE)
#keep copy for later comparison
data <- old.data
#data is now a deSolve object
plot(data)
#plot.DEB(data, scaled.length=FALSE, scale=100) # breaks with walb model
data<-make.obs(data, sds, params = c( shape = 0.721, gamma = 630), w.t=1, Tmax=Tmax, ode.pars=setparams.DEB.walb())#params are simply copied from daphnia model. doesn't make sense for walb. breaks for any w.t values != 1
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

p.start<-c(-0.001, 1.06) #what is this?
prop.sd<-c(f_slope=0.003, f_intercept=0.4)#what is this?

N<-5000

samps<-deb.mcmc(N, p.start, data, w.p, params, inits, sim=DEB.walb, sds, hyper, prop.sd, Tmax, cnt=10, burnin=0, plot=TRUE, sizestep=ss,which = 1)
##out<-mcmc(N=N, p.start=p.start, data, params, inits, sim=DEB1, sds, Tmax, burnin=0, cnt=50)

samps<-samps$samps
par(mfrow=c(2,2))
plot(samps$f_slope,type='l')
hist(samps$f_slope, freq=FALSE)
plot(density(samps$f_slope))
plot(samps$f_intercept,type='l')

#solve with estimated parameters
new.params <- params
new.params['f_slope'] <- min(samps$f_slope)
new.params['f_intercept'] <- min(samps$f_intercept)

new.data<-solve.DEB(DEB.walb, new.params, inits, Tmax=Tmax, numsteps=NULL,
                            which=1, sizestep=ss, verbose = FALSE)

plot(old.data)#plot debtool fit, which is the basis for the simulated data
plot(new.data)

get_Ww <- function(w_E = 23.9, d_v = 0.5, mu_E = 550000, sim.data=new.data, ode.pars=new.params ){
#wet weight. this is all hard coded now, should not be!
#w_E = 23.9 # molecular weight of reserve g mol^-1
#d_v = 0.5 # specific density of structure
#mu_E = 550000 # chemical potential of reserve J / mol
wdratio = -1.37e-3 * sim.data[, 'time'] + 2.09
omega = unname(ode.pars['p_Am'] * w_E / (ode.pars['v'] * d_v * mu_E))
Ww = sim.data[,'L']^3 * (1 + sim.data[,'f_n'] * omega) * d_v * wdratio
return(Ww)
}
new.Ww <- get_Ww(w_E = 23.9, d_v = 0.5, mu_E = 550000, sim.data=new.data, ode.pars=new.params )
old.Ww <- get_Ww(sim.data = old.data, ode.pars=params)
plot(old.data[,"time"],old.Ww,type='l')
lines(new.data[,"time"],new.Ww, col='red')
points(data$Ww, cex=0.1,col='grey')
