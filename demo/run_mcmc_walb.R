## this code estimates the f_n slope paramter (successfully) and f_n intercept (unsuccessfully, b/c the inits are not recalculated inside the macmc sampler) from simulated data
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
pdf("figs/chain_post_prior.pdf", width=12, height=8)
par(mfrow=c(2,2))
plot(samps$f_slope,type='l', main="mcmc trace", xlab = 'iteration')
#hist(samps$f_slope, freq=FALSE)
plot(density(samps$f_slope),xlim=c(-.01,0), main='f_slope')
lines(seq(-.01,0,length.out = 100),dnorm(seq(-0.01,0,length.out = 100), mean=hyper$f_slope[1], sd=hyper$f_slope[2]),col="red")
legend("topleft",legend=c("prior PDF", "posterior KDE"), lty=1, col=c('red','black'))
plot(samps$f_intercept,type='l', main="mcmc trace", xlab = 'iteration')
plot(density(samps$f_intercept), main='f_intercept')
lines(seq(0,4,length.out = 100),dnorm(seq(0,4,length.out = 100), mean=hyper$f_intercept[1], sd=hyper$f_intercept[2]),col="red")
dev.off()

#solve with estimated parameters
new.params <- params
new.params['f_slope'] <- mean(samps$f_slope[2000:5000])
new.params['f_intercept'] <- mean(samps$f_intercept[2000:5000])

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

pdf("figs/post_pred.pdf", width=12, height=8)
par(mfrow=c(1,1))
plot(old.data[,"time"],old.Ww,type='l',lwd=3,xlab = 'time (days since hatching)', ylab = 'wet weight (g)')
lines(new.data[,"time"],new.Ww, col='red',lty=3,lwd=4)
points(data$Ww, cex=1,col='grey')
legend("topleft",lwd=c(3,4,NA), lty=c(1,3,NA), pch=c(NA,NA,1), col=c('black','red','grey'),legend=c('"true" model', 'fitted model','simulated data'))
dev.off()
