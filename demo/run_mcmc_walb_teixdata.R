## this script runs mcmc inference on the actual data used by Teixeira et al in their paper
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



#set sizestep for simulator
ss<-0.5

#load data from Teixeira paper
lit.table <- read.csv('../walb_deb/data/teixeira2014tab02.csv', skip=4)
#remove last row, because the data contains a typo
lit.table <- lit.table[-nrow(lit.table),]
#make list fit for mcmc
lit.data <-list(t = lit.table$Age.d, L = lit.table$Culmen.length.cm, Ww = lit.table$Mass.kg*1000 )

#plot data
plot(lit.data$t,lit.data$L)
plot(lit.data$t,lit.data$Ww)

#extract Tmax
Tmax <- max(lit.data$t)
abline(v=Tmax)


#solve model for published pars
teix.fit<-solve.DEB(DEB.walb, params, inits, Tmax=Tmax, numsteps=NULL,
                    which=1, sizestep=ss, verbose = FALSE, data.times = lit.data$t)
plot(teix.fit)
plot(teix.fit, obs=data.frame(time=lit.data$t, L=lit.data$L*params['delta_M']))
## p.start is a list giving the initial guesses for the params we want
## to infer since we're not currently proposing X.h, the initial
## condition for e0 is set, so we don't need to reset the "inits"
## inits<-setinits.DEB()
hyper<-make.hypers() #make the prior slightly more vague than the std. error on the fit suggests
w.p<-c("f_slope") #name the parameters that are to be estimated??

p.start<-c(-0.0015) #initial values of parameters
prop.sd<-c(f_slope=0.00003)#what is this? Metropolis-Hastings Tuning parameter?!


sds<-list(L=0.5, Ww=250)

N<-5000

lit.samps<-deb.mcmc(N=N, p.start=p.start, data=lit.data, w.p=w.p, params=params, inits=inits, sim=DEB.walb, sds=sds, hyper=hyper, prop.sd=prop.sd, Tmax=Tmax, cnt=10, burnin=0, plot=TRUE, sizestep=ss, which = 1, data.times = lit.data$t)
##out<-mcmc(N=N, p.start=p.start, data, params, inits, sim=DEB1, sds, Tmax, burnin=0, cnt=50)

lit.samps<-lit.samps$samps
#pdf("figs/chain_post_prior.pdf", width=12, height=8)
par(mfrow=c(1,2))
plot(lit.samps$f_slope,type='l', main="mcmc trace", xlab = 'iteration')
#hist(samps$f_slope, freq=FALSE)
plot(density(lit.samps$f_slope),xlim=c(-.01,0), main='f_slope')
lines(seq(-.01,0,length.out = 100),dnorm(seq(-0.01,0,length.out = 100), mean=hyper$f_slope[1], sd=hyper$f_slope[2]),col="red")
#legend("topleft",legend=c("prior PDF", "posterior KDE"), lty=1, col=c('red','black'))
#plot(samps$f_intercept,type='l', main="mcmc trace", xlab = 'iteration')
#plot(density(samps$f_intercept), main='f_intercept')
#lines(seq(0,4,length.out = 100),dnorm(seq(0,4,length.out = 100), mean=hyper$f_intercept[1], sd=hyper$f_intercept[2]),col="red")
#dev.off()

#solve with estimated parameters
new.params <- params
new.params['f_slope'] <- mean(lit.samps$f_slope)
#new.params['f_intercept'] <- mean(samps$f_intercept[2000:5000])

new.data<-solve.DEB(DEB.walb, new.params, inits, Tmax=Tmax, numsteps=NULL,
                            which=1, sizestep=ss, verbose = FALSE)

#plot(old.data)#plot debtool fit, which is the basis for the simulated data
plot(new.data)
plot(new.data, obs=data.frame(time=lit.data$t, L=lit.data$L*params['delta_M']))

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
teix.Ww <- get_Ww(sim.data = teix.fit, ode.pars=params)

#pdf("figs/post_pred.pdf", width=12, height=8)
par(mfrow=c(1,2))
plot(lit.data$t, lit.data$Ww,type='p',lwd=3,xlab = 'time (days since hatching)', ylab = 'wet weight (g)')
lines(new.data[,"time"],new.Ww, col='red',lty=3,lwd=4)
lines(teix.fit[,"time"], teix.Ww)
legend("topleft",lwd=c(3,4,NA), lty=c(1,3,NA), pch=c(NA,NA,1), col=c('black','red','grey'),legend=c('Teixeira', 'MCMC fit','data'),bty='n')
plot(lit.data$t, lit.data$L,type='p',lwd=3,xlab = 'time (days since hatching)', ylab = 'wet weight (g)')
lines(new.data[,"time"],new.data[,"L"]/params["delta_M"], col='red',lty=3,lwd=4)
lines(teix.fit[,"time"], teix.fit[,"L"]/params["delta_M"])


#dev.off()
